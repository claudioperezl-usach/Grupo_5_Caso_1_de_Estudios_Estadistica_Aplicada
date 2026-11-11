# =============================================================================
# 16_intervalos_predictivos.R — Intervalos de predicción (cuantil empírico de residuos LOSO)
#                               y política de inspección
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/16_intervalos_predictivos.R")
#
# CALIBRACIÓN: residuos absolutos de las predicciones LOSO fuera de partición del
#   entrenamiento (paso 12, con los hiperparámetros finales). Cada residuo proviene de
#   una velocidad que el modelo no vio, igual que en el uso previsto. Los test NO se usan
#   para calibrar (versión anterior: calibraba con el test aleatorio; corregido).
# EVALUACIÓN: cobertura y ancho en los dos test, solo como verificación.
#   Cuantil: q = residuo ordenado en la posición ceil((n + 1) * nivel).
#   Relación con la materia: es un INTERVALO DE PREDICCIÓN (dónde puede caer y para un estado
#   nuevo), no un intervalo de confianza de la media. predict(lm, interval = "prediction")
#   no existe para Elastic-net, así que el ancho se estima con el cuantil empírico de los
#   errores de validación (LOSO), que ya incluyen el error de una velocidad no vista.
#   La garantía formal exige intercambiabilidad; entre velocidades es aproximada,
#   por eso la cobertura se verifica y se reporta.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto(character(0))
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 16: intervalos de predicción con el cuantil empírico de los residuos LOSO ===")

ruta_cal  <- file.path(DIR_TABLAS, "12_predicciones_loso.csv")
ruta_eval <- file.path(DIR_TABLAS, "12_predicciones_test.csv")
if (!file.exists(ruta_cal) || !file.exists(ruta_eval)) stop("Ejecute primero el paso 12.")
cal <- data.table::fread(ruta_cal)
eva <- data.table::fread(ruta_eval)
verificar(nrow(cal) > 0 && !anyNA(cal$pred) && !anyDuplicated(cal$id),
          sprintf("calibración con %d residuos LOSO del entrenamiento (uno por fila, sin faltantes)", nrow(cal)))
verificar(length(intersect(cal$id, eva$id)) == 0, "ningún residuo de calibración pertenece a los test")

res_cal <- abs(cal$y - cal$pred)
n <- length(res_cal)
rango <- c(0, 5)   # la pérdida y está acotada a [0, 5] puntos % por construcción
resumen <- do.call(rbind, lapply(cfg$niveles_intervalo, function(nivel) {
  k <- ceiling((n + 1) * nivel)
  q <- if (k > n) Inf else sort(res_cal)[k]
  do.call(rbind, lapply(split(eva, eva$conjunto), function(d) {
    lo <- pmax(d$pred - q, rango[1]); hi <- pmin(d$pred + q, rango[2])
    data.frame(nivel_nominal = nivel, n_calibracion = n, conjunto = d$conjunto[1], n_evaluacion = nrow(d),
               semiancho = q, cobertura = mean(d$y >= lo & d$y <= hi),
               ancho_promedio_recortado = mean(hi - lo))
  }))
}))
data.table::fwrite(resumen, file.path(DIR_TABLAS, "16_resumen_intervalos.csv"))
for (i in seq_len(nrow(resumen))) {
  log_msg(sprintf("Nivel %.0f%% | %-13s | semiancho = %.3f | cobertura = %.1f%% | ancho medio = %.3f",
                  100 * resumen$nivel_nominal[i], resumen$conjunto[i], resumen$semiancho[i],
                  100 * resumen$cobertura[i], resumen$ancho_promedio_recortado[i]))
}

# Intervalos individuales al 90 % (para revisión y figuras).
q90 <- resumen$semiancho[resumen$nivel_nominal == 0.90][1]
eva[, `:=`(limite_inferior = pmax(pred - q90, rango[1]), limite_superior = pmin(pred + q90, rango[2]))]
eva[, cubierto := y >= limite_inferior & y <= limite_superior]
data.table::fwrite(eva, file.path(DIR_TABLAS, "16_intervalos_90_test.csv"))

# --- Política de inspección: alerta si la pérdida estimada >= c ------------------------------
# c se elige con las predicciones LOSO del entrenamiento (regla fijada en 00_config.R) y
# luego solo se VERIFICA en los test.
critico <- function(y) y >= cfg$umbral_critico
evaluar_alerta <- function(y, p, c) {
  al <- p >= c; cr <- critico(y)
  data.frame(umbral = c, sensibilidad = mean(al[cr]), tasa_falsa_alarma = mean(al[!cr]),
             fraccion_alertada = mean(al), prevalencia_critica = mean(cr))
}
grid <- do.call(rbind, lapply(cfg$grilla_umbral_alerta, function(c) evaluar_alerta(cal$y, cal$pred, c)))
cand <- grid[grid$sensibilidad >= cfg$sensibilidad_objetivo, ]
verificar(nrow(cand) > 0, "existe un umbral con la sensibilidad objetivo en las predicciones LOSO")
c_star <- max(cand$umbral)
politica <- rbind(
  cbind(conjunto = "entrenamiento (LOSO)", evaluar_alerta(cal$y, cal$pred, c_star)),
  do.call(rbind, lapply(split(eva, eva$conjunto), function(d) cbind(conjunto = d$conjunto[1], evaluar_alerta(d$y, d$pred, c_star)))))
data.table::fwrite(grid, file.path(DIR_TABLAS, "16_grilla_umbral_loso.csv"))
data.table::fwrite(politica, file.path(DIR_TABLAS, "16_politica_alerta.csv"))
log_msg(sprintf("Umbral de alerta elegido con LOSO: %.1f puntos %% (sensibilidad objetivo %.0f %%)",
                c_star, 100 * cfg$sensibilidad_objetivo))
for (i in seq_len(nrow(politica))) {
  log_msg(sprintf("  %-21s sensibilidad = %.3f | falsas alarmas = %.3f | alertas = %.3f",
                  politica$conjunto[i], politica$sensibilidad[i], politica$tasa_falsa_alarma[i],
                  politica$fraccion_alertada[i]))
}
log_msg("Paso 16 terminado. Los test solo se usaron para verificar cobertura y política.")
