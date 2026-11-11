# =============================================================================
# 11_validacion_anidada_velocidad.R — Validación cruzada ANIDADA por velocidad
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/11_validacion_anidada_velocidad.R")
#
# Esquema (solo con el entrenamiento; NO abre el test aleatorio ni el de 15 nudos):
#   externo: Leave-One-Speed-Out (8 velocidades)  -> estima el error en una velocidad no vista
#   interno: Leave-One-Speed-Out (7 velocidades)  -> elige hiperparámetros (regla 1-SE)
# La receta (estandarización, nudos de splines) se ajusta dentro de cada partición.
# Es el esquema de clases (CV dentro del entrenamiento para elegir hiperparámetros) repetido
# dentro de cada velocidad externa: así el error externo no se usa para elegir lambda/alpha
# y estima honestamente el error del PROCEDIMIENTO completo en una velocidad no vista.
#
# Candidatos: Elastic-net, Ridge, Lasso y Splines naturales (+ baseline media global).
# Con cfg$anidar_ensembles = TRUE se agregan Árbol, Random Forest y XGBoost.
#
# Salidas (resultados_validacion_anidada/tablas):
#   11_detalle_anidado.csv        métricas por velocidad externa y modelo
#   11_parametros_internos.csv    hiperparámetros elegidos en cada partición externa
#   11_predicciones_externas.csv  predicciones fuera de partición (para incertidumbre)
#   11_ranking_anidado.csv        promedio entre velocidades
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_anidada)
log_msg("=== Paso 11: validación anidada por velocidad ===")

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
rm(part, datos)                                   # los test no se usan en este paso
predictores <- predictores_modelo(cfg)
velocidades <- sort(unique(entrena$v))
verificar(length(velocidades) == 8L && !(cfg$velocidad_reservada %in% velocidades),
          "el entrenamiento contiene 8 velocidades y excluye la velocidad reservada")

cfg$grilla_lambda <- cfg$grilla_lambda_loso       # grilla de lambda para LOSO (00_config.R)
claves <- c("enet", "ridge", "lasso", "splines")
if (isTRUE(cfg$anidar_ensembles)) claves <- c(claves, "arbol", "rf", "xgb")
log_msg("Candidatos anidados: ", paste(vapply(claves, function(k) FAMILIAS[[k]]$nombre, ""), collapse = ", "))

# crear_loso() y metricas_fila() están en 01_utilidades.R.

resultados <- list(); parametros <- list(); predicciones <- list()
for (v_externa in velocidades) {
  log_msg("Velocidad externa excluida: ", v_externa, " nudos")
  d_train <- entrena[entrena$v != v_externa, , drop = FALSE]
  d_test  <- entrena[entrena$v == v_externa, , drop = FALSE]
  verificar(!any(d_train$v == v_externa) && length(unique(d_train$v)) == 7L,
            paste0("partición externa ", v_externa, " nudos: 7 velocidades de ajuste y ninguna fila de la evaluada"))
  inner <- crear_loso(d_train)

  # Baseline: media del entrenamiento externo (sin la velocidad evaluada).
  p0 <- rep(mean(d_train$y), nrow(d_test))
  resultados[[length(resultados) + 1L]] <- metricas_fila(v_externa, "Baseline global", d_test$y, p0)
  predicciones[[length(predicciones) + 1L]] <- data.frame(velocidad = v_externa, modelo = "Baseline global",
                                                          id = d_test$id, y = d_test$y, pred = p0)
  for (clave in claves) {
    familia <- FAMILIAS[[clave]]
    wf <- construir_wf(clave, d_train, predictores, cfg)
    set.seed(cfg$semilla + as.integer(v_externa))
    ajuste_cv <- tune::tune_grid(wf, resamples = inner, grid = familia$grilla(cfg),
                                 metrics = yardstick::metric_set(yardstick::rmse),
                                 control = tune::control_grid(save_pred = FALSE, verbose = FALSE, allow_par = FALSE))
    # Regla 1-SE exclusivamente con las particiones internas.
    mejor <- tune::select_by_one_std_err(ajuste_cv, !!!familia$orden, metric = "rmse")
    # Solo los hiperparámetros de la grilla (select_by_one_std_err agrega columnas de métricas).
    pars <- dplyr::select(mejor, dplyr::all_of(names(familia$grilla(cfg))))
    set.seed(cfg$semilla + as.integer(v_externa))
    fit_externo <- generics::fit(tune::finalize_workflow(wf, pars), data = d_train)
    pred <- stats::predict(fit_externo, new_data = d_test)$.pred
    resultados[[length(resultados) + 1L]] <- metricas_fila(v_externa, familia$nombre, d_test$y, pred)
    predicciones[[length(predicciones) + 1L]] <- data.frame(velocidad = v_externa, modelo = familia$nombre,
                                                            id = d_test$id, y = d_test$y, pred = pred)
    parametros[[length(parametros) + 1L]] <- data.frame(
      velocidad = v_externa, modelo = familia$nombre,
      parametros = paste(names(pars), unlist(pars[1, ]), sep = "=", collapse = "; "))
    log_msg(sprintf("  %-20s RMSE externo = %.4f", familia$nombre, utils::tail(resultados, 1)[[1]]$rmse))
  }
  # Guardado incremental (permite revisar el avance si la ejecución se interrumpe).
  data.table::fwrite(data.table::rbindlist(resultados), file.path(DIR_TABLAS, "11_detalle_anidado.csv"))
  data.table::fwrite(data.table::rbindlist(parametros), file.path(DIR_TABLAS, "11_parametros_internos.csv"))
}

detalle <- data.table::rbindlist(resultados)
data.table::fwrite(data.table::rbindlist(predicciones), file.path(DIR_TABLAS, "11_predicciones_externas.csv"))
verificar(all(detalle[, .N, by = modelo]$N == 8L), "cada modelo tiene exactamente 8 evaluaciones externas")

ranking <- detalle[, .(rmse_promedio = mean(rmse), de_rmse_entre_velocidades = stats::sd(rmse),
                       mae_promedio = mean(mae), r2_promedio = mean(r2),
                       rmse_critico_promedio = mean(rmse_critico, na.rm = TRUE),
                       rmse_maximo = max(rmse)), by = modelo]
data.table::setorder(ranking, rmse_promedio)
data.table::fwrite(ranking, file.path(DIR_TABLAS, "11_ranking_anidado.csv"))
log_msg("Ranking anidado: ", paste(sprintf("%s = %.4f", ranking$modelo, ranking$rmse_promedio), collapse = " | "))

# --- Verificación automática de reproducibilidad (ver verificar_referencia_ranking) -------
# Primera ejecución con una configuración: guarda la referencia. Siguientes con la misma
# configuración (misma huella): exige diferencias <= cfg$tolerancia_reproducibilidad.
verificar_referencia_ranking(ranking, cfg)
log_msg("Paso 11 terminado. No se evaluaron los test.")
