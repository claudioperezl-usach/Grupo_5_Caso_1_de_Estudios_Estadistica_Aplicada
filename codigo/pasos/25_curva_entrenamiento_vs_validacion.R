# =============================================================================
# 25_curva_entrenamiento_vs_validacion.R — Error de entrenamiento vs validación (Elastic-net)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/25_curva_entrenamiento_vs_validacion.R")
# Requiere el paso 12 (tuning LOSO, parámetros elegidos y evaluación única en test).
#
# Curva de complejidad de clases ("Error train / Error validación", brecha train–valid):
#   * Entrenamiento: Elastic-net (alpha final) ajustado con TODO el entrenamiento para cada
#     lambda de la grilla y evaluado en esas mismas filas (dentro de muestra).
#   * Validación: RMSE LOSO medio ± 1 error estándar de cada lambda (paso 12, mismas particiones).
#   * Test: solo en el lambda elegido, leído del paso 12 (no se evalúa ningún otro lambda en test:
#     el test sigue siendo una evaluación única).
# Lambda chico = modelo flexible (riesgo de sobreajuste); lambda grande = modelo rígido (subajuste).
# Salidas (resultados_elastic_net_final): 25_curva_entrenamiento_validacion.csv y .png.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 25: curva de error de entrenamiento vs validación (Elastic-net) ===")

rutas <- file.path(DIR_TABLAS, c("12_tuning_loso.csv", "12_parametros_elegidos.csv",
                                 "12_resultados_test.csv", "17_convergencia_estricta.csv"))
if (any(!file.exists(rutas))) stop("Ejecute primero los pasos 12 y 17:\n", paste(rutas[!file.exists(rutas)], collapse = "\n"))
tuning <- data.table::fread(rutas[1])
pars   <- data.table::fread(rutas[2])
test   <- data.table::fread(rutas[3])
conv   <- data.table::fread(rutas[4])
alpha  <- pars$mixture[1]
lambda_final <- pars$penalty[1]

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
rm(part, datos)                                   # este paso no evalúa nada nuevo en test
predictores <- predictores_modelo(cfg)

# --- Error de entrenamiento (dentro de muestra) para cada lambda, alpha fijo -------------
cfg$grilla_lambda <- cfg$grilla_lambda_loso
curva <- data.table::rbindlist(lapply(cfg$grilla_lambda_loso, function(lam) {
  wf <- tune::finalize_workflow(construir_wf("enet", entrena, predictores, cfg),
                                tibble::tibble(penalty = lam, mixture = alpha))
  ajuste <- generics::fit(wf, data = entrena)
  p <- as.numeric(stats::predict(ajuste, new_data = entrena)$.pred)
  data.table::data.table(penalty = lam, rmse_entrenamiento = rmse_vec(entrena$y, p),
                         mae_entrenamiento = mean(abs(entrena$y - p)))
}))

# --- Error de validación LOSO del paso 12 (mismo alpha) ---------------------------------
val <- tuning[.metric == "rmse" & abs(mixture - alpha) < 1e-12,
              .(penalty, rmse_validacion = mean, ee_validacion = std_err)]
# Unión por lambda redondeado: el CSV no conserva todos los bits del valor de la grilla.
curva[, clave := signif(penalty, 6)]
val[, clave := signif(penalty, 6)]
curva <- merge(curva, val[, !"penalty"], by = "clave")[, clave := NULL]
curva[, brecha := rmse_validacion - rmse_entrenamiento]
data.table::setorder(curva, penalty)
verificar(nrow(curva) == length(cfg$grilla_lambda_loso), "curva: un punto de validación por cada lambda de la grilla")

i_fin <- which.min(abs(log(curva$penalty) - log(lambda_final)))
verificar(abs(curva$rmse_entrenamiento[i_fin] - conv$rmse_original[1]) < 1e-6,
          sprintf("RMSE de entrenamiento con el lambda elegido = paso 17 (%.3f)", conv$rmse_original[1]))
curva[, elegido := seq_len(.N) == i_fin]
guardar_tabla(curva, "25_curva_entrenamiento_validacion")
log_msg(sprintf("Lambda elegido %.4f: entrenamiento %.3f (MAE %.3f) | validación LOSO %.3f | test aleatorio %.3f | test 15 nudos %.3f",
                lambda_final, curva$rmse_entrenamiento[i_fin], curva$mae_entrenamiento[i_fin],
                curva$rmse_validacion[i_fin], test$rmse[test$conjunto == "interpolacion"],
                test$rmse[test$conjunto == "velocidad_15"]))
log_msg(sprintf("Lambda más chico (%.0e): entrenamiento %.3f | validación LOSO %.3f (brecha %.3f)",
                curva$penalty[1], curva$rmse_entrenamiento[1], curva$rmse_validacion[1], curva$brecha[1]))

# --- Figura -----------------------------------------------------------------------------
largo <- data.table::rbindlist(list(
  curva[, .(penalty, rmse = rmse_entrenamiento, lo = NA_real_, hi = NA_real_, serie = "Entrenamiento (dentro de muestra)")],
  curva[, .(penalty, rmse = rmse_validacion, lo = rmse_validacion - ee_validacion,
            hi = rmse_validacion + ee_validacion, serie = "Validación LOSO (media ± 1 EE)")]))
puntos_test <- data.table::data.table(
  penalty = lambda_final, rmse = test$rmse,
  serie = c(interpolacion = "Test aleatorio (λ elegido)", velocidad_15 = "Test 15 nudos (λ elegido)")[test$conjunto])
colores <- c("Entrenamiento (dentro de muestra)" = COL$secundario, "Validación LOSO (media ± 1 EE)" = COL$principal,
             "Test aleatorio (λ elegido)" = COL$tinta, "Test 15 nudos (λ elegido)" = COL$terciario)
p <- ggplot(largo, aes(penalty, rmse, colour = serie)) +
  geom_hline(yintercept = val[penalty == max(penalty), rmse_validacion], linetype = "dotted", colour = COL$baseline) +
  geom_vline(xintercept = lambda_final, linetype = "dashed", colour = COL$tinta_suave) +
  geom_ribbon(data = largo[!is.na(lo)], aes(ymin = lo, ymax = hi, fill = serie), alpha = 0.15, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  geom_point(data = puntos_test, size = 3.2, shape = 18) +
  annotate("text", x = min(curva$penalty) * 1.1, y = (curva$rmse_entrenamiento[1] + curva$rmse_validacion[1] - curva$ee_validacion[1]) / 2,
           hjust = 0, size = 3.2, colour = COL$tinta_suave,
           label = "λ chico: entrenamiento bajo,\nvalidación alta (sobreajuste)") +
  annotate("text", x = max(curva$penalty), y = min(largo$rmse, na.rm = TRUE), hjust = 1, vjust = 0, size = 3.2,
           colour = COL$tinta_suave, label = "λ grande: ambos ≈ media\n(subajuste)") +
  annotate("text", x = lambda_final * 1.15, y = max(largo$rmse, na.rm = TRUE) * 0.82, hjust = 0, size = 3.2,
           colour = COL$tinta_suave, label = sprintf("λ elegido (1-SE)\n= %s", format(signif(lambda_final, 2), decimal.mark = ","))) +
  scale_x_log10() +
  scale_colour_manual(values = colores, name = NULL) + scale_fill_manual(values = colores, guide = "none") +
  guides(colour = ggplot2::guide_legend(nrow = 2)) +
  labs(title = "Error de entrenamiento vs validación del Elastic-net",
       subtitle = sprintf("α = %s fijo; λ de la grilla de tuning. Punteada: baseline (media) en validación", format(alpha, decimal.mark = ",")),
       x = "λ (penalización, escala log)", y = "RMSE (puntos % de pérdida)",
       caption = "Test solo con el λ elegido (evaluación única del paso 12); no se usó para elegir λ.") +
  tema_proyecto()
guardar_figura(p, "25_curva_entrenamiento_validacion", ancho = 8, alto = 5)
log_msg("Paso 25 terminado. No se evaluó ningún λ nuevo en los test.")
