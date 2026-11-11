# =============================================================================
# 23_cv_aleatoria_vs_loso.R — ¿Por qué la partición se agrupa por velocidad?
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/23_cv_aleatoria_vs_loso.R")
#
# Respalda el §3 del informe SOLO con el entrenamiento (no abre ningún test):
#   1. Validación cruzada ALEATORIA de 5 particiones (vfold_cv, como en clases) para elegir
#      los hiperparámetros de cada candidato (menor RMSE).
#   2. La MISMA configuración se evalúa dejando una velocidad completa fuera (LOSO).
# Si el RMSE aleatorio es mucho menor que el LOSO, la CV aleatoria es optimista: en la grilla,
# filas vecinas comparten velocidad y niveles de degradación, y el modelo "reconoce" el régimen
# en vez de generalizar a uno nuevo (slides: "grupos repetidos condicionan la partición").
# Este paso es ilustrativo: NO elige el modelo final (eso lo hace la anidada del paso 11).
# Salida: resultados_validacion_anidada/tablas/23_cv_aleatoria_vs_loso.csv
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_anidada)
log_msg("=== Paso 23: validación aleatoria frente a LOSO (solo entrenamiento) ===")

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
rm(part, datos)                                   # los test no se usan en este paso
predictores <- predictores_modelo(cfg)
cfg$grilla_lambda <- cfg$grilla_lambda_loso

# Particiones: 5 aleatorias (como en clases) y una por velocidad (LOSO).
set.seed(cfg$semilla)
folds_aleatorios <- rsample::vfold_cv(entrena, v = 5)
folds_loso <- crear_loso(entrena)

metrica <- yardstick::metric_set(yardstick::rmse)
resumen <- list()

# Baseline: media del ajuste en cada partición.
rmse_baseline <- function(folds) mean(vapply(folds$splits, function(s) {
  rmse_vec(rsample::assessment(s)$y, mean(rsample::analysis(s)$y))
}, numeric(1)))
resumen[[1]] <- data.frame(modelo = "Baseline (media global)", configuracion = "—",
                           rmse_cv_aleatoria = rmse_baseline(folds_aleatorios),
                           rmse_loso = rmse_baseline(folds_loso))

for (clave in c("enet", "ridge", "lasso", "splines")) {
  familia <- FAMILIAS[[clave]]
  wf <- construir_wf(clave, entrena, predictores, cfg)
  set.seed(cfg$semilla)
  ajuste_cv <- tune::tune_grid(wf, resamples = folds_aleatorios, grid = familia$grilla(cfg), metrics = metrica,
                               control = tune::control_grid(save_pred = FALSE, allow_par = FALSE))
  mejor <- tune::select_best(ajuste_cv, metric = "rmse")       # lo que haría un análisis ingenuo
  pars <- dplyr::select(mejor, -dplyr::any_of(".config"))
  rmse_cv <- dplyr::filter(tune::collect_metrics(ajuste_cv), .config == mejor$.config)$mean
  set.seed(cfg$semilla)
  ev_loso <- tune::fit_resamples(tune::finalize_workflow(wf, pars), resamples = folds_loso, metrics = metrica,
                                 control = tune::control_resamples(save_pred = FALSE, allow_par = FALSE))
  resumen[[length(resumen) + 1L]] <- data.frame(
    modelo = familia$nombre,
    configuracion = paste(names(pars), signif(unlist(pars[1, ]), 3), sep = " = ", collapse = "; "),
    rmse_cv_aleatoria = rmse_cv,
    rmse_loso = tune::collect_metrics(ev_loso)$mean)
  log_msg(sprintf("  %-18s RMSE CV aleatoria = %.3f | RMSE LOSO = %.3f", familia$nombre,
                  rmse_cv, utils::tail(resumen, 1)[[1]]$rmse_loso))
}
tabla <- do.call(rbind, resumen)
tabla$cociente_loso_sobre_aleatoria <- tabla$rmse_loso / tabla$rmse_cv_aleatoria
guardar_tabla(tabla, "23_cv_aleatoria_vs_loso")
log_msg("Mejor por CV aleatoria: ", tabla$modelo[which.min(tabla$rmse_cv_aleatoria)],
        " | mejor por LOSO: ", tabla$modelo[which.min(tabla$rmse_loso)])
log_msg("Paso 23 terminado. Los test no se utilizaron.")
