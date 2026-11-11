# =============================================================================
# 12_entrenamiento_final_elastic_net.R — Sistema final congelado y evaluación en los test
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/12_entrenamiento_final_elastic_net.R")
#
# 1. Elige lambda y alpha de Elastic-net con LOSO (8 velocidades) y regla 1-SE, SOLO
#    con el entrenamiento.
# 2. Guarda las predicciones fuera de partición (LOSO) con esos hiperparámetros: son los
#    residuos que calibran los intervalos de predicción (paso 16) sin tocar los test.
# 3. Reajusta con todo el entrenamiento, guarda el modelo congelado y lo evalúa UNA vez
#    en el test aleatorio (interpolación) y en el test de 15 nudos.
# NOTA: ambos test se examinaron en un análisis exploratorio previo con splines (ya no incluido);
# esta evaluación es adicional y no totalmente independiente. Nada de lo que se calcula
# con los test se usa para elegir modelo, hiperparámetros ni intervalos.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 12: selección LOSO, modelo final congelado y evaluación en test ===")

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
predictores <- predictores_modelo(cfg)
n_resto <- nrow(entrena) + nrow(part$prueba_interp)
verificar(nrow(part$prueba_interp) == round(cfg$prop_test * n_resto) && nrow(part$prueba_vel) > 0,
          sprintf("partición coherente con la configuración: %d entrenamiento, %d test aleatorio, %d test de %s nudos",
                  nrow(entrena), nrow(part$prueba_interp), nrow(part$prueba_vel), cfg$velocidad_reservada))

folds <- crear_loso(entrena)

# --- 1. Selección de hiperparámetros (solo entrenamiento) -------------------------------
cfg$grilla_lambda <- cfg$grilla_lambda_loso
familia <- FAMILIAS[["enet"]]
wf <- construir_wf("enet", entrena, predictores, cfg)
set.seed(cfg$semilla)
tuning <- tune::tune_grid(wf, resamples = folds, grid = familia$grilla(cfg),
                          metrics = yardstick::metric_set(yardstick::rmse, yardstick::mae),
                          control = tune::control_grid(save_pred = TRUE, verbose = FALSE, allow_par = FALSE))
data.table::fwrite(tune::collect_metrics(tuning), file.path(DIR_TABLAS, "12_tuning_loso.csv"))
notas <- tune::collect_notes(tuning)
# Las notas de tune pueden contener columnas de tipo lista: se guardan en .rds.
if (nrow(notas) > 0) saveRDS(notas, file.path(DIR_TABLAS, "12_advertencias_tuning.rds"))
log_msg("Advertencias de tune durante el ajuste: ", nrow(notas))

mejor <- tune::select_by_one_std_err(tuning, !!!familia$orden, metric = "rmse")
# Solo los hiperparámetros de la grilla (select_by_one_std_err agrega columnas de métricas).
parametros <- dplyr::select(mejor, dplyr::all_of(names(familia$grilla(cfg))))
data.table::fwrite(parametros, file.path(DIR_TABLAS, "12_parametros_elegidos.csv"))
log_msg("Hiperparámetros elegidos (1-SE): penalty = ", signif(parametros$penalty, 4),
        "; mixture = ", parametros$mixture)

# --- 2. Predicciones LOSO fuera de partición con los hiperparámetros elegidos ------------
oof <- tune::collect_predictions(tuning, parameters = parametros)
oof <- data.frame(velocidad = entrena$v[oof$.row], id = entrena$id[oof$.row],
                  y = oof$y, pred = oof$.pred)
verificar(nrow(oof) == nrow(entrena) && !anyDuplicated(oof$id),
          "una predicción LOSO fuera de partición por fila de entrenamiento")
data.table::fwrite(oof, file.path(DIR_TABLAS, "12_predicciones_loso.csv"))

# --- 3. Sistema congelado y evaluación única en los test ---------------------------------
wf_final <- tune::finalize_workflow(wf, parametros)
set.seed(cfg$semilla)
ajuste <- generics::fit(wf_final, data = entrena)
saveRDS(ajuste, file.path(DIR_RESULTADOS, "elastic_net_final.rds"))
log_msg("Modelo congelado guardado en ", file.path(DIR_RESULTADOS, "elastic_net_final.rds"))

pruebas <- list(interpolacion = part$prueba_interp, velocidad_15 = part$prueba_vel)
for (nm in names(pruebas)) {
  verificar(length(intersect(entrena$id, pruebas[[nm]]$id)) == 0,
            paste0("ninguna fila del test '", nm, "' está en el entrenamiento"))
}
baseline <- mean(entrena$y)
met <- function(y, p) c(rmse = sqrt(mean((y - p)^2)), mae = mean(abs(y - p)),
                        r2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
resumen <- list(); predicciones <- list(); incert <- list()
for (nombre in names(pruebas)) {
  d <- pruebas[[nombre]]
  pred <- as.numeric(stats::predict(ajuste, new_data = d)$.pred)
  pred_ref <- rep(baseline, nrow(d))
  a <- met(d$y, pred); b <- met(d$y, pred_ref)
  resumen[[nombre]] <- data.frame(conjunto = nombre, n = nrow(d),
    rmse = unname(a["rmse"]), mae = unname(a["mae"]), r2 = unname(a["r2"]),
    rmse_baseline = unname(b["rmse"]), mae_baseline = unname(b["mae"]),
    mejora_rmse_pct = unname(100 * (1 - a["rmse"] / b["rmse"])))
  predicciones[[nombre]] <- data.frame(conjunto = nombre, id = d$id, velocidad = d$v,
                                       y = d$y, pred = pred, pred_baseline = pred_ref)
  # Incertidumbre: bootstrap de filas (B = cfg$B_bootstrap) y subgrupo crítico (y >= umbral).
  # Supuesto: filas intercambiables. En la grilla, las filas comparten velocidad y niveles de
  # kMc/kMt, así que este IC refleja solo la variación de muestreo de filas del test y es
  # optimista (más estrecho) frente a la incertidumbre entre regímenes; la variación entre
  # velocidades se reporta con la validación LOSO (pasos 11 y 21).
  set.seed(cfg$semilla)
  boot <- t(replicate(cfg$B_bootstrap, {
    i <- sample.int(nrow(d), replace = TRUE)
    e1 <- sqrt(mean((d$y[i] - pred[i])^2)); e0 <- sqrt(mean((d$y[i] - pred_ref[i])^2))
    c(e1, e1 - e0)
  }))
  crit <- d$y >= cfg$umbral_critico
  incert[[nombre]] <- data.frame(conjunto = nombre,
    rmse = unname(a["rmse"]),
    rmse_ic_inf = stats::quantile(boot[, 1], 0.025), rmse_ic_sup = stats::quantile(boot[, 1], 0.975),
    dif_vs_baseline = unname(a["rmse"] - b["rmse"]),
    dif_ic_inf = stats::quantile(boot[, 2], 0.025), dif_ic_sup = stats::quantile(boot[, 2], 0.975),
    n_critico = sum(crit),
    rmse_critico = sqrt(mean((d$y[crit] - pred[crit])^2)),
    rmse_critico_baseline = sqrt(mean((d$y[crit] - pred_ref[crit])^2)), row.names = NULL)
}
res <- do.call(rbind, resumen)
data.table::fwrite(res, file.path(DIR_TABLAS, "12_resultados_test.csv"))
data.table::fwrite(do.call(rbind, predicciones), file.path(DIR_TABLAS, "12_predicciones_test.csv"))
data.table::fwrite(do.call(rbind, incert), file.path(DIR_TABLAS, "12_test_incertidumbre_subgrupo.csv"))
for (i in seq_len(nrow(res))) {
  log_msg(sprintf("TEST %-13s RMSE = %.4f (baseline %.4f; mejora %.1f%%) | MAE = %.4f | R2 = %.4f",
                  res$conjunto[i], res$rmse[i], res$rmse_baseline[i], res$mejora_rmse_pct[i], res$mae[i], res$r2[i]))
}
log_msg("Paso 12 terminado. Los test no se usaron para elegir hiperparámetros.")
