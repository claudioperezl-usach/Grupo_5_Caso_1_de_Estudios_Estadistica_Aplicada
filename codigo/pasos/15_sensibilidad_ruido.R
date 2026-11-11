# =============================================================================
# 15_sensibilidad_ruido.R — Sensibilidad del Elastic-net al ruido de sensores (LOSO)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/15_sensibilidad_ruido.R")
# Requiere el paso 12 (hiperparámetros finales).
#
# Pregunta: si los sensores tuvieran un ruido de medición pequeño, ¿cuánto empeoraría el error?
# Se usa SOLO el entrenamiento, con la misma lógica de la validación LOSO:
#   * en cada velocidad excluida se reajusta el Elastic-net con los hiperparámetros finales
#     usando las otras 7 velocidades;
#   * a los sensores de la velocidad excluida se les suma ruido normal independiente con
#     DE = 1 % de la DE del sensor en las 7 velocidades de ajuste (cfg$B_ruido repeticiones);
#   * se compara el RMSE con y sin ruido.
# La velocidad (v) no se perturba: es la condición de operación. Los test no se usan.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 15: sensibilidad a ruido de sensores (LOSO, solo entrenamiento) ===")

ruta_par <- file.path(DIR_TABLAS, "12_parametros_elegidos.csv")
if (!file.exists(ruta_par)) stop("Ejecute primero el paso 12: falta ", ruta_par)
pars <- tibble::as_tibble(as.data.frame(data.table::fread(ruta_par))[1, c("penalty", "mixture"), drop = FALSE])

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
rm(part, datos)                                   # los test no se usan en este paso
predictores <- predictores_modelo(cfg)
sensores <- setdiff(predictores, "v")
nivel_ruido <- 0.01                               # DE del ruido = 1 % de la DE del sensor

set.seed(cfg$semilla_ruido)
detalle <- list()
for (vel in sort(unique(entrena$v))) {
  tr <- entrena[entrena$v != vel, , drop = FALSE]
  te <- entrena[entrena$v == vel, , drop = FALSE]
  verificar(!any(tr$v == vel), paste0("ruido, ", vel, " nudos: la velocidad evaluada no está en el ajuste"))
  ajuste <- generics::fit(tune::finalize_workflow(construir_wf("enet", tr, predictores, cfg), pars), data = tr)
  rmse_original <- rmse_vec(te$y, stats::predict(ajuste, new_data = te)$.pred)
  escalas <- vapply(sensores, function(nm) stats::sd(tr[[nm]]), numeric(1))   # solo ajuste
  rmse_ruido <- vapply(seq_len(cfg$B_ruido), function(b) {
    e2 <- te
    for (nm in sensores) e2[[nm]] <- te[[nm]] + stats::rnorm(nrow(te), 0, nivel_ruido * escalas[[nm]])
    rmse_vec(te$y, stats::predict(ajuste, new_data = e2)$.pred)
  }, numeric(1))
  detalle[[length(detalle) + 1L]] <- data.frame(
    velocidad = vel, rmse_original = rmse_original, rmse_ruido_promedio = mean(rmse_ruido),
    rmse_ruido_p05 = unname(stats::quantile(rmse_ruido, 0.05)),
    rmse_ruido_p95 = unname(stats::quantile(rmse_ruido, 0.95)),
    cambio_pct = 100 * (mean(rmse_ruido) / rmse_original - 1))
}
detalle <- do.call(rbind, detalle)
resumen <- data.frame(repeticiones = cfg$B_ruido, ruido_pct_de = 100 * nivel_ruido,
                      rmse_original_medio = mean(detalle$rmse_original),
                      rmse_ruido_medio = mean(detalle$rmse_ruido_promedio),
                      cambio_medio_pct = mean(detalle$cambio_pct),
                      cambio_max_pct = max(detalle$cambio_pct))
guardar_tabla(detalle, "15_sensibilidad_ruido_por_velocidad")
guardar_tabla(resumen, "15_resumen_sensibilidad_ruido")
log_msg(sprintf("Ruido de 1 %% de la DE en %d sensores: el RMSE LOSO sube %.2f %% en promedio (máximo %.2f %%)",
                length(sensores), resumen$cambio_medio_pct, resumen$cambio_max_pct))
log_msg("Paso 15 terminado. Simulación de ruido, no validación con sensores reales; los test no se utilizaron.")
