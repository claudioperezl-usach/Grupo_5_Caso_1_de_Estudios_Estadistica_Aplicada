# =============================================================================
# 13_graficos_finales.R — Figuras del sistema final (sin reentrenar)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/13_graficos_finales.R")
# Lee resultados de los pasos 11 y 12. Genera en resultados_elastic_net_final/figuras:
#   13_rmse_por_velocidad.png   (Elastic-net, Ridge y baseline en cada velocidad no vista)
#   13_test_vs_baseline.png
#   13_predicho_vs_observado.png

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto(character(0))
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 13: figuras finales ===")

archivos <- c(file.path(cfg$dir_anidada, "tablas", "11_detalle_anidado.csv"),
              file.path(DIR_TABLAS, "12_resultados_test.csv"),
              file.path(DIR_TABLAS, "12_predicciones_test.csv"))
if (any(!file.exists(archivos))) stop("Faltan resultados de los pasos 11 y 12:\n", paste(archivos[!file.exists(archivos)], collapse = "\n"))
detalle <- data.table::fread(archivos[1])
test    <- data.table::fread(archivos[2])
pred    <- data.table::fread(archivos[3])

# 1. Error en cada velocidad excluida (validación anidada).
d <- detalle[modelo %in% c("Elastic-net", "Ridge", "Baseline global")]
d[, modelo := factor(modelo, levels = c("Elastic-net", "Ridge", "Baseline global"))]
ancho <- data.table::dcast(d, velocidad ~ modelo, value.var = "rmse")
gana <- sum(ancho[["Elastic-net"]] < ancho[["Baseline global"]])
subtitulo <- sprintf("Elastic-net mejora al baseline en %d de %d velocidades", gana, nrow(ancho))
p1 <- ggplot(d, aes(factor(velocidad), rmse, colour = modelo, group = modelo)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("Elastic-net" = COL$principal, "Ridge" = COL$secundario,
                                 "Baseline global" = COL$baseline), name = NULL) +
  labs(title = "Error en velocidades no vistas (validación anidada LOSO)",
       subtitle = subtitulo,
       x = "Velocidad excluida (nudos)", y = "RMSE (puntos % de pérdida)",
       caption = "3 y 27 nudos son extrapolación de borde; el resto, interpolación entre velocidades.") +
  tema_proyecto()
guardar_figura(p1, "13_rmse_por_velocidad")

# 2. Test: Elastic-net frente al baseline.
tl <- data.table::melt(test[, .(conjunto, `Elastic-net` = rmse, Baseline = rmse_baseline)],
                       id.vars = "conjunto", variable.name = "modelo", value.name = "rmse")
p2 <- ggplot(tl, aes(conjunto, rmse, fill = modelo)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c(`Elastic-net` = COL$principal, Baseline = COL$baseline), name = NULL) +
  labs(title = "Evaluación final: Elastic-net frente al baseline",
       subtitle = "Ambos test se examinaron antes con el pipeline histórico: evaluación adicional",
       x = NULL, y = "RMSE (puntos %)") +
  tema_proyecto()
guardar_figura(p2, "13_test_vs_baseline")

# 3. Predicho vs observado en los test.
p3 <- ggplot(pred, aes(y, pred)) +
  geom_jitter(width = 0.03, height = 0, alpha = 0.15, size = 0.6, colour = COL$principal) +
  geom_abline(slope = 1, intercept = 0, colour = COL$secundario, linewidth = 0.8) +
  facet_wrap(~conjunto) +
  labs(title = "Pérdida observada vs estimada (Elastic-net)",
       subtitle = "La recta es la estimación perfecta; la dispersión muestra el error individual",
       x = "Pérdida observada (puntos %)", y = "Pérdida estimada (puntos %)") +
  tema_proyecto()
guardar_figura(p3, "13_predicho_vs_observado")
log_msg("Paso 13 terminado.")
