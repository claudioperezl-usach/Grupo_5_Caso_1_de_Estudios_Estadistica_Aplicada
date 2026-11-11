# =============================================================================
# 26_figuras_presentacion.R — Figuras de la presentación (Equipo 5)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/26_figuras_presentacion.R")
# Requiere los pasos 10 a 24. NO calcula resultados nuevos: solo grafica, en alta resolución
# (350 dpi) y con el tamaño de cada marco de la presentación, las tablas que ya generó el script:
# (número de lámina en el PowerPoint, contando la portada como 1)
#   lámina 4   relación T2-pérdida (entrenamiento; misma recta que el paso 10)
#   lámina 6   VIF                             <- 05_ols_vif.csv
#   lámina 7   RMSE por velocidad y Tabla 1    <- 11_detalle_anidado.csv, 21_tabla1_validacion.csv
#   lámina 8   PDP/ICE y SHAP                  <- 22_ice.csv, 22_pdp.csv, 24_shap_muestra.csv, 24_shap_importancia.csv
#   lámina 10  borde de 3 nudos                <- 11_detalle_anidado.csv
# Salida: resultados_elastic_net_final/figuras/presentacion/26_*.png. Los test no se usan.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto(c("02_datos.R", "03_particion.R"))
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 26: figuras de la presentación ===")
suppressPackageStartupMessages({ library(grid); library(gridExtra) })

T <- DIR_TABLAS; A <- file.path(cfg$dir_anidada, "tablas")
OUT <- file.path(DIR_FIGURAS, "presentacion")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
necesarios <- c(file.path(T, c("04_relacion_ic.csv", "05_ols_vif.csv", "21_tabla1_validacion.csv", "22_ice.csv", "22_pdp.csv",
                               "24_shap_muestra.csv", "24_shap_importancia.csv")), file.path(A, "11_detalle_anidado.csv"))
if (any(!file.exists(necesarios))) stop("Ejecute primero los pasos 10 a 24. Faltan:\n", paste(necesarios[!file.exists(necesarios)], collapse = "\n"))

DPI <- 350
# Fuente Inter si está instalada; si no, la sans del sistema (el contenido no cambia).
FUENTE <- tryCatch(if ("Inter" %in% systemfonts::system_fonts()$family) "Inter" else "sans", error = function(e) "sans")
usa_ragg <- requireNamespace("ragg", quietly = TRUE)
abrir_png <- function(ruta, w, h) {
  if (usa_ragg) ragg::agg_png(ruta, width = w, height = h, units = "in", res = DPI, background = "white")
  else grDevices::png(ruta, width = w, height = h, units = "in", res = DPI, bg = "white")
}
AZUL <- "#2A78D6"; NARANJA <- "#E8692F"; GRIS <- "#8E8E8A"; TINTA <- "#1F2A44"; SUAVE <- "#6B6F7B"; REJ <- "#E7E7E2"
coma <- function(x, d = 3) formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = ".")
tema <- function(base = 10) theme_minimal(base_size = base, base_family = FUENTE) +
  theme(plot.title = element_text(face = "bold", colour = TINTA, size = base * 1.15, margin = margin(b = 2)),
        plot.subtitle = element_text(colour = SUAVE, size = base * 0.85, margin = margin(b = 6)),
        plot.caption = element_text(colour = SUAVE, size = base * 0.72, hjust = 0),
        axis.title = element_text(colour = SUAVE, size = base * 0.85), axis.text = element_text(colour = TINTA, size = base * 0.8),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = REJ, linewidth = 0.35),
        legend.position = "top", legend.text = element_text(colour = TINTA, size = base * 0.8),
        plot.margin = margin(6, 8, 4, 6), plot.background = element_rect(fill = "white", colour = NA))
guardar <- function(p, nombre, w, h) {
  ruta <- file.path(OUT, paste0("26_", nombre)); abrir_png(ruta, w, h); print(p); invisible(grDevices::dev.off())
  log_msg("    figura: ", ruta)
}

# ---- Lámina 4: relación T2 - pérdida a igual velocidad (entrenamiento) ----------------
set.seed(cfg$semilla)
d <- data.table::as.data.table(particionar_datos(cargar_datos(cfg), cfg)$entrena)
d[, `:=`(x_c = t2 - mean(t2), y_c = y - mean(y)), by = v]
b <- stats::coef(stats::lm(y ~ t2 + factor(v), data = d))[["t2"]]
b_paso10 <- as.numeric(data.table::fread(file.path(T, "04_relacion_ic.csv"))$valor[1])
verificar(abs(b - b_paso10) < 1e-4, sprintf("pendiente de T2 igual a la del paso 10 (%.5f)", b_paso10))
p4 <- ggplot(d, aes(x_c, y_c)) +
  geom_point(colour = AZUL, alpha = 0.22, size = 0.8, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = NARANJA, linewidth = 1.1) +
  annotate("label", x = min(d$x_c) * 0.95, y = max(d$y_c) * 0.92, hjust = 0, vjust = 1, size = 3.1, family = FUENTE,
           colour = TINTA, fill = "white", label.size = 0, label = sprintf("Pendiente: +%s puntos por K\n(+%s puntos por cada +10 K)", coma(b, 3), coma(10 * b, 3))) +
  labs(title = "Pérdida de rendimiento vs T2, a igual velocidad",
       subtitle = "Ambas variables centradas en el promedio de su velocidad (8.486 filas de entrenamiento)",
       x = "T2 − promedio de su velocidad (K)", y = "Pérdida − promedio de su velocidad\n(puntos %)",
       caption = "Asociación, no causalidad. Datos simulados.") + tema(10)
guardar(p4, "lamina04_relacion_t2.png", 6.21, 3.49)

# ---- Lámina 6: VIF ---------------------------------------------------------------------
vif <- data.table::fread(file.path(T, "05_ols_vif.csv"))[order(VIF)]
vif[, variable := factor(toupper(variable), levels = toupper(variable))]
p6 <- ggplot(vif, aes(VIF, variable)) +
  geom_vline(xintercept = 10, colour = NARANJA, linetype = "dashed", linewidth = 0.6) +
  geom_segment(aes(x = 1, xend = VIF, yend = variable), colour = "#BFD3F2", linewidth = 1.6, lineend = "round") +
  geom_point(colour = AZUL, size = 2.6) +
  geom_text(aes(x = VIF * 1.45, label = trimws(format(round(VIF), big.mark = ".", decimal.mark = ","))), hjust = 0, size = 2.9, family = FUENTE, colour = TINTA) +
  annotate("text", x = 10, y = nrow(vif) + 0.75, label = "VIF = 10 (alerta)", hjust = -0.08, vjust = 0.5, size = 2.6, colour = NARANJA, family = FUENTE) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 1e5), labels = function(x) trimws(format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)),
                limits = c(1, 3e5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
  labs(title = "Inflación de la varianza en OLS (VIF)", subtitle = "Todos los predictores superan 10 · escala logarítmica",
       x = NULL, y = NULL) + tema(10) + theme(panel.grid.major.y = element_blank())
guardar(p6, "lamina06_vif.png", 4.87, 3.35)

# ---- Lámina 7: RMSE por velocidad + Tabla 1 ------------------------------------------------
det <- data.table::fread(file.path(A, "11_detalle_anidado.csv"))[modelo %in% c("Elastic-net", "Ridge", "Baseline global")]
det[, modelo := factor(modelo, levels = c("Elastic-net", "Ridge", "Baseline global"))]
t1 <- data.table::fread(file.path(T, "21_tabla1_validacion.csv"))
en <- t1[modelo == "Elastic-net"]
cols <- c("Elastic-net" = AZUL, "Ridge" = NARANJA, "Baseline global" = GRIS)
p7a <- ggplot(det, aes(factor(velocidad), rmse, colour = modelo, group = modelo)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "RMSE en velocidades no vistas",
       subtitle = sprintf("Validación anidada · EN − baseline: %s [%s; %s]\nmejora en %d de 8 velocidades",
                          coma(en$dif_vs_base), coma(en$dif_ic_inf), coma(en$dif_ic_sup), en$gana_baseline),
       x = "Velocidad excluida (nudos)", y = "RMSE (puntos %)") + tema(9.5) +
  theme(legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))
orden <- c("Media global (baseline)", "OLS", "Ridge", "Lasso", "Elastic-net", "Splines naturales", "GAM",
           "Árbol de regresión", "Random Forest", "XGBoost", "Promedio EN + Ridge")
tt <- t1[match(orden, modelo)]
tab <- data.frame(
  Modelo = ifelse(tt$modelo == "Elastic-net", "Elastic-net (final)", ifelse(tt$modelo == "Media global (baseline)", "Baseline (media)", tt$modelo)),
  `RMSE ± DE` = paste(coma(tt$rmse), "±", coma(tt$de_rmse)), MAE = coma(tt$mae), `R²` = coma(ifelse(abs(tt$r2) < 5e-4, 0, tt$r2)),
  `RMSE y≥4` = coma(tt$rmse_critico),
  `vs base` = ifelse(tt$modelo == "Media global (baseline)", "ref.", paste0(ifelse(tt$mejora_pct > 0, "+", ""), coma(tt$mejora_pct, 1), " %")),
  check.names = FALSE)
fila_en <- which(tt$modelo == "Elastic-net")
fondos <- rep(c("white", "#F4F6FA"), length.out = nrow(tab)); fondos[fila_en] <- "#DCE8FA"
pesos <- rep("plain", nrow(tab)); pesos[fila_en] <- "bold"
tema_tab <- ttheme_minimal(
  base_family = FUENTE, base_size = 7.6, padding = unit(c(3.2, 2.6), "mm"),
  core = list(fg_params = list(col = TINTA, hjust = rep(c(0, 1, 1, 1, 1, 1), each = nrow(tab)), x = rep(c(0.04, 0.95, 0.95, 0.95, 0.95, 0.95), each = nrow(tab)),
                               fontface = rep(pesos, 6)),
              bg_params = list(fill = fondos, col = NA)),
  colhead = list(fg_params = list(col = "white", fontface = "bold", hjust = c(0, 1, 1, 1, 1, 1), x = c(0.04, 0.95, 0.95, 0.95, 0.95, 0.95)),
                 bg_params = list(fill = "#1F3864", col = NA)))
g_tab <- tableGrob(tab, rows = NULL, theme = tema_tab)
titulo_tab <- textGrob("Tabla 1 · Validación (promedio de 8 velocidades)", x = 0.01, hjust = 0,
                       gp = gpar(fontfamily = FUENTE, fontface = "bold", col = TINTA, fontsize = 10.5))
nota_tab <- textGrob("RMSE y≥4: subgrupo crítico · vs base: reducción del RMSE frente a la media", x = 0.01, hjust = 0,
                     gp = gpar(fontfamily = FUENTE, col = SUAVE, fontsize = 6.8))
panel_tab <- arrangeGrob(nullGrob(), titulo_tab, g_tab, nota_tab, nullGrob(), ncol = 1,
                         heights = unit.c(unit(0.35, "in"), unit(2.2, "lines"), sum(g_tab$heights), unit(1.6, "lines"), unit(1, "null")))
ruta7 <- file.path(OUT, "26_lamina07_rmse_tabla1.png"); abrir_png(ruta7, 8.49, 4.12)
grid.arrange(ggplotGrob(p7a), panel_tab, ncol = 2, widths = c(0.43, 0.57))
invisible(grDevices::dev.off()); log_msg("    figura: ", ruta7)

# ---- Lámina 8: PDP/ICE y SHAP ---------------------------------------------------------------
ice <- data.table::fread(file.path(T, "22_ice.csv")); pdp <- data.table::fread(file.path(T, "22_pdp.csv"))
etq <- c(t2 = "T2 · temp. salida compresor (K)", pexh = "Pexh · presión de escape (bar)", p2 = "P2 · presión salida compresor (bar)", tic = "TIC · control de inyección (%)")
niv <- unique(pdp$variable)
ice[, panel := factor(etq[variable], levels = etq[niv])]; pdp[, panel := factor(etq[variable], levels = etq[niv])]
p8a <- ggplot() +
  geom_line(data = ice, aes(valor, pred, group = fila), colour = AZUL, alpha = 0.10, linewidth = 0.35) +
  geom_line(data = pdp, aes(valor, pred), colour = NARANJA, linewidth = 1.1) +
  facet_wrap(~panel, scales = "free_x", ncol = 2) +
  labs(title = "PDP (naranja) e ICE (azul)", subtitle = "4 variables más importantes · 150 filas de entrenamiento",
       x = NULL, y = "Pérdida predicha (puntos %)", caption = "Fuera de 0–5: combinaciones de sensores que no existen en la planta.") +
  tema(8.5) + theme(strip.text = element_text(colour = TINTA, size = 7, face = "bold", hjust = 0), panel.spacing = unit(6, "pt"))
guardar(p8a, "lamina08_pdp_ice.png", 4.12, 3.09)

sh <- data.table::fread(file.path(T, "24_shap_muestra.csv")); imp <- data.table::fread(file.path(T, "24_shap_importancia.csv"))
sh[, variable := factor(variable, levels = toupper(rev(imp$variable)))]
p8b <- ggplot(sh, aes(shap, variable, colour = valor)) +
  geom_vline(xintercept = 0, colour = SUAVE, linewidth = 0.4) +
  geom_jitter(height = 0.22, width = 0, size = 0.55, alpha = 0.55, stroke = 0) +
  scale_colour_gradient(low = NARANJA, high = AZUL, breaks = c(0, 1), labels = c("bajo", "alto"), name = "Valor del\nsensor") +
  labs(title = "Valores SHAP del Elastic-net", subtitle = "Aporte a cada predicción respecto de la media (2.000 filas)",
       x = "Valor SHAP (puntos % de pérdida)", y = NULL, caption = "SHAP = coeficiente × (valor − media). Aportes que se compensan: sensores correlacionados.") +
  tema(8.5) + theme(legend.position = "right", legend.key.height = unit(14, "pt"), legend.title = element_text(size = 7, colour = TINTA))
guardar(p8b, "lamina08_shap.png", 4.89, 2.99)

# ---- Lámina 10: borde de 3 nudos -------------------------------------------------------------
z <- det[velocidad %in% c(3, 6, 9)]
lab3 <- z[velocidad == 3]
p10 <- ggplot(z, aes(factor(velocidad), rmse, colour = modelo, group = modelo)) +
  annotate("rect", xmin = 0.6, xmax = 1.4, ymin = -Inf, ymax = Inf, fill = "#FDE7DC", alpha = 0.8) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  geom_text(data = lab3, aes(label = coma(rmse)), hjust = 1.35, size = 2.7, family = FUENTE, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.75, 0.3))) +
  labs(title = "Borde de 3 nudos", subtitle = "RMSE por velocidad no vista", x = "Velocidad (nudos)", y = "RMSE (puntos %)") +
  tema(8.5) + theme(legend.position = "none")
guardar(p10, "lamina10_borde_3_nudos.png", 2.43, 3.08)
log_msg("Paso 26 terminado. Solo grafica tablas ya generadas; los test no se utilizaron.")
