# =============================================================================
# 24_multicolinealidad_shap.R — Multicolinealidad de los sensores y valores SHAP del
#                               Elastic-net final
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/24_multicolinealidad_shap.R")
# Requiere los pasos 12 (modelo congelado) y 18 (VIF). Usa SOLO el entrenamiento.
#
# 1. Multicolinealidad: matriz de correlaciones entre predictores (cor(), como en clases) y
#    VIF de cada predictor (regresión auxiliar sobre los demás, paso 18).
# 2. Valores SHAP: reparten cada predicción entre los predictores. En un modelo lineal son
#    exactos y no requieren simulación:  phi_ij = beta_j * (z_ij - media(z_j)),
#    con z = predictores estandarizados por la receta del modelo final. Se verifica que
#    intercepto medio + suma de los phi = predicción del modelo (aditividad).
#    Con sensores muy correlacionados, los SHAP de un modelo lineal heredan la inestabilidad
#    de los coeficientes: describen cómo el modelo reparte la predicción, no efectos causales.
# Salidas (resultados_elastic_net_final): tablas 24_*.csv y figuras 24_*.png.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 24: multicolinealidad y valores SHAP (solo entrenamiento) ===")

ruta_modelo <- file.path(cfg$dir_final, "elastic_net_final.rds")
ruta_par    <- file.path(DIR_TABLAS, "12_parametros_elegidos.csv")
ruta_vif    <- file.path(DIR_TABLAS, "05_ols_vif.csv")
if (!all(file.exists(c(ruta_modelo, ruta_par, ruta_vif)))) {
  stop("Ejecute primero los pasos 12 y 18: faltan el modelo final, sus parámetros o el VIF.")
}

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
rm(part, datos)                                   # los test no se usan en este paso
predictores <- predictores_modelo(cfg)
etiqueta <- function(v) toupper(v)

# --- 1. Matriz de correlaciones entre predictores --------------------------------------
R <- stats::cor(as.matrix(entrena[, predictores]))
orden <- predictores[stats::hclust(stats::as.dist(1 - abs(R)))$order]   # agrupa sensores afines
cor_largo <- data.table::as.data.table(as.table(R))
data.table::setnames(cor_largo, c("var1", "var2", "correlacion"))
guardar_tabla(cor_largo, "24_correlaciones_predictores")
n_altas <- sum(abs(R[upper.tri(R)]) >= 0.9)
log_msg(sprintf("Pares de predictores con |correlación| >= 0,9: %d de %d", n_altas, sum(upper.tri(R))))

cor_largo[, `:=`(var1 = factor(etiqueta(var1), levels = etiqueta(orden)),
                 var2 = factor(etiqueta(var2), levels = rev(etiqueta(orden))))]
p_cor <- ggplot(cor_largo, aes(var1, var2, fill = correlacion)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", correlacion),
                colour = abs(correlacion) > 0.6), size = 3.3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = COL$tinta)) +
  scale_fill_gradient2(low = COL$secundario, mid = "#f4f4f1", high = COL$principal,
                       midpoint = 0, limits = c(-1, 1), name = "Correlación") +
  coord_fixed() +
  labs(title = "Los sensores se mueven casi juntos",
       subtitle = sprintf("Correlación de Pearson entre predictores (entrenamiento, %s filas)\n%d de %d pares con |r| ≥ 0,9",
                          format(nrow(entrena), big.mark = ".", decimal.mark = ","), n_altas, sum(upper.tri(R))),
       x = NULL, y = NULL,
       caption = paste("Correlación global, dominada por la velocidad: todos los sensores suben con ella.",
                       "\nCoeficientes individuales inestables: motiva Ridge, Lasso y Elastic-net.")) +
  tema_proyecto() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
guardar_figura(p_cor, "24_correlaciones_predictores", ancho = 6.8, alto = 6.4)

# --- 1b. VIF (paso 18, OLS con todo el entrenamiento) ---------------------------------
vif <- data.table::fread(ruta_vif)
verificar(setequal(vif$variable, predictores), "el VIF del paso 18 cubre los mismos predictores")
vif[, variable := factor(etiqueta(variable), levels = etiqueta(variable[order(VIF)]))]
p_vif <- ggplot(vif, aes(VIF, variable)) +
  geom_col(width = 0.7, fill = COL$principal) +
  geom_vline(xintercept = 10, colour = COL$secundario, linetype = "dashed") +
  geom_text(aes(label = format(round(VIF), big.mark = ".", decimal.mark = ",")), hjust = -0.15, size = 3, colour = COL$tinta) +
  scale_x_log10(labels = function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE),
                expand = ggplot2::expansion(mult = c(0, 0.18))) +
  labs(title = "Inflación de la varianza de OLS (VIF)",
       subtitle = "VIF = 1 / (1 - R²) de cada predictor sobre los demás; escala logarítmica",
       x = "VIF (log)", y = NULL,
       caption = "Línea discontinua: VIF = 10, nivel usual de alerta. Todos los predictores lo superan.") +
  tema_proyecto()
guardar_figura(p_vif, "24_vif", ancho = 7, alto = 4.5)

# --- 2. Valores SHAP exactos del Elastic-net final --------------------------------------
modelo <- readRDS(ruta_modelo)
pars <- data.table::fread(ruta_par)
receta <- workflows::extract_recipe(modelo)
Z <- as.matrix(recipes::bake(receta, new_data = entrena, recipes::all_predictors()))
beta_mat <- as.matrix(stats::coef(workflows::extract_fit_engine(modelo), s = pars$penalty[1]))
beta <- beta_mat[colnames(Z), 1]
centro <- colMeans(Z)
phi <- sweep(sweep(Z, 2, centro), 2, beta, `*`)                 # phi_ij = beta_j (z_ij - media_j)
base <- beta_mat["(Intercept)", 1] + sum(beta * centro)        # predicción media
pred <- as.numeric(stats::predict(modelo, new_data = entrena)$.pred)
verificar(max(abs(base + rowSums(phi) - pred)) < 1e-8,
          "aditividad SHAP: predicción media + suma de los SHAP = predicción del modelo")

importancia <- data.table::data.table(variable = colnames(phi),
                                      shap_abs_medio = colMeans(abs(phi)),
                                      coeficiente_estandarizado = beta)
data.table::setorder(importancia, -shap_abs_medio)
guardar_tabla(importancia, "24_shap_importancia")
log_msg("SHAP |medio| (puntos %): ", paste(sprintf("%s %.3f", importancia$variable, importancia$shap_abs_medio), collapse = " | "))
log_msg(sprintf("Predicción media (base de los SHAP) = %.3f puntos %%", base))

# Gráfico tipo enjambre sobre una muestra de filas de entrenamiento.
set.seed(cfg$semilla)
filas <- sample.int(nrow(phi), min(2000L, nrow(phi)))
niveles <- etiqueta(rev(importancia$variable))
puntos <- data.table::rbindlist(lapply(colnames(phi), function(v) data.table::data.table(
  variable = etiqueta(v), shap = phi[filas, v],
  valor = (Z[filas, v] - min(Z[, v])) / diff(range(Z[, v])))))      # 0 = mínimo, 1 = máximo
puntos[, variable := factor(variable, levels = niveles)]
guardar_tabla(puntos, "24_shap_muestra")
p_shap <- ggplot(puntos, aes(shap, variable, colour = valor)) +
  geom_vline(xintercept = 0, colour = COL$tinta_suave) +
  geom_jitter(height = 0.22, width = 0, size = 0.45, alpha = 0.5) +
  scale_colour_gradient(low = COL$secundario, high = COL$principal, breaks = c(0, 1),
                        labels = c("bajo", "alto"), name = "Valor del\npredictor") +
  labs(title = "Valores SHAP del Elastic-net final",
       subtitle = "Aporte de cada predictor a la predicción respecto de la media (2.000 filas de entrenamiento)",
       x = "Valor SHAP (puntos % de pérdida)", y = NULL,
       caption = paste("Modelo lineal: SHAP = coeficiente x (valor - media). Ordenados por |SHAP| medio.",
                       "\nAportes de hasta ±15 puntos que se compensan (la pérdida va de 0 a 5): sensores correlacionados.",
                       "\nPatrones del modelo, no efectos causales.")) +
  tema_proyecto() + theme(legend.position = "right")
guardar_figura(p_shap, "24_shap_enjambre", ancho = 8.5, alto = 5.2)

log_msg("Paso 24 terminado. Los test no se utilizaron.")
