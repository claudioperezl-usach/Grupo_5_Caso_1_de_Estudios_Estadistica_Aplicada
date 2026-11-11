# =============================================================================
# 22_interpretacion_pdp_ice.R — Importancia por permutación, PDP e ICE (Elastic-net final)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/22_interpretacion_pdp_ice.R")
# Requiere los pasos 11 y 12 (modelo congelado y parámetros elegidos).
#
# Todo lo que sigue describe PATRONES PREDICTIVOS del modelo, no efectos causales:
#   * Importancia por permutación: cuánto sube el RMSE al romper la relación de una variable
#     con la respuesta. Con sensores muy correlacionados la importancia se reparte.
#   * PDP/ICE fijan una variable en valores que pueden no coexistir con las demás
#     (p. ej., cambiar T2 sin cambiar la velocidad ni las presiones).
#
# NO usa los test:
#   * Permutación: dentro de la validación LOSO del entrenamiento. En cada velocidad
#     excluida se reajusta el Elastic-net con los hiperparámetros finales y se permuta cada
#     variable SOLO en las filas de esa velocidad. La velocidad (v) no se permuta: es
#     constante dentro de cada velocidad excluida.
#   * PDP/ICE: modelo final congelado (paso 12) sobre una muestra de filas de ENTRENAMIENTO,
#     para las variables más importantes. La rejilla va del percentil 2,5 al 97,5.
# Salidas (resultados_elastic_net_final): tablas 22_*.csv y figuras 22_*.png.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto()
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 22: importancia por permutación, PDP e ICE ===")

ruta_modelo <- file.path(cfg$dir_final, "elastic_net_final.rds")
ruta_par    <- file.path(DIR_TABLAS, "12_parametros_elegidos.csv")
if (!file.exists(ruta_modelo) || !file.exists(ruta_par)) {
  stop("Ejecute primero el paso 12: faltan el modelo final o sus parámetros.")
}

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
rm(part, datos)                                   # los test no se usan en este paso
predictores <- predictores_modelo(cfg)
modelo <- readRDS(ruta_modelo)
par_df <- as.data.frame(data.table::fread(ruta_par))
verificar(all(c("penalty", "mixture") %in% names(par_df)), "el archivo de parámetros tiene penalty y mixture")
pars <- tibble::as_tibble(par_df[1, c("penalty", "mixture"), drop = FALSE])
vels <- sort(unique(entrena$v))

# --- 1. Importancia por permutación en cada velocidad excluida (LOSO) ---------------------
variables_perm <- setdiff(predictores, "v")
perm <- list()
for (vel in vels) {
  tr <- entrena[entrena$v != vel, , drop = FALSE]
  te <- entrena[entrena$v == vel, , drop = FALSE]
  verificar(!any(tr$v == vel), paste0("permutación, ", vel, " nudos: la velocidad evaluada no está en el ajuste"))
  wf_v <- tune::finalize_workflow(construir_wf("enet", tr, predictores, cfg), pars)
  set.seed(cfg$semilla + as.integer(vel))
  ajuste_v <- generics::fit(wf_v, data = tr)
  rmse_base <- rmse_vec(te$y, stats::predict(ajuste_v, new_data = te)$.pred)
  for (variable in variables_perm) {
    aumento <- vapply(seq_len(cfg$n_perm), function(b) {
      e2 <- te
      e2[[variable]] <- sample(e2[[variable]])
      rmse_vec(e2$y, stats::predict(ajuste_v, new_data = e2)$.pred) - rmse_base
    }, numeric(1))
    perm[[length(perm) + 1L]] <- data.frame(velocidad = vel, variable = variable,
                                            rmse_base = rmse_base, aumento_rmse = mean(aumento))
  }
  log_msg("  permutación terminada: velocidad excluida ", vel, " nudos")
}
perm <- data.table::rbindlist(perm)
importancia <- perm[, .(aumento_rmse_medio = mean(aumento_rmse),
                        de_entre_velocidades = stats::sd(aumento_rmse),
                        velocidades_donde_empeora = sum(aumento_rmse > 0)), by = variable]
data.table::setorder(importancia, -aumento_rmse_medio)
guardar_tabla(perm, "22_importancia_permutacion_detalle")
guardar_tabla(importancia, "22_importancia_permutacion")
log_msg("Variables más importantes: ", paste(head(importancia$variable, cfg$n_pdp_variables), collapse = ", "))

imp_f <- data.table::copy(importancia)
imp_f[, variable := factor(variable, levels = rev(variable))]
p1 <- ggplot(imp_f, aes(aumento_rmse_medio, variable)) +
  geom_col(width = 0.7, fill = COL$principal) +
  geom_errorbar(aes(xmin = aumento_rmse_medio - de_entre_velocidades,
                    xmax = aumento_rmse_medio + de_entre_velocidades),
                width = 0.25, colour = COL$tinta_suave, orientation = "y") +
  geom_vline(xintercept = 0, colour = COL$tinta_suave) +
  labs(title = "Importancia por permutación del Elastic-net",
       subtitle = "Aumento del RMSE al permutar la variable en la velocidad no vista (LOSO)\nBarra = ± 1 DE entre las 8 velocidades",
       x = "Aumento del RMSE (puntos %)", y = NULL,
       caption = "Asociación predictiva, no efecto causal. Sensores correlacionados reparten su importancia.") +
  tema_proyecto()
guardar_figura(p1, "22_importancia_permutacion", alto = 5)

# --- 2. PDP e ICE con el modelo final congelado, sobre filas de ENTRENAMIENTO --------------
vars_pdp <- head(importancia$variable, cfg$n_pdp_variables)
set.seed(cfg$semilla)
filas <- entrena[sample.int(nrow(entrena), min(cfg$n_ice, nrow(entrena))), , drop = FALSE]
ice <- list()
for (variable in vars_pdp) {
  q <- stats::quantile(entrena[[variable]], c(0.025, 0.975), names = FALSE)
  rejilla <- seq(q[1], q[2], length.out = cfg$pdp_puntos)
  e2 <- filas[rep(seq_len(nrow(filas)), times = length(rejilla)), , drop = FALSE]
  e2[[variable]] <- rep(rejilla, each = nrow(filas))
  pred <- as.numeric(stats::predict(modelo, new_data = e2)$.pred)
  verificar(length(pred) == nrow(e2) && all(is.finite(pred)),
            paste0("PDP/ICE de ", variable, ": predicciones válidas"))
  ice[[variable]] <- data.table::data.table(variable = variable,
                                            fila = rep(seq_len(nrow(filas)), times = length(rejilla)),
                                            valor = e2[[variable]], pred = pred)
}
ice <- data.table::rbindlist(ice)
pdp <- ice[, .(pred = mean(pred)), by = .(variable, valor)]
guardar_tabla(ice, "22_ice")
guardar_tabla(pdp, "22_pdp")

# Etiquetas con nombre y unidad (diccionario de 04_variables.R).
source(file.path("codigo", "R", "04_variables.R"), encoding = "UTF-8")
etiqueta <- function(v) {
  i <- match(v, DICCIONARIO$variable)
  paste0(toupper(v), ": ", DICCIONARIO$descripcion[i], " (", DICCIONARIO$unidad[i], ")")
}
niveles <- etiqueta(vars_pdp)
ice[, variable := factor(etiqueta(as.character(variable)), levels = niveles)]
pdp[, variable := factor(etiqueta(as.character(variable)), levels = niveles)]
p2 <- ggplot() +
  geom_line(data = ice, aes(valor, pred, group = fila), alpha = 0.08, colour = COL$principal) +
  geom_line(data = pdp, aes(valor, pred), colour = COL$secundario, linewidth = 1.2) +
  facet_wrap(~variable, scales = "free_x", labeller = ggplot2::label_wrap_gen(width = 40)) +
  labs(title = "PDP (naranja) e ICE (azul) del Elastic-net final",
       subtitle = paste0("Variables más importantes; ", nrow(filas), " filas de entrenamiento; rejilla entre los percentiles 2,5 y 97,5"),
       x = "Valor de la variable (unidad original)", y = "Pérdida predicha (puntos %)",
       caption = paste("Asociación predictiva, no causal. Predicciones fuera de 0-5 puntos: al fijar un sensor sin mover",
                       "\nlos demás se crean combinaciones que no existen en la planta (sensores muy correlacionados).")) +
  tema_proyecto()
guardar_figura(p2, "22_pdp_ice", alto = 6)

log_msg("Paso 22 terminado. Los test no se utilizaron.")
