# =============================================================================
# 21_tabla1_final.R — Tabla 1 del informe: validación LOSO con incertidumbre
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/21_tabla1_final.R")
# Consolida las métricas POR VELOCIDAD de:
#   paso 11 (anidada: Elastic-net, Ridge, Lasso, Splines, baseline; y ensembles si se anidaron)
#   paso 19 (Árbol fijo y Promedio Elastic-net + Ridge)
#   paso 20 (OLS, GAM, Random Forest y XGBoost)
# y agrega la incertidumbre que pide la pauta: desviación estándar entre velocidades y
# diferencia PAREADA de RMSE frente al baseline con IC 95 % (t de Student, 7 g.l.).
# Supuesto del IC: las 8 diferencias por velocidad se tratan como independientes. Las
# particiones LOSO comparten 6 de 7 velocidades de ajuste, por lo que el IC es aproximado
# (tiende a ser estrecho); se usa para comparar modelos, no como garantía formal.
# No reentrena ni inventa métricas: si falta un archivo, se detiene.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto(character(0))
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 21: Tabla 1 consolidada ===")

rutas <- c(anidada = file.path(cfg$dir_anidada, "tablas", "11_detalle_anidado.csv"),
           arbol   = file.path(DIR_TABLAS, "19_detalle_arbol_promedio.csv"),
           otros   = file.path(DIR_TABLAS, "20_detalle_ols_rf_xgb.csv"))
if (any(!file.exists(rutas))) stop("Faltan archivos:\n", paste(rutas[!file.exists(rutas)], collapse = "\n"))
cols <- c("velocidad", "modelo", "rmse", "mae", "r2", "rmse_critico")
a <- data.table::fread(rutas["anidada"])[, ..cols][, evaluacion := "LOSO anidada (1-SE interno)"]
b <- data.table::fread(rutas["arbol"])[modelo %in% c("Arbol fijo", "Promedio ENet+Ridge"), ..cols]
b[, evaluacion := ifelse(modelo == "Arbol fijo", "LOSO, configuración fija", "LOSO (miembros con 1-SE interno)")]
o <- data.table::fread(rutas["otros"])[, ..cols]
o[, evaluacion := data.table::fcase(modelo == "OLS", "LOSO (sin hiperparámetros)",
                                    modelo == "GAM", "LOSO (suavidad por REML en cada partición)",
                                    default = "LOSO, configuración fija")]
d <- data.table::rbindlist(list(a, b, o))
# Si los ensembles se anidaron (paso 11), esa versión reemplaza a la de configuración fija.
nombres <- c("Baseline global" = "Media global (baseline)", "Arbol fijo" = "Árbol de regresión",
             "Árbol de regresión" = "Árbol de regresión", "Promedio ENet+Ridge" = "Promedio EN + Ridge",
             "Splines naturales" = "Splines naturales")
d[modelo %in% names(nombres), modelo := nombres[modelo]]
d <- d[order(modelo, evaluacion != "LOSO anidada (1-SE interno)")]
d <- d[, .SD[evaluacion == evaluacion[1]], by = modelo]
verificar(all(d[, .N, by = modelo]$N == 8L), "cada modelo de la Tabla 1 tiene 8 velocidades evaluadas")

familia <- c("Media global (baseline)" = "—", "OLS" = "Lineal", "Ridge" = "Regularizado",
             "Lasso" = "Regularizado", "Elastic-net" = "Regularizado", "Splines naturales" = "Flexible", "GAM" = "Flexible (GAM)",
             "Árbol de regresión" = "Árbol", "Random Forest" = "Ensemble (bagging)",
             "XGBoost" = "Ensemble (boosting)", "Promedio EN + Ridge" = "Ensemble (averaging)")
base <- d[modelo == "Media global (baseline)", .(velocidad, rmse_base = rmse)]
d <- merge(d, base, by = "velocidad")
tq <- stats::qt(0.975, df = 7)
tabla <- d[, .(familia = familia[[modelo[1]]], evaluacion = evaluacion[1],
               rmse = mean(rmse), de_rmse = stats::sd(rmse), mae = mean(mae), r2 = mean(r2),
               rmse_critico = mean(rmse_critico, na.rm = TRUE),
               mejora_pct = 100 * (1 - mean(rmse) / mean(rmse_base)),
               dif_vs_base = mean(rmse - rmse_base),
               dif_ic_inf = mean(rmse - rmse_base) - tq * stats::sd(rmse - rmse_base) / sqrt(8),
               dif_ic_sup = mean(rmse - rmse_base) + tq * stats::sd(rmse - rmse_base) / sqrt(8),
               gana_baseline = sum(rmse < rmse_base)), by = modelo]
tabla[modelo == "Media global (baseline)", c("dif_vs_base", "dif_ic_inf", "dif_ic_sup") := NA_real_]
data.table::setorder(tabla, rmse)
data.table::fwrite(tabla, file.path(DIR_TABLAS, "21_tabla1_validacion.csv"))

# Comparación pareada entre los dos mejores candidatos anidados (incertidumbre de la elección).
w <- data.table::dcast(d[modelo %in% c("Elastic-net", "Ridge")], velocidad ~ modelo, value.var = "rmse")
dif <- w[["Elastic-net"]] - w[["Ridge"]]
pareada <- data.frame(comparacion = "Elastic-net - Ridge", dif_media = mean(dif),
                      ic_inf = mean(dif) - tq * stats::sd(dif) / sqrt(8),
                      ic_sup = mean(dif) + tq * stats::sd(dif) / sqrt(8),
                      velocidades_en_que_gana_EN = sum(dif < 0))
data.table::fwrite(pareada, file.path(DIR_TABLAS, "21_comparacion_en_ridge.csv"))

# Versión para el informe (3 decimales).
pres <- tabla[, .(Modelo = modelo, Familia = familia, Validacion = evaluacion,
                  RMSE = sprintf("%.3f ± %.3f", rmse, de_rmse), MAE = sprintf("%.3f", mae),
                  R2 = sprintf("%.3f", r2), `RMSE subgrupo (y>=4)` = sprintf("%.3f", rmse_critico),
                  `Mejora vs baseline (%)` = sprintf("%.1f", mejora_pct),
                  `Delta RMSE vs baseline [IC 95%]` = ifelse(is.na(dif_vs_base), "referencia",
                       sprintf("%.3f [%.3f; %.3f]", dif_vs_base, dif_ic_inf, dif_ic_sup)))]
data.table::fwrite(pres, file.path(DIR_TABLAS, "21_tabla1_presentacion.csv"))
log_msg("Tabla 1: ", paste(sprintf("%s %.3f", tabla$modelo, tabla$rmse), collapse = " | "))
log_msg(sprintf("Elastic-net - Ridge: %.3f [%.3f; %.3f]", pareada$dif_media, pareada$ic_inf, pareada$ic_sup))
log_msg("Paso 21 terminado.")
