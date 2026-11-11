# =============================================================================
# 20_loso_ols_rf_xgb.R — Métricas LOSO de OLS, GAM, Random Forest y XGBoost
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/20_loso_ols_rf_xgb.R")
# Reemplaza a los antiguos scripts 20 y 21 (que reconstruían la Tabla 1 con "ND" y
# exigían un orden manual de ejecución). Calcula RMSE, MAE, R2 y RMSE del subgrupo
# crítico por velocidad excluida.
# OLS no tiene hiperparámetros: su LOSO es directamente comparable con la anidada.
# GAM (mgcv, como en clases): s() con splines cúbicos de regresión (bs = "cr") en los mismos
# 7 sensores que los splines naturales y efecto lineal en el resto. La suavidad de cada s()
# se estima por REML SOLO con las 7 velocidades de ajuste, así que no requiere CV interna
# ni usa la velocidad evaluada.
# RF y XGBoost con configuración fija (cfg$rf_arboles, cfg$xgb_fijo) SOLO si
# cfg$anidar_ensembles = FALSE; con TRUE, sus versiones ajustadas por LOSO anidado salen del paso 11.
# =============================================================================

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto(c("02_datos.R", "03_particion.R", "04_variables.R"))
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 20: métricas LOSO de OLS, GAM, Random Forest y XGBoost ===")
fijos <- !isTRUE(cfg$anidar_ensembles)    # RF y XGBoost fijos solo si no se anidan

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos, cfg)
train <- as.data.frame(part$entrena)
rm(part, datos)
predictoras <- predictores_modelo(cfg)
vels <- sort(unique(train$v))
f <- stats::reformulate(predictoras, response = "y")
vars_suaves <- intersect(c("gtt", "ggn", "t48", "t2", "p2", "mf", "tic"), predictoras)
f_gam <- stats::reformulate(c(sprintf('s(%s, bs = "cr")', vars_suaves), setdiff(predictoras, vars_suaves)),
                            response = "y")
log_msg("Fórmula GAM: ", deparse1(f_gam))

res <- list()
for (vel in vels) {
  tr <- train[train$v != vel, , drop = FALSE]
  te <- train[train$v == vel, , drop = FALSE]
  p_ols <- as.numeric(stats::predict(stats::lm(f, data = tr), newdata = te))
  gam_v <- mgcv::gam(f_gam, data = tr, method = "REML")
  p_gam <- as.numeric(stats::predict(gam_v, newdata = te))
  modelos <- list(list("OLS", p_ols), list("GAM", p_gam))
  if (fijos) {
    set.seed(cfg$semilla + as.integer(vel))
    rf <- ranger::ranger(f, data = tr, num.trees = cfg$rf_arboles, num.threads = cfg$n_hilos,
                         seed = cfg$semilla + as.integer(vel))
    p_rf <- as.numeric(predict(rf, data = te)$predictions)
    # xgb.train + xgb.DMatrix: misma interfaz en xgboost 1.7 y 3.x (xgboost() cambió en 3.0).
    # tree_method y base_score se fijan explícitamente: sus valores por defecto cambiaron entre
    # versiones (1.7: "exact" y 0.5; >= 2.0: "hist" y media de y), lo que altera el RMSE.
    # base_score = media de y en el AJUSTE de la partición (no usa la velocidad evaluada).
    xt <- xgboost::xgb.DMatrix(as.matrix(tr[, predictoras, drop = FALSE]), label = tr$y)
    xv <- as.matrix(te[, predictoras, drop = FALSE])
    set.seed(cfg$semilla + as.integer(vel))
    xgb <- xgboost::xgb.train(params = list(objective = "reg:squarederror", eta = cfg$xgb_fijo$eta,
                                            max_depth = cfg$xgb_fijo$profundidad,
                                            subsample = 1, colsample_bytree = 1, nthread = cfg$n_hilos,
                                            tree_method = "hist", max_bin = 256, base_score = mean(tr$y)),
                              data = xt, nrounds = cfg$xgb_fijo$rondas, verbose = 0)
    p_xgb <- as.numeric(predict(xgb, xgboost::xgb.DMatrix(xv)))
    modelos <- c(modelos, list(list("Random Forest", p_rf), list("XGBoost", p_xgb)))
  }
  for (m in modelos) {
    verificar(length(m[[2]]) == nrow(te) && all(is.finite(m[[2]])),
              paste0(m[[1]], ": predicciones válidas a ", vel, " nudos"))
    res[[length(res) + 1L]] <- metricas_fila(vel, m[[1]], te$y, m[[2]])
  }
}
det <- do.call(rbind, res)
data.table::fwrite(det, file.path(DIR_TABLAS, "20_detalle_ols_rf_xgb.csv"))

met <- data.table::as.data.table(det)[, .(rmse = mean(rmse), mae = mean(mae)), by = modelo]
log_msg("LOSO (media de 8 velocidades): ", paste(sprintf("%s RMSE = %.3f", met$modelo, met$rmse), collapse = " | "))
log_msg("Paso 20 terminado.")
