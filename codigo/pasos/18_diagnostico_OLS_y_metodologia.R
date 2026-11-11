# =============================================================================
# 18_diagnostico_OLS_y_metodologia.R — Diagnóstico de OLS: residuos, influencia y VIF
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/18_diagnostico_OLS_y_metodologia.R")
# 1. OLS con todo el entrenamiento: residuos vs ajustados, distancia de Cook y VIF
#    (diagnostico_ols, 06_modelos.R). Origen del "VIF máximo" del informe.
# 2. OLS en cada partición LOSO: error en la velocidad excluida y VIF por partición.
# Usa solo el entrenamiento; no abre los test.
# =============================================================================
if (!file.exists('codigo/R/00_config.R')) stop('Ejecute desde la raiz de equipo_regresion_naval')
source('codigo/R/00_config.R',encoding='UTF-8')
source('codigo/R/01_utilidades.R',encoding='UTF-8')
cargar_proyecto(c('02_datos.R','03_particion.R','04_variables.R','06_modelos.R'))
iniciar_salida(cfg$dir_final)   # define RUTA_LOG/DIR_TABLAS: ejecutable desde sesion limpia
set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos,cfg)
train <- as.data.frame(part$entrena)
predictoras <- predictores_modelo(cfg)
stopifnot(all(c('y','v',predictoras) %in% names(train)))
# 1. Diagnóstico con todo el entrenamiento (tablas 05_ols_* y figura de residuos).
diagnostico_ols(part$entrena,predictoras,cfg)
# VIF se calcula en el entrenamiento de cada fold, no en los datos de prueba.
# Regresión auxiliar, VIF infinito si R2 ~ 1.
calcular_vif <- function(d,vars) {
  sapply(vars,function(variable) {
    otros <- setdiff(vars,variable)
    m <- lm(reformulate(otros,response=variable),data=d)
    r2 <- summary(m)$r.squared
    if (!is.finite(r2) || r2 >= 1-1e-12) Inf else 1/(1-r2)
  })
}
# 2. OLS en cada partición LOSO (error en la velocidad excluida y VIF por partición).
vels <- sort(unique(train$v))
por_vel <- vector('list',length(vels))
vifs <- vector('list',length(vels))
for (i in seq_along(vels)) {
  v <- vels[i]
  tr <- train[train$v!=v,,drop=FALSE]
  te <- train[train$v==v,,drop=FALSE]
  formula <- reformulate(predictoras,response='y')
  m <- lm(formula,data=tr)
  p <- suppressWarnings(predict(m,newdata=te))
  co <- coef(m)
  vif <- calcular_vif(tr,predictoras)
  por_vel[[i]] <- data.frame(velocidad=v,n=nrow(te),rmse_OLS=rmse_vec(te$y,p),
    rmse_baseline=rmse_vec(te$y,rep(mean(tr$y),nrow(te))),
    rango_deficiente=anyNA(co),coeficientes_NA=sum(is.na(co)),
    rango_modelo=m$rank,columnas_diseno=ncol(model.matrix(m)),
    vif_max=max(vif),coeficiente_abs_max=max(abs(co),na.rm=TRUE))
  vifs[[i]] <- data.frame(velocidad=v,variable=names(vif),vif=as.numeric(vif))
}
por_vel <- do.call(rbind,por_vel)
vifs <- do.call(rbind,vifs)
guardar_tabla(por_vel,'18_diagnostico_ols_por_velocidad')
guardar_tabla(vifs,'18_vif_por_velocidad')
log_msg(sprintf('OLS LOSO: RMSE medio = %.3f; VIF máximo entre particiones = %.1f; particiones con coeficientes NA = %d',
                mean(por_vel$rmse_OLS),max(vifs$vif),sum(por_vel$coeficientes_NA>0)))
log_msg('Paso 18 terminado. Los test no se utilizaron.')
