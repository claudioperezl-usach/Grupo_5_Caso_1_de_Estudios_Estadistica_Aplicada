# =============================================================================
# 17_comprobaciones_adicionales.R — Convergencia y estabilidad de coeficientes (Elastic-net)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/17_comprobaciones_adicionales.R")
# 1. Reajuste estricto de glmnet (tolerancia 1e-10) con alpha/lambda finales: comprueba que
#    el workflow final convergió (sobre el ENTRENAMIENTO; los test no se usan).
# 2. Estabilidad de coeficientes: mismo alpha/lambda en las 8 particiones LOSO.
# NO reemplaza la validación anidada ni cambia el modelo elegido.
# =============================================================================
if (!file.exists('codigo/R/00_config.R')) stop('Ejecute desde la raiz del proyecto.')
source('codigo/R/00_config.R',encoding='UTF-8')
source('codigo/R/01_utilidades.R',encoding='UTF-8')
cargar_proyecto(c('02_datos.R','03_particion.R','04_variables.R'))
iniciar_salida(cfg$dir_final)   # define RUTA_LOG/DIR_TABLAS: ejecutable desde sesion limpia
for (pkg in c('data.table','recipes','glmnet','workflows')) {
  if (!requireNamespace(pkg,quietly=TRUE)) stop('Instale el paquete: ',pkg)
}
DIR <- DIR_TABLAS
dir.create(DIR,recursive=TRUE,showWarnings=FALSE)
set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos,cfg)
train <- as.data.frame(part$entrena)
# La comprobacion de convergencia se hace sobre el ENTRENAMIENTO: los test no se usan aqui.
tests <- list(entrenamiento=train)
predictoras <- predictores_modelo(cfg)
stopifnot(length(unique(train$v))==8,all(c('y','v',predictoras) %in% names(train)))
modelo <- readRDS(file.path(cfg$dir_final,'elastic_net_final.rds'))
param <- read.csv(file.path(DIR,'12_parametros_elegidos.csv'))
lam <- as.numeric(param$penalty[1]); alpha <- as.numeric(param$mixture[1])
stopifnot(is.finite(lam),lam>0,is.finite(alpha))
log_msg('Hiperparámetros finales: alpha = ',alpha,'; lambda = ',signif(lam,4))
borde_inf <- isTRUE(all.equal(lam,min(cfg$grilla_lambda_loso)))
log_msg('Grilla de lambda: ',format(min(cfg$grilla_lambda_loso)),' a ',format(max(cfg$grilla_lambda_loso)),
        '. Lambda en el borde inferior: ',borde_inf)
if (borde_inf) warning('Lambda en el borde inferior de la grilla: amplie cfg$grilla_lambda_loso hacia valores menores.')
# Escalado dentro de cada conjunto de entrenamiento; NO escalar antes de particionar.
escalar <- function(tr,te) {
  r <- recipes::recipe(stats::reformulate(predictoras,response='y'),data=tr)
  r <- recipes::step_normalize(r,recipes::all_predictors())
  r <- recipes::prep(r,training=tr,retain=TRUE)
  a <- recipes::bake(r,new_data=tr)
  b <- recipes::bake(r,new_data=te)
  list(x=as.matrix(a[,predictoras,drop=FALSE]),xt=as.matrix(b[,predictoras,drop=FALSE]))
}
# 1. Ajuste estricto con la misma alpha y lambda; comparación contra workflow original.
comp <- list()
for (nombre in names(tests)) {
  te <- tests[[nombre]]
  z <- escalar(train,te)
  fit2 <- glmnet::glmnet(z$x,train$y,alpha=alpha,lambda=lam,
                          thresh=1e-10,maxit=10000000L,standardize=TRUE)
  pred2 <- as.numeric(stats::predict(fit2,newx=z$xt,s=lam))
  pred1 <- as.numeric(stats::predict(modelo,new_data=te)$.pred)
  comp[[nombre]] <- data.frame(conjunto=nombre,n=nrow(te),
    rmse_original=rmse_vec(te$y,pred1),rmse_estricto=rmse_vec(te$y,pred2),
    diferencia_rmse=rmse_vec(te$y,pred2)-rmse_vec(te$y,pred1),
    max_diferencia_pred=max(abs(pred2-pred1)),jerr_estricto=fit2$jerr)
}
comp <- do.call(rbind,comp)
guardar_tabla(comp,'17_convergencia_estricta')
log_msg(sprintf('Convergencia estricta (tol. 1e-10): diferencia de RMSE = %.5f; diferencia máxima de predicción = %.4f',
                comp$diferencia_rmse[1],comp$max_diferencia_pred[1]))
# 2. Estabilidad de coeficientes: mismo alpha/lambda en los ocho folds externos.
# No son hiperparámetros elegidos dentro de cada fold: solo diagnostico descriptivo.
vels <- sort(unique(train$v))
coef_mat <- matrix(NA_real_,nrow=length(predictoras),ncol=length(vels),
                   dimnames=list(predictoras,paste0('vel_',vels)))
for (i in seq_along(vels)) {
  tr <- train[train$v!=vels[i],,drop=FALSE]
  z <- escalar(tr,tr)
  m <- glmnet::glmnet(z$x,tr$y,alpha=alpha,lambda=lam,
                      thresh=1e-10,maxit=10000000L,standardize=TRUE)
  if (m$jerr!=0) warning('Convergencia fold externo: ',vels[i],', jerr=',m$jerr)
  c0 <- as.matrix(stats::coef(m,s=lam))
  coef_mat[,i] <- as.numeric(c0[predictoras,1])
}
coef_det <- data.frame(variable=rownames(coef_mat),coef_mat,check.names=FALSE)
coef_res <- data.frame(variable=rownames(coef_mat),minimo=apply(coef_mat,1,min),
 maximo=apply(coef_mat,1,max),n_positivos=rowSums(coef_mat>1e-8),
 n_negativos=rowSums(coef_mat< -1e-8),n_ceros=rowSums(abs(coef_mat)<=1e-8))
coef_res$cambia_signo <- coef_res$n_positivos>0 & coef_res$n_negativos>0
guardar_tabla(coef_det,'17_coeficientes_por_fold')
guardar_tabla(coef_res,'17_estabilidad_coeficientes')
log_msg('Variables que cambian de signo entre particiones: ',
        if (any(coef_res$cambia_signo)) paste(coef_res$variable[coef_res$cambia_signo],collapse=', ') else 'ninguna')
log_msg('Paso 17 terminado. Las métricas LOSO de OLS, RF y XGBoost se calculan en el paso 20.')
