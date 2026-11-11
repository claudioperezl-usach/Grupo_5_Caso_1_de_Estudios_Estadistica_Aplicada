# =============================================================================
# 14_diagnosticos_elastic_net.R — Comparación por velocidad, coeficientes y error por banda
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/14_diagnosticos_elastic_net.R")
# Solo LEE resultados de los pasos 11 y 12; NO reentrena ni modifica el modelo.
# Los coeficientes están en predictores estandarizados por la receta: puntos % de pérdida
# por cada desviación estándar del sensor (origen de los coeficientes del informe, §6).
# =============================================================================

if (!file.exists('codigo/R/00_config.R')) stop('Ejecute desde la carpeta raiz del proyecto.')
source('codigo/R/00_config.R', encoding = 'UTF-8')
source('codigo/R/01_utilidades.R', encoding = 'UTF-8')
cargar_proyecto(character(0))
iniciar_salida(cfg$dir_final)
archivo_detalle <- file.path(cfg$dir_anidada, 'tablas', '11_detalle_anidado.csv')
archivo_modelo <- file.path(cfg$dir_final, 'elastic_net_final.rds')
archivo_parametros <- file.path(cfg$dir_final, 'tablas', '12_parametros_elegidos.csv')
archivo_predicciones <- file.path(cfg$dir_final, 'tablas', '12_predicciones_test.csv')
salida <- file.path(cfg$dir_final, 'tablas')
dir.create(salida, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(archivo_detalle)) stop('Primero ejecute el script 11: falta 11_detalle_anidado.csv')
if (!file.exists(archivo_modelo) || !file.exists(archivo_parametros))
  stop('Primero ejecute el script 12: faltan el modelo final o sus parametros.')

detalle <- read.csv(archivo_detalle, check.names=FALSE)
requeridas <- c('velocidad','modelo','rmse','mae')
if (!all(requeridas %in% names(detalle))) stop('Columnas inesperadas en 11_detalle_anidado.csv')
if (anyDuplicated(detalle[c('velocidad','modelo')])) stop('Hay duplicados velocidad-modelo en la validacion.')

modelos <- unique(detalle$modelo)
if (!all(c('Elastic-net','Ridge','Baseline global') %in% modelos))
  stop('Faltan Elastic-net, Ridge o Baseline global en el archivo 11.')

# 1. Tabla por velocidad: comparar errores PAREADOS del mismo fold externo.
sub <- detalle[detalle$modelo %in% c('Elastic-net','Ridge','Baseline global'),]
ancha <- reshape(sub[c('velocidad','modelo','rmse')],idvar='velocidad',timevar='modelo',direction='wide')
names(ancha) <- sub('^rmse\\.','rmse_',names(ancha))
ancha <- ancha[order(ancha$velocidad),]
ancha$tipo <- ifelse(ancha$velocidad %in% c(min(ancha$velocidad),max(ancha$velocidad)),
                     'Extrapolacion de borde','Interpolacion interior')
ancha$mejora_vs_baseline_pct <- 100*(1-ancha[['rmse_Elastic-net']]/ancha[['rmse_Baseline global']])
ancha$gana_al_baseline <- ancha[['rmse_Elastic-net']] < ancha[['rmse_Baseline global']]
ancha$gana_a_ridge <- ancha[['rmse_Elastic-net']] < ancha[['rmse_Ridge']]
guardar_tabla(ancha,'14_comparacion_por_velocidad')
log_msg('Elastic-net gana al baseline en ',sum(ancha$gana_al_baseline),' de ',nrow(ancha),
        ' velocidades y a Ridge en ',sum(ancha$gana_a_ridge),'.')
log_msg('Diferencia pareada media de RMSE (Elastic-net menos Ridge): ',
    round(mean(ancha[['rmse_Elastic-net']]-ancha[['rmse_Ridge']]),4))

# 2. Hiperparámetros y coeficientes: receta de normalización dentro del workflow.
pars <- read.csv(archivo_parametros,check.names=FALSE)
log_msg('Hiperparámetros elegidos: penalty = ',signif(pars$penalty[1],4),'; mixture = ',pars$mixture[1])
fit_wf <- readRDS(archivo_modelo)
if (!requireNamespace('workflows',quietly=TRUE)) stop('Instale el paquete workflows.')
if (!requireNamespace('glmnet',quietly=TRUE)) stop('Instale el paquete glmnet.')
fit <- workflows::extract_fit_parsnip(fit_wf)
# La receta estandariza predictores; coeficientes de glmnet se expresan
# en las variables ya procesadas por la receta.
coeficientes <- NULL
# En algunos parsnip, penalty es una expresión; usar el parámetro elegido del CSV.
if ('penalty' %in% names(pars)) {
  lambda <- as.numeric(pars$penalty[1]);
  coeficientes <- as.matrix(stats::coef(fit$fit,s=lambda))
} else {
  stop('No se encontro penalty en 12_parametros_elegidos.csv.')
}
coefs <- data.frame(variable=rownames(coeficientes),coeficiente=as.numeric(coeficientes[,1]))
coefs$seleccionada <- coefs$variable!='(Intercept)' & abs(coefs$coeficiente)>1e-8
coefs <- coefs[order(-abs(coefs$coeficiente)),]
guardar_tabla(coefs,'14_coeficientes_elastic_net')
log_msg('Coeficientes estandarizados (puntos % por DE): ',
        paste(sprintf('%s %+.2f',coefs$variable[coefs$seleccionada],coefs$coeficiente[coefs$seleccionada]),collapse=' | '))
log_msg('Variables con coeficiente no nulo: ',sum(coefs$seleccionada),'. Su estabilidad entre particiones está en el paso 17.')

# 3. Advertencias de convergencia, sin suponer que afectaron al lambda elegido.
notas_rds <- file.path(salida,'12_advertencias_tuning.rds')
notas_csv <- file.path(salida,'12_advertencias_tuning.csv')
if (file.exists(notas_rds)) {
  notas <- readRDS(notas_rds)
  log_msg('Advertencias del tuning guardadas: ',nrow(notas),' (convergencia de glmnet en lambdas pequeños; ver README).')
} else if (file.exists(notas_csv)) {
  notas <- read.csv(notas_csv)
  log_msg('Advertencias del tuning guardadas: ',nrow(notas))
} else log_msg('No hay advertencias del tuning guardadas.')

# 4. Test aleatorio y 15 nudos: resumen por bandas de degradación.
if (file.exists(archivo_predicciones)) {
  p <- read.csv(archivo_predicciones)
  if (all(c('conjunto','y','pred') %in% names(p))) {
    p$banda_y <- cut(p$y,breaks=c(-Inf,1,2,3,4,Inf),right=FALSE,
                     labels=c('[0,1)','[1,2)','[2,3)','[3,4)','[4,5+]'))
    grupos <- split(p,list(p$conjunto,p$banda_y),drop=TRUE)
    por_banda <- do.call(rbind,lapply(grupos,function(d) data.frame(
      conjunto=d$conjunto[1],banda_y=as.character(d$banda_y[1]),n=nrow(d),
      rmse=sqrt(mean((d$y-d$pred)^2)),mae=mean(abs(d$y-d$pred)))))
    rownames(por_banda) <- NULL
    guardar_tabla(por_banda,'14_error_por_degradacion')
  }
}
log_msg('Paso 14 terminado.')
