# =============================================================================
# 19_arbol_y_promedio_loso.R — Árbol de regresión y promedio Elastic-net + Ridge (LOSO)
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/19_arbol_y_promedio_loso.R")
# Árbol con configuración fija + promedio 50/50 (model averaging, como en clases) de
# Elastic-net y Ridge, reajustados con los hiperparámetros 1-SE elegidos dentro de cada
# partición externa por el paso 11. NO modifica el modelo final ni abre los test.
# =============================================================================

if (!file.exists('codigo/R/00_config.R')) stop('Abra el proyecto desde su carpeta raiz.')
source('codigo/R/00_config.R',encoding='UTF-8')
source('codigo/R/01_utilidades.R',encoding='UTF-8')
cargar_proyecto()
iniciar_salida(cfg$dir_final)   # define RUTA_LOG/DIR_TABLAS: ejecutable desde sesion limpia
for (pkg in c('rpart','tibble','workflows','tune','parsnip','recipes','generics')) {
  if (!requireNamespace(pkg,quietly=TRUE)) stop('Falta paquete: ',pkg)
}
archivo_par <- file.path(cfg$dir_anidada,'tablas','11_parametros_internos.csv')
if (!file.exists(archivo_par)) stop('Falta ',archivo_par, '. Ejecute primero el script 11.')
pars <- read.csv(archivo_par,check.names=FALSE,stringsAsFactors=FALSE)
if (!all(c('velocidad','modelo','parametros') %in% names(pars)))
  stop('El archivo de parametros internos tiene columnas diferentes a las esperadas.')

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
part <- particionar_datos(datos,cfg)
train <- as.data.frame(part$entrena)
predictores <- predictores_modelo(cfg)
vels <- sort(unique(train$v))
stopifnot(length(vels)==8,all(c('y','v',predictores) %in% names(train)))

# Se leen los hiperparámetros 1-SE seleccionados exclusivamente dentro de cada fold.
extraer_par <- function(vel,modelo) {
  p <- pars[pars$velocidad==vel & pars$modelo==modelo,,drop=FALSE]
  if (nrow(p)!=1L) stop('Falta parametro unico de ',modelo,' para velocidad ',vel)
  partes <- trimws(strsplit(p$parametros[1],';',fixed=TRUE)[[1]])
  valores <- list()
  for (parte in partes) {
    kv <- strsplit(parte,'=',fixed=TRUE)[[1]]
    if (length(kv)!=2L) stop('Formato inesperado de parametros: ',parte)
    valores[[trimws(kv[1])]] <- as.numeric(trimws(kv[2]))
  }
  if (!is.finite(valores$penalty)) stop('Penalty invalido para ',modelo,' en ',vel)
  if (modelo=='Elastic-net' && !is.finite(valores$mixture)) stop('Mixture invalido')
  valores
}

ajustar_regularizado <- function(clave,tr,te,vel) {
  nombre <- FAMILIAS[[clave]]$nombre
  p <- extraer_par(vel,nombre)
  # Receta escalada preparada SOLO con tr, no con el fold externo de prueba.
  rec <- receta_lineal_escalada(tr,predictores)
  modelo <- parsnip::linear_reg(penalty=p$penalty,
                               mixture=if (clave=='enet') p$mixture else 0) |>
    parsnip::set_engine('glmnet')
  wf <- workflows::workflow() |>
    workflows::add_recipe(rec) |>
    workflows::add_model(modelo)
  ajuste <- generics::fit(wf,data=tr)
  as.numeric(predict(ajuste,new_data=te)$.pred)
}

resultados <- list(); predicciones <- list()
for (vel in vels) {
  tr <- train[train$v!=vel,,drop=FALSE]
  te <- train[train$v==vel,,drop=FALSE]
  # Árbol de regresión con configuración FIJA (cfg$arbol_fijo). Si cfg$anidar_ensembles = TRUE,
  # la Tabla 1 usa en su lugar el árbol ajustado por LOSO anidado (paso 11).
  set.seed(cfg$semilla+as.integer(vel))
  f <- stats::reformulate(predictores,response='y')
  arbol <- rpart::rpart(f,data=tr,method='anova',
    control=rpart::rpart.control(cp=cfg$arbol_fijo$cp,maxdepth=cfg$arbol_fijo$profundidad,
                                 minsplit=cfg$arbol_fijo$min_division,xval=0))
  pred_arbol <- as.numeric(predict(arbol,newdata=te))
  # Promedio de modelos (averaging): combinación de predicciones de dos modelos con 1-SE interno.
  pred_enet <- ajustar_regularizado('enet',tr,te,vel)
  pred_ridge <- ajustar_regularizado('ridge',tr,te,vel)
  pred_avg <- (pred_enet+pred_ridge)/2
  pred_base <- rep(mean(tr$y),nrow(te))
  ps <- list('Arbol fijo'=pred_arbol,'Promedio ENet+Ridge'=pred_avg,
             'Elastic-net reajustado'=pred_enet,'Ridge reajustado'=pred_ridge,
             'Baseline global'=pred_base)
  for (nombre in names(ps)) {
    resultados[[length(resultados)+1L]] <- metricas_fila(vel,nombre,te$y,ps[[nombre]])
    predicciones[[length(predicciones)+1L]] <- data.frame(
      velocidad=vel,modelo=nombre,y=te$y,pred=ps[[nombre]])
  }
  log_msg(sprintf('  velocidad %d: árbol fijo RMSE = %.3f | promedio EN+Ridge RMSE = %.3f',
                  as.integer(vel),rmse_vec(te$y,pred_arbol),rmse_vec(te$y,pred_avg)))
}
res <- do.call(rbind,resultados)
preds <- do.call(rbind,predicciones)
salida <- DIR_TABLAS
dir.create(salida,recursive=TRUE,showWarnings=FALSE)
guardar_tabla(res,'19_detalle_arbol_promedio')
guardar_tabla(preds,'19_predicciones_arbol_promedio')

# RMSE/MAE/R2 medios POR VELOCIDAD: misma agregacion que el ranking del script 11.
ranking <- do.call(rbind,lapply(split(res,res$modelo),function(d) {
  data.frame(modelo=d$modelo[1],metodo=if (d$modelo[1]=='Arbol fijo')
      'LOSO, parametros fijos' else if (d$modelo[1]=='Baseline global')
      'LOSO, sin tuning' else 'LOSO, parametros 1-SE internos',
    rmse=mean(d$rmse),mae=mean(d$mae),r2=mean(d$r2,na.rm=TRUE),
    rmse_subgrupo_critico=mean(d$rmse_critico,na.rm=TRUE),
    n_velocidades= nrow(d))
}))
ref <- ranking$rmse[ranking$modelo=='Baseline global']
ranking$mejora_vs_baseline_pct <- 100*(1-ranking$rmse/ref)
ranking <- ranking[order(ranking$rmse),]
guardar_tabla(ranking,'19_ranking_arbol_promedio')
log_msg('Ranking LOSO: ',paste(sprintf('%s %.3f',ranking$modelo,ranking$rmse),collapse=' | '))
# Control: Elastic-net y Ridge reajustados aquí deben reproducir el paso 11.
if (file.exists(file.path(cfg$dir_anidada,'tablas','11_detalle_anidado.csv'))) {
  original <- read.csv(file.path(cfg$dir_anidada,'tablas','11_detalle_anidado.csv'))
  for (par in list(c('Elastic-net reajustado','Elastic-net'),c('Ridge reajustado','Ridge'))) {
    a <- res[res$modelo==par[1],c('velocidad','rmse')]
    b <- original[original$modelo==par[2],c('velocidad','rmse')]
    cmp <- merge(a,b,by='velocidad',suffixes=c('_19','_11'))
    dif <- max(abs(cmp$rmse_19-cmp$rmse_11))
    verificar(is.finite(dif) && dif <= 1e-4,
              paste0(par[2], ' reajustado en el paso 19 reproduce el paso 11 (dif. max. = ', signif(dif, 3), ')'))
  }
}
log_msg('Paso 19 terminado. Los test no se utilizaron.')
