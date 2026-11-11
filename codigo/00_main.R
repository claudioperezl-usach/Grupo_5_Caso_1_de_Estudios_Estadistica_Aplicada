# =============================================================================
# 00_main.R — PUNTO DE ENTRADA ÚNICO del proyecto (regresión naval, USACH) · Equipo 5
# =============================================================================
# Estimación de la degradación del compresor de una turbina a gas naval (kMc) a partir
# de sensores en régimen estacionario. Modelo final: Elastic-net (paso 12), cuya elección
# frente a los demás candidatos se verifica con la validación cruzada anidada dejando una
# velocidad fuera (LOSO, paso 11).
#
# USO (sesión limpia de R; directorio de trabajo = carpeta raíz equipo_regresion_naval/):
#   RStudio : abrir equipo_regresion_naval.Rproj, reiniciar R (Ctrl+Shift+F10) y ejecutar
#             source("codigo/00_main.R")
#   Terminal: Rscript codigo/00_main.R
# Ejecuta los pasos 10 a 26 (lista "pasos" más abajo).
# Tiempo aproximado: 60-100 min según el computador con anidar_ensembles = TRUE (60,8 min en la
# ejecución oficial; unos 5 min con FALSE); la validación
# anidada de Random Forest y XGBoost es la etapa más lenta.
# Los parámetros se cambian SOLO en codigo/R/00_config.R.
# Reproducibilidad: el paso 11 compara el ranking anidado con codigo/referencia/ (entregada
# con el proyecto). Si no existe, o si cambió la configuración, la crea: ejecute DOS veces.
# =============================================================================

rm(list = ls())
if (!file.exists(file.path("codigo", "R", "00_config.R"))) {
  stop("Fije el directorio de trabajo en la carpeta raíz equipo_regresion_naval/ ",
       "(debe existir codigo/R/00_config.R). Directorio actual: ", getwd(), call. = FALSE)
}
inicio <- Sys.time()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # generador por defecto de R >= 3.6
source(file.path("codigo", "R", "00_config.R"), encoding = "UTF-8")
set.seed(cfg$semilla)

# Resultados nuevos en cada ejecución: nada de una corrida anterior se reutiliza.
# (codigo/referencia/ no se borra: es la referencia de reproducibilidad.)
unlink(c(cfg$dir_anidada, cfg$dir_final), recursive = TRUE)

pasos <- c(
  "10_auditoria_inferencia.R",          # auditoría, partición, verificación de fuga, relación con IC
  "11_validacion_anidada_velocidad.R",  # selección del modelo: LOSO anidada + verificación de reproducibilidad
  "12_entrenamiento_final_elastic_net.R", # sistema congelado, test, incertidumbre y subgrupo
  "13_graficos_finales.R",              # figuras del sistema final
  "14_diagnosticos_elastic_net.R",      # comparación por velocidad, coeficientes, error por degradación
  "15_sensibilidad_ruido.R",            # robustez ante ruido de sensores (1 % de la DE)
  "16_intervalos_predictivos.R",        # intervalos de predicción: cuantil empírico de residuos LOSO
  "17_comprobaciones_adicionales.R",    # convergencia y estabilidad de coeficientes
  "18_diagnostico_OLS_y_metodologia.R", # OLS: residuos, Cook y VIF (entrenamiento y por velocidad)
  "19_arbol_y_promedio_loso.R",         # árbol y promedio Elastic-net + Ridge
  "20_loso_ols_rf_xgb.R",               # métricas LOSO de OLS, GAM, RF y XGBoost
  "21_tabla1_final.R",                  # Tabla 1 con incertidumbre frente al baseline
  "22_interpretacion_pdp_ice.R",        # importancia por permutación (LOSO), PDP e ICE
  "23_cv_aleatoria_vs_loso.R",          # CV aleatoria vs LOSO solo con entrenamiento (§3 del informe)
  "24_multicolinealidad_shap.R",        # correlaciones y VIF de los predictores; valores SHAP del Elastic-net
  "25_curva_entrenamiento_vs_validacion.R", # error de entrenamiento vs validación LOSO según lambda
  "26_figuras_presentacion.R"           # figuras de la presentación (solo grafica tablas ya generadas)
)
tiempos <- numeric(length(pasos))
for (i in seq_along(pasos)) {
  cat("\n", strrep("=", 65), "\n", "ETAPA ", i, "/", length(pasos), ": ", pasos[i], "\n",
      strrep("=", 65), "\n", sep = "")
  t0 <- Sys.time()
  tryCatch(
    source(file.path("codigo", "pasos", pasos[i]), encoding = "UTF-8", local = new.env()),
    error = function(e) stop("Falló ", pasos[i], ": ", conditionMessage(e),
                             "\nRevise el error, corríjalo y vuelva a ejecutar 00_main.R.", call. = FALSE)
  )
  tiempos[i] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
}

# --- Registro de la ejecución: configuración, versiones y tiempos ------------------------
source(file.path("codigo", "R", "00_config.R"), encoding = "UTF-8")
dir.create(cfg$dir_final, recursive = TRUE, showWarnings = FALSE)
writeLines(utils::capture.output(utils::str(cfg)), file.path(cfg$dir_final, "config_usada.txt"))
writeLines(utils::capture.output(utils::sessionInfo()), file.path(cfg$dir_final, "sessionInfo.txt"))
paquetes <- c("data.table", "dplyr", "tibble", "tidyr", "ggplot2", "rlang", "rsample", "recipes",
              "parsnip", "workflows", "tune", "yardstick", "dials", "generics", "glmnet",
              "ranger", "xgboost", "rpart", "sandwich", "mgcv", "gridExtra")
utils::write.csv(data.frame(paquete = c("R", paquetes),
                            version = c(paste(R.version$major, R.version$minor, sep = "."),
                                        vapply(paquetes, function(p) as.character(utils::packageVersion(p)), ""))),
                 file.path(cfg$dir_final, "versiones_paquetes.csv"), row.names = FALSE)

# El mejor modelo del ranking anidado se LEE del resultado (no se escribe a mano).
ranking <- data.table::fread(file.path(cfg$dir_anidada, "tablas", "11_ranking_anidado.csv"), encoding = "UTF-8")
mejor_anidado <- ranking$modelo[which.min(ranking$rmse_promedio)]
writeLines(c("EJECUCIÓN COMPLETA", paste("Fecha:", format(Sys.time())),
             paste("Semilla:", cfg$semilla),
             sprintf("%-40s %6.2f min", pasos, tiempos),
             sprintf("Tiempo total: %.1f min", as.numeric(difftime(Sys.time(), inicio, units = "mins"))),
             "Modelo final entrenado y evaluado en test (paso 12): Elastic-net.",
             paste0("Menor RMSE en la validación anidada LOSO (paso 11): ", mejor_anidado,
                    sprintf(" (RMSE medio %.4f).", min(ranking$rmse_promedio))),
             if (!identical(mejor_anidado, "Elastic-net"))
               "ATENCIÓN: el mejor modelo del ranking anidado NO es el Elastic-net; justifique la elección (regla 1-SE, simplicidad, diferencia pareada del paso 21)." else
               "El modelo final coincide con el mejor del ranking anidado.",
             "Los dos test se examinaron antes en un análisis exploratorio previo: evaluación confirmatoria, no totalmente independiente.",
             "Ninguna cifra calculada con los test se usa para elegir partición, modelo, hiperparámetros, intervalos ni umbral."),
           file.path(cfg$dir_final, "REGISTRO_EJECUCION.txt"))
cat("\nPROYECTO COMPLETADO. Resultados en '", cfg$dir_anidada, "/' y '", cfg$dir_final, "/'.\n", sep = "")
cat("Tabla 1 del informe: ", file.path(cfg$dir_final, "tablas", "21_tabla1_presentacion.csv"), "\n", sep = "")
