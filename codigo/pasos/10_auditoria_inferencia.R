# =============================================================================
# 10_auditoria_inferencia.R — Auditoría, partición, verificación de fuga e inferencia
# =============================================================================
# Ejecutar desde la raíz: source("codigo/pasos/10_auditoria_inferencia.R")
# Genera (en resultados_elastic_net_final/tablas):
#   00_diccionario_variables.csv, 01_auditoria_variables.csv, 01_estructura_de_grilla.csv,
#   02_particion.csv, 03_correlaciones_con_respuesta.csv, 04_relacion_ic.csv
# y las figuras 01_t2_vs_perdida_por_velocidad.png y 04_relacion_t2_perdida.png.
# Usa SOLO el entrenamiento para la inferencia y las verificaciones; no abre los test.

if (!file.exists("codigo/R/00_config.R")) stop("Ejecute desde la carpeta raíz del proyecto.")
source("codigo/R/00_config.R", encoding = "UTF-8")
source("codigo/R/01_utilidades.R", encoding = "UTF-8")
cargar_proyecto(c("02_datos.R", "03_particion.R", "04_variables.R", "05_inferencia.R"))
iniciar_salida(cfg$dir_final)
log_msg("=== Paso 10: auditoría, partición, verificación de fuga e inferencia ===")

set.seed(cfg$semilla)
datos <- cargar_datos(cfg)
auditar_datos(datos, cfg)
part <- particionar_datos(datos, cfg)
entrena <- part$entrena
predictores <- predictores_modelo(cfg)

guardar_tabla(DICCIONARIO, "00_diccionario_variables")
verificar_variables_sin_fuga(entrena, predictores, cfg)

# Relación estimada con intervalo: T2 y pérdida de rendimiento, a igual velocidad.
estimar_relacion(entrena, cfg)

rm(part, datos)
log_msg("Paso 10 terminado. Los test no se utilizaron.")
