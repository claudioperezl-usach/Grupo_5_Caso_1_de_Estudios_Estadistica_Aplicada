# =============================================================================
# 00_config.R — Parámetros del proyecto (ÚNICO lugar donde se cambian)
# Proyecto Integrador: Regresión · Caso 1 · Equipo 5 · Estadística Aplicada (USACH)
# Problema: estimar la degradación del compresor de una turbina a gas naval
# =============================================================================
# Respuesta: y = 100 * (1 - kMc)  = pérdida de rendimiento del compresor (puntos %).
#   kMc = 1 es compresor sano; kMc = 0.95 es la mayor degradación del dataset (y = 5).
#   Es una transformación lineal de kMc: solo sirve para leer los errores en
#   "puntos porcentuales" en vez de milésimas.
# Si cambia CUALQUIER valor de este archivo, la huella de configuración cambia y el paso 11
# crea una nueva referencia de reproducibilidad (hay que ejecutar dos veces para verificarla).

cfg <- list(

  # --- Reproducibilidad -------------------------------------------------------
  semilla       = 20261016L,
  semilla_ruido = 20261031L,         # paso 15 (sensibilidad a ruido de sensores, solo entrenamiento)
  # Un solo hilo: ranger/xgboost dan resultados distintos según el número de hilos.
  n_hilos       = 1L,

  # --- Rutas (relativas a la carpeta raíz equipo_regresion_naval/) -----------
  dir_datos      = "datos",
  archivo_datos  = "data.txt",
  filas_datos    = 11934L,                            # UCI CBM (Coraddu et al.): 9 x 51 x 26
  firma_datos    = 474761108.888,                     # suma de todas las celdas (ver 02_datos.R)
  dir_anidada    = "resultados_validacion_anidada",   # pasos 11 y 23
  dir_final      = "resultados_elastic_net_final",    # pasos 10 y 12-22

  # --- Variables excluidas (decididas por propiedades de los datos) ------------
  excluir_constantes    = c("t1", "p1"),   # no varían: no informan
  excluir_duplicadas    = c("tp", "lp"),   # tp = ts (idénticas); lp equivale a v (cor = 0,99992)
  excluir_no_disponible = "kmt",           # degradación de la turbina: no se conoce al decidir

  # --- Partición -----------------------------------------------------------------
  # Dos conjuntos de test, apartados del entrenamiento:
  #   (a) test aleatorio: 20% de las filas        -> mide interpolación (optimista: en una
  #       grilla, filas del test comparten velocidad y kMc con filas de entrenamiento)
  #   (b) velocidad no vista: una velocidad completa fuera del entrenamiento
  # NOTA: ambos test se examinaron antes en un análisis exploratorio previo (ya no incluido),
  # por lo que son una evaluación confirmatoria y no totalmente independiente. La elección
  # del modelo se hace SOLO con la validación anidada LOSO sobre el entrenamiento.
  velocidad_reservada = 15,          # nudos
  prop_test           = 0.20,

  # --- Modelos ----------------------------------------------------------------------
  rf_arboles         = 500L,         # árboles de Random Forest (paso 11 anidado y paso 20 fijo)
  grilla_lambda_loso = 10^seq(-4, 0, length.out = 10),  # grilla de lambda con LOSO (pasos 11, 12, 17 y 23)
  # Configuraciones FIJAS (solo si anidar_ensembles = FALSE; pasos 19 y 20). Con TRUE, la
  # Tabla 1 usa las versiones ajustadas por LOSO anidado del paso 11.
  arbol_fijo = list(cp = 0.001, profundidad = 10L, min_division = 20L),
  xgb_fijo   = list(rondas = 200L, eta = 0.05, profundidad = 4L),
  # TRUE : Árbol, Random Forest y XGBoost ajustan sus hiperparámetros por LOSO anidado
  #        (paso 11, regla 1-SE), igual que los modelos regularizados y los splines (60-100 min).
  # FALSE: se evalúan con LOSO y la configuración fija de arriba (más rápido).
  anidar_ensembles = TRUE,

  # --- Inferencia (relación estimada con intervalo) -----------------------------------
  inf_variable = "t2",               # temperatura de salida del compresor
  inf_delta    = 10,                 # se reporta el efecto de +10 K (= +10 °C)

  # --- Evaluación final ----------------------------------------------------------------
  B_bootstrap       = 500L,
  B_ruido           = 100L,          # repeticiones de la sensibilidad a ruido (paso 15)
  umbral_critico    = 4,             # subgrupo crítico: pérdida >= 4 puntos % (degradación severa)
  niveles_intervalo = c(0.80, 0.90, 0.95),   # intervalos de predicción (paso 16)
  # Política de inspección (paso 16): se alerta si la pérdida estimada >= c. El umbral c se
  # elige SOLO con las predicciones LOSO del entrenamiento: el mayor c (menos alarmas) que
  # detecta al menos esta fracción de los casos críticos (y >= umbral_critico).
  sensibilidad_objetivo = 0.90,
  grilla_umbral_alerta  = seq(1, 4, by = 0.1),

  # --- Interpretación (paso 22) ---------------------------------------------------------
  n_perm          = 10L,             # repeticiones de la permutación por variable y velocidad
  n_ice           = 150L,            # filas de entrenamiento para las curvas ICE
  n_pdp_variables = 4L,              # variables más importantes con PDP/ICE
  pdp_puntos      = 20L,             # puntos de la rejilla de cada PDP/ICE

  # --- Verificación de reproducibilidad -----------------------------------------------------
  # La primera ejecución guarda el ranking anidado como referencia (con la huella de esta
  # configuración); las siguientes con la MISMA configuración deben reproducirlo (paso 11).
  # Entregue codigo/referencia/ junto con el proyecto.
  ruta_referencia = file.path("codigo", "referencia", "ranking_anidado_referencia.csv"),
  tolerancia_reproducibilidad = 1e-3
)
