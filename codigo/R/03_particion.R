# =============================================================================
# 03_particion.R — Entrenamiento y dos conjuntos de test
# =============================================================================
#   TODAS las filas
#   |-- velocidad reservada ............ test_vel    (velocidad no vista)
#   `-- resto
#       |-- 20% aleatorio .............. test_interp (interpolación; optimista en una grilla)
#       `-- 80% ........................ entrenamiento -> validación LOSO (pasos 11 y 12)
# En el pipeline actual ninguna decisión (modelo, hiperparámetros, variables, intervalos)
# usa los dos test. Aviso: ambos se examinaron antes en un análisis exploratorio previo,
# así que son evaluación confirmatoria y no totalmente independiente.

particionar_datos <- function(dt, cfg) {
  es_vel <- abs(dt$v - cfg$velocidad_reservada) < 1e-8
  verificar(any(es_vel), paste0("existe la velocidad reservada (", cfg$velocidad_reservada, " nudos)"))
  resto <- dt[!es_vel]

  set.seed(cfg$semilla)
  i_test <- sample.int(nrow(resto), round(cfg$prop_test * nrow(resto)))
  prueba_interp <- resto[i_test]
  entrena       <- resto[-i_test]
  prueba_vel    <- dt[es_vel]

  verificar(length(intersect(entrena$id, prueba_interp$id)) == 0, "sin filas compartidas: entrenamiento / test aleatorio")
  verificar(length(intersect(entrena$id, prueba_vel$id)) == 0, "sin filas compartidas: entrenamiento / test de velocidad")
  verificar(length(intersect(prueba_interp$id, prueba_vel$id)) == 0, "sin filas compartidas entre los dos test")
  verificar(!(cfg$velocidad_reservada %in% entrena$v), "la velocidad reservada no aparece en el entrenamiento")
  verificar(nrow(entrena) + nrow(prueba_interp) + nrow(prueba_vel) == nrow(dt), "las tres particiones suman todas las filas")

  guardar_tabla(data.table::data.table(
    conjunto = c("Entrenamiento (con CV)", "Test aleatorio (interpolación)", "Test velocidad no vista"),
    filas = c(nrow(entrena), nrow(prueba_interp), nrow(prueba_vel))), "02_particion")
  log_msg("Entrenamiento: ", nrow(entrena), " | test aleatorio: ", nrow(prueba_interp),
          " | test velocidad ", cfg$velocidad_reservada, ": ", nrow(prueba_vel))
  list(entrena = tibble::as_tibble(entrena), prueba_interp = tibble::as_tibble(prueba_interp),
       prueba_vel = tibble::as_tibble(prueba_vel))
}
