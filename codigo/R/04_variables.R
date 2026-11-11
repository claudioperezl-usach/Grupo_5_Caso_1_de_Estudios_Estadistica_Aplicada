# =============================================================================
# 04_variables.R — Diccionario, predictores y verificaciones anti-fuga
# =============================================================================
# Instante de predicción: cuando los sensores registran el régimen estacionario.
# Se estima el estado ACTUAL del compresor (no es un pronóstico futuro).

DICCIONARIO <- data.table::data.table(
  variable = c("y", "v", "gtt", "gtn", "ggn", "ts", "t48", "t2", "p48", "p2", "pexh", "tic", "mf",
               "lp", "tp", "t1", "p1", "kmt"),
  descripcion = c("Pérdida de rendimiento del compresor = 100*(1 - kMc)", "Velocidad del buque",
                  "Torque del eje de la turbina a gas", "Revoluciones de la turbina a gas",
                  "Revoluciones del generador de gas", "Torque de hélice (estribor)",
                  "Temperatura de salida de la turbina de alta presión", "Temperatura de salida del compresor",
                  "Presión de salida de la turbina de alta presión", "Presión de salida del compresor",
                  "Presión de los gases de escape", "Control de inyección de la turbina",
                  "Flujo de combustible", "Posición de palanca", "Torque de hélice (babor)",
                  "Temperatura de entrada del compresor", "Presión de entrada del compresor",
                  "Coeficiente de degradación de la turbina"),
  # Temperaturas en K: UCI las rotula en °C, pero T1 = 288 constante solo es físico en Kelvin
  # (15 °C, condición ISO). Un cambio de +10 K equivale a +10 °C.
  unidad = c("puntos %", "nudos", "kN m", "rpm", "rpm", "kN", "K", "K", "bar", "bar", "bar", "%", "kg/s",
             "-", "kN", "K", "bar", "-"),
  rol = c("respuesta", rep("predictor (sensor en régimen estacionario)", 12),
          "excluida: equivale a v (correlación 0,9999)", "excluida: idéntica a ts",
          "excluida: constante", "excluida: constante",
          "excluida: no se conoce al decidir (fuga de información)")
)

predictores_modelo <- function(cfg) {
  todos <- setdiff(NOMBRES_COLUMNAS, "kmc")           # kMc es la respuesta
  setdiff(todos, c(cfg$excluir_constantes, cfg$excluir_duplicadas, cfg$excluir_no_disponible))
}

# Verificaciones calculadas SOLO con el entrenamiento.
verificar_variables_sin_fuga <- function(entrena, predictores, cfg) {
  d <- data.table::as.data.table(entrena)

  verificar(!("kmt" %in% predictores), "kMt (degradación de la turbina) no es predictor: no se conoce al decidir")
  verificar(!any(c("kmc", "y") %in% predictores), "ni la respuesta ni kMc están entre los predictores")
  verificar(all(abs(d$y - 100 * (1 - d$kmc)) < 1e-9), "y = 100 * (1 - kMc)")
  verificar(length(predictores) >= 8L, paste0(length(predictores), " predictores informativos (>= 8)"))

  const_ok <- vapply(cfg$excluir_constantes, function(v) data.table::uniqueN(d[[v]]) == 1L, logical(1))
  verificar(all(const_ok), paste0("variables excluidas por constantes: ", paste(cfg$excluir_constantes, collapse = ", ")))

  M <- as.matrix(d[, ..predictores])
  dup_ok <- vapply(cfg$excluir_duplicadas, function(v) max(abs(stats::cor(d[[v]], M))) >= 0.9999, logical(1))
  verificar(all(dup_ok), paste0("variables excluidas por duplicadas: ", paste(cfg$excluir_duplicadas, collapse = ", ")))

  # Ningún predictor reproduce la respuesta (global ni dentro de cada velocidad)
  centrar <- function(x, g) x - stats::ave(x, g)
  cor_global <- vapply(predictores, function(v) stats::cor(d[[v]], d$y), numeric(1))
  sin_v <- setdiff(predictores, "v")
  cor_intra <- vapply(sin_v, function(v) stats::cor(centrar(d[[v]], d$v), centrar(d$y, d$v)), numeric(1))
  verificar(all(abs(cor_global) < 0.98), paste0("ninguna |correlación global| con y >= 0.98 (máx: ",
                                              names(which.max(abs(cor_global))), " = ", round(max(abs(cor_global)), 3), ")"))
  verificar(all(abs(cor_intra) < 0.98), paste0("ninguna |correlación intra-velocidad| con y >= 0.98 (máx: ",
                                             names(which.max(abs(cor_intra))), " = ", round(max(abs(cor_intra)), 3), ")"))

  guardar_tabla(data.table::data.table(
    variable = predictores, cor_global = round(cor_global, 3),
    cor_intra_velocidad = round(unname(cor_intra[predictores]), 3)), "03_correlaciones_con_respuesta")
  invisible(TRUE)
}
