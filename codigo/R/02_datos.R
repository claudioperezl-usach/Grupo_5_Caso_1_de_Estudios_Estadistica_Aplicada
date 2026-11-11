# =============================================================================
# 02_datos.R — Carga, definición de la respuesta y auditoría
# =============================================================================
# Cada fila = un estado de operación en régimen estacionario de la planta
# (velocidad x degradación del compresor x degradación de la turbina).

NOMBRES_COLUMNAS <- c("lp", "v", "gtt", "gtn", "ggn", "ts", "tp", "t48", "t1", "t2",
                      "p48", "p1", "p2", "pexh", "tic", "mf", "kmc", "kmt")

cargar_datos <- function(cfg) {
  ruta <- file.path(cfg$dir_datos, cfg$archivo_datos)
  if (!file.exists(ruta)) {
    stop("No se encuentra ", ruta, ". Descargue el dataset (ver README.md: UCI, DOI 10.24432/C5K31K) ",
         "y copie 'UCI CBM Dataset/data.txt' en la carpeta datos/.", call. = FALSE)
  }
  log_msg("Leyendo ", ruta)
  d <- utils::read.table(ruta, header = FALSE, col.names = NOMBRES_COLUMNAS)
  d <- data.table::as.data.table(d)

  verificar(ncol(d) == 18L && nrow(d) >= 1000L, paste0("dataset con 18 columnas y ", nrow(d), " filas (>= 1000)"))
  verificar(!anyNA(d), "sin valores faltantes")
  verificar(!anyDuplicated(d), "sin filas duplicadas")
  # Huella del CONTENIDO (no del archivo): no depende de saltos de línea ni del formato
  # numérico, así que detecta si se usa otra versión del dataset.
  # Tolerancia relativa: la suma puede variar en los últimos dígitos según la plataforma.
  firma <- sum(as.matrix(d))
  verificar(nrow(d) == cfg$filas_datos && abs(firma / cfg$firma_datos - 1) < 1e-9,
            paste0("los datos coinciden con la versión usada en el informe (", nrow(d), " filas; firma ", format(firma, nsmall = 3), ")"))

  d[, id := .I]
  d[, y := 100 * (1 - kmc)]
  verificar(all(d$y >= -1e-9 & d$y <= 5 + 1e-9), "respuesta y dentro de [0, 5] puntos %")
  d[]
}

auditar_datos <- function(dt, cfg) {
  vars <- setdiff(NOMBRES_COLUMNAS, character(0))
  aud <- data.table::rbindlist(lapply(vars, function(v) {
    x <- dt[[v]]
    data.table::data.table(variable = v, n_valores_distintos = data.table::uniqueN(x),
                           media = mean(x), de = stats::sd(x), minimo = min(x), maximo = max(x),
                           constante = stats::sd(x) == 0)
  }))
  guardar_tabla(aud, "01_auditoria_variables")
  log_msg("Variables constantes: ", paste(aud[constante == TRUE, variable], collapse = ", "))

  grilla <- dt[, .(velocidades = data.table::uniqueN(v), niveles_kmc = data.table::uniqueN(kmc),
                   niveles_kmt = data.table::uniqueN(kmt), filas = .N)]
  grilla[, producto := velocidades * niveles_kmc * niveles_kmt]
  guardar_tabla(grilla, "01_estructura_de_grilla")
  verificar(grilla$producto == grilla$filas,
            "el dataset es una grilla completa: velocidades x niveles kMc x niveles kMt = filas")

  p <- ggplot(dt, aes(t2, y, colour = factor(v))) +
    geom_point(alpha = 0.3, size = 0.5) +
    scale_colour_viridis_d(name = "Velocidad (nudos)") +
    labs(title = "La velocidad domina a los sensores; la degradación es una señal pequeña",
         subtitle = "Temperatura de salida del compresor (T2) vs pérdida de rendimiento, por velocidad",
         x = "T2 (K)", y = "Pérdida de rendimiento del compresor (puntos %)",
         caption = "Datos simulados (UCI, Coraddu et al. 2014).") +
    tema_proyecto()
  guardar_figura(p, "01_t2_vs_perdida_por_velocidad")
  invisible(aud)
}
