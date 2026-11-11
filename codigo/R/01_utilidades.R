# =============================================================================
# 01_utilidades.R — Registro, verificaciones, métricas y estilo de figuras
# =============================================================================

# --- Paquetes ----------------------------------------------------------------
PAQUETES <- c(
  "data.table", "dplyr", "tibble", "tidyr", "ggplot2", "rlang",
  "rsample", "recipes", "parsnip", "workflows", "tune", "yardstick", "dials",
  "generics", "glmnet", "ranger", "xgboost", "rpart", "sandwich", "mgcv", "gridExtra"
)

verificar_paquetes <- function(paquetes = PAQUETES) {
  faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
  if (length(faltan) > 0) {
    stop("Faltan paquetes: ", paste(faltan, collapse = ", "),
         "\nInstálelos con:\ninstall.packages(c(",
         paste0('"', faltan, '"', collapse = ", "), "))", call. = FALSE)
  }
  invisible(TRUE)
}

cargar_paquetes <- function() {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)       # se carga después de data.table: usamos sintaxis data.table explícita
    library(tibble)
    library(ggplot2)
    library(rsample)
    library(recipes)
    library(parsnip)
    library(workflows)
    library(tune)
    library(yardstick)
    library(dials)
  })
}

# --- Carpeta de salida de cada paso --------------------------------------------------
# Define DIR_RESULTADOS, DIR_TABLAS, DIR_FIGURAS y RUTA_LOG (variables globales que usan
# log_msg(), guardar_tabla() y guardar_figura()). Cada paso la llama al inicio, de modo
# que también puede ejecutarse por separado desde una sesión limpia.
iniciar_salida <- function(dir) {
  assign("DIR_RESULTADOS", dir, envir = .GlobalEnv)
  assign("DIR_TABLAS", file.path(dir, "tablas"), envir = .GlobalEnv)
  assign("DIR_FIGURAS", file.path(dir, "figuras"), envir = .GlobalEnv)
  dir.create(file.path(dir, "tablas"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "figuras"), recursive = TRUE, showWarnings = FALSE)
  assign("RUTA_LOG", file.path(dir, "log.txt"), envir = .GlobalEnv)
  invisible(dir)
}

# Carga configuración, utilidades, módulos y paquetes (inicio común de cada paso).
cargar_proyecto <- function(modulos = c("02_datos.R", "03_particion.R", "04_variables.R", "06_modelos.R")) {
  if (!file.exists(file.path("codigo", "R", "00_config.R"))) {
    stop("Ejecute desde la carpeta raíz 'equipo_regresion_naval/' (debe existir codigo/R/00_config.R).",
         call. = FALSE)
  }
  for (m in modulos) source(file.path("codigo", "R", m), encoding = "UTF-8")
  verificar_paquetes()
  cargar_paquetes()
  invisible(TRUE)
}

# --- Registro de la ejecución -----------------------------------------------------
# RUTA_LOG se define con iniciar_salida()
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  message(msg)
  cat(msg, "\n", file = RUTA_LOG, append = TRUE, sep = "")
}

# Detiene la ejecución si una condición de integridad no se cumple.
# Cada verificación superada queda registrada en el log (evidencia para la rúbrica).
verificar <- function(condicion, descripcion) {
  if (!isTRUE(condicion)) {
    stop("VERIFICACIÓN FALLIDA: ", descripcion, call. = FALSE)
  }
  log_msg("  [OK] ", descripcion)
  invisible(TRUE)
}

# --- Métricas -------------------------------------------------------------------
rmse_vec <- function(y, pred) sqrt(mean((y - as.numeric(pred))^2))

# Métricas de una velocidad evaluada: RMSE, MAE, R2 (1 - SSE/SST, puede ser negativo) y RMSE
# del subgrupo crítico (y >= umbral). Usada por los pasos 11, 19 y 20.
metricas_fila <- function(vel, modelo, y, p, umbral = cfg$umbral_critico) {
  p <- as.numeric(p)
  crit <- y >= umbral
  data.frame(velocidad = vel, modelo = modelo, n = length(y),
             rmse = rmse_vec(y, p), mae = mean(abs(y - p)),
             r2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2),
             n_critico = sum(crit),
             rmse_critico = if (any(crit)) rmse_vec(y[crit], p[crit]) else NA_real_)
}

# --- Particiones dejando UNA velocidad completa fuera (LOSO) ---------------------------
crear_loso <- function(d) {
  vv <- sort(unique(d$v))
  idx <- lapply(vv, function(v) list(analysis = which(d$v != v), assessment = which(d$v == v)))
  splits <- lapply(idx, function(i) rsample::make_splits(i, data = d))
  rsample::manual_rset(splits, ids = paste0("Vel_", vv))
}

# --- Guardado de tablas y figuras -------------------------------------------------
guardar_tabla <- function(x, nombre) {
  ruta <- file.path(DIR_TABLAS, paste0(nombre, ".csv"))
  data.table::fwrite(as.data.frame(x), ruta)
  log_msg("    tabla: ", ruta)
  invisible(ruta)
}

guardar_figura <- function(p, nombre, ancho = 8, alto = 4.5) {
  ruta <- file.path(DIR_FIGURAS, paste0(nombre, ".png"))
  ggplot2::ggsave(ruta, p, width = ancho, height = alto, dpi = 200, bg = "white")
  log_msg("    figura: ", ruta)
  invisible(ruta)
}

# Paleta: un tono por rol, con orden fijo (el baseline siempre en gris).
COL <- list(
  principal = "#2a78d6",   # sistema / modelo de interés
  secundario = "#eb6834",  # segunda serie
  terciario = "#1baf7a",   # tercera serie
  baseline = "#8a8a85",    # referencia
  tinta = "#2b2b29",
  tinta_suave = "#6b6b66"
)

tema_proyecto <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", colour = COL$tinta),
      plot.subtitle    = ggplot2::element_text(colour = COL$tinta_suave),
      plot.caption     = ggplot2::element_text(colour = COL$tinta_suave, hjust = 0),
      axis.title       = ggplot2::element_text(colour = COL$tinta_suave),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e5e0", linewidth = 0.3),
      legend.position  = "top"
    )
}

# --- Verificación automática de reproducibilidad --------------------------------------
# La primera ejecución con una configuración guarda el ranking como referencia (junto con
# una huella de la configuración). Las ejecuciones siguientes con la MISMA configuración
# deben reproducirlo dentro de la tolerancia.
huella_cfg <- function(cfg) {
  tmp <- tempfile(fileext = ".txt")
  # Se excluyen n_hilos (no define el análisis) y grilla_lambda, que los pasos 11, 12 y 23
  # copian de grilla_lambda_loso en tiempo de ejecución: la huella depende solo de 00_config.R.
  cfg_sin_hilos <- cfg[setdiff(names(cfg), c("n_hilos", "grilla_lambda"))]
  # Texto construido con los VALORES (números con 12 cifras significativas), no con str():
  # así la huella no cambia entre versiones de R ni sistemas operativos.
  valores <- rapply(cfg_sin_hilos, function(z) if (is.numeric(z)) sprintf("%.12g", z) else as.character(z),
                    how = "unlist")
  writeLines(paste(names(valores), valores, sep = "=", collapse = "|"), tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}

# Compara el ranking anidado (paso 11) con la referencia guardada, SOLO si la configuración
# es la misma (misma huella). Si la configuración cambió, o la referencia no tiene huella,
# se crea una referencia nueva y se avisa en el log; ejecute otra vez para verificarla.
verificar_referencia_ranking <- function(ranking, cfg) {
  actual <- data.table::copy(ranking[, .(modelo, rmse_promedio, mae_promedio)])
  actual[, huella := huella_cfg(cfg)]
  ruta <- cfg$ruta_referencia
  if (file.exists(ruta)) {
    ref <- data.table::fread(ruta, encoding = "UTF-8")
    if (!("huella" %in% names(ref))) {
      log_msg("La referencia existente no tiene huella de configuración: se reemplaza por una nueva.")
    } else if (!identical(as.character(ref$huella[1]), actual$huella[1])) {
      log_msg("La configuración cambió respecto de la referencia: se crea una nueva referencia (no se compara).")
    } else {
      comp <- merge(ref, actual, by = "modelo", suffixes = c("_ref", ""))
      verificar(nrow(comp) == nrow(actual) && nrow(ref) == nrow(actual),
                "la referencia y la ejecución actual comparan los mismos modelos")
      dif <- max(abs(comp$rmse_promedio - comp$rmse_promedio_ref),
                 abs(comp$mae_promedio - comp$mae_promedio_ref))
      verificar(dif <= cfg$tolerancia_reproducibilidad,
                paste0("el ranking anidado reproduce la ejecución de referencia (dif. máx. = ", signif(dif, 3), ")"))
      return(invisible(TRUE))
    }
  }
  dir.create(dirname(ruta), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(actual, ruta)
  log_msg("Referencia de reproducibilidad guardada en ", ruta,
          ". Vuelva a ejecutar desde una sesión limpia para verificarla.")
  invisible(FALSE)
}

# Formato "estimación [IC]" para las tablas (sin exceso de decimales).
fmt_ic <- function(m, lo, hi, d = 2) sprintf(paste0("%.", d, "f [%.", d, "f; %.", d, "f]"), m, lo, hi)
