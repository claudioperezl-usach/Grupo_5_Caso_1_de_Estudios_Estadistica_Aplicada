# =============================================================================
# 05_inferencia.R — Relación estimada con intervalo: T2 y pérdida de rendimiento
# =============================================================================
# Pregunta: a igual velocidad, ¿cuántos puntos % de pérdida de rendimiento se asocian
# a +10 K (equivalente a +10 °C) en la temperatura de salida del compresor (T2)?
#   y = b * T2 + efecto fijo de velocidad + error
# Decisiones de rigor:
#   * Solo entrenamiento. Se controla la velocidad (efecto fijo), que domina a los sensores.
#   * IC clásico de lm()/confint() (visto en clases) y, como sensibilidad, IC con EE
#     robustos a heterocedasticidad (HC3, paquete sandwich): los residuos cambian de
#     dispersión entre velocidades, y HC3 no supone varianza constante. Si ambos IC
#     coinciden en lo sustantivo, la conclusión no depende de ese supuesto.
#     (HC3 corrige heterocedasticidad, NO la dependencia entre filas de la grilla.)
#   * Los datos son una grilla simulada sin ruido: los p-valores son muy pequeños casi
#     por construcción. Por eso se juzga la RELEVANCIA con el tamaño del efecto.
#   * Es una asociación, no causalidad (kMt, no observada en la práctica, también afecta T2;
#     en la grilla kMc y kMt son ortogonales, por lo que no sesga b).

estimar_relacion <- function(entrena, cfg) {
  d <- data.table::as.data.table(entrena)
  xv <- cfg$inf_variable
  delta <- cfg$inf_delta
  d[, v_f := factor(v)]

  m  <- stats::lm(stats::as.formula(paste0("y ~ ", xv, " + v_f")), data = d)
  m0 <- stats::lm(y ~ v_f, data = d)

  b      <- stats::coef(m)[[xv]]
  ee_rob <- sqrt(sandwich::vcovHC(m, type = "HC3")[xv, xv])
  ee_cla <- sqrt(stats::vcov(m)[xv, xv])
  gl     <- stats::df.residual(m)
  z      <- stats::qt(0.975, df = gl)
  ic     <- b + c(-1, 1) * z * ee_rob
  ic_cla <- stats::confint(m, xv, level = 0.95)[1, ]          # IC clásico (MCO)

  xc_vec   <- stats::resid(stats::lm(stats::as.formula(paste0(xv, " ~ v_f")), data = d))
  sd_intra <- stats::sd(xc_vec)
  r2_parcial <- 1 - sum(stats::resid(m)^2) / sum(stats::resid(m0)^2)

  res <- data.table::data.table(
    cantidad = c(paste0("Coeficiente b (puntos % por 1 unidad de ", xv, ")"),
                 "EE robusto (HC3)", "EE clásico (MCO)", "IC 95% de b (HC3)", "IC 95% de b (clásico, confint)",
                 paste0("Efecto de +", delta, " unidades de ", xv, " [IC 95% HC3] (puntos %)"),
                 paste0("Efecto de +", delta, " unidades de ", xv, " [IC 95% clásico] (puntos %)"),
                 paste0("Efecto de +1 DE de ", xv, " dentro de cada velocidad (puntos %)"),
                 "DE de la respuesta (puntos %), para dimensionar el efecto",
                 "R2 parcial de la variable (sobre el modelo solo con velocidad)",
                 "Valor p (robusto)", "Observaciones"),
    valor = c(sprintf("%.5f", b), sprintf("%.5f", ee_rob), sprintf("%.5f", ee_cla),
              sprintf("[%.5f; %.5f]", ic[1], ic[2]),
              sprintf("[%.5f; %.5f]", ic_cla[1], ic_cla[2]),
              fmt_ic(delta * b, delta * ic[1], delta * ic[2], 3),
              fmt_ic(delta * b, delta * ic_cla[1], delta * ic_cla[2], 3),
              sprintf("%.3f", b * sd_intra), sprintf("%.3f", stats::sd(d$y)),
              sprintf("%.3f", r2_parcial),
              format.pval(2 * stats::pt(-abs(b / ee_rob), df = gl), digits = 3),
              as.character(nrow(d))))
  guardar_tabla(res, "04_relacion_ic")
  log_msg("Relación ", xv, ": +", delta, " unidades -> ", fmt_ic(delta * b, delta * ic[1], delta * ic[2], 3),
          " puntos % de pérdida (IC 95%, EE robusto); R2 parcial = ", round(r2_parcial, 3))

  d[, x_c := xc_vec]
  d[, y_c := y - mean(y), by = v]
  p <- ggplot(d, aes(x_c, y_c)) +
    geom_point(alpha = 0.15, size = 0.6, colour = COL$principal) +
    geom_smooth(method = "lm", formula = y ~ x, colour = COL$secundario, linewidth = 1, se = FALSE) +
    labs(title = paste0("Pérdida de rendimiento vs ", toupper(xv), " (descontada la velocidad)"),
         subtitle = "Ambas variables centradas por velocidad; la recta es el ajuste lineal",
         x = paste0(toupper(xv), " - promedio de su velocidad"), y = "Pérdida - promedio de su velocidad (puntos %)",
         caption = "Asociación, no causalidad. Datos simulados.") +
    tema_proyecto()
  guardar_figura(p, "04_relacion_t2_perdida")
  invisible(list(modelo = m, tabla = res))
}
