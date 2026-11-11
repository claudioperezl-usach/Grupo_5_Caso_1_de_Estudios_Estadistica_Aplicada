# =============================================================================
# 06_modelos.R — Recetas, modelos, hiperparámetros y validación cruzada
# =============================================================================
# Principio anti-fuga: toda transformación aprendida de los datos (imputación,
# escalamiento, nudos de splines) es un paso de la receta. tune ajusta la receta
# SOLO con la parte de análisis de cada fold y la aplica a la de evaluación.

# --- Recetas ---------------------------------------------------------------------
receta_base <- function(datos, predictores) {
  recipes::recipe(stats::reformulate(predictores, response = "y"), data = datos) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>   # mediana del fold (no hay NA; queda como salvaguarda)
    recipes::step_zv(recipes::all_predictors())
}

receta_lineal_escalada <- function(datos, predictores) {
  receta_base(datos, predictores) |>
    recipes::step_normalize(recipes::all_numeric_predictors())          # media/de del fold
}

receta_splines <- function(datos, predictores) {
  vars_spl <- intersect(c("gtt", "ggn", "t48", "t2", "p2", "mf", "tic"), predictores)
  receta_base(datos, predictores) |>
    recipes::step_ns(dplyr::all_of(vars_spl), deg_free = tune::tune())  # nudos en cuantiles del fold
}

# --- Familias de modelos ------------------------------------------------------------
# 'orden' define qué configuración es "más simple" para la regla de 1 error estándar
# (equivale a lambda.1se de cv.glmnet: el modelo más penalizado a 1 EE del mínimo).
# Árbol, Random Forest y XGBoost solo se usan si cfg$anidar_ensembles = TRUE (paso 11).
FAMILIAS <- list(
  ols = list(nombre = "OLS", familia = "Lineal", receta = "base",
             modelo = function(cfg) parsnip::linear_reg() |> parsnip::set_engine("lm"),
             grilla = function(cfg) NULL, orden = NULL),

  ridge = list(nombre = "Ridge", familia = "Regularizado", receta = "escalada",
               modelo = function(cfg) parsnip::linear_reg(penalty = tune::tune(), mixture = 0) |>
                 parsnip::set_engine("glmnet"),
               grilla = function(cfg) tibble::tibble(penalty = cfg$grilla_lambda),
               orden = rlang::exprs(dplyr::desc(penalty))),

  lasso = list(nombre = "Lasso", familia = "Regularizado", receta = "escalada",
               modelo = function(cfg) parsnip::linear_reg(penalty = tune::tune(), mixture = 1) |>
                 parsnip::set_engine("glmnet"),
               grilla = function(cfg) tibble::tibble(penalty = cfg$grilla_lambda),
               orden = rlang::exprs(dplyr::desc(penalty))),

  enet = list(nombre = "Elastic-net", familia = "Regularizado", receta = "escalada",
              modelo = function(cfg) parsnip::linear_reg(penalty = tune::tune(), mixture = tune::tune()) |>
                parsnip::set_engine("glmnet"),
              grilla = function(cfg) tidyr::expand_grid(penalty = cfg$grilla_lambda,
                                                         mixture = c(0.25, 0.5, 0.75)),
              orden = rlang::exprs(dplyr::desc(penalty), dplyr::desc(mixture))),

  splines = list(nombre = "Splines naturales", familia = "Flexible", receta = "splines",
                 modelo = function(cfg) parsnip::linear_reg() |> parsnip::set_engine("lm"),
                 grilla = function(cfg) tibble::tibble(deg_free = c(2L, 3L, 4L, 6L, 8L)),
                 orden = rlang::exprs(deg_free)),

  arbol = list(nombre = "Árbol de regresión", familia = "Árbol", receta = "base",
               modelo = function(cfg) parsnip::decision_tree(cost_complexity = tune::tune(),
                                                            tree_depth = tune::tune(), min_n = 5L) |>
                 parsnip::set_engine("rpart") |> parsnip::set_mode("regression"),
               grilla = function(cfg) tidyr::expand_grid(cost_complexity = 10^c(-6, -5, -4, -3),
                                                          tree_depth = c(6L, 10L, 15L)),
               orden = rlang::exprs(dplyr::desc(cost_complexity), tree_depth)),

  rf = list(nombre = "Random Forest", familia = "Ensemble", receta = "base",
            modelo = function(cfg) parsnip::rand_forest(mtry = tune::tune(), min_n = tune::tune(),
                                                        trees = cfg$rf_arboles) |>
              parsnip::set_engine("ranger", num.threads = cfg$n_hilos, seed = cfg$semilla) |>
              parsnip::set_mode("regression"),
            grilla = function(cfg) tidyr::expand_grid(mtry = c(3L, 6L, 10L), min_n = c(2L, 5L, 20L)),
            orden = rlang::exprs(dplyr::desc(min_n), mtry)),

  xgb = list(nombre = "XGBoost", familia = "Ensemble", receta = "base",
             modelo = function(cfg) parsnip::boost_tree(trees = tune::tune(), learn_rate = tune::tune(),
                                                       tree_depth = tune::tune(), min_n = tune::tune(),
                                                       sample_size = 0.8) |>
               parsnip::set_engine("xgboost", nthread = cfg$n_hilos, tree_method = "hist") |>  # "hist": igual en xgboost 1.7 y 3.x
               parsnip::set_mode("regression"),
             grilla = function(cfg) tidyr::expand_grid(trees = c(300L, 800L), learn_rate = c(0.05, 0.15),
                                                        tree_depth = c(4L, 8L), min_n = c(2L, 10L)),
             orden = rlang::exprs(tree_depth, trees, dplyr::desc(min_n)))
)

construir_wf <- function(clave, datos, predictores, cfg) {
  f <- FAMILIAS[[clave]]
  rec <- switch(f$receta,
                base     = receta_base(datos, predictores),
                escalada = receta_lineal_escalada(datos, predictores),
                splines  = receta_splines(datos, predictores))
  workflows::workflow() |>
    workflows::add_recipe(rec) |>
    workflows::add_model(f$modelo(cfg))
}

# --- Diagnóstico de OLS: residuos, multicolinealidad e influencia ---------------------------
diagnostico_ols <- function(entrena, predictores, cfg) {
  ajuste <- generics::fit(construir_wf("ols", entrena, predictores, cfg), data = entrena)
  m <- workflows::extract_fit_engine(ajuste)

  # Multicolinealidad: VIF_j = 1 / (1 - R2_j), con los predictores ya procesados
  X <- as.matrix(recipes::bake(workflows::extract_recipe(ajuste), new_data = entrena,
                               recipes::all_predictors()))
  vif <- vapply(seq_len(ncol(X)), function(j) {
    fit_j <- stats::lm.fit(cbind(1, X[, -j, drop = FALSE]), X[, j])
    r2 <- 1 - sum(fit_j$residuals^2) / sum((X[, j] - mean(X[, j]))^2)
    1 / (1 - r2)
  }, numeric(1))
  tab_vif <- data.table::data.table(variable = colnames(X), VIF = round(vif, 1))[order(-VIF)]
  guardar_tabla(tab_vif, "05_ols_vif")
  log_msg("VIF máximo: ", tab_vif$variable[1], " = ", tab_vif$VIF[1])

  # Influencia: distancia de Cook
  cook <- stats::cooks.distance(m)
  n_inf <- sum(cook > 4 / length(cook), na.rm = TRUE)
  tab_inf <- data.table::data.table(id = entrena$id, v = entrena$v, y = entrena$y, cook = cook)[order(-cook)][1:10]
  guardar_tabla(tab_inf, "05_ols_observaciones_influyentes")
  log_msg("Observaciones con Cook > 4/n: ", n_inf, " (", round(100 * n_inf / length(cook), 2), "%)")

  # Residuos vs ajustados
  set.seed(cfg$semilla)
  i <- sample.int(length(cook), min(5000L, length(cook)))
  d <- data.frame(ajustado = stats::fitted(m)[i], residuo = stats::residuals(m)[i])
  p <- ggplot(d, aes(ajustado, residuo)) +
    geom_hline(yintercept = 0, colour = COL$baseline) +
    geom_point(alpha = 0.2, size = 0.8, colour = COL$principal) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = COL$secundario, linewidth = 0.8) +
    labs(title = "OLS: residuos vs valores ajustados",
         subtitle = "Una curva en los residuos indica que la relación no es lineal",
         x = "Pérdida ajustada (puntos %)", y = "Residuo (puntos %)",
         caption = "Muestra aleatoria de hasta 5.000 filas de entrenamiento.") +
    tema_proyecto()
  guardar_figura(p, "05_ols_residuos_vs_ajustados")

  invisible(list(vif = tab_vif, influyentes = tab_inf))
}
