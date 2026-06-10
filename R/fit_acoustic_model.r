#' Fit an acoustic movement model
#'
#' @param formula Formula describing movement covariates. Covariates in `formula`
#'   must match columns in `covariate_grid`.
#' @param covariate_grid Gridded covariate surface with `x`, `y`, and the
#'   movement covariates referenced by `formula`. Defines the spatial state
#'   space over which movement is modelled.
#' @param detectors Detector locations with `detector_id`, `x`, and `y` columns.
#' @param detections Detection data with `animal_id`, `time`, and `detector_id`
#'   columns.
#' @param emission_rate Known tag emission rate per unit time.
#' @param study_end End of the study in the same units as detection times.
#' @param releases Optional release locations and times.
#' @param known_fates Optional known final locations or fates.
#' @param absorbing_states Optional locations representing absorbing states.
#' @param attraction Whether to estimate a centre of attraction.
#' @param mortality Whether to estimate a mortality rate.
#' @param start Optional named list of starting values for `alpha`, `beta`, `q`,
#'   `gamma`, and `mu`.
#' @param optimizer_control Optional control list passed to [stats::nlminb()].
#' @param model_control Optional control list passed to
#'   `acousticModel$makeADFun()`.
#'
#' @return A list containing parameter estimates, the optimizer result, and the
#'   fitted `acousticModel` object.
#' @importFrom assertthat assert_that
#' @examples
#' library(acousticMove)
#'
#' # Simulate a small example dataset:
#' set.seed(125)
#' grid <- simulate_gmrf(
#'   x = seq(0, 1, length.out = 16),
#'   y = seq(0, 1, length.out = 16),
#'   nhabitat = 1,
#'   kappa = 0.95
#' )
#' names(grid)[names(grid) == "habitat_1"] <- "habitat"
#' grid$habitat <- as.numeric(scale(grid$habitat))
#' 
#' detectors <- expand.grid(
#'   x = seq(0.1, 0.9, length.out = 6),
#'   y = seq(0.1, 0.9, length.out = 6)
#' )
#' detectors$detector_id <- seq_len(nrow(detectors))
#' 
#' # Simulation params and values:
#' alpha <- 0.2 # diffusion parameter
#' beta <- 0.02 # advection coefficients
#' q <- 0.1 # prob detecting ping while in a detector state
#' emission_rate <- 30 # known pings by tag per model time unit
#' study_end <- 10 # end of study in time units
#' 
#' set.seed(1)
#' simulation <- acousticModel$new(grid = grid, detectors = detectors)
#' simulation$modelSetUp(~ habitat)
#' simulation$simulate(
#'   N = 5,
#'   alpha = alpha,
#'   beta = beta,
#'   q = q,
#'   emissionrate = emission_rate,
#'   studyperiod = study_end
#' )
#' # Input data frames for fitting:
#' detections <- subset(
#'   simulation$sim_obs,
#'   !is.na(detector_id),
#'   select = c(animal_id, time, detector_id)
#' )
#' 
#' # Visualize what we simulated:
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(simulation$statespace, aes(x, y)) +
#'     geom_tile(aes(fill = habitat)) +
#'     geom_path(
#'       data = subset(simulation$sim_move, animal_id == 1),
#'       aes(x, y)
#'     ) +
#'     geom_point(
#'       data = simulation$detectors,
#'       aes(x, y),
#'       shape = 3,
#'       colour = "red"
#'     ) +
#'     scale_fill_viridis_c()
#' }
#' 
#' # Our model inputs:
#' head(grid)
#' head(detections)
#' head(detectors)
#' emission_rate
#' max(detections$time)
#' study_end
#' 
#' # Fit a model to the simulated data:
#' fit <- fit_acoustic_model(
#'   formula = ~ habitat,
#'   covariate_grid = grid,
#'   detectors = detectors,
#'   detections = detections,
#'   emission_rate = emission_rate,
#'   study_end = study_end
#' )
#' 
#' fit$estimates
#' fit$optimizer
#'
#' @export
fit_acoustic_model <- function(formula, covariate_grid, detectors, detections,
                               emission_rate, study_end, releases = NULL,
                               known_fates = NULL, absorbing_states = NULL,
                               attraction = FALSE, mortality = FALSE,
                               start = NULL, optimizer_control = NULL,
                               model_control = NULL) {
  assert_that(inherits(formula, "formula"))
  assert_that(
    is.data.frame(detections) &&
      all(c("animal_id", "time", "detector_id") %in% names(detections)) &&
      nrow(detections) > 0 && is.numeric(detections$time),
    msg = "`detections` must contain `animal_id`, numeric `time`, and `detector_id`."
  )
  assert_that(
    is.data.frame(covariate_grid) &&
      all(c("x", "y", all.vars(formula)) %in% names(covariate_grid)),
    msg = "`covariate_grid` must contain `x`, `y`, and all variables in `formula`."
  )
  assert_that(
    is.data.frame(detectors) &&
      all(c("detector_id", "x", "y") %in% names(detectors)) &&
      !anyDuplicated(detectors$detector_id),
    msg = "`detectors` must contain unique `detector_id` values and `x`, `y`."
  )
  assert_that(
    assertthat::is.number(emission_rate) && emission_rate > 0,
    assertthat::is.number(study_end) && study_end > 0,
    all(is.finite(detections$time)) && max(detections$time) <= study_end,
    msg = "Rates and times must be finite, positive, and no later than `study_end`."
  )

  detections <- detections[, c("animal_id", "time", "detector_id")]
  detector_match <- match(detections$detector_id, detectors$detector_id)
  assert_that(
    !anyNA(detector_match),
    msg = "Every detection `detector_id` must occur in `detectors`."
  )
  detections$x <- detectors$x[detector_match]
  detections$y <- detectors$y[detector_match]

  model <- acousticModel$new(grid = covariate_grid, detectors = detectors)
  formula <- stats::update.formula(formula, ~ . - 1) # always drop intercept
  model$modelSetUp(formula)
  model$processData(
    detections = detections,
    releases = releases,
    known_fates = known_fates,
    absorbing_states = absorbing_states
  )

  alpha <- if (is.null(start$alpha)) model$resolution[1] else start$alpha
  beta <- if (is.null(start$beta)) {
    rep(0, ncol(model$designmatrix) + as.integer(attraction))
  } else {
    start$beta
  }
  q <- if (is.null(start$q)) 0.1 else start$q
  gamma <- if (!attraction) {
    NULL
  } else if (is.null(start$gamma)) {
    c(mean(model$statespace$x), mean(model$statespace$y))
  } else {
    start$gamma
  }
  mu <- if (!mortality) NULL else if (is.null(start$mu)) 0.001 else start$mu

  model$makeADFun(
    alpha = alpha,
    beta = beta,
    q = q,
    mu = mu,
    gamma = gamma,
    emissionrate = emission_rate,
    studyperiod = study_end,
    control = if (is.null(model_control)) list() else model_control
  )

  optimizer <- suppressWarnings(stats::nlminb(
    start = model$negll$par,
    objective = model$negll$fn,
    gradient = model$negll$gr,
    control = if (is.null(optimizer_control)) list() else optimizer_control
  ))
  estimates <- reList(optimizer$par)
  model$estimated_pars <- estimates

  list(estimates = estimates, optimizer = optimizer, model = model)
}
