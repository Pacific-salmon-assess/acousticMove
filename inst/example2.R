# library(acousticMove)
library(ggplot2)

set.seed(125)
grid <- simulate_gmrf(
  x = seq(0, 1, length.out = 16),
  y = seq(0, 1, length.out = 16),
  nhabitat = 1,
  kappa = 0.95
)
names(grid)[names(grid) == "habitat_1"] <- "habitat"
grid$habitat <- as.numeric(scale(grid$habitat))

detectors <- expand.grid(
  x = seq(0.1, 0.9, length.out = 6),
  y = seq(0.1, 0.9, length.out = 6)
)
detectors$detector_id <- seq_len(nrow(detectors))

alpha <- 0.2 # diffusion param
beta <- 0.02 # advection coefs
q <- 0.1 # prob detecting ping while in a detector state
emission_rate <- 30 # known pings by tag per model time unit
study_end <- 10 # end of study in time units

# sim a tiny dataset
set.seed(1)
simulation <- acousticModel$new(grid = grid, detectors = detectors)
simulation$modelSetUp(~ habitat)
simulation$simulate(
  N = 5,
  alpha = alpha,
  beta = beta,
  q = q,
  emissionrate = emission_rate,
  studyperiod = study_end
)

ggplot(simulation$statespace, aes(x, y)) +
  geom_tile(aes(fill = habitat)) +
  geom_path(
    data = subset(simulation$sim_move, animal_id == 1),
    aes(x, y)
  ) +
  geom_point(
    data = simulation$detectors,
    aes(x, y),
    shape = 3,
    colour = "red"
  ) +
  scale_fill_viridis_c()

detections <- subset(
  simulation$sim_obs,
  !is.na(detector_id),
  select = c(animal_id, time, detector_id)
)

head(detectors)
head(detections)
head(grid)

fit <- fit_acoustic_model(
  detections = detections,
  covariate_grid = grid,
  detectors = detectors,
  formula = ~ habitat,
  emission_rate = emission_rate,
  study_end = study_end
)

fit$estimates
fit$optimizer
