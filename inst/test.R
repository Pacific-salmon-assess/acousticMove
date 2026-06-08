library(RTMB)
library(R6)
library(Matrix)
library(ggplot2)
# library(acousticMove)

# remotes::install_github("Pacific-salmon-assess/acousticMove")

# files <- paste0("R/", dir("R"))
# invisible(lapply(files, source))

## Create a shape file and simulation surface:
set.seed(125)
grid <- simulate_gmrf(x = seq(-1, 1, 0.05), y = seq(-1, 1, 0.05), nhabitat = 1, kappa = 1)
grid$habitat_1 <- scale(grid$habitat_1)
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat_1)) +
  scale_fill_viridis_c("Habitat") +
  theme_bw()

seed = sample(10000, 20)
for( i in 1:20 ){
set.seed(seed[i])
grid <- simulate_gmrf(x = seq(0, 1, length = 41), y = seq(0, 1, length = 41), nhabitat = 1, kappa = 1)
grid$habitat_1 <- scale(grid$habitat_1)
  p <- ggplot(data = grid, aes(x=x, y=y)) + 
    geom_tile(aes(fill = habitat_1)) +
    scale_fill_viridis_c("Habitat") +
    theme_bw() +
    ggtitle(paste(seed[i]))
  dev.new()
  print(p)
}

set.seed(7263)
grid <- simulate_gmrf(x = seq(0, 1, length = 41), y = seq(0, 1, length = 41), nhabitat = 1, kappa = 1)
names(grid) <- c("x", "y", "habitat")
xlim <- range(grid$x) + c(1, -1)*0.05*diff(range(grid$x))
ylim <- range(grid$y) + c(1, -1)*0.05*diff(range(grid$x))
recs <- expand.grid(x = seq(xlim[1], xlim[2], length = 10), y = seq(ylim[1], ylim[2], length = 10))
recs$id <- 1:nrow(recs)

## Sim1 Example:
sim_1 <- list(statespace = grid, detectors = recs)
usethis::use_data(sim_1)


grid$habitat_1 <- scale(grid$habitat_1)
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat_1)) +
  scale_fill_viridis_c("Habitat") +
  geom_point(data = recs, aes(x = x, y = y), col = 'red', shape = 3, size = 2)
  theme_bw()

alpha <- 0.15
beta <- c(0.03, 0.10)
q <- 0.03
mu <- NULL
N <- 10
gamma <- c(0.5, 0.5)
# gamma <- cbind(runif(N,-1,1), runif(N,-1,1))
emissionrate <- 720
studyperiod <- 50
startbbox <- NULL

## Fast version of simulation:
obj <- acousticModel$new(grid = grid, detectors = recs)
obj$modelSetUp(formula = ~ 0 + habitat_1)
obj$simulate(N=3, alpha = alpha, beta=beta, q=q, gamma = gamma, emissionrate=emissionrate, studyperiod=studyperiod)
id <- 1
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat_1)) +
  geom_path(data = obj$sim_move |> subset(animal_id == id & time > 5), aes(x = x, y = y), col = 'black', linetype = 2) + 
  geom_point(data = recs, aes(x=x, y=y), shape = 3, col = 'red', size = 3) +
  scale_fill_viridis_c("Habitat") +
  # geom_point(data = data.frame(x=gamma[id,1], y=gamma[id,2]), col = 'red', size = 3, pch = 16) +
  theme_bw() +
  geom_point(data = obj$detectors, aes(x = x, y = y), col = 'red', shape = 3, size = 2)

obj$observations

obj$makeADFun(alpha, beta, q, mu, gamma, emissionrate = emissionrate, studyperiod = studyperiod, control=list())
fit <- nlminb(obj$negll$par, obj$negll$fn, obj$negll$gr)

make_ad_fun_mmpp(obj, alpha, beta, q, gamma = gamma, mu = mu, control = list())

self <- obj
control <- list()
fit <- nlminb(self$negll$par, self$negll$fn, self$negll$gr)
fit$par

# fit <- nlminb(obj$negll$par, obj$negll$fn, obj$negll$gr)
# reList(fit$par)

process_data_movement(obj, delta_t = 0.5)
make_ad_fun_move(obj, alpha, beta, gamma = gamma, mu = mu, control = list())
fit_move <- nlminb(obj$negll$par, obj$negll$fn, obj$negll$gr)
reList(fit_move$par)

# Q <- obj$calculateQ(alpha, beta, mu, gamma = gamma)
# newmat <- expm::expm(as.matrix(Q)*100)
# ggplot(data = grid, aes(x = x, y = y)) + geom_tile(aes(fill = diag(newmat))) + 
  # geom_point(data = data.frame(x = gamma[1], y = gamma[2]), aes(x = x, y = y), col = 'red')

path <- obj$solveExpectedPath(alpha = alpha, beta = beta, q = q, mu = mu, gamma = gamma, delta_t = 0.1)
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat_1)) +
  geom_path(data = obj$sim_move, aes(x = x, y = y), colour = "red") + 
  geom_path(data = path$path, aes(x = x, y = y), colour = "green") + 
  theme_bw()

Q <- obj$calculateQ(alpha, beta, mu, gamma)
expQ <- expm::expm(as.matrix(Q)*50)

ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = diag(expQ))) +
  theme_bw()

ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = path$forwardprob[1,,50])) +
  theme_bw()

p_total <- rowSums(path$forwardprob[1,,])*0.1
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(x = x, y = y, fill = p_total)) + 
  geom_path(data = obj$sim_move, aes(x = x, y = y), colour = "red", lty = 3) + 
  theme_bw()
