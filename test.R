library(RTMB)
library(R6)
library(Matrix)
library(ggplot2)

files <- paste0("R/", dir("R"))
invisible(lapply(files, source))


## Create a shape file and simulation surface:
set.seed(1234)
grid <- simulate_gmrf(x = seq(-1, 1, 0.05), y = seq(-1, 1, 0.05), nhabitat = 1, kappa = 0.95)
grid$habitat_1 <- scale(grid$habitat_1)

recs <- expand.grid(x = seq(-0.9, 0.9, 0.25), y = seq(-0.9, 0.9, 0.25))
recs$id <- 1:nrow(recs)

alpha <- 0.15
beta <- c(0.01, 0.15)
q <- 0.03
mu <- NULL
gamma <- c(0, 0)
emissionrate <- 720
studyperiod <- 50

## Fast version of simulation:
obj <- acousticModel$new(grid = grid, detectors = recs)
obj$modelSetUp(formula = ~ 0 + habitat_1)
obj$simulate(N=15, alpha = alpha, beta=beta, q=q, gamma = gamma, emissionrate=emissionrate, studyperiod=studyperiod)
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat_1)) +
  geom_path(data = obj$sim_move |> subset(animal_id == 1), aes(x = x, y = y, colour = factor(animal_id))) + 
  theme_bw()

obj$makeADFun(alpha, beta, q, mu, gamma, emission_rate = emissionrate, study_period = studyperiod, control=list())
fit <- nlminb(obj$negll$par, obj$negll$fn, obj$negll$gr)
reList(fit$par)

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

p_total <- rowSums(path$forwardprob[1,,])*0.1
ggplot(data = grid, aes(x=x, y=y)) + 
  geom_tile(aes(x = x, y = y, fill = p_total)) + 
  geom_path(data = obj$sim_move, aes(x = x, y = y), colour = "red", lty = 3) + 
  theme_bw()
