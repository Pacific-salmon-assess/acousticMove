library(RTMB)
library(R6)
library(Matrix)
library(ggplot2)
# library(acousticMove)

# remotes::install_github("Pacific-salmon-assess/acousticMove")

# sourceCpp("src/expAv.cpp")

## New Version:
# remotes::install_github("Pacific-salmon-assess/acousticMove")
library(acousticMove)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
data("sim_1")

alpha <- 0.15
beta <- c(0.03, 0.1)
q <- 0.03
mu <- NULL
gamma <- c(0.5, 0.5)
emissionrate <- 720
studyperiod <- 50

obj <- acousticModel$new(grid = sim_1$statespace, detectors = sim_1$detectors)
obj$modelSetUp(formula = ~ 0 + habitat)
Q <- obj$calculateQ(alpha, beta, mu, gamma)
v <- numeric(obj$nstates)
v[obj$detectors$state_id[1]] <- 1
ev <- acousticMove:::expAv_cpp(Q*10, v, 1e-8, 5, FALSE)
evR <- RTMB::expAv(Q*10, v)
plot(ev, evR)

ev2 <- expAv_cpp(Q*10, v, 1e-8, 5, TRUE)
evR2 <- expAv(t(Q)*10, v, control = list(tranpose = TRUE, rescale_freq = 5))
plot(ev2 - evR2)
plot(evR2, ev3[,1])
plot(ev2, ev3[,1])

eQ <- Matrix::expm(Q*10)
ev3 <- t(eQ) %*% v

obj <- acousticModel$new(grid = sim_1$statespace, detectors = sim_1$detectors)
obj$modelSetUp(formula = ~ 0 + habitat)
obj$simulate(N=15, alpha = alpha, beta=beta, q=q, gamma = gamma, emissionrate=emissionrate, studyperiod=studyperiod)
id <- 1
ggplot(data = obj$statespace, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat)) +
  geom_path(data = obj$sim_move |> subset(animal_id == id & time > 5), aes(x = x, y = y), col = 'black', linetype = 2) + 
  geom_point(data = obj$detectors, aes(x=x, y=y), shape = 3, col = 'red', size = 3) +
  scale_fill_viridis_c("Habitat") +
  geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y),  col = 'purple', size = 4, shape = 20) +
  theme_bw()

plot_limit(obj, alpha = alpha, beta = beta, gamma = gamma, mu = NULL, "habitat", "Habitat")
v <- numeric(obj$nstates)
v[obj$statespace$x == 0.5 & obj$statespace$y == 0.5] <- 1
plot_path(obj, deltat = 0.1, tstart = 0, tend = 0.5, s_init = v, alpha = alpha, beta = beta, gamma = gamma, mu = NULL, expected = TRUE)
plot_path(obj, deltat = 0.1, tstart = 0, tend = 0.5, s_init = v, alpha = alpha, beta = beta, gamma = gamma, mu = NULL, expected = FALSE)

## Fit a model:
obj$makeADFun(alpha, beta, q, mu, gamma, emissionrate = emissionrate, studyperiod = studyperiod)
obj$negll$gr(obj$negll$par)
fit <- nlminb(obj$negll$par, obj$negll$fn, obj$negll$gr)
reList(fit$par)

Q <- obj$calculateQ(alpha, beta, mu, gamma)
expQ <- expm::expm(as.matrix(Q)*100)
ggplot(data = obj$statespace, aes(x=x, y=y)) + 
  geom_tile(aes(fill = diag(expQ))) +
  theme_bw() + 
  scale_fill_viridis_c("Utility") +
  geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y),  col = 'purple', size = 4, shape = 20)

## USE Lapack to get Eigen Decomp for limiting dist.

ggplot(data = obj$statespace, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat)) +
  theme_bw() + 
  scale_fill_viridis_c("Utility")

M <- expm::expm(Q)
eigQ <- eigen(t(Q))
mine <- which.min(abs(Re(eigQ$values)))
P <- abs(Re(eigQ$vectors[, 1]))
P <- P/sum(P)
ggplot(data = obj$statespace, aes(x=x, y=y)) + 
  geom_tile(aes(fill = prob )) +
  theme_bw() + 
  scale_fill_viridis_c("Utility")


calcLimit <- function(alpha, beta, gamma, mu){
  Q <- obj$calculateQ(alpha, beta, mu, gamma)

  ## Linear Algebra solve: pi*Q = 0 and sum(pi) == 0.
  A <- rbind(t(Q), rep(1, nrow(Q)))
  b <- c(numeric(nrow(Q)), 1)
  prob <- solve(A[-1,], b[-1])
  return(prob)

  ## Eigen Version: Slower
  if(FALSE){
    eP <- eigen(t(Q))
    prob <- as.numeric(eP$vectors[,which.min(abs(Re(eP$values)))])
    prob <- prob / sum( prob )
  }
}
p <- calcLimit(alpha, beta, gamma, NULL)
contours <- quantile(p, 0.95) 
xy <- obj$statespace[which(p > contours - 0.001 & p < contours + 0.001),]
ggplot(data = obj$statespace, aes(x=x, y=y)) + 
  geom_tile(aes(fill = habitat)) +
  scale_fill_viridis_c("Habitat") +
  theme_bw() + 
  geom_contour(aes(x = x, y = y, z = p, linetype = factor(after_stat(level))), 
               breaks = quantile(p, c(0.95, 0.8, 0.5)), colour = "black", linewidth = 1) + 
  coord_fixed() +
  scale_linetype("Quantile", labels = c("95%", "80%", "50%")) +
  geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y),  col = 'red', size = 3, shape = 4) +
  xlab("X") + ylab()




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
