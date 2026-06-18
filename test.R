library(RTMB)
library(R6)
library(Matrix)
library(ggplot2)
library(acousticMove)

## New Version:
# remotes::install_github("Pacific-salmon-assess/acousticMove")
library(acousticMove)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(RTMB)
data("sim_1")
sourceCpp("src/expAv.cpp")

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

vid <- calc_states(obj, data.frame(x=0.5, y=0.5))
v <- numeric(obj$nstates)
v[vid] <- 1
lambda <- numeric(obj$nstates)
for( i in 1:nrow(obj$detectors) ) lambda[obj$detectors$state_id[i]] <- lambda[obj$detectors$state_id[i]] + emissionrate
tmp <- expAv_approx_cpp(Q*5, v, lambda, 1e-8)
Q2 <- Q-diag(lambda)
tmp2 <- expAv(Q2*5, v)
plot(tmp[,1], tmp2)
abline(0, 1, col = 'red')
idx <- which(abs(tmp2-tmp[,1]) > 1e-4)
points(tmp[obj$detectors$state_id,1], tmp2[obj$detectors$state_id], col = 'red')

niter <- qpois(1e-8,-min(diag(Q)),lower.tail = FALSE)
# niter2 <- qpois(1e-8,-min(diag(Q2)),lower.tail = FALSE)
# Pij <- Matrix::expm(Q/niter)
rho_nested <- -min(diag(Q))/niter
nested <- qpois(1e-8, rho_nested, lower.tail = FALSE)
Pij_approx <- Q/niter + rho_nested*diag(obj$nstates)
# Pij_approx <- Q/niter/2 + diag(obj$nstates)
term <- v
for( i in 1:niter ){
  term <- exp(-lambda/niter/2)*term
  ans <- term
  for(j in 1:nested) {
    term <- Pij_approx %*% term/j
    ans <- ans + term
  }
  term <- ans*exp(-rho_nested)
  term <- exp(-lambda/niter/2)*term  
}
term2 <- expAv(Q2, v)
plot(term, term2)
abline(0,1,col='red')
points(term[obj$detectors$state_id,1], term2[obj$detectors$state_id], col = 'red')
plot(term[obj$detectors$state_id,1], term2[obj$detectors$state_id])
abline(0,1,col='red')

## Make PMat:
maxrate <- -min(diag(Q))
deltat <- (pracma::lambertWp(-(1-1e-8)/exp(1)) + 1)/maxrate
sum(dpois(0:2, maxrate*0.001))
sum(dpois(0:1, maxrate*0.001))
niters <- qpois(1e-8, -min(diag(Q)), lower.tail = FALSE)
deltat <- 1/niters/4
P <- Matrix::Matrix(0, nrow = obj$nstates, ncol = obj$nstates)
for( i in 1:obj$nstates ){
  ## No events:
  P[i,i] <- exp(-lambda[i]*deltat + Q[i,i]*deltat)

  ## One event:
  posi <- which(Q[i,] > 0)
  eq <- posi[which(abs(Q[cbind(posi, posi)] - Q[i,i]) < 1e-8)]
  pos <- posi
  if(length(eq) > 0) pos <- pos[!pos %in% eq]
  
  P[i,pos] <- P[i, pos] + Q[i, pos]/(-lambda[i] + lambda[pos] + Q[i,i] - Q[cbind(pos, pos)])*(exp((-lambda[i] + Q[i,i])*deltat) - exp((-lambda[pos] + Q[cbind(pos, pos)])*deltat))
  if(length(eq) > 0) P[i, eq] <- P[i, eq] + Q[i, eq]/(-lambda[i] + lambda[eq] + Q[i,i]) * (exp(-lambda[i]*deltat + Q[i,i]*deltat) - exp(-lambda[eq]*deltat))

  # Two events:
  # for( j in seq_along(posi) ){
    # posj <- which(Q[posi[j],] > 0)    
    # idj <- posi[j]
    
    # tmp <- Q[i, idj]*Q[idj, posj]/(-lambda[i] + lambda[idj] + Q[i,i] - Q[idj,idj])
    # tmp <- tmp*(lambda[posj] - Q[cbind(posj,posj)] - lambda[i] + Q[i,i])*(lambda[posj] - Q[cbind(posj,posj)] - lambda[idj] + Q[idj,idj])
    # tmp <- tmp*(exp((-lambda[i]+Q[i,i])*deltat) - exp((-lambda[posj]+Q[cbind(posj, posj)])*deltat) - exp((-lambda[idj]+Q[idj,idj])*deltat) + exp((-lambda[idj]+Q[cbind(idj, idj)])*deltat))
    # P[idj,posj] <- tmp
    # P[idj,pos] <- P[i,pos]*(exp(-lambda[i]*deltat + Q[i,i]*deltat) - exp(-lambda[pos]*deltat + Q[cbind(pos, pos)]*deltat))
    # if(length(eq) > 0) P[i, eq] <- Q[i, eq]/(-lambda[i] + lambda[eq] + Q[i,i]) * (exp(-lambda[i]*deltat + Q[i,i]*deltat) - exp(-lambda[eq]*deltat))
    # Q[pos[j]]
  # }
}
P2 <- P%*%P%*%P%*%P

term <- v
for( i in 1:niters ){
  term <- P2 %*% term
}
plot(term)
term2 <- expAv(Q2, v)
points(term2, col = 'red')
plot(term[lambda>0], term2[lambda>0])
abline(0, 1, col='red')

qpois(1e-8, -min(diag(Q2)), lower.tail = FALSE)

id <- 85
rates <- Q[id, Q[id,] > 0]
sum_rates <- sum(rates)
ni <- rpois(100000, sum_rates*deltat)
hist(ni)
abline(v = 1, col = 'red')
sum(P[id, -id][P[id,-id] > 0])
P[id,id] - mean(ni < 1)*exp(-lambda[id]*deltat)

qpois(1e-8, -min(diag(Q2*deltat)), lower.tail = FALSE)
term2 <- expAv(Q2*deltat, v)
plot(term2, P[v==1,])
abline(0,1,col = 'red')





library(pske)
library(tictoc)
tic()
eP <- pske::skeletoid_expm(Q, eps = 1e-8)
toc()
tic()
eP2 <- pske::unif_expm(Q, eps = 1e-8)
toc()
tmp <- skeletoid_vtexpm(Q, t_pow = 1, v=Matrix::Matrix(v, nrow = obj$nstates, ncol = 1), eps = 1e-8)
tmp2 <- unif_vtexpm(Q, t_pow = 1, v=Matrix::Matrix(v, nrow = obj$nstates, ncol = 1), eps = 1e-8)
K = pske:::skeletoid_auto_tune(Q = Q,t_pow = 1,eps = 1e-8)$K
part_K = pske:::get_opt_partition(K=K,N=obj$nstates,n3=pske:::glovars$N3,nr=1)
sp_cost = 1*Matrix::nnzero(Q)*(2^K)


qpois(1e-8, -min(diag(Q)), lower.tail = FALSE)
P <- Matrix::Matrix(0, nrow = obj$nstates, ncol = obj$nstates)
niter <- 300
deltat <- 1/niter
for( i in 1:obj$nstates ){
  pos <- which(Q[i,] > 0)
  eq <- pos[which(abs(Q[cbind(pos, pos)] - Q[i,i]) < 1e-8)]
  P[i,pos] <- Q[i, pos]/(Q[cbind(pos, pos)] - Q[i,i]) * (exp(Q[cbind(pos,pos)]*deltat) - exp(Q[i,i]*deltat))
  if(length(eq) > 0) P[i,eq] <- Q[i, eq]*deltat*exp(Q[i,i]*deltat)
  P[i,i] <- exp(Q[i,i]*deltat)
}

term <- v
for( i in 1:niter ){
  term <- P %*% term
}
term2 <- expAv(Q2, v)

ggplot(obj$statespace, aes(x=x,y=y)) + geom_tile(aes(fill = term[,]))
ggplot(obj$statespace, aes(x=x,y=y)) + geom_tile(aes(fill = term2))

plot(term[,1], term2)

plot(term)
term2 <- expAv(Q2, v)
plot(term, term2)
abline(0,1,col='red')

ggplot(data = obj$statespace) + geom_tile(aes(x=x,y=y,fill=habitat)) + 
  geom_point(data = obj$statespace[idx,], aes(x=x,y=y), col = 'red')

min(diag(Q*20))
qpois(1e-8, 1440, lower.tail = TRUE)
ppois(1, 0.5)

obj <- acousticModel$new(grid = sim_1$statespace, detectors = sim_1$detectors)
obj$modelSetUp(formula = ~ 0 + habitat)
obj$simulate(N=10, alpha = alpha, beta=beta, q=q, gamma = gamma, emissionrate=emissionrate, studyperiod=studyperiod)
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
obj$negll$fn(obj$negll$par)

tictoc::tic()
fit <- fit_negll(obj, alpha, beta, q, gamma = gamma, mu = mu, control = list(trace = 1))
tictoc::toc()

tictoc::tic()
fit2 <- fit_negll(obj, alpha, beta, q, gamma = gamma, mu = mu, control = list(trace = 1, jump_min = 1/24*0.25))
tictoc::toc()
## 371.17 seconds

tictoc::tic()
start <- initValues(obj, obj$negll$par)
fit2 <- nlminb(start, obj$negll$fn, obj$negll$gr, control = list(trace = 1))
tictoc::toc()

tictoc::tic()
obj$negll$fn(obj$negll$par)
obj$negll$gr(obj$negll$par)
tictoc::tic()

tictoc::tic()
negll_gr(x)
tictoc::toc()

fit <- nlminb(obj$negll$par, obj$negll$fn, obj$negll$gr, control = list(trace = 1))
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
