##' Plot surfaces
##'
##' @param self R6 object containing statespace and receivers.
##'
# plot <- function(self, layer, label){
  # p1 <- ggplot(data = self$statespace, aes(x=x, y=y)) +
    # geom_tile(aes_string(fill = layer)) +
    # coord_equal() +
    # scale_fill_viridis_c(label) +
    # xlab("Easting (m)") + ylab("Northing (m)") +
    # ggsidekick::theme_sleek()
  # invisible(print(p1))
# }

#' Plot Limiting Distribution
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion rate.
#' @param beta Advection rates.
#' @param mu Mortality rate.
#' @param gamma Centre of attraction.
#' @param Layer Habitat layer name.
#' @param Label to name the habitat in the plot.
#' @param quantiles Quantiles to plot for the limting distribution (default = c(0.95, 0.8, 0.5))
#' 
#' @details Plot habitat variable and the limiting distribution as a contour.
#' 
#' @export
plot_limit <- function(self, alpha, beta, gamma, mu, layer, label, quantiles = NULL){
  p <- calcLimit(self, alpha, beta, gamma, mu)
  p[p < 0] <- 0
  p <- p/sum(p)
  pord <- order(p)
  p <- p[pord]
  statespace <- self$statespace[pord,]
  cdf <- cumsum(p)
  if(is.null(quantiles)) quantiles <- c(0.95, 0.8, 0.5)
  breaks <- NULL
  for( i in seq_along(quantiles) ) breaks <- c(breaks, min(p[cdf >= 1 - quantiles[i]]))
  plot_1 <- ggplot(data = statespace, aes(x=x, y=y)) + 
    geom_tile(aes(fill = .data[[layer]])) +
    scale_fill_viridis_c(label) +
    theme_bw() + 
    geom_contour(aes(x = x, y = y, z = p, linetype = factor(after_stat(level))), 
                 breaks = breaks, colour = "black", linewidth = 1) + 
    coord_fixed() +
    scale_linetype("Quantile", labels = paste0(100*sort(quantiles, decreasing = TRUE), "%")) +
    xlab("X") + ylab("Y")
  if(!is.null(gamma))  plot_1 <- plot_1 + geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y), col = 'red', size = 1, shape = 4, stroke = 2)
  invisible(print(plot_1))
  return(plot_1)
}

#' Plot Expected Path
#'
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion rate.
#' @param beta Advection rates.
#' @param mu Mortality rate.
#' @param gamma Centre of attraction.
#' 
#' @details Solve for limiting distribution pi, as pi %*% Q = 0, and sum(pi) = 1.
#' 
#' @return Vector of limiting probabilities, assuming that a steady state solution exists.
#' @export
plot_path <- function(self, deltat = 0.1, tstart = 0, tend = 1, s_init = NULL, alpha, beta, gamma, mu, expected = TRUE){
  if(is.null(s_init)){
    s_init <- numeric(self$nstates)
    s_init[sample(self$nstates, 1)] <- 1
  }
  Q <- self$calculateQ(alpha, beta, mu, gamma)
  tx <- seq(tstart, tend, deltat)
  psum <- numeric(nrow(Q))
  pnew <- s_init
  for( i in seq_along(tx) ) {
    pnew <- acousticMove:::expAv_cpp(Q*deltat, pnew, 1e-8, 25, TRUE)
    psum <- psum + pnew*deltat
  }
  if(expected){
    plot_1 <- ggplot(data = self$statespace, aes(x=x, y=y)) + 
      geom_tile(aes(fill = as.numeric(psum))) +
      scale_fill_viridis_c("Expected Time") +
      theme_bw() + 
      coord_fixed() +
      xlab("X") + ylab("Y")
    if(!is.null(gamma))  plot_1 <- plot_1 + geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y), col = 'red', size = 1, shape = 4, stroke = 2)
    invisible(print(plot_1))
  }else{
    plot_1 <- ggplot(data = self$statespace, aes(x=x, y=y)) + 
      geom_tile(aes(fill = as.numeric(pnew))) +
      scale_fill_viridis_c(paste0("State Probability at time ", tend - tstart)) +
      theme_bw() + 
      coord_fixed() +
      xlab("X") + ylab("Y")
    if(!is.null(gamma))  plot_1 <- plot_1 + geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y), col = 'red', size = 1, shape = 4, stroke = 2)
    invisible(print(plot_1))
  }
}

#' Plot Expected Time along with Path
#'
#'
#' @param self R6 object containing statespace and receivers.
#' @param deltat Time steps to discretize
#' @param tstart Time to start computation.
#' @param tend Time to end computation.
#' @param s_init Initial location of the animal as a vector e.g. (0,0,0,1).
#' @param s_end Final location of the animal as a vector e.g. (0,0,0,1).
#' @param alpha Diffusion rate.
#' @param beta Advection rates.
#' @param mu Mortality rate.
#' @param gamma Centre of attraction.
#' 
#' @details Solve for limiting distribution pi, as pi %*% Q = 0, and sum(pi) = 1.
#' 
#' @return Vector of limiting probabilities, assuming that a steady state solution exists.
#' @export
plot_expected_time <- function(self, deltat = 0.1, tstart = 0, tend = 1, s_init = NULL, s_end = NULL, alpha, beta, gamma, mu, q, emissionrate){
  Q <- self$calculateQ(alpha, beta, mu, gamma)
  dt <- (tend - tstart)
  ndiv <- dt %/% deltat
  deltat <- dt/ndiv
  ## Midpoint Rule to integrate v1*exp(Qs)exp(Q(t-s))v2ds to then compute expected time.
  tx <- seq(tstart + deltat/2, tend - deltat/2, deltat)
  
  lambda <- numeric(self$nstates)
  for( i in 1:nrow(self$detectors) ) lambda[self$detectors$state_id[i]] <- lambda[self$detectors$state_id[i]] + q*emissionrate
  A <- Q - diag(lambda)
  
  nsteps <- length(tx)  
  pforward <- matrix(acousticMove:::expAv_cpp(A*(tx[1]-tstart), s_init, 1e-8, 25, TRUE), ncol=1)
  tf <- tx[1]
  for( i in 1:(nsteps-1) ){
    pforward <- cbind(pforward, acousticMove:::expAv_cpp(A*(tx[i+1]-tx[i]), pforward[,i], 1e-8, 25, TRUE))
    tf <- c(tf, tx[i+1])
  }

  preverse <- acousticMove:::expAv_cpp(A*(tend-tx[nsteps]), s_end, 1e-8, 25, FALSE)
  expected_time <- deltat*pforward[,nsteps]*preverse
  for( i in 1:(nsteps-1)) {
    preverse <- acousticMove:::expAv_cpp(A*(tx[nsteps-i+1]-tx[nsteps-i]), preverse, 1e-8, 25, FALSE)    
    expected_time <- expected_time + deltat*pforward[,nsteps-i]*preverse
  }
  
  pall <- acousticMove:::expAv_cpp(A*dt, s_init, 1e-8, 25, TRUE)[which.max(s_end)]
  pall <- acousticMove:::expAv_cpp(A*(tend-tx[nsteps]), pforward[,nsteps], 1e-8, 25, TRUE)[which.max(s_end)]
  expected_time <- expected_time/pall
  # sum(expected_time)
  
  ## Make path:
  v <- s_end
  path <- numeric(nsteps + 2)  
  path[nsteps+1] <- which.max(s_end)
  path[1] <- which.max(s_init)
  for(i in 1:(nsteps-1) ){
    if(i == 1 ){  d <- tend - tx[nsteps]
    }else{ d <- tx[i+1]-tx[i] }
    tmp <- pforward[,nsteps-i+1]*acousticMove:::expAv_cpp(A*d, v, 1e-8, 25, FALSE)
    maxk <- which.max(tmp)
    v <- numeric(self$nstates)
    v[maxk] <- 1
    path[nsteps-i+1] <- maxk
  }
  locs_path <- data.frame(self$statespace[path,])
  locs <- data.frame(self$statespace[which(s_init > 0 | s_end > 0),])
  
  plot_1 <- ggplot(data = self$statespace, aes(x=x, y=y)) + 
    geom_tile(aes(fill = expected_time)) +
    scale_fill_viridis_c("Expected Time") +
    theme_bw() + 
    coord_fixed() +
    xlab("X") + ylab("Y") + 
    geom_point(data = self$statespace[self$detectors$state_id,], aes(x = x, y = y), shape = 3, col = 'grey', size = 1) +
    geom_point(data = locs, aes(x = x, y = y), col = 'red', shape = 1, size = 4) +
    geom_path(data = locs_path, aes(x = x, y = y), col = 'red')# +
    # geom_point(data = locs_path, aes(x = x, y = y), col = 'red')    
  invisible(print(plot_1))
  return(plot_1)
}
