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
  if(is.null(quantiles)) quantiles <- c(0.95, 0.8, 0.5)
  plot_1 <- ggplot(data = obj$statespace, aes(x=x, y=y)) + 
    geom_tile(aes(fill = .data[[layer]])) +
    scale_fill_viridis_c(label) +
    theme_bw() + 
    geom_contour(aes(x = x, y = y, z = p, linetype = factor(after_stat(level))), 
                 breaks = quantile(p, quantiles), colour = "black", linewidth = 1) + 
    coord_fixed() +
    scale_linetype("Quantile", labels = paste0(100*sort(quantiles, decreasing = TRUE), "%")) +
    xlab("X") + ylab("Y")
  if(!is.null(gamma))  plot_1 <- plot_1 + geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y), col = 'red', size = 1, shape = 4, stroke = 2)
  invisible(print(plot_1))
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
plot_path <- function(self, deltat = 0.1, tstart = 0, tend = 1, s_init = NULL, alpha, beta, gamma, mu, layer, label){
  if(is.null(s_init)){
    s_init <- numeric(self$nstates)
    s_init[sample(self$nstates, 1)] <- 1
  }
  Q <- self$calculateQ(alpha, beta, mu, gamma)
  Pt <- expm::expm(Q*deltat)
  tx <- seq(tstart, tend, deltat)
  psum <- numeric(nrow(Q))
  pnew <- s_init
  for( i in seq_along(tx) ) {
    pnew <- pnew %*% Pt
    psum <- psum + pnew*deltat
  }
  plot_1 <- ggplot(data = obj$statespace, aes(x=x, y=y)) + 
    geom_tile(aes(fill = drop(psum))) +
    scale_fill_viridis_c("Expected Time") +
    theme_bw() + 
    coord_fixed() +
    xlab("X") + ylab("Y")
  if(!is.null(gamma))  plot_1 <- plot_1 + geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y), col = 'red', size = 1, shape = 4, stroke = 2)
  invisible(print(plot_1))
  # plot_1 <- ggplot(data = obj$statespace, aes(x=x, y=y)) + 
    # geom_tile(aes(fill = drop(pnew))) +
    # scale_fill_viridis_c(paste0("State Probability at time ", tend - tstart)) +
    # theme_bw() + 
    # coord_fixed() +
    # xlab("X") + ylab("Y")
  # if(!is.null(gamma))  plot_1 <- plot_1 + geom_point(data = data.frame(x=gamma[1], y=gamma[2]), aes(x=x,y=y), col = 'red', size = 1, shape = 4, stroke = 2)
  # invisible(print(plot_1))
}
