#' Make CTMC Generator
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion parameter
#' @param beta Advection parameters, final value relates to OU process if gamma values provided.
#' @param mu Mortality rate per unit time (default = NULL).
#' @param gamma Two values relating to the point of attraction in space (default = NULL).
#' @param nstates Number of states in CTMC.
#' @param dx Grid cell width in each state.
#' 
#' @details Build the generator matrix for the animal movement CTMC.
#'
#' @return Generator Matrix.
#'
#' @export
make_generator <- function(alpha, beta, mu = NULL, gamma = NULL, nstates, dx, absorbingstates, itoj, designmatrix, design_ou){
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
    ndesign <- ncol(designmatrix)

    Q <- AD(Matrix(0, nrow = nstates, ncol = nstates))
    xbeta <- designmatrix %*% beta[1:ndesign]
    if(!is.null(gamma)){
      xbeta <- xbeta + beta[ndesign+1]*design_ou[,1]*(gamma[self$design_ou[,2]] - design_ou[,3])
    }

    dx2 <- dx*dx
    Q[itoj] <- (alpha^2 + dx*xbeta)/(2*dx2)

    if(!is.null(mu)) {
      Q[-nstates, nstates] <- mu
    }
    if(!is.null(absorbingstates)) Q[absorbingstates,] <- 0
    
    diag(Q) <- -rowSums(Q)
    return(Q)
}

#' Make expAv Atomic for memory efficiency
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion parameter
#' @param beta Advection parameters, final value relates to OU process if gamma values provided.
#' @param q Detection rate within a state.
#' @param mu Mortality rate per unit time (default = NULL).
#' @param gamma Two values relating to the point of attraction in space (default = NULL).
#' @param control List of values to input in the `expAv` function (default = list()).
#' 
#' @details Builds the atomic function to compute expAv with memory efficiency. If \code{control$tranpose} is TRUE then we compute vexp(A) or FALSE we compute exp(A)v. 
#' \code{control$mmpp} will include the detection process or will just compute movement process.
#'
#' @return Atomic function of expAv.
#'
#' @export
make_expav_atomic <- function(self, alpha, beta, q, mu = NULL, gamma = NULL, control = list()){

  rescale_freq <- extractControls(control$rescale_freq, 25)
  tolerance <- extractControls(control$tolerance, 1e-5)
  Nmax <- extractControls(control$Nmax, 1e9)
  uniformization <- extractControls(control$uniformization, TRUE)
  trace <- extractControls(control$trace, FALSE)
  mmpp <- extractControls(control$mmpp, TRUE)
  transpose <- extractControls(control$transpose, TRUE)

  ouprocess <- !is.null(gamma)
  includemortality <- !is.null(mu)

  nstates <- self$nstates + includemortality

  theta <- c(numeric(nstates), 0, alpha, beta)
  if(mmpp) theta <- c(theta, q)
  if(includemortality) theta <- c(theta, mu)
  if(ouprocess) theta <- c(theta, gamma)
  
  delta_x <- self$resolution[1]
  if(self$resolution[1] != self$resolution[2]) stop("Currently must assume that delta_x = delta_y (square grid).")
  delta_xy <- prod(self$resolution)
  
  nalpha <- length(alpha)
  nbeta <- length(beta)
  ntheta <- length(theta)

  absorbingstates = self$absorbingstates
  itoj <- self$itoj
  designmatrix <- self$designmatrix
  design_ou <- self$design_ou

  expAv_atomic <- function(theta){
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
  
    ## Process theta vector: 'v, tdiff, alpha, beta, q, mu'
    iter <- 1
    v <- theta[iter:(iter + nstates - 1)]
    iter <- iter + nstates
    deltat <- theta[iter]
    iter <- iter + 1
    alpha <- theta[iter:(iter + nalpha - 1)]
    iter <- iter + nalpha
    beta <- theta[iter:(iter + nbeta - 1)]
    iter <- iter + nbeta
    if(mmpp){
      q <- theta[iter]
      iter <- iter + 1
    }
    if(includemortality){
      mu <- theta[iter]
      iter <- iter + 1
    }else{
      mu <- NULL
    }
    if(ouprocess){
      gamma <- theta[iter:(iter+1)]
      iter <- iter + 2
    }else{
      gamma <- NULL
    }
    
    ## Build Generator:
    Q <- make_generator(alpha, beta, mu, gamma, nstates, delta_x, absorbingstates, itoj, designmatrix, design_ou)
    
    ## Remove Detection rate from generator if MMPP.
    if(mmpp){
      detRate <- self$emissionrate*q
      for(i in seq_along(self$detectors$state_id)) Q[self$detectors$state_id[i], self$detectors$state_id[i]] <- Q[self$detectors$state_id[i], self$detectors$state_id[i]] - detRate
    }
    
    pstate <- expAv(A = Q*deltat, v = v, rescale_freq = rescale_freq, transpose = transpose, tol = tolerance, uniformization = uniformization, Nmax = Nmax, trace = trace)
    return(pstate)
  }
  F <- MakeTape(expAv_atomic, theta)
  expAv_forward <- F$atomic()
  return(expAv_forward)
}


make_expav_atomic_approx <- function(self, n = 1, delta = 0.1, alpha, beta, q, mu = NULL, gamma = NULL, control = list()){

  rescale_freq <- extractControls(control$rescale_freq, 25)
  tolerance <- extractControls(control$tolerance, 1e-5)
  Nmax <- extractControls(control$Nmax, 1e9)
  uniformization <- extractControls(control$uniformization, TRUE)
  trace <- extractControls(control$trace, FALSE)
  mmpp <- extractControls(control$mmpp, TRUE)
  transpose <- extractControls(control$transpose, TRUE)

  ouprocess <- !is.null(gamma)
  includemortality <- !is.null(mu)

  nstates <- self$nstates + includemortality

  theta <- c(numeric(nstates), 0, alpha, beta)
  if(mmpp) theta <- c(theta, q)
  if(includemortality) theta <- c(theta, mu)
  if(ouprocess) theta <- c(theta, gamma)
  
  delta_x <- self$resolution[1]
  if(self$resolution[1] != self$resolution[2]) stop("Currently must assume that delta_x = delta_y (square grid).")
  delta_xy <- prod(self$resolution)
  
  nalpha <- length(alpha)
  nbeta <- length(beta)
  ntheta <- length(theta)

  absorbingstates = self$absorbingstates
  itoj <- self$itoj
  designmatrix <- self$designmatrix
  design_ou <- self$design_ou

  expAv_atomic <- function(theta){
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")
  
    ## Process theta vector: 'v, tdiff, alpha, beta, q, mu'
    iter <- 1
    v <- theta[iter:(iter + nstates - 1)]
    iter <- iter + nstates
    deltat <- theta[iter]
    iter <- iter + 1
    alpha <- theta[iter:(iter + nalpha - 1)]
    iter <- iter + nalpha
    beta <- theta[iter:(iter + nbeta - 1)]
    iter <- iter + nbeta
    if(mmpp){
      q <- theta[iter]
      iter <- iter + 1
    }
    if(includemortality){
      mu <- theta[iter]
      iter <- iter + 1
    }else{
      mu <- NULL
    }
    if(ouprocess){
      gamma <- theta[iter:(iter+1)]
      iter <- iter + 2
    }else{
      gamma <- NULL
    }
    
    ## Build Generator:
    Q <- make_generator(alpha, beta, mu, gamma, nstates, delta_x, absorbingstates, itoj, designmatrix, design_ou)
    
    ## Remove Detection rate from generator if MMPP.
    detRate <- self$emissionrate*q

    pstate <- v
    for( i in 1:n ){
      pstate <- expAv(A = Q*delta, v = pstate, rescale_freq = rescale_freq, transpose = transpose, tol = tolerance, uniformization = uniformization, Nmax = Nmax, trace = trace)
      for(i in seq_along(self$detectors$state_id)) pstate[self$detectors$state_id[i]] <- pstate[self$detectors$state_id[i]] * exp(-detRate*delta)
    }
    for(i in seq_along(self$detectors$state_id)) Q[self$detectors$state_id[i], self$detectors$state_id[i]] <- Q[self$detectors$state_id[i], self$detectors$state_id[i]] - detRate
    pstate <- expAv(A = Q*(deltat - delta*n), v = pstate, rescale_freq = rescale_freq, transpose = transpose, tol = tolerance, uniformization = uniformization, Nmax = Nmax, trace = trace)
    return(pstate)
  }
  F <- MakeTape(expAv_atomic, theta)
  expAv_forward <- F$atomic()
  return(expAv_forward)
}