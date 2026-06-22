#' Simulate a Marokov Modulate Poisson process for acoustic telemetry.
#'
#' @param self R6 object containing statespace and receivers.
#' @param N Number of animals to simulate.
#' @param alpha Difussion parameter in movement process.
#' @param beta Advection parameter relating to covariates in movement process. If OU process is included, the final term relates to the OU process.
#' @param q Detection rate per unit time within a state (default = 0.3).
#' @param gamma Vector of length 2, of a point of attraction for the animals in an OU process (default = NULL). Not included when NULL.
#' @param startbbox bounding box as a vector with values \code{c('xmin', 'ymin', 'xmax', 'ymax')}
#' @param mu Mortality rate per unit time (default = NULL). No mortality assumed when NULL.
#' 
#' @export
simulate_mmpp <- function(self, N, alpha, beta, q = 0.3, gamma = NULL, startbbox = NULL, mu = NULL){
  studyperiod <- self$studyperiod

  if(is.vector(gamma)){
    if(N > 1){ gamma <- t(replicate(N, gamma))
    }else{gamma <- matrix(gamma, nrow = 1, ncol = 2)} 
  }

  m <- nrow(self$statespace)
  if(!is.null(startbbox)){
    x0poss <- self$statespace$x >= startbbox['xmin'] & self$statespace$x <= startbbox['xmax'] & 
      self$statespace$y >= startbbox['ymin'] & self$statespace$y <= startbbox['ymax']
    x0 <- sample(m, N, replace = TRUE, prob = x0poss*1)
  }else{
    x0 <- sample(m, N, replace = TRUE)
  }
  nstates <- self$nstates + !is.null(mu)

  ## check by making gamma the centroid:
  if(!is.null(gamma)){ check_gamma <- c(mean(self$statespace$x), mean(self$statespace$y))
  }else{ check_gamma <- NULL }

  test <- initCheck(self, alpha, beta, check_gamma, verbose = FALSE)
  if(test <= 0) stop("alpha value is too small relative to beta. Minimum is ", sqrt(alpha^2 + abs(test)),".\n")
  
  rookstates <- cbind(self$statespace[,grep("Rook", names(self$statespace))])  ## Add mortality state
  if(!is.null(mu)) rookstates <- cbind(rookstates, m+1)  

  haz <- numeric(m+1)
  ndets <- table(self$detectors$state_id)
  haz[as.numeric(names(ndets))] <- ndets*self$emissionrate*q
  J <- nrow(self$detectors)

  ## Quick subset of detectors in states:
  posdetectors <- lapply(self$statespace$state_id, FUN = function(i){self$detectors |> subset(state_id == i)})

  true.movement <- NULL
  observations <- list()
  for( i in 1:N ){
    Q <- make_generator(alpha = alpha, beta = beta, mu = mu, gamma = gamma[i,], nstates = nstates, dx = self$delta_xy, 
                        self$absorbingstates, self$itoj, self$designmatrix, self$design_ou)
    Qd <- -diag(Q)  ## Exit rate:
    if(any(Qd < 0)) stop("Negative exit rate (positive diagonal) for some state in the generator Q.")
    
    statei <- x0[i]
    timei <- 0
    obsi <- data.frame(animal_id = i, state_id = x0[i], detector_id = NA, time = 0)
    ti <- 0
    j <- 1
    while(ti < studyperiod){
      if(statei[j] > m) break;  ## Death occurs.
      if(statei[j] %in% self$absorbingstates){ ## Stays here the rest of the time.
        z <- studyperiod + 1
      }else{
        z <- rexp(1, Qd[statei[j]]) ## Residence time.
      }
      ## Check if stay over the study period.
      zextra <- 0
      if(ti+z > studyperiod) {
        zextra <- ti + z - studyperiod
        z <- z - zextra
      }
      ## Check obs:
      nobs <- rpois(1, z*haz[statei[j]])
      if(nobs > 0){
        detz <- posdetectors[[statei[j]]]
        smp <- sample(nrow(detz), nobs, replace = TRUE)
        obsij <- data.frame(animal_id = i, state_id = statei[j], detector_id = detz[smp, "detector_id"], time = ti + sort(runif(nobs, 0, z)))
        obsi <- rbind(obsi, obsij)
      }
      if(statei[j] %in% c(self$absorbingstates, m+1)){
        statei <- c(statei, statei[j])
      }else{
        toStates <- as.numeric(rookstates[statei[j],])
        toStates <- toStates[!is.na(toStates)]
        pij  <- Q[statei[j], toStates]
        
        # fix: when only one neighbour, no need to sample
        if (length(toStates) == 1) {
          nextstate <- toStates[1]
        } else {
          nextstate = sample(toStates, 1, prob = pij)
        }        
        statei <- c(statei, nextstate)
      }
      timei <- c(timei, ti+z)
      j <- j+1
      ti <- ti + z + zextra
    }
    true.movement <- rbind(true.movement, data.frame(animal_id = i, state_id = statei, time = timei))
    observations <- rbind(observations, obsi)
  }
  self$sim_obs <- observations
  self$sim_move <- cbind(true.movement, self$statespace[true.movement$state_id, c("x", "y")])
}

#' Simulate a Gaussian Markov Random Field
#'
#' @param x Vector of centroids of states along x axis.
#' @param y Vector of centroids of states along y axis.
#' @param nhabitat Number of habitat variables to simulate.
#' @param kappa Vector of kappa values to simulate the amount of habitat.
#' 
#' @export
simulate_gmrf <- function(x, y, nhabitat = 1, kappa = 0.95){

  nx <- length(x)
  ny <- length(y)
  df <- expand.grid( y = sort(y), x = sort(x) )
  grid <- expand.grid(y = y, x = x)

  Px <- bandSparse(n = nx, k = c(-1,1), diag = list(rep(0.5, nx), rep(0.5, nx)))
  Ix <- Diagonal( n = nx)
  Py <- bandSparse(n = ny, k = c(-1,1), diag = list(rep(0.5, ny), rep(0.5, ny)))
  Iy <- Diagonal( n = ny)
  I <- Diagonal(n = nx*ny)

  P <- kronecker(Py, Ix) + kronecker(Iy, Px)
  P <- Diagonal(n = nx*ny, x = 1/rowSums(P) ) %*% P

  if(length(kappa) == 1) kappa <- rep(kappa, nhabitat)
  for(i in 1:nhabitat) {

    IminusP <- I - kappa[i] * P
    Q <- t(IminusP) %*% (IminusP)
    habi <- RTMB:::rgmrf0( n = 1, Q = Q )
    df[[paste0("habitat_", i)]] <- habi
  }

  return(df)
}
