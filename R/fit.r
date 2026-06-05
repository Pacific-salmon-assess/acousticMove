#' Vertibi Algorithm
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion parameter
#' @param beta Advection parameters, final value relates to OU process if gamma values provided.
#' @param q Detection rate within a state.
#' @param gamma Two values relating to the point of attraction in space (default = NULL).
#' @param mu Mortality rate per unit time (default = NULL).
#' @param tstart Start time of algorithm (default = 0).
#' @param tend End time of algoirthm (default = NULL).
#' @param deltat Time step to use (default = NULL).
#' @param control List of control values (default = list()).
#' 
#' @details Run the forward algorithm for each animal and detection, then the backward algorithm and then finally the Vertibi algorithm to compute the expected path.
#'
#' @return List of `forwardprob` an array for each animal each time step the probability in each location,  `backwardprob` an array for each animal each step the probability in each location,
#' `times` vector of each time step, and `path` matrix of the expected path for each animal.
#'
#' @export
vertibi_algorithm <- function(self, alpha, beta, q, gamma = NULL, mu = NULL, tstart = 0, tend = NULL, deltat = NULL, control = list()){

  N <- length(self$observations)
  if(is.vector(gamma)){
    if(N > 1){ gamma <- t(replicate(N, gamma))
    }else{gamma <- matrix(gamma, nrow = 1, ncol = 2)} 
  }

  if(is.null(tend)) tend <- self$studyperiod
  if(is.null(deltat)) deltat <- 1

  tsteps = seq(tstart, tend, by = deltat)
  nsteps <- length(tsteps)
  nstates <- self$nstates + !is.null(mu)
  
  control$tranpose <- FALSE
  control$tolerance <- 1e-8
  expAv_atomic <- make_expav_atomic(self, alpha, beta, q, mu, gamma[1,1:2], control = control)

  control$tranpose <- TRUE
  expAv_tranpose <- make_expav_atomic(self, alpha, beta, q, mu, gamma[1,1:2], control = control)

  theta <- c(numeric(nstates), 0, alpha, beta, q)
  if(!is.null(mu)) theta <- c(theta, mu)
  if(!is.null(gamma)){
    indx_gamma <- length(theta) + 1:2
    theta <- c(theta, gamma[1,1:2])
  }
  N <- length(self$observations)
  forward_all <- array(0, c(N, nstates, nsteps))
  reverse_all <- array(0, c(N, nstates, nsteps))
  path_all <- NULL

  # Q <- make_generator(self, alpha, beta, nstates, delta_x)
  # detRate <- self$emissionrate*q
  # for(i in seq_along(self$detectors$state_id)) Q[self$detectors$state_id[i], self$detectors$state_id[i]] <- Q[self$detectors$state_id[i], self$detectors$state_id[i]] - detRate

  for( k in 1:N ){
    timesk <- self$observations[[k]]$detections[,"time"]
    obstatek <- self$observations[[k]]$detections[,"state_id"]
  
    forwardprob <- matrix(0, nrow = nstates, ncol = nsteps)
    forwardprob[obstatek[1],1] <- 1

    theta[indx_gamma] <- gamma[i,1:2]

    ## Forward algorithm for MMPP.
    for(i in 2:nsteps ){
      t1 <- tsteps[i]
      t0 <- tsteps[i-1]
      indxi <- which(timesk <= t1 & timesk > t0)
      nobsi <- length(indxi)
      detsi <- obstatek[indxi]
      timesi <- timesk[indxi]
      forwardprob[,i] <- forwardprob[,i-1]
      if(nobsi > 0){
        for(j in 1:nobsi){ 
          theta[1:nstates] <- forwardprob[,i]
          theta[nstates+1] <- (timesi[j]-t0)
          forwardprob[,i] <- expAv_tranpose(theta)
          t0 <- timesi[j]
          forwardprob[detsi[j],i] <- forwardprob[detsi[j],i]*self$emissionrate*q
          forwardprob[-detsi[j],i] <- 0
        }
      }
      theta[1:nstates] <- forwardprob[,i]
      theta[nstates+1] <- (t1-t0)
      forwardprob[,i] <- expAv_tranpose(theta)
      lfp <- log(forwardprob[,i])
      maxlfp <- max(lfp)
      forwardprob[,i] <- exp(lfp - (log(sum(exp(lfp-maxlfp))) + maxlfp))
    }

    ## Backward algorithm for MMPP.
    backwardprob <- matrix(0, nrow = nstates, ncol = nsteps)
    backwardprob[,nsteps] <- 1  
    # if(!is.null(self$observations[[k]]$known_fate)){
      # backwardprob[,nsteps] <- observations[[i]]$known_fate["time"]-timesi[nobsi]
    # }
    for(i in 1:(nsteps-1) ){
      t0 <- tsteps[nsteps-i]
      t1 <- tsteps[nsteps-i+1]
      indxi <- which(timesk < t1 & timesk >= t0)
      nobsi <- length(indxi)
      detsi <- obstatek[indxi]
      timesi <- timesk[indxi]
      backwardprob[,nsteps-i] <- backwardprob[,nsteps-i+1]
      if(nobsi > 0){
        for(j in 1:nobsi){
          theta[1:nstates] <- backwardprob[,nsteps-i]
          theta[nstates+1] <- (t1 - timesi[nobsi-j+1])
          backwardprob[,nsteps-i] <- expAv_atomic(theta)
          t1 <- t1 - (t1 - timesi[nobsi-j+1])
          backwardprob[detsi[nobsi-j+1],nsteps-i] <- backwardprob[detsi[nobsi-j+1],nsteps-i]*self$emissionrate*q
          backwardprob[-detsi[nobsi-j+1],nsteps-i] <- 0
        }
      }
      theta[1:nstates] <- backwardprob[,nsteps-i]
      theta[nstates+1] <- (t1-t0)  
      backwardprob[,nsteps-i] <- expAv_atomic(theta)
      lfp <- log(backwardprob[,nsteps-i])
      maxlfp <- max(lfp)
      backwardprob[,nsteps-i] <- exp(lfp - (log(sum(exp(lfp-maxlfp))) + maxlfp))
    }
  
    ## Vertibi algorithm:
    path <- numeric(nsteps)
    path[nsteps] <- which.max(forwardprob[,nsteps])
    for(i in 1:(nsteps-1) ){
      t0 <- tsteps[nsteps-i]
      t1 <- tsteps[nsteps-i+1]
      theta[1:nstates] <- forwardprob[,nsteps-i+1]
      theta[nstates+1] <- (t1-t0)
      path[nsteps-i] <- which.max(expAv_atomic(theta))
    }
    forward_all[k,,] <- forwardprob
    reverse_all[k,,] <- backwardprob
    path_all <- rbind(path_all, cbind(data.frame(self$statespace[path,c("x", "y", "state_id")]), time = tsteps, animal_id = k))
  }
  return(list(forwardprob = forward_all, backwardprob = reverse_all, times = tsteps, path = path_all))
}

#' Build Likelihood for Markov modulated Poisson process
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion parameter
#' @param beta Advection parameters, final value relates to OU process if gamma values provided.
#' @param q Detection rate within a state.
#' @param gamma Two values relating to the point of attraction in space (default = NULL).
#' @param mu Mortality rate per unit time (default = NULL).
#' @param control List of control values (default = list()).
#' 
#' @details Compute the negative log likelihood for the Markov modulated Poisson process. Defines an AD model object as part of the `self` R6 object.
#'
#' @export
make_ad_fun_mmpp <- function(self, alpha, beta, q, gamma = NULL, mu = NULL, control = list()){

  silent <- extractControls(control$silent, TRUE)
  nbeta <- ncol(self$designmatrix)
  
  par <- list()
  if(missing(beta)) beta <- rep(0, nbeta)
  par$beta <- beta
  if(missing(alpha)){
    logalpha <- 0
  }else{ 
    logalpha <- log(alpha)
  }
  par$logalpha <- logalpha
  if(missing(q)){
    logitq <- 0
  }else{ 
    logitq <- log(q/(1-q))
  }
  par$logitq <- logitq
  if(!is.null(mu)){
    logmu <- log(mu)
    par$logmu <- logmu
  }
  
  N <- length(self$observations)

  if(is.vector(gamma)){
    gamma <- matrix(gamma, nrow = 1, ncol = 2)
  }
  if(!is.null(gamma)) par$gamma <- gamma
  
  nstates <- self$nstates + !is.null(mu)
  observations <- self$observations

  ## nstates + time + alpha + beta + q + mu + gamma
  ntheta <- nstates + 1 + 1 + length(beta) + 1 + length(mu) + length(gamma[1,])
  alpha_indx <- nstates + 1 + 1
  beta_indx <- nstates + 1 + 1 + 1:length(beta)
  q_indx <- max(beta_indx) + 1
  mu_indx <- max(q_indx) + 1
  gamma_indx <-  max(q_indx) + length(mu) + 1:2

  m <- nrow(self$statespace)
  xmax <- max(self$designmatrix)

  emissionrate <- as.numeric(self$emissionrate)
  studyperiod <- as.numeric(self$studyperiod)

  ## Build e^(Q-Lambda)*t * v as an atomic.
  control$tranpose <- TRUE
  control$mmpp <- TRUE
  gamma_check <- NULL
  if(!is.null(gamma)) gamma_check <- c(0,0)
  expAv_forward <- make_expav_atomic(self, alpha = alpha, beta = beta, q = q, mu = mu, gamma = gamma_check, control = control)
  
  mmpp_negll <- function(par){
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")

    getAll(par)
    alpha <- exp(logalpha)
    nbeta <- length(beta)
    q <- 1/(1+exp(-logitq))

    ## theta = v, time, pars:
    theta <- AD(numeric(ntheta))
    theta[alpha_indx] <- alpha
    theta[beta_indx] <- beta
    theta[q_indx] <- q

    ## Add in mortality if there is any.
    if(!is.null(mu)){
      mu <- exp(logmu)
      theta[mu_indx] <- mu
    }

    detRate <- emissionrate*q
    logdetRate <- log(detRate)
        
    N <- length(observations)
    ll <- 0
    for( i in 1:N ){
      if(!is.null(gamma)){
        if(nrow(gamma) == 1){ theta[gamma_indx] <- drop(gamma[1, 1:2])
        }else{ theta[gamma_indx] <- drop(gamma[i, 1:2])}
      }    
      timesi <- observations[[i]]$detections[,"time"]
      obstatei <- observations[[i]]$detections[,"state_id"]
      v <- AD(numeric(nstates))
      v[obstatei[1]] <- 1
      nobsi <- observations[[i]]$ndets  ## Ensure we stop at the end of study:
      theta[1:nstates] <- v
      nobsi <- observations[[i]]$ndets
      if(nobsi > 1){
        for(j in 2:nobsi){
          theta[nstates+1] <- (timesi[j]-timesi[j-1])
          pstate <- expAv_forward(theta)
          v <- AD(numeric(nstates))
          v[obstatei[j]] <- 1
          theta[1:nstates] <- v          
          ll <- ll + log(pstate[obstatei[j]]) + logdetRate
        }
      }
      ## If known state of the animal then include it here:
      if(!is.null(observations[[i]]$known_fate)){
        theta[nstates+1] <- observations[[i]]$known_fate["time"]-timesi[nobsi]
        pstate <- expAv_forward(theta)
        ll <- ll + log(pstate[observations[[i]]$known_fate["state_id"]])
      }else{
        ## Right censoring
        theta[nstates+1] <- studyperiod - timesi[nobsi]      
        pstate <- expAv_forward(theta)
        logpstate <- log(pstate)
        maxlp <- max(logpstate)
        ll <- ll + log(sum(exp(logpstate - maxlp))) + maxlp
      }
    }
    ADREPORT(q)
    ## Return negative log likelihood:
    -ll
  }
  nll <- RTMB::MakeADFun(func = mmpp_negll,
                         parameters = par,
                         silent=silent)                         
  self$negll <- nll
}

#' Build Likelihood for CTMC animal movement.
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion parameter
#' @param beta Advection parameters, final value relates to OU process if gamma values provided.
#' @param gamma Two values relating to the point of attraction in space (default = NULL).
#' @param mu Mortality rate per unit time (default = NULL).
#' @param control List of control values (default = list()).
#' 
#' @details Compute the negative log likelihood for the CTMC animal movement. Defines an AD model object as part of the `self` R6 object.
#'
#' @export
make_ad_fun_move <- function(self, alpha, beta, gamma = NULL, mu = NULL, control = list()){

  silent <- extractControls(control$silent, TRUE)
  nbeta <- ncol(self$designmatrix)
  
  par <- list()
  if(missing(beta)) beta <- rep(0, nbeta)
  par$beta <- beta
  if(missing(alpha)){
    logalpha <- 0
  }else{ 
    logalpha <- log(alpha)
  }
  par$logalpha <- logalpha
  if(!is.null(mu)){
    logmu <- log(mu)
    par$logmu <- logmu
  }
  if(is.vector(gamma)){
    gamma <- matrix(gamma, nrow = 1, ncol = 2)
  }
  if(!is.null(gamma)) par$gamma <- gamma
  
  nstates <- self$nstates + !is.null(mu)
  observations <- self$observations
  
  ## nstates + time + alpha + beta +  mu + gamma
  ntheta <- nstates + 1 + 1 + length(beta) + length(mu) + length(gamma[1,])
  alpha_indx <- nstates + 1 + 1
  beta_indx <- nstates + 1 + 1 + 1:length(beta)
  mu_indx <- max(beta_indx) + 1
  gamma_indx <-  max(beta_indx) + length(mu) + 1:2

  m <- nrow(self$statespace)
  xmax <- max(self$designmatrix)

  ## Build e^(Q-Lambda)*t * v as an atomic.
  control$tranpose <- TRUE
  control$mmpp <- FALSE
  gamma_check <- NULL
  if(!is.null(gamma)) gamma_check <- c(0,0)
  expAv_forward <- make_expav_atomic(self, alpha = alpha, beta = beta, mu = mu, gamma = gamma_check, control = control)

  mmpp_negll <- function(par){
    ####### Sometimes necessary to avoid rtmb errors #######
    # "[<-" <- RTMB::ADoverload("[<-")
    # "diag<-" <- RTMB::ADoverload("diag<-")
    # "c" <- RTMB::ADoverload("c")

    getAll(par)
    alpha <- exp(logalpha)
    nbeta <- length(beta)

    ## theta = v, time, pars:
    theta <- AD(numeric(ntheta))
    theta[alpha_indx] <- alpha
    theta[beta_indx] <- beta

    ## Add in mortality if there is any.
    if(!is.null(mu)){
      mu <- exp(logmu)
      theta[mu_indx] <- mu
    }
    
    N <- length(observations)
    ll <- 0
    for( i in 1:N ){
      if(!is.null(gamma)){
        if(nrow(gamma) == 1){ theta[gamma_indx] <- drop(gamma[1, 1:2])
        }else{ theta[gamma_indx] <- drop(gamma[i, 1:2])}
      }    
      timesi <- observations[[i]]$detections[,"time"]
      obstatei <- observations[[i]]$detections[,"state_id"]
      v <- AD(numeric(nstates))
      v[obstatei[1]] <- 1
      nobsi <- observations[[i]]$ndets  ## Ensure we stop at the end of study:
      theta[1:nstates] <- v
      nobsi <- observations[[i]]$ndets
      if(nobsi > 1){
        for(j in 2:nobsi){
          theta[nstates+1] <- (timesi[j]-timesi[j-1])
          pstate <- expAv_forward(theta)
          v <- AD(numeric(nstates))
          v[obstatei[j]] <- 1
          theta[1:nstates] <- v          
          ll <- ll + log(pstate[obstatei[j]])
        }
      }
      ## If known state of the animal then include it here:
      if(!is.null(observations[[i]]$known_fate)){
        theta[nstates+1] <- observations[[i]]$known_fate["time"]-timesi[nobsi]
        pstate <- expAv_forward(theta)
        ll <- ll + log(pstate[observations[[i]]$known_fate["state_id"]])
      }else{
        ## Right censoring
        theta[nstates+1] <- self$studyperiod - timesi[nobsi]      
        pstate <- expAv_forward(theta)
        logpstate <- log(pstate)
        maxlp <- max(logpstate)
        ll <- ll + log(sum(exp(logpstate - maxlp))) + maxlp
      }
    }
    ## Return negative log likelihood:
    -ll
  }
  
  nll <- RTMB::MakeADFun(func = mmpp_negll,
                         parameters = par,
                         silent=silent)                         
  self$negll <- nll
}

