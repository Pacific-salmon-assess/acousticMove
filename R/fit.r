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

#' Fit Acoustic Model with joint loglikelihood and gradient evaluation
#'
#' @param self R6 object containing statespace and receivers.
#' @param alpha Diffusion parameter
#' @param beta Advection parameters, final value relates to OU process if gamma values provided.
#' @param q Detection rate within a state.
#' @param gamma Two values relating to the point of attraction in space (default = NULL).
#' @param mu Mortality rate per unit time (default = NULL).
#' @param control List of control values (default = list()).
#' 
#' @details Compute the negative log likelihood for the Markov modulated Poisson process and fit it using `nlminb`.
#' @return nlminb fitted object.
#'
#' @export
fit_negll <- function(self, alpha, beta, q, gamma = NULL, mu = NULL, control = list()){

  silent <- extractControls(control$silent, TRUE)
  trace <- extractControls(control$trace, 0)
  rescale_freq <- extractControls(control$rescale_freq, 25)
  tol <- extractControls(control$tolerance, 1e-8)

  nbeta <- ncol(self$designmatrix)
  
  par <- list()
  if(missing(alpha)){
    logalpha <- 0
  }else{ 
    logalpha <- log(unname(alpha))
  }
  par$logalpha <- logalpha
  if(missing(beta)) beta <- rep(0, nbeta)
  par$beta <- unname(beta)
  if(missing(q)){
    logitq <- 0
  }else{ 
    logitq <- log(q/(1-q))
  }
  par$logitq <- unname(logitq)
  if(!is.null(mu)){
    logmu <- log(mu)
    par$logmu <- unname(logmu)
  }
  
  N <- length(self$observations)

  if(is.vector(gamma)){
    gamma <- matrix(gamma, nrow = 1, ncol = 2)
  }
  if(!is.null(gamma)) par$gamma <- unname(gamma)
  
  nstates <- self$nstates + !is.null(mu)
  observations <- self$observations

  m <- nrow(self$statespace)
  xmax <- max(self$designmatrix)

  emissionrate <- as.numeric(self$emissionrate)
  studyperiod <- as.numeric(self$studyperiod)

  ## Build e^(t(Q)-Lambda)*t * v as an atomic.
  control$tranpose <- TRUE
  control$mmpp <- TRUE
  gamma_check <- NULL
  if(!is.null(gamma)) gamma_check <- c(0,0)
  calc_Q <- make_Q_rtmb(self, alpha = alpha, beta = beta, q = q, mu = mu, gamma = gamma_check, control = control) 
  theta <- as.numeric(do.call("c", par))
  Q_gr <- calc_Q$jacfun()

  x <- do.call("c", par)
  cached_env <- new.env()
  cached_env$grad <- rep(Inf, length(x))
  cached_env$pars <-  NULL
  cached_env$negll <-  NULL

  log_det_fn <- MakeTape(\(x){log(emissionrate * plogis(x[grep("logitq", names(x))]))}, x)
  log_det_gr <- log_det_fn$jacfun()

  negll_gr <- function(x){
    Q <- calc_Q(x)
    if(any(min(Q - diag(diag(Q))) < 0)) return(NaN) ## Not valid pars.
    Qgr <- Q_gr(x)    
    
    logdetrate <- log_det_fn(x)
    gr_logdetrate <- log_det_gr(x)

    N <- length(observations)
    ll <- 0
    gr <- numeric(length(x))

    for( i in 1:N ){
      timesi <- observations[[i]]$detections[,"time"]
      obstatei <- observations[[i]]$detections[,"state_id"]
      v <- numeric(nstates)
      v[obstatei[1]] <- 1
      nobsi <- observations[[i]]$ndets  ## Ensure we stop at the end of study:
      if(nobsi > 1){
        for(j in 2:nobsi){
          p_gr <- acousticMove:::expAv_gr_cpp(Q*(timesi[j]-timesi[j-1]), Qgr*(timesi[j]-timesi[j-1]), v, tol, rescale_freq, Q@i, Q@p)
          v <- numeric(nstates)
          v[obstatei[j]] <- 1
          ll <- ll + log(p_gr[["ans"]][obstatei[j]]) + logdetrate
          gr <- gr + p_gr[["gr_ans"]][obstatei[j],]/p_gr[["ans"]][obstatei[j]] + gr_logdetrate
        }
      }
      ## If known state of the animal then include it here:
      if(!is.null(observations[[i]]$known_fate)){
        deltat <- observations[[i]]$known_fate["time"]-timesi[nobsi]
        p_gr <- acousticMove:::expAv_gr_cpp(Q*deltat, Qgr*deltat, v, tol, rescale_freq, Q@i, Q@p)
        ll <- ll + log(p_gr[["ans"]][observations[[i]]$known_fate["state_id"]])
        gr <- gr + p_gr[["gr_ans"]][obstatei[j],]/p_gr[["ans"]][[observations[[i]]$known_fate["state_id"]]] 
      }else{
        ## Right censoring
        deltat <- studyperiod - timesi[nobsi]      
        p_gr <- acousticMove:::expAv_gr_cpp(Q*deltat, Qgr*deltat, v, tol, rescale_freq, Q@i, Q@p)
        logpstate <- log(p_gr[["ans"]])
        maxlp <- max(logpstate)
        log_p_marg <- log(sum(exp(logpstate - maxlp))) + maxlp
        ll <- ll + log_p_marg
        gr <- gr + colSums(p_gr[["gr_ans"]])/exp(log_p_marg)
      }
    }
    nll <- -ll

    ## Return negative log likelihood:
    assign("grad", -gr, envir = cached_env)
    if(!is.nan(nll)){
      assign("pars", x, envir = cached_env)
      assign("negll", nll, envir = cached_env) 
    }else{
      assign("pars_NaN", x, envir = cached_env) 
    }
    return(-ll)
  }
  x <- do.call('c', par)
  x <- initValues(self, x)
  ## Make sure beta is not negative for OU process... Best to confirm that should be true.
  if(!is.null(gamma) & x[grep("beta", names(x))][length(beta)] < 0) x[grep("beta", names(x))][length(beta)] <- 0.0001
  fit <- nlminb(x, negll_gr, \(x){cached_env$grad}, control = list(trace = trace))
  self$estimated_pars <- reList(fit$par)
  return(fit)
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

