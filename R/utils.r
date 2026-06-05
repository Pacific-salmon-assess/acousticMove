#' Process Acoustic Telemetry Data for Modelling
#'
#' @param self R6 object containing statespace and receivers.
#' @param detections Detection history dataframe that includes x,y coordinates of detectors (default = NULL). Default implies simulated data.
#' @param releases Dataframe of where the animals were released (Default NULL). 
#' @param known_fates Dataframe of where the animals ended up (Default NULL).
#' @param absorbing_states Vector of which states are considered absorbing, such as entering a river system (default = NULL).
#'
#' @export
process_data <- function(self, detections = NULL, releases = NULL, known_fates = NULL, absorbing_states = NULL){
  if(is.null(detections)){
    detections <- self$sim_obs
  }
  
  lookup_animal_id <- unique(detections$animal_id)
  animal_id <- as.numeric(factor(lookup_animal_id))
  names(animal_id) <- lookup_animal_id
  self$lookup_animal_id <- lookup_animal_id
  
  detections$animal_id <- animal_id[detections$animal_id]
  N <- max(animal_id)
  
  ## Add state_id to detections:
  if(!"state_id" %in% names(detections)){
    ## Remove all detections that are not within the statespace:
    detections[, c("x", "y")] 
    detections$state_id <- calc_states(self, locs = detections[, c("x", "y")])
    if(any(is.na(detections$state_id))){
      cat("Removing detections outside of the statespace.\n")
      detections <- detections |> subset(!is.na(state_id))
    }
  }
  if(any(is.na(detections$state_id))){
    navals <- which(is.na(detections$state_id))
    cat("Removing", length(navals), "that had NA state_id.\n")
  }
  
  if(!is.null(absorbing_states)){
    self$absorbingstates <- unique(calc_states(self, absorbing_states[, c("x", "y")]))
  }else{
    self$absorbingstates <- NULL
  }
    
  ## Do same thing with known_fates and releases:
  if(!is.null(known_fates)){
    known_fates$animal_id <- animal_id[known_fates$animal_id]
    known_fates$state_id <- calc_states(self, known_fates[, c("x", "y")])
    if(any(known_fates$fate == "failure")){
      known_fate$state_id <- self$nstates + 1 ## Mortality state for a known loss.
    }
  }
  if(!is.null(releases)){
    releases$animal_id <- animal_id[releases$animal_id]
    releases$state_id <-  calc_states(self, releases[, c("x", "y")])
    if(!"time" %in% releases){
      cat("Assuming all releases occur at time 0.\n")
      releases$time <- 0
    }
  }
  
  observations <- list()
  for( i in 1:N ){
    deti <- detections |> subset(animal_id == i)
    deti <- deti[order(deti$time),] 
    ## Add release state if it's provided separately:
    if(!is.null(releases)){
      reli <- releases |> subset(animal_id == i)
      reli$state_id <- NA
      deti <- rbind(reli, deti)
    }   
    observations[[i]] <- list()
    observations[[i]]$detections <- as.matrix(deti[, c("time", "state_id")])
    observations[[i]]$ndets <- nrow(deti)
    ## Add known fates if provided
    if(!is.null(known_fates)){
      fatei <- known_fates |> subset(animal_id == i)
      if(nrow(fatei) == 0) observations[[i]]$known_fate  <- NULL
      else observations[[i]]$known_fate <- c("time" = fatei$time,  "state_id" = fatei$state_id)
    }
  }
  self$observations <- observations
}

#' Process Telemetry Data for Modelling
#'
#' @param self R6 object containing statespace and receivers.
#' @param detections Detection history dataframe that includes x,y coordinates of detectors (default = NULL). Default implies simulated data.
#' @param releases Dataframe of where the animals were released (Default NULL). 
#' @param known_fates Dataframe of where the animals ended up (Default NULL).
#' @param absorbing_states Vector of which states are considered absorbing, such as entering a river system (default = NULL).
#'
#' @export
process_data_movement <- function(self, delta_t = 1, absorbing_states = NULL){
  detections <- self$sim_move
  
  lookup_animal_id <- unique(detections$animal_id)
  animal_id <- as.numeric(factor(lookup_animal_id))
  names(animal_id) <- lookup_animal_id
  self$lookup_animal_id <- lookup_animal_id
  
  detections$animal_id <- animal_id[detections$animal_id]
  N <- max(animal_id)
    
  ## Add state_id to detections:
  if(!"state_id" %in% names(detections)){
    ## Remove all detections that are not within the statespace:
    detections[, c("x", "y")] 
    detections$state_id <- calc_states(self, locs = detections[, c("x", "y")])
    if(any(is.na(detections$state_id))){
      cat("Removing detections outside of the statespace.\n")
      detections <- detections |> subset(!is.na(state_id))
    }
  }
  if(any(is.na(detections$state_id))){
    navals <- which(is.na(detections$state_id))
    cat("Removing", length(navals), "that had NA state_id.\n")
  }
  
  if(!is.null(absorbing_states)){
    self$absorbingstates <- unique(calc_states(self, absorbing_states[, c("x", "y")]))
  }else{
    self$absorbingstates <- NULL
  }
      
  observations <- list()
  for( i in 1:N ){
    deti <- detections |> subset(animal_id == i)
    deti <- deti[order(deti$time),]
    ## Thin the times to the detection window.
    subset <- deti$time %/% delta_t
    deti <- deti[which(!duplicated(subset)),]
    deti <- deti |> within(time <- floor(time/delta_t)*delta_t)
    ## Add release state if it's provided separately:
    observations[[i]] <- list()
    observations[[i]]$detections <- as.matrix(deti[, c("time", "state_id")])
    observations[[i]]$ndets <- nrow(deti)
  }
  self$observations <- observations
}


#' Process State Space
#'
#' @param self R6 object containing statespace and receivers.
#' @param statespace Dataframe of grid cells with x,y values that define the discrete space and includes any covariates.
#' @param detectors Dataframe of detector locations.
#'
#' @details Process the locations in both state space and receivers to discrete space. Annotate each with which state is left, right, up and down.
#'
#' @export
process_states_detectors <- function(self, statespace, detectors){
  statespace$state_id <- 1:nrow(statespace)
  if(!all( c("x", "y") %in% names(statespace)) ){
    stop("Coordinates must be named 'x' and 'y' in statespace data frame.")
  }
  if(!all( c("x", "y") %in% names(detectors)) ){
    stop("Coordinates must be named 'x' and 'y' in detectors data frame.")
  }

  ## Find resolution:
  ydiff <- statespace$y |> diff() |> abs() |> round(7)
  yres <- min(ydiff[ydiff != 0])
  xdiff <- statespace$x |> diff() |> abs() |> round(7)
  xres <- min(xdiff[xdiff != 0])
  
  ## Add an extra state to allow for indexing out of the statespace for convenience of the +/-.
  xmin <- min(statespace$x) - 3*xres/2
  ymin <- min(statespace$y) - 3*yres/2
  self$xymin <- c(xmin, ymin)

  ## Identities:
  idx_state <- floor((statespace$x - xmin) / xres) + 1
  idy_state <- floor((statespace$y - ymin) / yres) + 1
  mat <- matrix(NA, nrow = max(idx_state)+1, ncol = max(idy_state)+1)
  mat[cbind(idx_state, idy_state)] <- statespace$state_id
  dirs <- data.frame(
    "RookL" = mat[cbind(idx_state-1, idy_state)], 
    "RookR" = mat[cbind(idx_state+1, idy_state)], 
    "RookD" =  mat[cbind(idx_state, idy_state-1)],
    "RookU" =  mat[cbind(idx_state, idy_state+1)]
  )
  
  drop <- rowSums(is.na(dirs)) == 4
  if(any(drop)){
    cat("Warning: Removing isolated states that have no rook adjacency.\n")
    statespace <- statespace[!drop,]
    ## Reorder:
    statespace$state_id <- 1:nrow(statespace)
    idx_state <- idx_state[!drop]
    idy_state <- idy_state[!drop]
    mat <- matrix(NA, nrow = max(idx_state)+1, ncol = max(idy_state)+1)
    mat[cbind(idx_state, idy_state)] <- statespace$state_id
    dirs <- data.frame(
      "RookL" = mat[cbind(idx_state-1, idy_state)], 
      "RookR" = mat[cbind(idx_state+1, idy_state)], 
      "RookD" =  mat[cbind(idx_state, idy_state-1)],
      "RookU" =  mat[cbind(idx_state, idy_state+1)]
    )
  }
  statespace <- cbind(statespace, dirs)
  idxt <- floor((detectors$x - xmin) / xres) + 1
  idyt <- floor((detectors$y - ymin) / yres) + 1
  keep <- which(idxt > 0 & idxt <= max(idx_state) & idyt > 0 & idyt <= max(idy_state))
  if(length(keep) != nrow(detectors)) cat("Warning: Removing detectors that are not contained with the statespace.\n")

  self$detectors <- cbind(detectors[keep,], state_id = mat[cbind(idxt, idyt)[keep,]]) |> subset(!is.na(state_id))
  detector_names <- unique(self$detectors$detector_id)
  lookup_detector_id <- as.numeric(factor(detector_names))
  names(lookup_detector_id) <- detector_names
  self$lookup_detector_id <- lookup_detector_id
  self$detectors$detector_id <- lookup_detector_id[self$detectors$detector_id]
  self$detectors <- self$detectors[order(self$detectors$detector_id),]
  self$resolution <- c("x" = xres, "y" = yres)
  self$statespace <- statespace
  self$lookup_statespace <- mat
  self$nstates <- nrow(statespace)
}

#' Calculate state ID
#'
#' @param self R6 object containing statespace and receivers.
#' @param Locs Dataframe of locations x,y.
#' 
#' @details Compute which state a location is within.
#'
#' @export
calc_states <- function(self, locs){
  idx1 <- floor((locs$x - self$xymin[1]) / self$resolution["x"]) + 1
  idy1 <- floor((locs$y - self$xymin[2]) / self$resolution["y"]) + 1
  good <- which(idx1 %in% 1:nrow(self$lookup_statespace) & idy1 %in% 1:ncol(self$lookup_statespace))
  states <- rep(NA, nrow(locs))
  states[good] <- self$lookup_statespace[cbind(idx1[good], idy1[good])]
  states
}

#' Compute the advection-diffusion covariate space
#'
#' @param self R6 object containing statespace and receivers.
#' @param formula Formula (must exclude intercept) for which coviarate to include.
#' 
#' @details Generate a design matrix that is computed as the central difference of the covariates for each direction Left/Right and Up/Down.
#'
#' @export
generator_design_gr <- function(self, formula){
  vars <- attr(terms.formula(formula), "term.labels")
  if(attr(terms.formula(formula), "intercept") == 1){
    cat("[Warning] An intercept term is not allowed and is being removed.\n")
    formula <- update(formula, ~ . - 1)
  }

  ## Make sure everything is ordered:
  self$statespace <- self$statespace[order(self$statespace$state_id),]
  ids <- self$statespace$state_id
  delta_x <- self$resolution[1]
  DF <- data.frame()
  ij <- data.frame()
  
  L <- self$statespace[self$statespace[, c("RookL")],]
  R <- self$statespace[self$statespace[, c("RookR")],]
  U <- self$statespace[self$statespace[, c("RookU")],]
  D <- self$statespace[self$statespace[, c("RookD")],]
  
  ## Central Difference (f(x+h) - f(x-h))/(2*delta)
  grad <- list()
  grad[["LR"]] <- (R[, vars, drop = FALSE] - L[, vars])/(2*delta_x)
  grad[["UD"]] <- (U[, vars, drop = FALSE] - D[, vars])/(2*delta_x)

  
  ## IF R/L doesn't exist and only one direction, then use an approximation a single direction.
  Lb <- which(is.na(self$statespace[, c("RookL")]))
  grad[["LR"]][Lb,] <- (R[Lb, vars, drop = FALSE] - self$statespace[Lb, vars])/(delta_x)
  Rb <- which(is.na(self$statespace[, c("RookR")]))
  grad[["LR"]][Rb,] <- (self$statespace[Rb, vars, drop = FALSE] - L[Rb, vars, drop = FALSE])/(delta_x)
  Db <- which(is.na(self$statespace[, c("RookD")]))
  grad[["UD"]][Db,] <- (U[Db, vars, drop = FALSE] - self$statespace[Db, vars, drop = FALSE])/(delta_x)    
  Ub <- which(is.na(self$statespace[, c("RookU")]))
  grad[["UD"]][Ub,] <- (self$statespace[Ub, vars, drop = FALSE] - D[Ub, vars, drop = FALSE])/(delta_x)

  ij <- data.frame()  
  jnames <- grep("Rook", names(self$statespace), value = TRUE)
  dir <- c("RookL" = -1, "RookR" = 1, "RookD" = -1, "RookU" = 1)
  coord <- c("RookL" = 1, "RookR" = 1, "RookD" = 2, "RookU" = 2)
  df_dir <- c("RookL" = "LR", "RookR" = "LR", "RookD" = "UD", "RookU" = "UD")
  jxy <- c("RookL" = "x", "RookR" = "x", "RookD" = "y", "RookU" = "y")
  
  for( i in 1:4 ){
    j <- self$statespace[, jnames[i]]
    deltaf <- data.frame(dir[jnames[i]]*grad[[df_dir[i]]])
    xy <- as.numeric(self$statespace[[jxy[i]]][ids])
    names(deltaf) <- vars
    ij <- rbind(ij, data.frame(i = ids, j = j, dir = as.numeric(dir[jnames[i]]), coord = as.numeric(coord[jnames[i]]), xy = xy))
    DF <- rbind(DF, deltaf)
  }
  naj <- which(is.na(ij$j)) ## Drop the ones that do not connect:
  DF <- DF[-naj, , drop = FALSE]
  ij <- ij[-naj, , drop = FALSE]

  self$design_ou <- ij[,3:5]
  ij <- ij[,1:2]
  
  X <- model.matrix(formula, data = DF)
  self$designmatrix <- as.matrix(X)
  self$itoj <- as.matrix(ij)
}

#' Extract Default Values
#'
#' @param controlValue Any value which includes NULL.
#' @param defaultValue A set, non NULL value to return if controlValue is NULL.
#' 
#' @return controlValue if it is not NULL, or defaultValue.
#'
#' @export
extractControls <- function(controlValue, defaultValue){
  if(!is.null(controlValue))  
    return(controlValue)
  else 
    return(defaultValue)
}

#' Re-List parameters
#'
#' @param pars vector from a fitted object, or RTMB object parameters.
#' 
#' @return List of parameters back on the real scale.
#'
#' @export
reList <- function(pars){
  out <- list()
  out$alpha <- as.numeric(exp(pars[names(pars) == "logalpha"]))
  out$beta <- as.numeric(pars[names(pars) == "beta"])
  if("logitq" %in% names(pars)) out$q <- as.numeric(plogis(pars[names(pars) == "logitq"]))
  if("gamma" %in% names(pars)) out$gamma <- matrix(as.numeric(pars[names(pars) == "gamma"]), ncol = 2)
  if("logmu" %in% names(pars)) out$mu <- as.numeric(exp(pars[names(pars) == "logmu"]))
  return(out)
}

#' Check initial values
#'
#' @param self R6 object for fitting Acoustic telemetry data.
#' @param alpha Difussion parameter in movement process.
#' @param beta Advection parameter relating to covariates in movement process. If OU process is included, the final term relates to the OU process.
#' @param gamma Vector of length 2, of a point of attraction for the animals in an OU process (default = NULL). Not included when NULL.
#' @param verbose Whether or not to print a warning if values are not valid.
#'
#' @details Test if alpha^2 - dx*(max(abs(beta %*% X))) is negative.
#'  
#' @export
initCheck <- function(self, alpha, beta, gamma = NULL, verbose = FALSE){
  dx <- self$resolution[1]
  ndesign <- ncol(self$designmatrix)
  xbeta <- self$designmatrix %*% beta[1:ndesign]
  if(!is.null(gamma)){
    xbeta <- xbeta + beta[ndesign+1]*self$design_ou[,1]*(gamma[self$design_ou[,2]] - self$design_ou[,3])
  }
  test <- alpha^2 - dx * max(abs(xbeta))
  if(test <= 0 & verbose){
    cat("[Warning] alpha value is too small relative to beta. Minimum is ", sqrt(alpha^2 + abs(test)),".\n")
  }
  return(test)
}

#' Initialize random starting values
#'
#' @param self R6 object for fitting Acoustic telemetry data.
#' @param pars vector from a fitted object, or RTMB object parameters.
#' 
#' @return pars vector ready for fitting.
#'
#' @export
initValues <- function(self, pars){
  pars["logalpha"] <- pars["logalpha"] + rnorm(1, 0, 0.2)
  pars[names(pars) == "beta"] <- pars[names(pars) == "beta"] + rnorm(length(pars[names(pars) == "beta"]), 0, 0.2)
  pars[names(pars) == "logitq"] <- pars[names(pars) == "logitq"] + rnorm(1, 0, 0.5)
  pars[names(pars) == "logmu"] <- pars[names(pars) == "logmu"] + rnorm(1, 0, 0.5)
  pars[names(pars) == "gamma"] <- pars[names(pars) == "gamma"] + rnorm(sum(names(pars) == "gamma"), 0, 3)*self$resolution[1]
  pars_norm <- reList(pars)
  test <- initCheck(self, pars_norm$alpha, pars_norm$beta, drop(pars_norm$gamma[1,]), FALSE)
  if(test <= 0) pars["logalpha"] <- log(sqrt(pars_norm$alpha^2 + abs(test))) + 0.01
  return(pars)
}