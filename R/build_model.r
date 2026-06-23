#' @import RTMB
#' @importFrom R6 R6Class 
#' @import Matrix
NULL

#' Base class for Acoustic Telemetry Models
#'
#' @description An R6 base class that is used to fit a Markov modulated Poisson process to acoustic telemetry data.
#'
#' @field detectors Dataframe of detector locations and states.
#' @field studyperiod Study period stored for convenience.
#' @field statespace Dataframe of state space centroids, indices, and their neighbourhood structure.
#' @field nstates Number of physical location state (excluding mortality state).
#' @field resolution Width of grid cell in a state.
#' @field xymin Lower corner of the statespace.
#' @field itoj State i to j indices for building the generator matrix.
#' @field design_ou Matrix that shows direction coordinates and x,y values at same scale as itoj.
#' @field designmatrix Design matrix for advection on itoj for building generator matrix.
#' @field absorbingstates Holds which states are physically holding if they are supplied.
#' @field sim_obs Simulated observations when running a simulation.
#' @field sim_move Simulated true movement when running a simulation.
#' @field observations Processed observations ready for analysis.
#' @field emissionrate Rate that the tags phyiscally ping.
#' @field pars initial values for fitting the model.
#' @field negll Negative log likelihood RTMB object.
#' @field gr_negll Gradient of the negative log likelihood.
#' @field lookup_animal_id Lookup matrix for matching animal id.
#' @field lookup_detector_id Lookup matrix for matching detectors.
#' @field lookup_statespace Lookup matrix for matching state space and locations.
#' @field estimated_pars List of estimated parameters.
#' @field delta_xy Distance between states pairwise.
#'
#' @export
acousticModel <- R6::R6Class("acousticModel",
  public = list(
    # --- Fields ---  
    detectors = NULL,
    studyperiod = NULL,
    statespace = NULL,
    nstates = NULL,
    resolution = NULL,
    xymin = NULL,
    itoj = NULL,
    design_ou = NULL,
    designmatrix = NULL,
    absorbingstates = NULL,
    sim_obs = NULL,
    sim_move = NULL,
    observations = NULL,
    emissionrate = NULL,
    pars = NULL,
    negll = NULL,
    gr_negll = NULL,
    lookup_animal_id = NULL,
    lookup_detector_id = NULL,
    lookup_statespace = NULL,
    estimated_pars = NULL,
    delta_xy = NULL,

    #' @description Initialize the R6 object
    #' @param grid State space that represents where the animal can move and any covariates of interest, must contain x,y columns.
    #' @param detectors Dataframe of detector locations, must contain columns x,y.
    #' @param formula Formula for covariates to relate statespace to advection (optional).
    initialize = function(grid, detectors, formula = NULL){
      if(missing(detectors)) stop("Need to provide a detector dataframe with columns x, y for locations.")
      if(missing(grid)) stop("Need to provide the state space via 'grid', which contains x, y for the centroid of each state.")
      if(!all(c('x', 'y') %in% names(grid))) stop("Add locations 'x', 'y' to the grid.")
      if(!"detector_id" %in% names(detectors)) detectors$detector_id <- 1:nrow(detectors)
      detectors <- detectors[order(detectors$detector_id),]
      
      process_states_detectors(self, grid, detectors)
      
      ## Allow user to pass formula in on initialize.
      if(!is.null(formula)) generator_design_gr(self, formula)
    },
    #' @description Process Statespace and Detectors
    #' @param statespace Dataframe that represents where the animal can move and any covariates of interest, must contain x,y columns.
    #' @param detectors Dataframe of detector locations, must contain columns x,y.
    #' @details Create indices to the state space and match with detectors. Also computes match Left, Right, Up, Down moves. Output is stored in R6 object.
    processStatesDetectors = function(statespace, detectors){
      process_states_detectors(self, statespace, detectors)
    },
    #' @description Set up advection covariates
    #' @param formula Formula to relate glm like relationship of advection to the covariates.
    #' @details Uses central differences to compute the local gradient of each covariate to include in advection model.
    modelSetUp = function(formula = NULL){
      generator_design_gr(self, formula)      
    },
    #' @description Simulate
    #' @param N Number of animals.
    #' @param alpha Diffusion rate.
    #' @param beta Advection rates.
    #' @param q Detection probability.
    #' @param gamma Centre of attraction.
    #' @param formula Formula to relate covariates to advection (optional).
    #' @param emissionrate Rate that the tags send a signal for detection.
    #' @param studyperiod Time period of the study.
    #' @param absorbingstates Vector of states indicating which locations the animal will not leave if enters.
    #' @param startbbox Vector indicating an area that the animal taggging is taking place.
    #' @param mu Rate of mortality per unit time.
    #' @details Simulate from the Markov modulated Poisson process data that is equivalent to the acoustic telemetry study using diffusion-advection for movement. See \code{?simulate_mmpp} for more details.
    simulate = function(N  = 1, alpha = 0, beta = NULL, q = NULL, gamma = NULL, formula = NULL, emissionrate = 30, studyperiod = NULL, absorbingstates = NULL, startbbox = NULL, mu = NULL){
      ## Update formula if needed.
      if(!is.null(formula)) generator_design(self, formula)
      if(!is.null(studyperiod)) self$studyperiod <- studyperiod
      if(!is.null(emissionrate)) self$emissionrate <- emissionrate

      test <- initCheck(self, alpha, beta, gamma, verbose = FALSE)
      if(test <= 0) stop("alpha value is too small relative to beta. Minimum is ", sqrt(alpha^2 + abs(test)),".\n")

      self$absorbingstates <- absorbingstates
      simulate_mmpp(self, N=N, alpha=alpha, beta=beta, q=q, gamma=gamma, startbbox = startbbox, mu = mu)
      if(!is.null(self$sim_obs)) process_data(self)
    },
    #' @description Make negative log likelihood RTMB object.
    #' @param alpha Initial diffusion rate.
    #' @param beta Initial advection rates.
    #' @param q Initial detection probability.
    #' @param mu Initial mortality rate.
    #' @param gamma Initial centre of attraction.
    #' @param emissionrate Rate that the tags send a signal for detection.
    #' @param studyperiod Time period of the study.
    #' @param control List of object controls which include tolerance, rescale_freq, Nmax, uniformization, trace that controls the RTMB function \code{expAv}.
    #' @details Simulate from the Markov modulated Poisson process data that is equivalent to the acoustic telemetry study using diffusion-advection for movement. See \code{?simulate_mmpp} for more details.
    makeADFun = function(alpha, beta, q, mu = NULL, gamma = NULL, emissionrate = 720, studyperiod = NULL, control=list()){
      self$emissionrate <- emissionrate
      if(!is.null(studyperiod)) self$studyperiod <- studyperiod
      make_ad_fun_mmpp(self, alpha, beta, q, gamma, mu, control)
    },
    #' @description Solve Expected Paths for Observed Detections
    #' @param alpha Initial diffusion rate.
    #' @param beta Initial advection rates.
    #' @param q Initial detection probability.
    #' @param mu Initial mortality rate.
    #' @param gamma Initial centre of attraction.
    #' @param tstart Start time of solving.
    #' @param tend Time end of solving.
    #' @param delta_t Change in time to step through.
    #' @param control List of object controls which include tolerance, rescale_freq, Nmax, uniformization, trace that controls the RTMB function \code{expAv}.
    #' @details Run forward and backward algorithm to compute probability of where the animals were seen. Then use the Vertibi algorithm to estimate the most likely path the animal takes.
    solveExpectedPath = function(alpha, beta, q, mu = NULL, gamma = NULL, tstart = 0, tend = NULL, delta_t = NULL, control = list()){
      res <- vertibi_algorithm(self, alpha, beta, q, gamma, mu, tstart, tend, delta_t, control = list())
      return(res)
    },
    #' @description Compute the Movement Generator
    #' @param alpha Diffusion rate.
    #' @param beta Advection rates.
    #' @param mu Mortality rate.
    #' @param gamma Centre of attraction.
    calculateQ = function(alpha, beta, mu = NULL, gamma = NULL){
      nstates <- self$nstates + !is.null(mu)
      Q <- make_generator(alpha = alpha, beta = beta, mu = mu, gamma = gamma, nstates = nstates, dx = self$delta_xy, 
                          self$absorbingstates, self$itoj, self$designmatrix, self$design_ou)
      return(Q)
    },
    #' @description Process either data from real observations or from the simulation.
    #' @param detections Dataframe for when and where each animal detection occurs.
    #' @param releases Where and when the animals were released.
    #' @param known_fates Time and states of known fates of animals.
    #' @param absorbing_states States that are absorbing if the animal arrives.
    #' @param studyperiod How long was the study run for.
    #' @details See \code{?process_data} for details.
    processData = function(detections = NULL, releases = NULL, known_fates = NULL, absorbing_states = NULL, studyperiod = NULL){
      process_data(self, detections, releases, known_fates, absorbing_states, studyperiod)
    },
    #' @description Simulate and Fit Acoustic Telemetry Data
    #' @param N Number of animals.
    #' @param alpha Diffusion rate.
    #' @param beta Advection rates.
    #' @param q Detection probability.
    #' @param gamma Centre of attraction.
    #' @param emissionrate Rate that the tags send a signal for detection.
    #' @param studyperiod Time period of the study.
    #' @param mu Rate of mortality per unit time.
    #' @param absorbingstates Vector of states indicating which locations the animal will not leave if enters.
    #' @param startbbox Vector indicating an area that the animal taggging is taking place.
    #' @param formula Formula to relate covariates to advection (optional).
    #' @param control List of object controls which include tolerance, rescale_freq, Nmax, uniformization, trace that controls the RTMB function \code{expAv}.
    #' @details Run simulation and then fit using `nlminb`.
    simFit = function(N, alpha, beta, q, gamma = NULL, emissionrate = 720, studyperiod, mu = NULL, absorbingstates = NULL, startbbox = NULL, formula = NULL, control = list()){
      self$simulate(N, alpha, beta, q, gamma, formula, emissionrate, studyperiod, absorbingstates, startbbox, mu)
      # self$makeADFun(alpha = alpha, beta = beta, q = q, mu, gamma, control = control)
      # pars <- initValues(self, self$negll$par)
      trace <- extractControls(control$trace, 0)
      self$emissionrate <- emissionrate
      self$studyperiod <- studyperiod
      make_negll_gr(self, alpha, beta, q, gamma, mu, control = control)
      fit <- suppressWarnings(nlminb(self$pars, self$negll, self$gr_negll, control = list(trace = trace)))
      self$estimated_pars <- reList(fit$par)
      return(fit)
    },
    #' @description Simulate and Fit Movement Using Telemetry Data on a fixed Interval.
    #' @param N Number of animals.
    #' @param delta_t Time between detections of the telemetry data.
    #' @param alpha Diffusion rate.
    #' @param beta Advection rates.
    #' @param gamma Centre of attraction.
    #' @param emissionrate Rate that the tags send a signal for detection.
    #' @param studyperiod Time period of the study.
    #' @param mu Rate of mortality per unit time.
    #' @param absorbingstates Vector of states indicating which locations the animal will not leave if enters.
    #' @param startbbox Vector indicating an area that the animal taggging is taking place.
    #' @param formula Formula to relate covariates to advection (optional).
    #' @param control List of object controls which include tolerance, rescale_freq, Nmax, uniformization, trace that controls the RTMB function \code{expAv}.
    #' @details Run simulation only if simulation hasn't been run. Then process data as fixed interval telemetry and then fit using `nlminb`.
    fitMove = function(N, delta_t, alpha, beta, gamma = NULL, emissionrate = 720, studyperiod, mu = NULL, absorbingstates = NULL, startbbox = NULL, formula = NULL, control = list()){
      if(is.null(sim_move)) self$simulate(N, alpha, beta, q = 0.1, gamma, formula, emissionrate, studyperiod, absorbingstates, startbbox, mu)
      process_data_movement(self, delta_t)
      make_ad_fun_move(self, alpha, beta, gamma, mu, control)
      pars <- initValues(self, self$negll$par)
      fit <- nlminb(pars, self$nll$fn, self$nll$gr) #suppressWarnings(
      self$estimated_pars <- reList(pars = fit$par)
      return(fit)
    },
    #' @description Build Model with Joint Likelihood Gradient
    #' @param alpha Diffusion rate.
    #' @param beta Advection rates.
    #' @param q Detection probability.
    #' @param gamma Centre of attraction.
    #' @param emissionrate Rate that the tags send a signal for detection.
    #' @param studyperiod Time period of the study.
    #' @param mu Rate of mortality per unit time.
    #' @param control List of object controls which include tolerance, rescale_freq, Nmax, uniformization, trace that controls the RTMB function \code{expAv}.
    #' @details Build negative log likelihood and gradient to be used for optimization. 
    buildModel = function(alpha, beta, q, gamma = NULL, emissionrate = 720, studyperiod, mu = NULL, control = list()){
      self$emissionrate <- emissionrate
      self$studyperiod <- studyperiod
      make_negll_gr(self, alpha, beta, q, gamma, mu, control = control)
    }
  )
)