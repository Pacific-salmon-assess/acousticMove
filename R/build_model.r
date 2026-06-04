acousticModel <- R6Class("acousticModel",
  public = list(
    ## Fields
    detectors = NULL,
    detections = NULL,
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
    negll = NULL,
    lookup_animal_id = NULL,
    lookup_detector_id = NULL,
    lookup_statespace = NULL,
    estimated_pars = NULL,

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
    processStatesDetectors = function(statespace, detectors){
      process_states_detectors(self, statespace, detectors)
    },
    modelSetUp = function(formula = NULL){
      generator_design_gr(self, formula)      
    },
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
    makeADFun = function(alpha, beta, q, mu = NULL, gamma = NULL, emission_rate = 720, study_period = NULL, control=list()){
      self$emissionrate <- emission_rate
      if(!is.null(study_period)) self$studyperiod <- study_period
      make_ad_fun_mmpp(self, alpha, beta, q, gamma, mu, control)
    },
    solveExpectedPath = function(alpha, beta, q, mu = NULL, gamma = NULL, tstart = 0, tend = NULL, delta_t = NULL, control = list()){
      res <- vertibi_algorithm(self, alpha, beta, q, gamma, mu, tstart, tend, delta_t, control = list())
      return(res)
    },
    calculateQ = function(alpha, beta, mu = NULL, gamma = NULL){
      nstates <- self$nstates + !is.null(mu)
      Q <- make_generator(self, alpha = alpha, beta = beta, mu = mu, gamma = gamma, nstates = nstates, dx = self$resolution[1])
    },
    processData = function(detections = NULL, releases = NULL, known_fates = NULL, absorbing_states = NULL){
      process_data(self, detections, releases, known_fates, absorbing_states)
    },
    simFit = function(N, alpha, beta, q, gamma = NULL, emissionrate = 720, studyperiod, mu = NULL, absorbingstates = NULL, startbbox = NULL, formula = NULL, control = list()){
      self$simulate(N, alpha, beta, q, gamma, formula, emissionrate, studyperiod, absorbingstates, startbbox, mu)
      self$makeADFun(alpha = alpha, beta = beta, q = q, mu, gamma, control = control)
      pars <- initValues(self, self$negll$par)
      fit <- nlminb(pars, self$nll$fn, self$nll$gr)
      self$estimated_pars <- reList(pars = fit$par)
      return(fit)
    },
    fitMove = function(N, delta_t, alpha, beta, gamma = NULL, emissionrate = 720, studyperiod, mu = NULL, absorbingstates = NULL, startbbox = NULL, formula = NULL, control = list()){
      if(is.null(sim_move)) self$simulate(N, alpha, beta, q = 0.1, gamma, formula, emissionrate, studyperiod, absorbingstates, startbbox, mu)
      process_data_movement(self, delta_t)
      make_ad_fun_move(self, alpha, beta, gamma, mu, control)
      pars <- initValues(self, self$negll$par)
      fit <- nlminb(pars, self$nll$fn, self$nll$gr)
      self$estimated_pars <- reList(pars = fit$par)
      return(fit)
    }
  )
)