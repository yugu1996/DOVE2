# Cox Regression with time-varying hazard ratio
# Time is measured from the beginning of the study 
# Use piece-wise linear approximations for log hazard ratio.

# @param dt A data.frame containing all time variables
# @param X A matrix with covariates
# @param changePts A vector of change points of the piecewise linear function
#   for the log hazard ratio
# @param s.vec A vector of dates of vaccination
# @param t0 A scalar, set to 4 weeks by default to reflect the ramping vaccine effect after dose 1
# @param tau A scalar specifying the study end time

# @return A list of 9 data frames: 
#         (1) cumulative incidence (CI) with vaccination at dates specified in 's.vec'
#         (2) VE_CI estimate over time period (t0, t] with vaccination at dates specified in 's.vec'
#         (3) VE_HR estimate and 95% confidence interval

VE_newCox <- function(dt, X, changePts, s.vec, t0, tau) {
  
  n <- nrow(x = dt)
  
  if (is.null(x = X)) {
    p <- 0L
    varname <- NULL
    SD <- NULL
    X <- matrix(data = 0.0, nrow = n, ncol = 1L)
  } else {
    p <- ncol(x = X)
    varname <- colnames(x = X)
    
    # standardize covariates X (centered at median and divided by SD)
    SD <- apply(X = X, MARGIN = 2L, FUN = sd, na.rm = TRUE)
    med = apply(X = X, MARGIN = 2L, FUN = median, na.rm = TRUE)
    X <- scale(x = X, center = med, scale = SD)
  }
  
  npc <- length(x = changePts) + 1L
  
  YD <- dt$event_time
  delta = dt$event_status
  
  # sort the data by YD
  sorted_index <- order(YD)
  YD <- YD[sorted_index]
  delta <- delta[sorted_index]
  dt <- dt[sorted_index,]
  X <- X[sorted_index,,drop=FALSE]
  
  if (p == 0L) {
    funcCox <- "newCox_noX"
    args <- list("gamma" = rep(x = 0.0, times = npc), 
                 "time" = YD,
                 "delta" = delta,
                 "W" = dt$entry_time,
                 "S" = dt$vaccination_time,
                 "knots" = changePts)
  } else {
    funcCox <- "newCox"
    args <- list("beta"= rep(x = 0.0, times = p),
                 "gamma" = rep(x = 0.0, times = npc), 
                 "time" = YD,
                 "delta" = delta,
                 "X" = X,
                 "W" = dt$entry_time,
                 "S" = dt$vaccination_time,
                 "knots" = changePts)
  }
  
  # fit the standard Cox model
  
  fit.cox <- tryCatch(expr = do.call(what = funcCox, args = args),
                      error = function(e){message(e$message); return( NULL )})
  
  if (is.null(x = fit.cox)) stop("calculation aborted", call. = FALSE)
  
  if (p != 0L) {
    beta = args[[ "beta" ]]
    gamma = args[[ "gamma" ]]
    theta <- c(beta, gamma)
  } else {
    theta = gamma = args[[ "gamma" ]]
  }
  
  jump = fit.cox[[1L]]
  covar = fit.cox[[ 2L ]]
  covgamma = as.matrix(covar[p+(1L:npc), p+(1L:npc)])
  
  if (anyNA(x = theta)) {
    stop("NA values were produced in the standard Cox regression")
  } 
  
  
  # summarize the results for regression coefficients
  
  if (p != 0L) {
    # scale back beta and its SE
    beta <- beta/SD
    covbeta <- as.matrix(covar[1L:p, 1L:p])
    sebeta <- sqrt(x = diag(x = covbeta))/SD
    zbeta <- beta/sebeta
    
    beta.output <- cbind(beta, 
                         sebeta,
                         zbeta,
                         2*pnorm(q = abs(x = zbeta), lower.tail = FALSE),
                         exp(x = beta),
                         exp(x = beta-1.96*sebeta),
                         exp(x = beta+1.96*sebeta))
    
    colnames(x = beta.output) <- c("coef",
                                   "se(coef)", "z", "Pr(>|z|)",
                                   "exp(coef)",
                                   "lower .95", "upper .95")
    
    rownames(x = beta.output) <- varname
  } else {
    beta.output <- NA
  }
  

  # Calculate VE estimates 
  
  uniqt = unique(YD)
  
  getLambda = function(t) {
    ind = which(uniqt<=t)
    Lambda = sum(jump[ind])
    return(Lambda)
  }
  
  getInt = function(t, s) {
    ind = which(uniqt>=s & uniqt<=s+t)
    if(length(ind)) {
      expeta = sapply(uniqt[ind], function(x) exp(BS(x-s, changePts) %*% gamma))
      intLambda = jump[ind] %*% expeta
      return(intLambda)
    } else {
      return(0)
    }
  }
  
  # calculate cumulative incidence based on mean of X
  getCI = function(t, s) {
    CI_vaccine = 1-exp(-getInt(t, s))
    CI_placebo = 1-exp(getLambda(s)-getLambda(s+t))
    return(c(CI_placebo, CI_vaccine))
  }
  
  ns = length(s.vec)
  
  # cumulative incidence from time 0
  
  for(i in 1:ns) {
    uniqti = unique(YD[YD>s.vec[i] & delta==1])
    tempCI = t(sapply(uniqti-s.vec[i], function(t) getCI(t, s.vec[i])))
    assign(paste0("CI_", i), 
           data.frame("time" = c(0, uniqti-s.vec[i], tau),
                      "placebo" = c(0, tempCI[,1], tail(tempCI[,1], n=1)),
                      "vaccine" = c(0, tempCI[,2], tail(tempCI[,2], n=1))))
  }
  
  # VE_CI excluding events in the first t0 days
  for(i in 1:ns) {
    Sur0 = getCI(t0, s.vec[i])
    uniqti = unique(YD[YD>s.vec[i]+t0 & delta==1])
    tempCI = t(sapply(uniqti-s.vec[i], function(t) getCI(t, s.vec[i])-Sur0))
    tempVE = 1-tempCI[,2]/tempCI[,1]
    assign(paste0("VE_CI_", i), 
           data.frame("time" = c(uniqti-s.vec[i], tau),
                      "VE_CI" = c(tempVE, tail(tempVE, n=1))))
  }
  
  getVh = function(t) {
    tempZ = BS(t, changePts)
    tempeta = tempZ %*% gamma
    sdeta = sqrt(tempZ %*% covgamma %*% t(tempZ))
    est = 1-exp(tempeta)
    sdest = (1-est)*sdeta
    lower = 1-exp(tempeta+1.96*sdeta)
    upper = 1-exp(tempeta-1.96*sdeta)
    return(c(est, sdest, lower, upper))
  }
  
  t.vec = seq(0, tau, 0.1)
  Vh = t(sapply(t.vec, getVh))
  VE_HR = as.data.frame(cbind(t.vec, Vh))
  colnames(VE_HR) = c("time", "VE", "sdest", "lower", "upper")
  
  res <- list("CI_1" = CI_1, "VE_CI_1" = VE_CI_1,
              "CI_2" = CI_2, "VE_CI_2" = VE_CI_2,
              "CI_3" = CI_3, "VE_CI_3" = VE_CI_3,
              "CI_4" = CI_4, "VE_CI_4" = VE_CI_4,
              "VE_HR" = VE_HR,
              "covariates" = beta.output)
  
  
  return( res )
}
