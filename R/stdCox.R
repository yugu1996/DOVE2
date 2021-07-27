# Standard Cox model
# @param dt A data.frame containing all time variables
# @param X A matrix with covariates
# @param t0 A scalar, set to 4 weeks by default to reflect the ramping vaccine effect after dose 1
# @return A vector of VE estimator and its standard error
#' @importFrom survival Surv coxph

VE_stdCox = function(dt, X, t0) {
  # discard all the events that occur within t0 days after entry
  ind = (dt$event_time-dt$entry_time > t0) 
  dt = dt[ind, ]; X = X[ind,,drop=FALSE]
  # measure event.time from study entry + 4 weeks
  event_time = dt$event_time - dt$entry_time - t0
  event_status = dt$event_status
  group = ifelse(dt$vaccination_time==dt$entry_time, 1, 0)
  data = as.data.frame(cbind(event_time, event_status, group, X))
  
  # fit standard Cox model
  fit = summary(survival::coxph(survival::Surv(event_time, event_status) ~ ., 
                                data = data, robust = T, ties = "breslow"))
  HR = fit$coef[1,2]; VE = 1-HR
  seVE = fit$coef[1,4]*HR
  
  return(c(VE, seVE))
}
