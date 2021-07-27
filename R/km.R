# Kaplan-Meier Analysis
# time 0 is study entry
# @param data dataset generated from simulate()
# @return A list of two data frames containing jump time points and KM estimators,
#         one for placebo group, one for vaccine group.
#' @importFrom survival Surv survfit 

KM = function(data) {
  
  Y = data$event_time - data$entry_time
  Delta = data$event_status
  
  group = ifelse(data$vaccination_time==data$entry_time, 1, 0)
  message("Number of subjects: ", length(group)) 
  message("Number of subjects in the placebo group: ", length(which(group==0)))
  message("Number of subjects in the vaccine group: ", length(which(group==1)))
  
  # vaccine group
  ind1 = which(group == 1)
  fit1 = summary(survival::survfit(survival::Surv(Y[ind1], Delta[ind1]) ~ 1))
  time1 = c(0, fit1$time)
  S1 = c(1, fit1$surv)
  
  # placebo group
  ind0 = which(group == 0)
  fit0 = summary(survival::survfit(survival::Surv(Y[ind0], Delta[ind0]) ~ 1))
  time0 = c(0, fit0$time)
  S0 = c(1, fit0$surv)
  
  df1 = data.frame("time" = time1, "Sur" = S1)
  df0 = data.frame("time" = time0, "Sur" = S0)
  
  return(list("placebo" = df0,
              "vaccine" = df1))
}


##################################################################################################################
# Calculate VE_CI based on KM estimators (time 0 is study entry)
# @param t0: if greater than 0, the program will discard all the events 
#            that occurred before entry.time + t0
# @param df0 KM estimators for the placebo group, obtained from KM()
# @param df1 KM estimators for the vaccine group, obtained from KM()

# @return A data frame containing jump time points and VE_CI estimators

VE_km = function(t0, df0, df1) {
  # obtain KM estimate at t0 for vaccine group
  S10 = tail(df1$Sur[df1$time<=t0], n=1)
  
  # obtain KM estimate at t0 for placebo group
  S00 = tail(df0$Sur[df0$time<=t0], n=1)
  
  # extract jump points after t0 for vaccine group
  ind1 = which(df1$time>t0)
  time1 = df1$time[ind1]
  
  # extract jump points after t0 for placebo group
  ind0 = which(df0$time>t0)
  time0 = df0$time[ind0]
  
  # combine all jump points after t0 in two groups
  # they are jump points of the VE_CI estimator
  time = sort(unique(time1, time0))
  
  # obtain KM estimators at these jump points
  S1 = sapply(time, function(t) tail(df1$Sur[df1$time<=t], n=1))
  S0 = sapply(time, function(t) tail(df0$Sur[df0$time<=t], n=1))
  
  # compute VE_CI estimate
  VE = 1-(S10-S1)/(S00-S0)

  return(df = data.frame("time" = time,
                         "VE" = VE))
}


##################################################################################################################
### An Example of Vaccine Efficacy Estimation Based on Kaplan-Meier Analysis
# obtain KM estimators
# Sur_km = KM(data)

# obtain estimated VE_CI
# t0 = 1
# VE_CI = VE_km(t0, Sur_km[[1]], Sur_km[[2]]) 
###
##################################################################################################################