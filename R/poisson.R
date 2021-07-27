# Standard Poisson model
# @param data dataset generated from simulate()
# @param t0 discard all events that occur within t0 days after entry, set to 28 days by default
# @return VE estimator 

VE_poisson = function(data, t0) {
  # discard all the events that occur within t0 days after entry
  ind = (data$event_time-data$entry_time > t0) 
  data = data[ind, ]
  
  # vaccine group (1) vs. placebo group (0)
  group = ifelse(data$vaccination_time==data$entry_time, 1, 0)
  
  # the numbers of events in the placebo and vaccine groups between month 1 and crossover
  n0 = sum(data$event_status[group==0])
  n1 = sum(data$event_status[group==1])
  
  # VE = 1 - incidence ratio
  VE = 1-n1/n0
  
  return(VE)
}


##################################################################################################################
### An Example of Vaccine Efficacy Estimation Based on Standard Poisson model
# t0 = 1
# VE = VE_poisson(t0, data)
###
##################################################################################################################