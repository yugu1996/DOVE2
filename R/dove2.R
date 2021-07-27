#' Estimation of the Vaccine Efficacy
#' 
#' Estimate vaccine efficacy (VE) with Kaplan-Meier (KM) estimator, standard Cox and Poisson models, 
#' and the new Cox model. 
#' Event times are potentially right-censored. 
#' Curves of the estimated VE in reducing the cumulative incidence (VE_CI) 
#' and the estimated VE in reducing the hazard rate (VE_HR) can be plotted, 
#' based on the KM analysis or the new Cox model.
#' 
#' The information required for an analysis is 
#'   \describe{
#'     \item{Entry Time:}{Calendar time when the participant entered the trial (in days).}
#'     \item{Event Time:}{Calendar time when the 
#'        participant experienced symptomatic COVID-19 or when the follow-up ended, 
#'        whichever occurred first (in days).}
#'     \item{Event Status:}{Binary indicator on whether symptomatic COVID-19 occurred 
#'        before the end of follow-up.}
#'     \item{Vaccination Status:}{Binary indicator taking value 1 if
#'       vaccination occurred before the end of follow-up and 0 otherwise.}
#'     \item{Vaccination Time:}{Calendar time when vaccination took place,
#'       with an arbitrary value if the participant was not vaccinated.}
#'     \item{Covariates:}{Baseline covariates (e.g., priority group, age, sex, ethnicity).}
#'    }
#'    
#' Note that all the time variables are measured from the beginning of the 
#' clinical trial and are specified in days. For each 
#' individual, the entry_time and event_time must satisfy
#' entry_time \eqn{\le} event_time. If entry_time > event_time, the case will be 
#' removed from the analysis and a message will be generated.  
#' 
#' The general structure of the formula input is
#'   \preformatted{
#'   Surv(event_time, event_status) ~ covariates + 
#'     vaccine(entry_time, vaccination_status, vaccination_time)
#'   }
#' 
#' The response variable must be a survival object as returned by the
#' 'Surv()' function of package \pkg{survival}, where event_time is the
#' follow up time (formal argument 'time') and event_status is the status
#' indicator input (formal argument 'event'). Specifically, 
#' \preformatted{Surv(time = event_time, event = event_status)}
#' 
#' The covariates can be either numerical or categorical.
#' A model without covariates is also allowed.
#'
#' The vaccination and entry_time information must be specified through function 
#' 'vaccine()'. Specifically, 
#' \preformatted{vaccine(entry_time, vaccination_status, vaccination_time)}
#' 
#' @rdname dove2
#' @name dove2 
#' 
#' @param formula A formula object, with the response on the left hand side of a
#'   '~' operator, and the covariates and vaccine() function on the right.  
#'   The response must be a survival object as returned by the 'Surv'
#'   function of the \pkg{survival} package. See Details for further information.
#'   The vaccine() function must be used to specify the entry time and
#'   vaccination information. See ?vaccine and Details for further
#'   information. 
#'  
#' @param data A data.frame object. The data.frame in which to interpret the
#'   variable names in formula. Must contain the entry time, the event time,
#'   the event status, the vaccination status, the vaccination time, and the covariates. 
#'   See Details.
#'   
#' @param changePts A numerical vector object. The potential change
#'   points (in days) of the piece-wise log-linear hazard ratio. The default
#'   change points are 30 and 60 days.
#'   
#' @param plots A logical object. If TRUE (default), plots of the estimated 
#'   VE_CI (based on KM estimator and the new Cox model), 
#'   the estimated VE_HR and the 95\% confidence interval 
#'   (based on the new Cox model) will be automatically generated 
#'   and saved to a PDF file named ``figure.pdf'' in the current working directory. 
#'   If FALSE, plots will not be generated.
#'   
#' @return A list containing four elements, one for each method. Specifically, 
#'   \itemize{
#'     \item{KM method: }{A list containing two elements: `Sur_km' and `VE_CI_km'. 
#'       `Sur_km' is a list containing the time points and survival probabilities of the KM estimates
#'       for both the placebo and the vaccine groups.
#'       `VE_CI_km' is a data.frame containing the estimated VE_CI at all jump points, 
#'       with all events within the first 4 weeks after the 1st dose excluded.}
#'     \item{Standard Cox model: }{A vector of the estimated vaccine efficacy and its standard error.}
#'     \item{Standard Poisson model: }{A scalar of the estimated vaccine efficacy.}
#'     \item{New Cox model: }{A list containing 10 data frames. Four of them contain the estimated
#'        cumulative incidence for individuals who were randomized on the following four dates:
#'        (1) first entry time; (2) medium entry time; (3) last entry time; (4) last entry time + 2 months. 
#'        Four of them contain the estimated VE_CI for individuals who were vaccinated at the above dates, 
#'        excluding the first 4 weeks of events. One of them contains the estimated VE_HR and the 95\% 
#'        confidence interval. 
#'        The last data.frame contains the estimation results on the hazard ratio of each covariate
#'        (NA if there's no covariate in the model).}
#'   }
#' 
#' @export
#' 
#' @include km.R newCox.R stdCox.R poisson.R plot.fig.R
#' @import Rcpp ggplot2
#' @importFrom stats model.response update.formula complete.cases terms
#' @importFrom survival Surv is.Surv survfit coxph 
#' @importFrom ggpubr ggarrange
#' 
#' @examples
#' data(dove2Data)
#' 
#' set.seed(1234)
#' smp <- sample(1L:nrow(x = dove2Data), size = 500L)
#' 
#' # Include priority and sex as covariates
#' res <- dove2(formula = Surv(event.time, event.status) ~ priority + sex + 
#'                        vaccine(entry.time, vaccine.status, vaccine.time),
#'              data = dove2Data[smp,])
#' 

dove2 = function(formula, data, changePts = c(30L, 60L), plots = TRUE) {
  
  if (missing(x = formula)) {
    stop("a formula argument must be provided", call. = FALSE)
  }
  
  if (missing(x = data)) {
    stop("a data argument must be provided", call. = FALSE)
  }
  
  # reset options to allow for keeping na values
  opt <- options()
  options(na.action = 'na.pass')
  on.exit(options(opt))
  
  # add intercept from model if not provided to ensure that factors are handled properly
  if (attr(x = stats::terms(x = formula), which = "intercept") == 0L) {
    formula = update.formula(old = formula, new = .~. +1)
  }
  
  # try to obtain the model.frame
  mf <- tryCatch(expr = stats::model.frame(formula = formula, data = data),
                 error = function(e) {
                   message("unable to obtain model.frame")
                   stop(e$message, call. = FALSE)
                 })
  
  # extract covariates
  X <- suppressMessages(stats::model.matrix(object = mf, data = data))
  # remove intercept
  int <- attr(x = X, which = "assign") != 0L
  X <- X[,int, drop = FALSE]
  
  # identify the columns that correspond to the returns returns by 
  # vaccine()
  
  lbl <- attr(x = stats::terms(x = formula), which = "term.labels")
  
  lbl1 <- paste0(lbl, "vaccination_time")
  lbl2 <- paste0(lbl, "vaccination_status")
  lbl3 <- paste0(lbl, "entry_time")
  
  if (!any(lbl1 %in% colnames(x = X)) || 
      !any(lbl2 %in% colnames(x = X)) || 
      !any(lbl3 %in% colnames(x = X))) {
    stop("the RHS of formula did not contain an appropriate vaccine() object",
         call. = FALSE)
  }
  
  i1 <- colnames(x = X) %in% lbl1
  i2 <- colnames(x = X) %in% lbl2
  i3 <- colnames(x = X) %in% lbl3
  
  vacTime <- X[,i1]
  vacStatus <- X[,i2]
  entryTime <- X[,i3]
  # added drop=FALSE to properly handle case when only 1 covariate is in the model.
  X <- X[,-c(which(i1),which(i2),which(i3)), drop = FALSE]
  
  if (ncol(x = X) == 0L) X <- NULL
  
  # dt will be a Surv object with columns "time" and "status"
  dt <- suppressMessages(stats::model.response(data = mf))
  
  if (!survival::is.Surv(x = dt)) {
    stop("the LHS of formula must be a Surv object", call. = FALSE)
  }
  
  if (ncol(x = dt) != 2L) {
    stop("the Surv object must include event time, and", 
         " event status", call. = FALSE)
  }
  
  eventTime <- dt[,1L]
  eventStatus <- dt[,2L]
  
  # remove any cases that have NA in the 
  # entry_time, event_time, event_status, covariates, or vacStatus
  use <- stats::complete.cases(cbind(X,entryTime,eventTime,eventStatus,vacStatus))
  
  if (all(!use, na.rm = TRUE)) {
    stop("input checks result in all NA -- verify inputs",
         call. = FALSE)
  }
  
  if (any(!use, na.rm = TRUE)) {
    entryTime <- entryTime[use]
    eventTime <- eventTime[use]
    eventStatus <- eventStatus[use]
    if (!is.null(x = X)) X <- X[use,, drop = FALSE]
    vacTime <- vacTime[use]
    vacStatus <- vacStatus[use]
  } 
  
  # ensure that event times are after entry times and 
  # vaccine() already tested to ensure that vacTime >= entryTime
  tst <- {eventTime < entryTime}
  
  if (any(tst, na.rm = TRUE)) {
    entryTime <- entryTime[!tst]
    eventTime <- eventTime[!tst]
    eventStatus <- eventStatus[!tst]
    if (!is.null(x = X)) X <- X[!tst,, drop = FALSE]
    vacTime <- vacTime[!tst]
    vacStatus <- vacStatus[!tst]
  } 
  
  if (sum(!use, na.rm = TRUE) > 0L || sum(tst, na.rm = TRUE) > 0L) {
    message(sum(!use | tst, na.rm = TRUE), 
            " cases removed from the analysis due to NA values")
  }
  
  # set tau to 1 plus the maximum event time
  tau <- round(max(eventTime, na.rm = TRUE)) + 1L
  message("tau set to Day ", tau)
  
  # remove NA vaccination times for convenience
  vacTime[vacStatus == 0L] <- Inf
  
  # create a data.frame containing all time variables
  dt <- data.frame("entry_time" = entryTime,
                  "event_time" = eventTime,
                  "event_status" = as.integer(x = round(x = eventStatus, digits = 0L)),
                  "vaccination_time" = vacTime)


  ### main analysis
  
  # dates of vaccination
  s = c(min(dt$entry_time), median(dt$entry_time), max(dt$entry_time), max(dt$entry_time)+61L)
  
  # potentially exclude events in the first 4 weeks
  t0 = 28L
  
  ## KM analysis
  # obtain KM estimator
  Sur_km = KM(dt)
  # VE_CI based on KM estimator
  VE_CI_km = VE_km(t0, Sur_km[[1]], Sur_km[[2]]) 
  
  # Standard Cox model
  VE.stdCox = VE_stdCox(dt, X, t0)
  
  # Standard Poisson model
  VE.poisson = VE_poisson(dt, t0)
  
  # new Cox model
  res = VE_newCox(dt, X, changePts, s, t0, tau)
  
  # generate figures
  if(plots) {
    plot.fig(Sur_km, VE_CI_km, res, tau, s)
  }
  
  return(list("KM" = list("Sur_km" = Sur_km, "VE_CI_km" = VE_CI_km), 
              "stdCox" = VE.stdCox, 
              "Poisson" = VE.poisson,
              "newCox" = res))
}