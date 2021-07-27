## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
opt <- options()
options(continue="  ", width=70, prompt=" ")
on.exit(options(opt))
library(DOVE2, quietly=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  vaccine(entry_time, vaccination_status, vaccination_time)

## ----eval=FALSE---------------------------------------------------------------
#  dove2(formula, data, changePts = c(30L, 60L), plots = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  Surv(event_time, event_status) ~ covariates +
#    vaccine(entry_time, vaccination_status, vaccination_time)

## -----------------------------------------------------------------------------
data(dove2Data)

## -----------------------------------------------------------------------------
head(dove2Data)

## -----------------------------------------------------------------------------
result <- dove2(formula = Surv(event.time, event.status) ~ priority + sex + 
                  vaccine(entry.time, vaccine.status, vaccine.time), 
                data = dove2Data)

## ----pressure, echo=FALSE, fig.align='center', fig.cap="\\label{fig: dove2Figs}Estimation of vaccine efficacy in a clinical trial: A. Kaplan-Meier estimates of the cumulative incidence of disease for the vaccine and placebo groups and the cumulative incidence curves for participants who were vaccinated at various dates under the new Cox model; B. Estimates of $VE_{CI}$ based on the Kaplan-Meier method and the new Cox model for participants who were vaccinated at various dates; C. Estimates and $95\\%$ confidence intervals for $VE_{HR}$ under the new Cox model.", out.width = '100%'----
knitr::include_graphics("figure.pdf")

## -----------------------------------------------------------------------------
str(result, max=2)

## -----------------------------------------------------------------------------
result$newCox$covariates

