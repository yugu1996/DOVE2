# Plot the figure 
#' @import ggplot2
#' @importFrom ggpubr ggarrange

plot.fig = function(Sur_km, VE_CI_km, res, tau, s.vec) {
 
  ex=0.03 
  martb = 0
  marlr = 20
  mar0 = 5
  curvesize = 0.625
  
  # Figure A: cumulative incidence based on KM estimator and new Cox model
  
  # KM analysis
  df.km_placebo = Sur_km$placebo 
  df.km_placebo$CI = (1-df.km_placebo$Sur)*100
  df.km_placebo = rbind(df.km_placebo, c(tau, min(df.km_placebo$Sur), max(df.km_placebo$CI)))
  
  df.km_vaccine = Sur_km$vaccine
  df.km_vaccine$CI = (1-df.km_vaccine$Sur)*100
  df.km_vaccine = rbind(df.km_vaccine, c(tau, min(df.km_vaccine$Sur), max(df.km_vaccine$CI)))
  
  # new Cox model
  S = c(rep(1, nrow(res$CI_1)),
        rep(2, nrow(res$CI_2)),
        rep(3, nrow(res$CI_3)),
        rep(4, nrow(res$CI_4)))
  df.newCox = as.data.frame(rbind(res$CI_1, res$CI_2, res$CI_3, res$CI_4)) 
  df.newCox[,2:3] = 100*df.newCox[,2:3]
  df.newCox$S = factor(S, levels = c(0:4))
  colnames(df.newCox) = c("time", "placebo", "vaccine", "S")
  
  ymax = max(c(max(df.km_placebo$CI), max(df.newCox$placebo)))
  ymax = ymax - ymax %% 10L + 10.0
  
  fig.A = ggplot2::ggplot()+
    ggplot2::geom_step(data = df.newCox, ggplot2::aes(x=time, y=placebo, group = S, color = S), size = curvesize)+
    ggplot2::geom_step(data = df.newCox, ggplot2::aes(x=time, y=vaccine, group = S, color = S), size = curvesize)+
    ggplot2::geom_step(data = df.km_placebo, ggplot2::aes(x=time, y=CI), size = curvesize)+
    ggplot2::geom_step(data = df.km_vaccine, ggplot2::aes(x=time, y=CI), size = curvesize)+
    ggplot2::scale_x_continuous(name="Days Since Dose 1", limits=c(0, tau+15), 
                       breaks = seq(0, tau+15, 50), expand = c(0,0))+
    ggplot2::scale_y_continuous(name="Cumulative Incidence (%)", limits = c(0, ymax+2), 
                       breaks = seq(0, ymax+2, 5), expand=c(0,ex))+
    ggplot2::scale_color_manual(name = NULL,
                       values = c("black", "red", "green", "blue", "orange"),
                       labels = c("Kaplan-Meier estimate", paste0("Vaccinated on Day ", s.vec)),
                       drop = FALSE) +
    ggplot2::theme(axis.line.x = ggplot2::element_line(arrow=ggplot2::arrow(length = ggplot2::unit(3, 'mm'))), 
          axis.line.y=ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(3, 'mm'))),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t=6)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r=8)),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
          legend.key = ggplot2::element_blank(),
          legend.position = c(0.32, 0.81),
          legend.box = "horizontal",
          legend.margin = ggplot2::margin(t=-5),
          legend.spacing = ggplot2::unit(0.15, "cm"),
          plot.margin = ggplot2::margin(l = mar0*4/3, r = marlr))+
    ggplot2::annotate("text", x = tau*0.9, y = max(c(max(df.km_placebo$CI), max(df.newCox$placebo)))+0.4, label = "Placebo")+
    ggplot2::annotate("text", x = tau*0.9, y = max(c(max(df.km_vaccine$CI), max(df.newCox$vaccine)))+0.4, label = "Vaccine")+
    ggplot2::coord_fixed(ratio = (tau+15)/(ymax+2))
  
  
  # Figure B: VE_CI based on KM estimator and new Cox model
  # KM analysis
  VE_CI_km = VE_CI_km 
  VE_CI_km$VE = VE_CI_km$VE*100
  VE_CI_km = rbind(VE_CI_km, c(tau, tail(VE_CI_km$VE,n=1)))
  
  # new Cox model
  S = c(rep(1, nrow(res$VE_CI_1)),
        rep(2, nrow(res$VE_CI_2)),
        rep(3, nrow(res$VE_CI_3)),
        rep(4, nrow(res$VE_CI_4)))
  
  VE_CI_newCox = as.data.frame(rbind(res$VE_CI_1, res$VE_CI_2, res$VE_CI_3, res$VE_CI_4)) 
  VE_CI_newCox[,2] = 100*VE_CI_newCox[,2]
  VE_CI_newCox$S = factor(S, levels = c(0:4))
  colnames(VE_CI_newCox) = c("time", "VE_CI", "S")
  
  ymax = max(c(max(VE_CI_km$VE), max(VE_CI_newCox$VE_CI)))
  ymax = min(ymax - ymax %% 10L + 20.0, 100.0)
  
  ymin = min(c(min(VE_CI_km$VE), min(VE_CI_newCox$VE_CI)))
  ymin = max(ymin - ymin %% 10L - 10.0, 0.0)
  
  fig.B = ggplot2::ggplot()+
    ggplot2::geom_step(data = VE_CI_km, ggplot2::aes(x=time, y=VE), size = curvesize)+
    ggplot2::geom_step(data = VE_CI_newCox, ggplot2::aes(x=time, y=VE_CI, group = S, color = S), size = curvesize)+
    ggplot2::scale_x_continuous(name="Days Since Dose 1", limits=c(28, tau+15), 
                       breaks = c(28, seq(50, tau+15, 50)), expand = c(0,0))+
    ggplot2::scale_y_continuous(name="Vaccine Efficacy on Cumulative Incidence (%)", 
                       limits = c(ymin, ymax+5), breaks = seq(ymin, ymax+5, 10), expand=c(0,ex))+
    ggplot2::scale_color_manual(name = NULL, 
                       values = c("black", "red", "green", "blue", "orange"),
                       labels = c("Kaplan-Meier estimate", paste0("Vaccinated on Day ", s.vec)),
                       drop = FALSE)+
    ggplot2::theme(axis.line.x = ggplot2::element_line(arrow=ggplot2::arrow(length = ggplot2::unit(3, 'mm'))), 
          axis.line.y=ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(3, 'mm'))),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t=6)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r=4)),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
          legend.key = ggplot2::element_blank(),
          legend.position = c(0.32, 0.2),
          legend.box = "horizontal",
          legend.margin = ggplot2::margin(t=-5),
          legend.spacing = ggplot2::unit(0.15, "cm"),
          plot.margin = ggplot2::margin(l = marlr/2, r = marlr/2))+
    ggplot2::guides(color = ggplot2::guide_legend(order=1,ncol=1)) +
    ggplot2::coord_fixed(ratio = (tau+15-28)/(ymax+5-ymin))
  
  
  # Figure C: VE_HR estimators and confidence intervals based on new Cox model
  
  df.VE_HR = data.frame("time" = rep(res$VE_HR$time, 3),
                        "VE" = c(res$VE_HR$VE, res$VE_HR$lower, res$VE_HR$upper),
                        "type" = factor(c(rep(1, nrow(res$VE_HR)), rep(2, 2*nrow(res$VE_HR)))),
                        "group" = factor(rep(1:3, each = nrow(res$VE_HR))))
  df.VE_HR$VE = 100*df.VE_HR$VE
  
  fig.C = ggplot2::ggplot()+
    ggplot2::geom_line(data = df.VE_HR, ggplot2::aes(x=time, y=VE, group = group, color = type), size = curvesize)+
    ggplot2::scale_x_continuous(name="Days Since Dose 1", limits=c(0, tau+15), 
                       breaks = seq(0, tau+15, 50), expand = c(0,0))+
    ggplot2::scale_y_continuous(name="Vaccine Efficacy on Hazard Rate (%)", limits=c(0,111), 
                       breaks = seq(0,100,20), expand=c(0,ex))+
    ggplot2::scale_color_manual(name = NULL, 
                       values = c("black", "green"),
                       labels = c("Estimate", "95% CI"))+
    ggplot2::theme(axis.line.x = ggplot2::element_line(arrow=ggplot2::arrow(length = ggplot2::unit(3, 'mm'))), 
          axis.line.y=ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(3, 'mm'))),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t=6)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r=4)),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
          legend.key = ggplot2::element_blank(),
          legend.position = c(0.2, 0.1),
          legend.box = "horizontal",
          legend.margin = ggplot2::margin(t=-5),
          legend.spacing = ggplot2::unit(0.15, "cm"),
          plot.margin = ggplot2::margin(l = marlr*3/4, r = mar0)) +
    ggplot2::guides(color = ggplot2::guide_legend(order=1,ncol=1)) +
    ggplot2::coord_fixed(ratio = (tau+15)/111)
  
  ggpubr::ggarrange(fig.A, fig.B, fig.C, 
            ncol = 3, nrow = 1,
            labels = LETTERS[1:3],
            hjust = c(-1,-1.5,-2.2),
            vjust = c(4,4,4))
  
  ggplot2::ggsave("figure.pdf", width = 12, height = 5)
  
  return()
}
