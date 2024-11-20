

###---###---###---###---###---###---###---###---###---###---###---###---###---##
#
# CORRECT SCENARIO
#
###---###---###---###---###---###---###---###---###---###---###---###---###---##

load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24_n100.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Miss-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Miss-Treatment_11.15.24.RData")

#           DR-oc     DR-all     IPW      OR-all    OR-oc
grey  <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442" )
#             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
lines <- c("dotted","dotted",  "dotdash",   "solid",    "solid")


data_plot <- data_plot %>% filter(method=="Gcomp-OC-WE" |
                                    method=="Gcomp-All" |
                                    method=="IPW-OC-WE" |
                                    method=="AIPW-OC-WE" |
                                    method=="AIPW-All")


efficiency_gains <- c(ratio_gcomp_oc_w_vs_all,
                      ratio_gcomp_oc_we_vs_all,
                      ratio_aipw_oc_w_vs_all,
                      ratio_aipw_oc_we_vs_all)


method <- c(rep("Gcomp-OC-W/Gcomp-All",length(seq_th)),
            rep("Gcomp-OC-WE/Gcomp-All",length(seq_th)),
            rep("AIPW-OC-W/AIPW-All",length(seq_th)),
            rep("AIPW-OC-WE/AIPW-All",length(seq_th)))

threshold <- rep(seq_th,4)

data_plot_ratio <- data.frame(efficiency_gains,method,threshold)

data_plot_ratio <- data_plot_ratio %>% filter(  method=="Gcomp-OC-WE/Gcomp-All" |
                                                  method=="AIPW-OC-WE/AIPW-All" )

######


data_plot$method[data_plot$method == "AIPW-OC-WE"] <- "DR-oc"
data_plot$method[data_plot$method == "AIPW-All"] <- "DR-ac"
data_plot$method[data_plot$method == "Gcomp-OC-WE"] <- "OR-oc"
data_plot$method[data_plot$method == "Gcomp-All"] <- "OR-ac"
data_plot$method[data_plot$method == "IPW-OC-WE"] <- "IPW"


data_plot_ratio$method[data_plot_ratio$method=="Gcomp-OC-WE/Gcomp-All"] <- 
  "OR-oc/OR-ac"
data_plot_ratio$method[data_plot_ratio$method=="AIPW-OC-WE/AIPW-All"] <- 
  "DR-oc/DR-ac"


SIZET <- 16
SIZETi <- 16
LENGTHLe <- 2
YLIM <- 0.6

bias2 <- ggplot(data_plot, aes(x = threshold, y = (bias^2), col = method)) + 
  #geom_point(size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Bias^2")+
  ylim(0,YLIM) + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")




var <- ggplot(data_plot, aes(x = threshold, y = (se_sv^2), col = method) )+ 
  #geom_point(size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Var - Sampling variability") +
  ylim(0,YLIM) + 
  theme_bw() + 
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")


mse <- ggplot(data_plot, aes(x = threshold, y = mse, col = method) )+ 
  #geom_point(size = 4) +
  geom_line(aes(linetype=method), size = 2) +  
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("MSE")+
  ylim(0,YLIM) + 
  theme_bw() + 
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")



cove <- ggplot(data_plot, aes(x = threshold, y = coverage, col = method) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Coverage of the 95% CI")+
  ylim(0.9,1) + 
  theme_bw() + 
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")  + # scale_colour_grey() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") 




noprint <- ggplot(data_plot, aes(x = threshold, y = coverage, col = method) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Coverage of the 95% CI")+
  ylim(0,1) + 
  theme_bw() + 
  theme(legend.position="bottom", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"),
        legend.text = element_text(size=SIZET))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")  + # scale_colour_grey() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") 




#########################################################################################################


# Extract the legend. Returns a gtable
leg_b <- get_legend(noprint)

grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c( 0.5, #size of rows
                                                                   5,
                                                                   5, 
                                                                   0.5), "null"))))  

grid.text("Correct models",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=14))

print(bias2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

print(var, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# grid.text("Continuous Treatment", 
#           vp = viewport(layout.pos.row = 1, layout.pos.col = 2),gp=gpar(fontsize=14))

print(mse, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))

print(cove, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

print(as_ggplot(leg_b), vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2),gp=gpar(fontsize=20))





###---###---###---###---###---###---###---###---###---###---###---###---###---##
#
# MISSPECIFIED OUTCOME MODEL
#
###---###---###---###---###---###---###---###---###---###---###---###---###---##

# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24.RData")
load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24_n100.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Miss-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Miss-Treatment_11.15.24.RData")

#           DR-oc     DR-all     IPW      OR-all    OR-oc
grey  <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442" )
#             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
lines <- c("dotted","dotted",  "dotdash",   "solid",    "solid")

YLIM <- 5

data_plot <- data_plot %>% filter(method=="Gcomp-OC-WE" |
                                    method=="Gcomp-All" |
                                    method=="IPW-OC-WE" |
                                    method=="AIPW-OC-WE" |
                                    method=="AIPW-All")


efficiency_gains <- c(ratio_gcomp_oc_w_vs_all,
                      ratio_gcomp_oc_we_vs_all,
                      ratio_aipw_oc_w_vs_all,
                      ratio_aipw_oc_we_vs_all)


method <- c(rep("Gcomp-OC-W/Gcomp-All",length(seq_th)),
            rep("Gcomp-OC-WE/Gcomp-All",length(seq_th)),
            rep("AIPW-OC-W/AIPW-All",length(seq_th)),
            rep("AIPW-OC-WE/AIPW-All",length(seq_th)))

threshold <- rep(seq_th,4)

data_plot_ratio <- data.frame(efficiency_gains,method,threshold)

data_plot_ratio <- data_plot_ratio %>% filter(  method=="Gcomp-OC-WE/Gcomp-All" |
                                                  method=="AIPW-OC-WE/AIPW-All" )

######


data_plot$method[data_plot$method == "AIPW-OC-WE"] <- "DR-oc"
data_plot$method[data_plot$method == "AIPW-All"] <- "DR-ac"
data_plot$method[data_plot$method == "Gcomp-OC-WE"] <- "OR-oc"
data_plot$method[data_plot$method == "Gcomp-All"] <- "OR-ac"
data_plot$method[data_plot$method == "IPW-OC-WE"] <- "IPW"

data_plot_ratio$method[data_plot_ratio$method=="Gcomp-OC-WE/Gcomp-All"] <- 
  "OR-oc/OR-ac"
data_plot_ratio$method[data_plot_ratio$method=="AIPW-OC-WE/AIPW-All"] <- 
  "DR-oc/DR-ac"




SIZET <- 16
SIZETi <- 16
LENGTHLe <- 1

bias2 <- ggplot(data_plot, aes(x = threshold, y = (bias^2), col = method)) + 
  #geom_point(size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Bias^2")+
  ylim(0,YLIM) + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")




var <- ggplot(data_plot, aes(x = threshold, y = (se_sv^2), col = method) )+ 
  #geom_point(size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Var - Sampling variability") +
  ylim(0,1) + 
  theme_bw() + 
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")


mse <- ggplot(data_plot, aes(x = threshold, y = mse, col = method) )+ 
  #geom_point(size = 4) +
  geom_line(aes(linetype=method), size = 2) +  
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("MSE")+
  ylim(0,YLIM) + 
  theme_bw() + 
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")



cove <- ggplot(data_plot, aes(x = threshold, y = coverage, col = method) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Coverage of the 95% CI")+
  ylim(0,1) + 
  theme_bw() + 
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")  + # scale_colour_grey() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") 



######


ratio <- ggplot(data_plot_ratio, aes(x = threshold, y = efficiency_gains) )+ 
  geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(colour = factor(method)), size = 1)  + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Ratios of SEs")+
  # ylim(0.95,1.25) + 
  geom_hline(yintercept=1, linetype="dashed", color = "red") + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(1,"cm"))  +
  labs(color='Method')  #+ scale_colour_grey() 




noprint <- ggplot(data_plot, aes(x = threshold, y = coverage, col = method) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(linetype=method), size = 2) + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Coverage of the 95% CI")+
  ylim(0,1) + 
  theme_bw() + 
  theme(legend.position="bottom", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(LENGTHLe,"cm"),
        legend.text = element_text(size=SIZET))  +
  #labs(color='Method') +
  scale_color_manual(values=grey) +
  scale_linetype_manual(values = lines) +
  scale_fill_manual(values=grey,
                    name="Method")  + # scale_colour_grey() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") 





#########################################################################################################


# Extract the legend. Returns a gtable
leg_b <- get_legend(noprint)

grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c( 0.5, #size of rows
                                                                   5,
                                                                   5, 
                                                                   0.5), "null"))))  

grid.text("Misspecified outcome - Correct treatment",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=14))

print(bias2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

print(var, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# grid.text("Continuous Treatment", 
#           vp = viewport(layout.pos.row = 1, layout.pos.col = 2),gp=gpar(fontsize=14))

print(mse, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))

print(cove, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

print(as_ggplot(leg_b), vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2),gp=gpar(fontsize=20))







###############################################################################
###############################################################################
###############################################################################
#
#
# RATIO
#
#
################################################################################
###############################################################################
###############################################################################


SIZET <- 16
SIZETi <- 16
LENGTHLe <- 2
YLIM <- 0.04

load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24_n100.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Miss-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Miss-Treatment_11.15.24.RData")

######

efficiency_gains <- c(ratio_gcomp_oc_w_vs_all,
                      ratio_gcomp_oc_we_vs_all,
                      ratio_aipw_oc_w_vs_all,
                      ratio_aipw_oc_we_vs_all)


method <- c(rep("Gcomp-OC-W/Gcomp-All",length(seq_th)),
            rep("Gcomp-OC-WE/Gcomp-All",length(seq_th)),
            rep("AIPW-OC-W/AIPW-All",length(seq_th)),
            rep("AIPW-OC-WE/AIPW-All",length(seq_th)))

threshold <- rep(seq_th,4)

data_plot_ratio <- data.frame(efficiency_gains,method,threshold)

data_plot_ratio <- data_plot_ratio %>% filter(  method=="Gcomp-OC-WE/Gcomp-All" |
                                                  method=="AIPW-OC-WE/AIPW-All" )

######


data_plot_ratio$method[data_plot_ratio$method=="Gcomp-OC-WE/Gcomp-All"] <- 
  "OR-oc/OR-ac"
data_plot_ratio$method[data_plot_ratio$method=="AIPW-OC-WE/AIPW-All"] <- 
  "DR-oc/DR-ac"


ratiocc <- ggplot(data_plot_ratio, aes(x = threshold, y = efficiency_gains) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(colour = factor(method)), size = 2)  + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Ratios of SEs")+
  ylim(0.8,1.35) + 
  geom_hline(yintercept=1, linetype="dashed", color = "red") + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(1,"cm"))  +
  labs(color='Method')  #+ scale_colour_grey() 

####----------------------------------------------------------------------------

# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24.RData")
load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24_n100.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Miss-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Miss-Treatment_11.15.24.RData")

######

efficiency_gains <- c(ratio_gcomp_oc_w_vs_all,
                      ratio_gcomp_oc_we_vs_all,
                      ratio_aipw_oc_w_vs_all,
                      ratio_aipw_oc_we_vs_all)


method <- c(rep("Gcomp-OC-W/Gcomp-All",length(seq_th)),
            rep("Gcomp-OC-WE/Gcomp-All",length(seq_th)),
            rep("AIPW-OC-W/AIPW-All",length(seq_th)),
            rep("AIPW-OC-WE/AIPW-All",length(seq_th)))

threshold <- rep(seq_th,4)

data_plot_ratio <- data.frame(efficiency_gains,method,threshold)

data_plot_ratio <- data_plot_ratio %>% filter(  method=="Gcomp-OC-WE/Gcomp-All" |
                                                  method=="AIPW-OC-WE/AIPW-All" )

######


data_plot_ratio$method[data_plot_ratio$method=="Gcomp-OC-WE/Gcomp-All"] <- 
  "OR-oc/OR-ac"
data_plot_ratio$method[data_plot_ratio$method=="AIPW-OC-WE/AIPW-All"] <- 
  "DR-oc/DR-ac"


ratioMout_Ctre <- ggplot(data_plot_ratio, aes(x = threshold, y = efficiency_gains) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(colour = factor(method)), size = 2)  + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Ratios of SEs")+
  ylim(0.8,1.35) + 
  geom_hline(yintercept=1, linetype="dashed", color = "red") + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(1,"cm"))  +
  labs(color='Method')  #+ scale_colour_grey()    



####----------------------------------------------------------------------------

# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24.RData")
load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24 - Stochastic_n100.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Miss-Treatment_11.15.24.RData")

######

efficiency_gains <- c(ratio_gcomp_oc_w_vs_all,
                      ratio_gcomp_oc_we_vs_all,
                      ratio_aipw_oc_w_vs_all,
                      ratio_aipw_oc_we_vs_all_S)


method <- c(rep("Gcomp-OC-W/Gcomp-All",length(seq_th)),
            rep("Gcomp-OC-WE/Gcomp-All",length(seq_th)),
            rep("AIPW-OC-W/AIPW-All",length(seq_th)),
            rep("AIPW-OC-WE/AIPW-All",length(seq_th)))

threshold <- rep(seq_th,4)

data_plot_ratio <- data.frame(efficiency_gains,method,threshold)

data_plot_ratio <- data_plot_ratio %>% filter(  method=="Gcomp-OC-WE/Gcomp-All" |
                                                  method=="AIPW-OC-WE/AIPW-All" )

######


data_plot_ratio$method[data_plot_ratio$method=="Gcomp-OC-WE/Gcomp-All"] <- 
  "OR-oc/OR-ac"
data_plot_ratio$method[data_plot_ratio$method=="AIPW-OC-WE/AIPW-All"] <- 
  "DR-oc/DR-ac"


ratiocc_Stoch <- ggplot(data_plot_ratio, aes(x = threshold, y = efficiency_gains) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(colour = factor(method)), size = 2)  + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Ratios of SEs")+
  ylim(0.8,1.35) + 
  geom_hline(yintercept=1, linetype="dashed", color = "red") + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(1,"cm"))  +
  labs(color='Method')  #+ scale_colour_grey()   


####----------------------------------------------------------------------------

# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Correct-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24.RData")
# load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Correct-Outcome Miss-Treatment_11.15.24.RData")
load("/Users/santam13/Documents/NYU/Non concurrent/Ident-Estima-Parametric/rcode/non-concurrent/data_simu/Miss-Outcome Correct-Treatment_11.15.24 - Stochastic_n100.RData")

######

efficiency_gains <- c(ratio_gcomp_oc_w_vs_all,
                      ratio_gcomp_oc_we_vs_all,
                      ratio_aipw_oc_w_vs_all,
                      ratio_aipw_oc_we_vs_all_S)


method <- c(rep("Gcomp-OC-W/Gcomp-All",length(seq_th)),
            rep("Gcomp-OC-WE/Gcomp-All",length(seq_th)),
            rep("AIPW-OC-W/AIPW-All",length(seq_th)),
            rep("AIPW-OC-WE/AIPW-All",length(seq_th)))

threshold <- rep(seq_th,4)

data_plot_ratio <- data.frame(efficiency_gains,method,threshold)

data_plot_ratio <- data_plot_ratio %>% filter(  method=="Gcomp-OC-WE/Gcomp-All" |
                                                  method=="AIPW-OC-WE/AIPW-All" )

######


data_plot_ratio$method[data_plot_ratio$method=="Gcomp-OC-WE/Gcomp-All"] <- 
  "OR-oc/OR-ac"
data_plot_ratio$method[data_plot_ratio$method=="AIPW-OC-WE/AIPW-All"] <- 
  "DR-oc/DR-ac"


ratioMout_Ctre_Stoch <- ggplot(data_plot_ratio, aes(x = threshold, y = efficiency_gains) )+ 
  #geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(colour = factor(method)), size = 2)  + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Ratios of SEs")+
  ylim(0.8,1.35) +  
  geom_hline(yintercept=1, linetype="dashed", color = "red") + 
  theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(1,"cm"))  +
  labs(color='Method')  #+ scale_colour_grey()   




noprint <- ggplot(data_plot_ratio, aes(x = threshold, y = efficiency_gains) )+ 
  geom_point(aes(colour = factor(method)), size = 4) +
  geom_line(aes(colour = factor(method)), size = 2)  + 
  scale_x_continuous(breaks=seq_th,
                     labels = paste(100*(1-seq_th),"%",sep="")) +
  xlab("% of concurrent controls") +
  ylab("Ratios of SEs")+
  # ylim(0.95,1.25) + 
  geom_hline(yintercept=1, linetype="dashed", color = "red") + 
  theme_bw() +
  theme(legend.position="bottom", 
        axis.title = element_text(size = SIZET),
        axis.text = element_text(size = SIZET),
        plot.title = element_text(hjust = 0.5, size = SIZETi),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size=SIZET))  +
  labs(color='Method')  #+ scale_colour_grey() 

#########################################################################################################


# Extract the legend. Returns a gtable
leg_b <- get_legend(noprint)

grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 2, heights = unit(c( 0.5, #size of rows
                                                                   5,
                                                                   0.5,
                                                                   5,
                                                                   0.5), "null"))))

grid.text("Correct models",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=14))

print(ratiocc, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

grid.text("Misspecified outcome",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 2),gp=gpar(fontsize=14))

print(ratioMout_Ctre, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

grid.text("Correct models - Vk stochastic",
          vp = viewport(layout.pos.row = 3, layout.pos.col = 1),gp=gpar(fontsize=14))

print(ratiocc_Stoch, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))

grid.text("Misspecified outcome - Vk stochastic",
          vp = viewport(layout.pos.row = 3, layout.pos.col = 2),gp=gpar(fontsize=14))

print(ratioMout_Ctre_Stoch, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))

print(as_ggplot(leg_b), vp = viewport(layout.pos.row = 5, layout.pos.col = 1:2),gp=gpar(fontsize=20))




