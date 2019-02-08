############################################
# genpwr examples for manuscript
############################################

# load libraries
library(ggplot2)
library(gridExtra)
library(genpwr)

# Calculate power to detect and OR of 1.5 for a GWAS with N=5,000 and N=20,000
# Assume 1 control per case and 5 x 10^-8 significance level
pw <- genpwr.calc(calc = "pow", model = "logistic", ge.interaction = NULL,
  N=c(5000, 20000), Case.Rate = 0.5, OR=1.5, Alpha=0.00000005, MAF=seq(0.01, 0.8, 0.005))
x <- "RAF"


pw$Power <- pw$`Power_at_Alpha_5e-08`
pw$Test.Model <- factor(pw$Test.Model, levels = c('Dominant', "Additive", "Recessive", "2df"))
pw$True.Model <- factor(pw$True.Model, levels = c('Dominant', "Additive", "Recessive", "2df"))
pw$RAF <- pw$MAF
panel.by <- "True.Model"

plot_obj_5000 <- ggplot(data = pw[pw$N_total==5000,], aes(x = pw[pw$N_total==5000, x], 
                                          y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw[pw$N_total==5000,][, panel.by])


plot_obj_20000 <- ggplot(data = pw[pw$N_total==20000,], aes(x = pw[pw$N_total==20000, x], 
                                           y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw[pw$N_total==20000,][, panel.by])

plot_obj_5000 + ggtitle("N = 5,000") 

plot_obj_20000 + ggtitle("N = 20,000")


pdf('genpwr_plots_OR_1.25_alpha_genome.pdf', height = 3, width=9)

# plot_obj_2000 + ggtitle("N=2,000")
# plot_obj_3000 + ggtitle("N=3,000")
# plot_obj_4000 + ggtitle("N=4,000")
plot_obj_5000 + ggtitle("N=5,000")

# plot_obj_10000 + ggtitle("N=10,000")
# plot_obj_10000 + xlim(c(0,0.15))+ ggtitle("N=10,000, Zoom 1")
# plot_obj_10000 + xlim(c(0,0.15)) + ylim(c(0.7,1))+ ggtitle("N=10,000, Zoom 2")

# plot_obj_15000+ ggtitle("N=15,000")
# plot_obj_15000 + xlim(c(0,0.1))+ ggtitle("N=15,000, Zoom 1")
# plot_obj_15000 + xlim(c(0,0.1))+ ylim(c(0.7,1))+ ggtitle("N=15,000, Zoom 2")

plot_obj_20000+ ggtitle("N=20,000")
plot_obj_20000 + xlim(c(0,0.1))+ ggtitle("N=20,000, Zoom 1")
plot_obj_20000 + xlim(c(0,0.1))+ ylim(c(0.7,1))+ ggtitle("N=20,000, Zoom 2")
dev.off()

# Calculate detectable odds ratios
or <- odds_ratio_function(N=c(5000, 20000), Case.Rate = 0.5, power=0.8, Alpha=0.00000005, MAF=c( 0.1, 0.3) )


##############################################################################
# Power for G x E Interaction in IPF
##############################################################################
# Parameters
Ncase <- 3624
Ncontrol <- 4442
Ntotal <- Ncase+Ncontrol
cr <- Ncase/Ntotal
p_e <- p_ever <- (0.72*Ncase + 0.63*Ncontrol)/Ntotal

################################################
# Calculate power for a range of scenarios:
# OR_G = 1.25, 1.5, 5
# OR_GE = 1.5, 1.75, 2
################################################
# OR_G = 5
pw_5 <- genpwr.calc(calc = "pow", model = "logistic", ge.interaction = "binary", 
                           N = Ntotal, Case.Rate = cr, MAF = seq(0.05, 0.8, 0.05),
                           OR_G = c(5), OR_E = 1.6, OR_GE = c(1.5,1.75,2), P_e = p_e,
                           Alpha = 0.00005, True.Model = "All", Test.Model = "All")
pw_5$RAF <- pw_5$MAF
pw_5$Power <- pw_5$`Power_at_Alpha_5e.05`
pw_5$True.Model <- factor(pw_5$True.Model, levels = c('Dominant', 'Additive', 'Recessive'))
pw_5 <-pw_5[order(pw_5$True.Model),]
pw_5$Test.Model <- factor(pw_5$Test.Model, levels = c('Dominant', 'Additive', 'Recessive', '2df'))

panel.by <- "True.Model"
x <- 'RAF'

plot_obj_5_1.5 <- ggplot(data = pw_5[pw_5$OR_GE==1.5,], aes(x = pw_5[pw_5$OR_GE==1.5,x], 
                                             y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_5[pw_5$OR_GE==1.5,][, panel.by])

plot_obj_5_1.75 <- ggplot(data = pw_5[pw_5$OR_GE==1.75,], aes(x = pw_5[pw_5$OR_GE==1.75,x], 
                                              y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_5[pw_5$OR_GE==1.75,][, panel.by])

plot_obj_5_2 <- ggplot(data = pw_5[pw_5$OR_GE==2,], aes(x =pw_5[pw_5$OR_GE==2,x], 
                                           y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_5[pw_5$OR_GE==2,][, panel.by])


# OR_G = 1.5
pw_1.5 <- genpwr.calc(calc = "pow", model = "logistic", ge.interaction = "binary", 
                            N = Ntotal, Case.Rate = cr, MAF = seq(0.05, 0.8, 0.05),
                            OR_G = c(1.5), OR_E = 1.6, OR_GE = c(1.5,1.75,2), P_e = p_ever,
                            Alpha = 0.00005, True.Model = "All", Test.Model = "All")
pw_1.5$RAF <- pw_1.5$MAF
pw_1.5$Power <- pw_1.5$`Power_at_Alpha_5e.05`
pw_1.5$True.Model <- factor(pw_1.5$True.Model, levels = c('Dominant', 'Additive', 'Recessive'))
pw_1.5 <-pw_1.5[order(pw_1.5$True.Model),]
pw_1.5$Test.Model <- factor(pw_1.5$Test.Model, levels = c('Dominant', 'Additive', 'Recessive', '2df'))

plot_obj_1.5_1.5 <- ggplot(data = pw_1.5[pw_1.5$OR_GE==1.5,], aes(x = pw_1.5[pw_1.5$OR_GE==1.5,x], 
                                              y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_1.5[pw_1.5$OR_GE==1.5,][, panel.by])

plot_obj_1.5_1.75 <- ggplot(data = pw_1.5[pw_1.5$OR_GE==1.75,], aes(x = pw_1.5[pw_1.5$OR_GE==1.75,x], 
                                               y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_1.5[pw_1.5$OR_GE==1.75,][, panel.by])

plot_obj_1.5_2 <- ggplot(data = pw_1.5[pw_1.5$OR_GE==2,], aes(x = pw_1.5[pw_1.5$OR_GE==2,x], 
                                            y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_1.5[pw_1.5$OR_GE==2,][, panel.by])


# OR_G = 1.25
pw_1.25 <- genpwr.calc(calc = "pow", model = "logistic", ge.interaction = "binary", 
                           N = Ntotal, Case.Rate = cr, MAF = seq(0.05, 0.8, 0.05),
                           OR_G = c(1.25), OR_E = 1.6, OR_GE = c(1.5,1.75,2), P_e = p_ever,
                           Alpha = 0.00005, True.Model = "All", Test.Model = "All")
pw_1.25$RAF <- pw_1.25$MAF
pw_1.25$Power <- pw_1.25$`Power_at_Alpha_5e.05`
pw_1.25$True.Model <- factor(pw_1.25$True.Model, levels = c('Dominant', 'Additive', 'Recessive'))
pw_1.25$Test.Model <- factor(pw_1.25$Test.Model, levels = c('Dominant', 'Additive', 'Recessive', '2df'))
pw_1.25 <-pw_1.25[order(pw_1.25$True.Model),]

plot_obj_1.25_1.5 <- ggplot(data = pw_1.25[pw_1.25$OR_GE==1.5,], aes(x = pw_1.25[pw_1.25$OR_GE==1.5,x], 
                                               y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_1.25[pw_1.25$OR_GE==1.5,][, panel.by])

plot_obj_1.25_1.75 <- ggplot(data = pw_1.25[pw_1.25$OR_GE==1.75,], aes(x = pw_1.25[pw_1.25$OR_GE==1.75,x], 
                                              y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_1.25[pw_1.25$OR_GE==1.75,][, panel.by])

plot_obj_1.25_2 <- ggplot(data = pw_1.25[pw_1.25$OR_GE==2,], aes(x = pw_1.25[pw_1.25$OR_GE==2,x], 
                                               y = Power, group = Test.Model, colour = Test.Model)) + 
  geom_line() + 
  xlab(x) + geom_hline(aes(yintercept=0.8), linetype="dashed",color='black')+ facet_wrap(~pw_1.25[pw_1.25$OR_GE==2,][, panel.by])

#################################
# Make Plots
#################################
# Plots for ORge = 1.5
p1<-plot_obj_1.25_1.5 + ggtitle(expression('OR'[g]*" = 1.25")) + ylim(c(0,1))
p2<-plot_obj_1.5_1.5 + ggtitle(expression('OR'[g]*" = 1.5")) + ylim(c(0,1))
p3<-plot_obj_5_1.5 + ggtitle(expression('OR'[g]*" = 5")) + ylim(c(0,1))
grid.arrange(p1, p2, p3, nrow = 3)

# Plots for ORge = 1.75
p1<-plot_obj_1.25_1.75 + ggtitle(expression('OR'[g]*" = 1.25")) + ylim(c(0,1))
p2<-plot_obj_1.5_1.75 + ggtitle(expression('OR'[g]*" = 1.5")) + ylim(c(0,1))
p3<-plot_obj_5_1.75 + ggtitle(expression('OR'[g]*" = 5")) + ylim(c(0,1))
grid.arrange(p1, p2, p3, nrow = 3)

# Plots for ORge = 2
p1<-plot_obj_1.25_2 + ggtitle(expression('OR'[g]*" = 1.25")) + ylim(c(0,1))
p2<-plot_obj_1.5_2 + ggtitle(expression('OR'[g]*" = 1.5")) + ylim(c(0,1))
p3<-plot_obj_5_2 + ggtitle(expression('OR'[g]*" = 5")) + ylim(c(0,1))
grid.arrange(p1, p2, p3, nrow = 3)

