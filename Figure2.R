##Figures for PNAS manuscript
##Figure 2 - k and diagnostic sensitivity
require(scales)
source("funcs.R")

#Read in model output
m <- readRDS("posteriors/wm-pl-nbh-2p.rds")

############################## 
## Panel A - estimates of k ##
##############################

k_mu <- mean(m$k_mu)
k_sd <- mean(m$k_sd)
k_mat <- rbind(m$k1, m$k2, m$k3, m$k4, m$k5)
k_upper = apply(k_mat, 1, quantile, probs=0.025)
k_lower = apply(k_mat, 1, quantile, probs=0.975)
k_pred <- qnorm(p=c(0.05, 0.95), mean=k_mu, sd=k_sd)

k.df <- data.frame(
  mean = c(apply(k_mat, 1, mean), k_mu),
  study = factor(c("TH1", "TH2", "TH3", "TH4", "LAO1", "mean"),
  levels = c("TH1", "TH2", "TH3", "TH4", "LAO1", "mean")),
  x = c(1:5, 6.5),
  lower = c(k_lower, k_pred[1]),
  upper = c(k_upper, k_pred[2]))

val.df <- data.frame(
  study = factor(c("LAO2", "LAO3"), levels=c("LAO2", "LAO3")),
  mean = c(0.08511544, 0.04751100),
  x= c(7, 8))

col_vec <- c("grey20", "darkred", "darkgreen", "orange", "darkblue", "black")
val_cols <- c("darkolivegreen1",  "darkolivegreen4")

#pdf("Fig2A-k-surveys.pdf", height=5, width=7)
par(mar=c(2,5,1,1), xpd=F)

plot(k.df$mean~k.df$x, ylim=c(0,1), axes=F,
     cex.lab=1.4, xlab="", xlim=c(1,9), type="n",
     ylab=expression(paste("Worm burden dispersion ", italic("k"))))

arrows(angle=90, x0=k.df$x, x1=k.df$x, y0=k.df$mean, y1=k.df$upper, length = 0.08,
       col=alpha(col_vec, 0.75), lwd=2.4)

arrows(angle=90, x0=k.df$x, x1=k.df$x, y1=k.df$lower, y0=k.df$mean, length = 0.08,
       col=alpha(col_vec, 0.75), lwd=2.4)

points(k.df$mean~k.df$x, cex=1.7, pch=16, col=alpha(col_vec, 0.8))

points(k.df$mean~k.df$x, cex=1.7, lwd=1.7)

axis(2, las=2, cex.axis=1.3, lwd=1.2)

text(c("TH1", "TH2", "TH3", "TH4"), 
     x=c(k.df$x[1:4])+0.45, 
     y=c(k.df$mean[1:4]),
     family="sans", cex=1.2)

text(c("LAO1", "LAO2", "LAO3"), 
     x=c(k.df$x[5], val.df$x)+0.6, 
     y=c(k.df$mean[5], val.df$mean)+0.012,
     family="sans", cex=1.2)  

text("Prediction \ninterval",
     x=7.5, y=0.5,
     family="sans", cex=1.3) 

points(val.df$mean~val.df$x, col=val_cols, pch=15, cex=1.9)
points(val.df$mean~val.df$x, pch=0, cex=1.9, lwd=1.7)

#dev.off()

############################################
## Panel B - relationship between k and M ##
############################################

val <- data.frame(
  survey=c("TH5", "LAO2", "LAO3"),
  worms=c(2588, 11.1, 2.6),
  epg=c(32712, 768, 31),
  prev=c(0.999, 0.384, 0.205))

k.all <- data.frame(
  k = c(k.df$mean[1:5], val.df$mean),
  M = c(mean(m$M1), mean(m$M2), mean(m$M3), mean(m$M4), mean(m$M5), 25.1, 5.9)
)

k_cols <- c(alpha(col_vec[1:5],0.65), val_cols)

## Plot
#pdf("Fig2B-k-vs-M.pdf", height=5, width=7)
par(mar=c(5,5,1,1), xpd=F)
plot(k.all$k~ log(k.all$M+1), ylim=c(0,0.6), 
     xlab="", ylab="", axes=F, cex=1.9, xlim=log(c(0, 500)+1),
     col=k_cols, pch=c(rep(16, 5), rep(15, 3)))
points(k.all$k~ log(k.all$M+1), pch=c(rep(1, 5), rep(0, 3)), cex=1.9, lwd=1.65)
axis(1, cex.axis=1.3, lwd=1.2, at=log(c(0, 1, 5, 20, 100, 500)+1),
     labels=c(0, 1, 5, 20, 100, 500))
axis(2, las=2, cex.axis=1.3, lwd=1.2)
title(ylab=expression(paste("Worm burden dispersion ", italic("k"))), cex.lab=1.4, line=3.4)
title(xlab="Population mean worm burden", line=3.0, cex.lab=1.4)
#dev.off()

###################################### 
## Panel C - individual sensitivity ##
######################################

b <- mean(m$b)
b_lower <- quantile(m$b, probs=0.025)
b_upper <- quantile(m$b, probs=0.975)
worms <- c(0, 1, 5, 20, 100, 500)
sens_worms <- nbh_sens(worms, b=b)
sens_interval <- t(cbind(nbh_sens(worms, b=b_upper), nbh_sens(worms, b=b_lower)))

#pdf("Fig2C-indiv-sens.pdf",height=5, width=7)
par(mar=c(5,5,1,1), xpd=F)
plot(log(worms+1), sens_worms, xlab="Individual worm burden", ylab="", 
     ylim=c(0,1), axes=F, cex=1.5, cex.lab=1.4,
     col="grey10", pch=16)
title(ylab="Probability of detecting fecal eggs", line=3.6, cex.lab=1.4)
axis(1, at = log(c(0, 1, 5, 20, 100, 500)+1), 
     labels = c(0, 1, 5, 20, 100, 500),
     cex.axis=1.3, lwd=1.2)
axis(2, las=2, cex.axis=1.3, lwd=1.2)
shade(sens_interval, log(worms+1), col=alpha("grey30", 0.35))
#dev.off()

###############################################
## Panel D - pop sensitivity and worm burden ##
###############################################

mu <- c(0.0001, seq(0.001, 2, by=0.01), seq(2, 50, by=0.5), seq(50, 100, by=5), seq(100,500,by=20)) #Population mean worm burden
pop_sens_kmu <- sapply(mu, FUN=pop_sensitivity, k=k_mu, b=b) #Population sensitivity
pop_sens_k10 <- sapply(mu, FUN=pop_sensitivity, k=0.1, b=b) #Population sensitivity
pop_sens_k70 <- sapply(mu, FUN=pop_sensitivity, k=0.7, b=b) #Population sensitivity

leg_text <- c(expression(paste(italic("k")," = 0.10")), 
                     expression(paste(italic("k")," = 0.36")), 
                     expression(paste(italic("k")," = 0.70")))

#pdf("Fig2D-pop-sens-worm.pdf",height=5, width=7)
par(mar=c(5,5,1,1), xpd=F)
plot(log(mu+1), pop_sens_kmu, type="l", xlab="", lwd=2, col="grey4",
     ylab="", axes=F, ylim=c(0.2, 1), xlim=log(c(1,501)))
lines(log(mu+1), pop_sens_k10, lwd=2, col="blue3")
lines(log(mu+1), pop_sens_k70, lwd=2, col="red3")
axis(1, cex.axis=1.25, at=log(c(501, 101, 21, 6, 2, 1)),
     labels=c(500, 100, 20, 5, 1, 0), lwd=1.3)
axis(2, las=2, cex.axis=1.25, lwd=1.3)
title(xlab="Population mean worm burden", cex.lab=1.4)
title(ylab="Diagnostic sensitivity", line=3.8, cex.lab=1.4)

legend(title="Worm burden dispersion", legend=leg_text, 
       x = log(1), y=1.04, lty=c(1,1,1), lwd=2.5, 
       col=c("blue3","grey5","red3"), cex=1.25, bty = "n")
#dev.off()

###########################################
## Panel E - True vs observed prevalence ##
###########################################

true_prev_kmu <- prev(x=mu, k = k_mu)
true_prev_k10 <- prev(x=mu, k = 0.1)
true_prev_k70 <- prev(x=mu, k = 0.7)

observed_prev_kmu <- obs_prev_func(prev=true_prev_kmu, se = pop_sens_kmu, sp = 1)
observed_prev_k10 <- obs_prev_func(prev=true_prev_k10, se = pop_sens_k10, sp = 1)
observed_prev_k70 <- obs_prev_func(prev=true_prev_k70, se = pop_sens_k70, sp = 1)

#pdf("Fig2E-prev-true-vs-obs.pdf", height=5, width=7)
par(mar=c(5,5,1,1), xpd=F)
plot(observed_prev_kmu, true_prev_kmu , type="l", ylim=c(0,0.49), xlim=c(0,0.5), 
     xlab="", ylab="", axes=F, lwd=2.2, col="grey5")
axis(1, cex.axis=1.3, lwd=1.3, at = c(0, 0.1, 0.2,0.3,0.4,0.5), labels=c(0,10,20,30,40,50))
axis(2, las=2, cex.axis=1.3, lwd=1.3, at = c(0, 0.1, 0.2,0.3,0.4,0.5), labels=c(0,10,20,30,40,50))
lines(observed_prev_k10, true_prev_k10, col="blue3", lwd=2.2)
lines(observed_prev_k70, true_prev_k70, col="red3", lwd=2.2)
abline(0, 1, lty=5, lwd=2)
title(xlab="Observed prevalence from fecal eggs (%)", cex.lab=1.4)
title(ylab="True prevalence (%)", line=3.7, cex.lab=1.4)

lines(x=c(-0.06, 0.2), y=rep(true_prev_k10[median(which(round(observed_prev_k10, digits = 2)==0.2))],2),
      lty=2, lwd=2, col="blue3")
lines(x=c(-0.06, 0.2), y=rep(true_prev_kmu[median(which(round(observed_prev_kmu, digits = 2)==0.2))],2),
      lty=2, lwd=2, col="grey5")
lines(x=c(-0.06, 0.2), y=rep(true_prev_k70[median(which(round(observed_prev_k70, digits = 2)==0.2))],2),
      lty=2, lwd=2, col="red3")
lines(x=rep(0.2,2), y=c(-0.06, true_prev_k70[median(which(round(observed_prev_k70, digits = 2)==0.2))]),
      lty=2, lwd=2, col="grey5")
#dev.off()

###########################################
## Panel E - Impact of reduced specificity ##
###########################################

obs_p_kmu_sp100 <- obs_prev_func(prev=true_prev_kmu, se = pop_sens_kmu, sp = 1)
obs_p_kmu_sp95 <- obs_prev_func(prev=true_prev_kmu, se = pop_sens_kmu, sp = 0.95)
obs_p_kmu_sp90 <- obs_prev_func(prev=true_prev_kmu, se = pop_sens_kmu, sp = 0.90)

sp_leg_text <- c(expression(paste(italic("sp")," = 1")), 
              expression(paste(italic("sp")," = 0.95")), 
              expression(paste(italic("sp")," = 0.90")))

#pdf("Fig2F-specificity-prev.pdf", height=5, width=7)
par(mar=c(5,5,1,1), xpd=F)
plot(obs_p_kmu_sp100, true_prev_kmu , type="l", ylim=c(0,0.49), xlim=c(0,0.5), 
     xlab="", ylab="", axes=F, lwd=2.2, col="grey5")
axis(1, cex.axis=1.3, lwd=1.3, at = c(0, 0.1, 0.2,0.3,0.4,0.5), labels=c(0,10,20,30,40,50))
axis(2, las=2, cex.axis=1.3, lwd=1.3, at = c(0, 0.1, 0.2,0.3,0.4,0.5), labels=c(0,10,20,30,40,50))
lines(obs_p_kmu_sp95, true_prev_kmu, col="darkorchid3", lwd=2.2)
lines(obs_p_kmu_sp90, true_prev_kmu, col="goldenrod2", lwd=2.2)
abline(0, 1, lty=5, lwd=2)
title(xlab="Observed prevalence from fecal eggs (%)", cex.lab=1.4)
title(ylab="True prevalence (%)", line=3.7, cex.lab=1.4)

legend(title="Diagnostic specificity", legend=sp_leg_text, 
       x = 0.001, y=0.52, lty=c(1,1,1), lwd=2.5, 
       col=c("grey5","darkorchid3","goldenrod2"), cex=1.25, bty = "n")
#dev.off()
