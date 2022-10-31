##Figures for PNAS manuscript
##Figure 1 - relationship between worm burden and egg counts

library(scales)
require(lamW)

#Round function
round2 = function(x, n=0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

#Estimate k from prevalence and mean worm burden
k_estim <- function(p,m){
  numer = -m*log(1-p)
  denom1 = log(1-p)
  denom2 = (1-p)^(1/m)*log(1-p)/m
  lambert = lambertWm1(denom2)
  k = numer/(denom1-m*lambert)
  return(k)
}

#Power Law function
PL_func <- function(x, y1=53, gamma=0.96) y1*x^gamma #power law function

#Read in datasets
d  <- read.csv("datasets.csv")

#Read in model
m <- readRDS("posteriors/wm-pl-nbh-2p.rds")

#Params from model
stoll <- mean(m$stoll_factor)
y1_TH <- mean(m$y1_TH)
gamma_TH <- mean(m$gamma_TH)
y1_LAO <- mean(m$y1_LAO)
gamma_LAO <- mean(m$gamma_LAO)
h_TH <- mean(m$h_TH)
h_LAO <- mean(m$h_LAO)

#Study indices
TH1.indx <- which(d$survey=="TH1")
TH2.indx <- which(d$survey=="TH2")
TH3.indx <- which(d$survey=="TH3")
TH4.indx <- which(d$survey=="TH4")
LAO1.indx <- which(d$survey=="LAO1")

#Log scale
d$worms_log <- log(d$worms + 1)
d$epg_log <- log(d$epg + 1)

#############################
## Panel A - Original data ##
#############################

#pdf("Fig1A-observed-data.pdf", width=7.5, height=6)

par(mar=c(5, 6.5 ,7.5, 2), xpd=T)

#TH1 - Autopsy Study
plot(d$worms_log[TH1.indx], d$epg_log[TH1.indx], 
     pch=16, col=alpha("grey20", 0.75),
     xlim=log(c(1,5000)), ylim=log(c(1,60000)),
     axes=F, yaxt="n", ylab="",
     xaxt="n", xlab="", cex=1.2)
points(d$worms_log[TH1.indx], d$epg_log[TH1.indx],
       col=alpha("black", 0.75), cex=1.2)
#TH2 - Expulsion study 
points(d$worms_log[TH2.indx], log(d$epg[TH2.indx]),
       pch=16, col=alpha("darkred", 0.75), cex=1.2)
points(d$worms_log[TH2.indx], log(d$epg[TH2.indx]),
       col=alpha("black", 0.75), cex=1.2)
#TH3 - Expulsion study 
points(d$worms_log[TH3.indx], d$epg_log[TH3.indx],
       pch=16, col=alpha("darkgreen", 0.75), cex=1.2)
points(d$worms_log[TH3.indx], d$epg_log[TH3.indx],
       col=alpha("black", 0.75), cex=1.2)
#TH4 - Expulsion study 
points(d$worms_log[TH4.indx], d$epg_log[TH4.indx],
       pch=16, col=alpha("orange", 0.75), cex=1.2)
points(d$worms_log[TH4.indx], d$epg_log[TH4.indx],
       col=alpha("black", 0.75), cex=1.2)
#LAO1 - Expulsion study 
points(d$worms_log[LAO1.indx], d$epg_log[LAO1.indx],
       pch=16, col=alpha("darkblue", 0.75), cex=1.2)
points(d$worms_log[LAO1.indx], d$epg_log[LAO1.indx],
       col=alpha("black", 0.75), cex=1.2)

axis(1, at=log(c(1, 2, 10, 100, 1000, 5000)), 
     labels=c("0", "1" ,"10", "100", "1000", "5000"), 
     cex.axis=1.2)
axis(2, las=2, at=log(c(1, 10, 100, 1000, 10000, 60000)),
     labels=c("0", "10", "100", "1000", "10000", "60000"),
     cex.axis=1.2)
title(ylab = "Eggs per gram of stool", line=4.5, cex.lab=1.3)  
title(xlab = "Worm burden recovered", line=2.5, cex.lab=1.3)
legend("topright", 
       legend=c("TH1 (Autopsy)",
                "TH2 (Expulsion)",
                "TH3 (Expulsion)",
                "TH4 (Expulsion)",
                "LAO1 (Expulsion)"),
       pch=c(rep(16,5), rep(1,5)), inset=c(0.05,-0.4), xpd=T, horiz=F,
       pt.cex=1.5, cex=1.2,
       col=alpha(c("grey20","darkred","darkgreen","orange","darkblue"),0.9),
       bty="n")

legend("topright",
       legend=c("",
                "",
                "",
                "",
                ""),
       inset=c(0.311,-0.4), xpd=T, horiz=F,
       pt.cex=1.5, cex=1.2,
       pch=1,
       col=alpha("black", 0.75),
       bty="n")
#dev.off()

######################################
## Panel 2 - inferred worms burdens ##
######################################

autopsy <- d[which(d$study_type=="autopsy"),]
expulsion <- d[which(d$study_type=="expulsion"),]
expulsion$worms_inferred <- expulsion$worms + (m$expulsion_state-1)
expulsion$worms_inferred_log <- log(expulsion$worms + (m$expulsion_state-1)+1)
autopsy$epg_inferred_log <- log(autopsy$epg/m$autopsy_state+1)

all_worms <- c(autopsy$worms, expulsion$worms_inferred)
summary(all_worms)

#update survey indices
TH2_idx <- which(expulsion$survey=="TH2")
TH3_idx <- which(expulsion$survey=="TH3")
TH4_idx <- which(expulsion$survey=="TH4")
LAO1_idx <- which(expulsion$survey=="LAO1")

#pdf("Fig1B-inferred-relationship.pdf", width=7.5, height=6)

par(mar=c(5, 6.5 ,7.5, 2), xpd=T)

#TH1 - Autopsy Study
plot(autopsy$worms_log, autopsy$epg_inferred_log, 
     pch=16, col=alpha("grey20", 0.75),
     xlim=log(c(1,5000)), ylim=log(c(1,60000)),
     axes=F, yaxt="n", ylab="",
     xaxt="n", xlab="", cex=1.2)
points(autopsy$worms_log, autopsy$epg_inferred_log,
       col=alpha("black", 0.75), cex=1.2)
#TH2 - Expulsion study 
points(expulsion$worms_inferred_log[TH2_idx], log(d$epg[TH2.indx]/stoll),
       pch=16, col=alpha("darkred", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[TH2_idx], log(d$epg[TH2.indx]/stoll),
       col=alpha("black", 0.75), cex=1.2)
#TH3 - Expulsion study 
points(expulsion$worms_inferred_log[TH3_idx], d$epg_log[TH3.indx],
       pch=16, col=alpha("darkgreen", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[TH3_idx], d$epg_log[TH3.indx],
       col=alpha("black", 0.75), cex=1.2)
#TH4 - Expulsion study 
points(expulsion$worms_inferred_log[TH4_idx], d$epg_log[TH4.indx],
       pch=16, col=alpha("orange", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[TH4_idx], d$epg_log[TH4.indx],
       col=alpha("black", 0.75), cex=1.2)
#LAO1 - Expulsion study 
points(expulsion$worms_inferred_log[LAO1_idx], d$epg_log[LAO1.indx],
       pch=16, col=alpha("darkblue", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[LAO1_idx], d$epg_log[LAO1.indx],
       col=alpha("black", 0.75), cex=1.2)

axis(1, at=log(c(1, 2, 10, 100, 1000, 5000)), 
     labels=c("0", "1" ,"10", "100", "1000", "5000"), 
     cex.axis=1.2)
axis(2, las=2, at=log(c(1, 10, 100, 1000, 10000, 60000)),
     labels=c("0", "10", "100", "1000", "10000", "60000"),
     cex.axis=1.2)
title(ylab = "Eggs per gram of stool", line=4.5, cex.lab=1.3)  
title(xlab = "Worm burden inferred", line=2.5, cex.lab=1.3)

#dev.off()

##############################################
## Panel 3 - Fitted relationship by country ##
##############################################
x_worms <- c(0:10, seq(10, 100, by=10), seq(100, 1000, by=100), seq(1000,5000, by=1000))
x_log <- log(x_worms+1)

y_PL_TH <- sapply(x_worms, FUN=PL_func, y1=y1_TH, gamma=gamma_TH)
y_PL_LAO <- sapply(x_worms, FUN=PL_func, y1=y1_LAO, gamma=gamma_LAO)
                  
y_PL_TH_log <- log(y_PL_TH+1)
y_PL_LAO_log <- log(y_PL_LAO+1)

TH_pred_matrix <- LAO_pred_matrix <- matrix(nrow=length(x_worms), ncol=1000, data=0)

for(i in 1:1000){
  TH_pred_matrix[,i] <- sapply(x_worms, FUN=PL_func, y1=m$y1_TH[i], gamma=m$gamma_TH[i])
  LAO_pred_matrix[,i] <- sapply(x_worms, FUN=PL_func, y1=m$y1_LAO[i], gamma=m$gamma_LAO[i])
}

TH_CrI <- apply(TH_pred_matrix, 1, FUN=quantile, probs=c(0.025, 0.975))
LAO_CrI <- apply(LAO_pred_matrix, 1, FUN=quantile, probs=c(0.025, 0.975))
TH_CrI_log <- log(TH_CrI+1)
LAO_CrI_log <- log(LAO_CrI+1)

#TH1 - Autopsy Study
#pdf("Fig1C-PL-curve.pdf", width=7.5, height=6)

par(mar=c(5, 6.5 ,7.5, 2), xpd=T)

plot(autopsy$worms_log, autopsy$epg_inferred_log, 
     pch=16,
     col=alpha("tomato3", 0.75),
     xlim=log(c(1,5000)), ylim=log(c(1,60000)),
     axes=F, yaxt="n", ylab="",
     xaxt="n", xlab="", cex=1.2)
points(autopsy$worms_log, autopsy$epg_inferred_log,
       col=alpha("black", 0.75), cex=1.2)
#TH2 - Expulsion study 
points(expulsion$worms_inferred_log[TH2_idx], log(d$epg[TH2.indx]/stoll),
       pch=16, col=alpha("tomato3", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[TH2_idx], log(d$epg[TH2.indx]/stoll),
       col=alpha("black", 0.75), cex=1.2)
#TH3 - Expulsion study 
points(expulsion$worms_inferred_log[TH3_idx], d$epg_log[TH3.indx],
       pch=16, col=alpha("tomato3", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[TH3_idx], d$epg_log[TH3.indx],
       col=alpha("black", 0.75), cex=1.2)
#TH4 - Expulsion study 
points(expulsion$worms_inferred_log[TH4_idx], d$epg_log[TH4.indx],
       pch=16, col=alpha("tomato3", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[TH4_idx], d$epg_log[TH4.indx],
       col=alpha("black", 0.75), cex=1.2)
#LAO1 - Expulsion study 
points(expulsion$worms_inferred_log[LAO1_idx], d$epg_log[LAO1.indx],
       pch=16, col=alpha("darkblue", 0.75), cex=1.2)
points(expulsion$worms_inferred_log[LAO1_idx], d$epg_log[LAO1.indx],
       col=alpha("black", 0.75), cex=1.2)

axis(1, at=log(c(1, 2, 10, 100, 1000, 5000)), 
     labels=c("0", "1" ,"10", "100", "1000", "5000"), 
     cex.axis=1.2)
axis(2, las=2, at=log(c(1, 10, 100, 1000, 10000, 60000)),
     labels=c("0", "10", "100", "1000", "10000", "60000"),
     cex.axis=1.2)
title(ylab = "Eggs per gram of stool", line=4.5, cex.lab=1.3)  
title(xlab = "Worm burden inferred", line=2.5, cex.lab=1.3)

legend("topright", 
       legend=c("Thailand (TH1-4)",
                "Laos (LAO1)"),
       pch=c(rep(16,5), rep(1,5)), inset=c(0.05,-0.36), xpd=T, horiz=F,
       pt.cex=1.5, cex=1.2,
       col=alpha(c("tomato3","darkblue"),0.9),
       bty="n")

legend("topright",
       legend=c("",
                ""),
       inset=c(0.312,-0.36), xpd=T, horiz=F,
       pt.cex=1.5, cex=1.2,
       pch=1,
       col=alpha("black", 0.75),
       bty="n")

lines(y_PL_TH_log~x_log, lty=2, col="tomato3", lwd=4)
lines(y_PL_LAO_log~x_log, lty=2, col="darkblue", lwd=4)
rethinking::shade(TH_CrI_log, x_log, col=alpha("tomato3", 0.4))
rethinking::shade(LAO_CrI_log, x_log, col=alpha("darkblue", 0.4))
#dev.off()

####################################
## Panel 4 - Prediction intervals ##
####################################

#Validation data from other studies
val <- data.frame(
  survey=c("TH5", "LAO2", "LAO3"),
  worms=c(2588, 11.1, 2.6),
  epg=c(32712, 768, 31),
  prev=c(0.999, 0.384, 0.205))

#H. taichui egg output
h_taichui <- data.frame(
  worms = c(20.3, 1.6, 9.5),
  eggs = c(45.25, 8.33, 18.75)
)

Ht_epg_per_worm <- mean(h_taichui$worms/h_taichui$eggs)
Ht_epg_2008 <- PL_func(351, y1=Ht_epg_per_worm, gamma=0.80)
Ht_epg_2011 <- PL_func(265, y1=Ht_epg_per_worm, gamma=0.80)

Ht_Ov_epg_2008_corrected <- 768-Ht_epg_2008
ov_epg_correction_factor <- Ht_Ov_epg_2008_corrected/768

ov_worms_corrected <- val$worms[2:3]*(1/mean(m$pr_recovery))

val$worms_corrected <- c(2588, ov_worms_corrected)
val$epg_corrected <- c(32712, Ht_Ov_epg_2008_corrected, 31*ov_epg_correction_factor)
val$k <- k_estim(p=val$prev, m=val$worms_corrected)

PL_TH_pred <- sapply(y_PL_TH, function(x) qnbinom(p=c(0.05, 0.95), mu=x, size= h_TH))
PL_LAO_pred <- sapply(y_PL_LAO, function(x) qnbinom(p=c(0.05, 0.95), mu=x, size= h_LAO))
PL_TH_pred_log <- log(PL_TH_pred + 1)
PL_LAO_pred_log <- log(PL_LAO_pred + 1)

val_cols <- c("goldenrod2", "darkolivegreen1",  "darkolivegreen")

#Plot
#pdf("Fig1D-pred-interval-90.pdf", width=7.5, height=6)

par(mar=c(5, 6.5 ,7.5, 2), xpd=T)
#TH1 - Autopsy Study
plot(autopsy$worms_log, autopsy$epg_inferred_log, 
     pch=16,
     col=alpha("tomato3", 0.5),
     xlim=log(c(1,5000)), ylim=log(c(1,60000)),
     axes=F, yaxt="n", ylab="",
     xaxt="n", xlab="", cex=1.2)
points(autopsy$worms_log, autopsy$epg_inferred_log,
       col=alpha("black", 0.5), cex=1.2)
#TH2 - Expulsion study 
points(expulsion$worms_inferred_log[TH2_idx], log(d$epg[TH2.indx]/stoll),
       pch=16, col=alpha("tomato3", 0.5), cex=1.2)
points(expulsion$worms_inferred_log[TH2_idx], log(d$epg[TH2.indx]/stoll),
       col=alpha("black", 0.5), cex=1.2)
#TH3 - Expulsion study 
points(expulsion$worms_inferred_log[TH3_idx], d$epg_log[TH3.indx],
       pch=16, col=alpha("tomato3", 0.5), cex=1.2)
points(expulsion$worms_inferred_log[TH3_idx], d$epg_log[TH3.indx],
       col=alpha("black", 0.5), cex=1.2)
#TH4 - Expulsion study 
points(expulsion$worms_inferred_log[TH4_idx], d$epg_log[TH4.indx],
       pch=16, col=alpha("tomato3", 0.5), cex=1.2)
points(expulsion$worms_inferred_log[TH4_idx], d$epg_log[TH4.indx],
       col=alpha("black", 0.5), cex=1.2)
#LAO1 - Expulsion study 
points(expulsion$worms_inferred_log[LAO1_idx], d$epg_log[LAO1.indx],
       pch=16, col=alpha("darkblue", 0.4), cex=1.2)
points(expulsion$worms_inferred_log[LAO1_idx], d$epg_log[LAO1.indx],
       col=alpha("black", 0.4), cex=1.2)

axis(1, at=log(c(1, 2, 10, 100, 1000, 5000)), 
     labels=c("0", "1" ,"10", "100", "1000", "5000"), 
     cex.axis=1.2)
axis(2, las=2, at=log(c(1, 10, 100, 1000, 10000, 60000)),
     labels=c("0", "10", "100", "1000", "10000", "60000"),
     cex.axis=1.2)
title(ylab = "Eggs per gram of stool", line=4.5, cex.lab=1.3)  
title(xlab = "Worm burden inferred", line=2.5, cex.lab=1.3)

lines(y_PL_TH_log~x_log, lty=2, col="tomato3", lwd=4)
lines(y_PL_LAO_log~x_log, lty=2, col="darkblue", lwd=4)
rethinking::shade(PL_TH_pred_log, x_log, col=alpha("tomato3", 0.3))
rethinking::shade(PL_LAO_pred_log, x_log, col=alpha("darkblue", 0.3))

points(log(val$epg_corrected+1)~log(val$worms_corrected+1), pch=15, cex=2.2, 
       col= val_cols)
points(log(val$epg_corrected+1)~log(val$worms_corrected+1), pch=0, cex=2.2,
       lwd=1.8)

legend("topright", 
       legend=c("TH5 (Autopsy)",
                "LAO2 (Expulsion)",
                "LAO3 (Expulsion)"),
       pch=rep(15,3), inset=c(0.15,-0.36), xpd=T, horiz=F,
       pt.cex=2.2, cex=1.2,
       col=val_cols,
       bty="n")

legend("topright",
       legend=c("",
                "",
                ""),
       inset=c(0.412,-0.36), xpd=T, horiz=F,
       pt.cex=2.2, cex=1.2,
       pch=0,
       bty="n")

#dev.off()
