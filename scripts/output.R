## Analysis for R.kanmonensis swarming
## Nobuaki Mizumoto

#------------------------------------------------------------------------------#
{
  library(MASS)
  library(car)
  library(lme4)
  library(fitdistrplus)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
Idelta <- function(x){
  return(length(x)*sum(x*(x-1)) / ( sum(x) * (sum(x)-1) ))
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
d_lab <- read.csv("data/raw/lab_swarm.csv",header=T)
d_field <- read.csv("data/raw/semifield_swarm.csv",header=T)

d_lab_5 <- d_lab[d_lab$temp==5,]
d_lab_20 <- d_lab[d_lab$temp==20,]

library(ggplot2)
ggplot(d20, aes(x = day, y=total, col=ID))  +
  #scale_y_log10() +
  geom_path()

ggplot(d_field, aes(x = day, y=alates, col=as.factor(colony_id)))  +
  geom_path() +
  scale_y_log10()

# i-delta to examine how the occurence is concentrated
L <- length(swarm_data[,1])
flight_term <- swarm_data[1:L,]
for(i in 1:length(swarm_data[1,])){
  fd <- flight_term[,i]
  end_flight <- max( (1:length(fd))[fd>0] )
  fd <- fd[1:end_flight]
  print(Idelta(fd))
}

# Exponential: time-invarient occurence
# Weibull: shape parameter > 1 the probability increase according to time
par(mfrow=c(3,2), pin=c(3,1))
for(i in 1:length(swarm_data[1,])){
  fd <- flight_term[,i]
  end_flight <- max( (1:length(fd))[fd>0] )
  fd <- fd[1:end_flight]
  
  L2 <- length(fd)
  fit_data <- rep(1:L2, flight_term[1:L2,i])
  id <- colnames(flight_term)[i]
  matplot(cuml_swarm_data[1:L,i]/max(cuml_swarm_data[,i]), type="o", pch=19,
          xlab="day", ylab="Cumulative proportion", main=id, axes=F)
  axis(1)
  axis(2, at=c(0,0.5,1), las=1)
  box()
  
  #weibull_fit <- fitdistr(fit_data, densfun="weibull")
  weibull_fit <- fitdist(fit_data, distr="weibull", method="mle")
  r <- summary(weibull_fit)
  #print(r$aic)
  #print(r$estimate[2])
  points(1:L2, pweibull(1:L2, weibull_fit$estimate[1], weibull_fit$estimate[2]), type="l", col=2 )
  
  #exp_fit <- fitdistr(fit_data, densfun="exponential")
  exp_fit <- fitdist(fit_data, distr="exp", method="mle")
  r <- summary(exp_fit)
  print(r$aic)
  points(1:L2, pexp(1:L2, exp_fit$estimate[1]), type="l", col=4 )
  
}


## Analysis for R.kanmonensis swarming (field)
## 170817 N.Mizumoto

## packages
library(rEDM)
library(car)
library(lme4)

## functions

## main

## data
#d <- read.delim("clipboard",header=T)
#dt <- read.delim("clipboard",header=T)
d # data for swarm
dt # data for temperature

swarm_data <- tapply(d$total, d[,c("Exp_day","Colony_num")], sum) # make table
sum_swarm_data <- apply(swarm_data,1,sum, na.rm=T) # colony pool
num_swarm_col <- apply(swarm_data>0, 1, sum, na.rm=T)
cuml_swarm_data <- tapply(d$Culm, d[,c("Exp_day","Colony_num")], sum) # make table

dt <- data.frame(ave = tapply(d$ave, d$Exp_day, mean),
                 max = tapply(d$max, d$Exp_day, mean),
                 min = tapply(d$min, d$Exp_day, mean))


#########################
## plot for the first time 

# each colony
par(mfrow=c(4,2), pin=c(3,1))
for(i in 1:length(swarm_data[1,])){
  matplot(swarm_data[,i], type="o", pch=19, axes=F, ann=F)
  axis(1)
  axis(2, las=1)
  par(new=T)
  matplot(dt[,1:3], type="l", pch=19, col=c(3,2,4), axes=F, ann=F, lty=2)
  axis(4, las=1)
  box()
  mtext("day", side=1, padj=3)
  mtext("Num. of Individuals", side=2, padj=-5)
  mtext("temperature", side=4, padj=3)
  mtext(paste("colony",i), side=3, padj=-1)
  legend(0,45, c("max","ave","min"), col=c(2,3,4), lty=2, cex=0.8)
}

# colony pool
par(mfrow=c(1,1),pin=c(5,2))
plot(1:40, sum_swarm_data, type="o", pch=19, axes=F, ann=F)
axis(1)
axis(2, las=1)
par(new=T)
matplot(dt[,1:3], type="l", pch=19, col=c(3,2,4), axes=F, ann=F, lty=2)
axis(4, las=1)
box()
mtext("day", side=1, padj=3)
mtext("Num. of Individuals", side=2, padj=-5)
mtext("temperature", side=4, padj=3)
legend(0,45, c("max","ave","min"), col=c(2,3,4), lty=2, cex=0.8)


# num flight col
par(mfrow=c(1,1),pin=c(5,2))
plot(1:40, num_swarm_col, type="o", pch=19, axes=F, ann=F)
axis(1)
axis(2, las=1)
par(new=T)
matplot(dt[,1:3], type="l", pch=19, col=c(3,2,4), axes=F, ann=F, lty=2)
axis(4, las=1)
box()
mtext("day", side=1, padj=3)
mtext("Num. of flight colony", side=2, padj=-5)
mtext("temperature", side=4, padj=3)
legend(0,45, c("max","ave","min"), col=c(2,3,4), lty=2, cex=0.8)


#############
## GLMM 

d <- read.table("clipboard",header=T)
head(d)

r <- glmer(cbind(total,remain)~max+(1|Colony_num), family=binomial(link="logit"), data=d[d$end ==0,])
summary(r)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: cbind(total, remain) ~ max + (1 | Colony_num)
# AIC      BIC   logLik deviance df.resid 
# 23715.0  23724.5 -11854.5  23709.0      173 

# Scaled residuals: 
#  Min     1Q Median     3Q    Max 
# -42.95  -2.09  -0.37  -0.04 381.44 

# Random effects:
#  Groups     Name        Variance Std.Dev.
# Colony_num (Intercept) 6.699    2.588   
# Number of obs: 176, groups:  Colony_num, 7

# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -19.377967   0.997053  -19.44   <2e-16 ***
#  max           0.598848   0.006189   96.76   <2e-16 ***

Anova(r)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: cbind(total, remain)
# Chisq Df Pr(>Chisq)    
# max 9363.1  1  < 2.2e-16 ***

par(pin=c(4,3))

touka <- 0.3
d$Col[d$Colony_num == 7] <- rgb(1, 0, 0, alpha=touka)
d$Col[d$Colony_num == 2] <- rgb(0, 1, 0, alpha=touka)
d$Col[d$Colony_num == 3] <- rgb(0, 0, 1, alpha=touka)
d$Col[d$Colony_num == 4] <- rgb(1, 1, 0, alpha=touka)
d$Col[d$Colony_num == 1] <- rgb(1, 0, 1, alpha=touka)
d$Col[d$Colony_num == 6] <- rgb(0, 1, 1, alpha=touka)
d$Col[d$Colony_num == 5] <- rgb(0, 0, 0, alpha=touka)

MinMax <- min(d[d$end==0,]$max)
MaxMax <- max(d[d$end==0,]$max)
plot(at_flight~max, data=d[d$end ==0,], col=Col, bg=Col, pch=21, cex=log10(total+remain),
     xlab="maximum temperature (degree)", ylab="ratio of flight individual", las=1,
     ylim=c(-0.05, 1.05), xlim=c(8,38))

logistic <- function(x){ 1/(1+exp(-x)) }#???ʂ̃??W?X?e?B?b?N?Ȑ???????
Est <- summary(r)$coefficient[,1]#?W???̗\???l???????o??
xdata <- seq(MinMax, MaxMax,length.out=120) #x?͈̔͂??o??
ydata <- logistic(Est[1]+Est[2]*xdata) #y?̒l???v?Z
lines(xdata, ydata)#?????㏑??


r <- glmer(cbind(total,remain)~ave+(1|Colony_num), family=binomial(link="logit"), data=d[d$end ==0,])
summary(r)
# Fixed effects:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -12.158487   0.672050  -18.09   <2e-16 ***
#  ave           0.971652   0.008791  110.52   <2e-16 ***

Anova(r)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: cbind(total, remain)
# Chisq Df Pr(>Chisq)    
# ave 12216  1  < 2.2e-16 ***



####################
## ???????ȕ??z?W???x?w???ŏW?????ČQ?򂪋N???Ă????Ƃ??��Ă݂?

apply(swarm_data>0, 2, sum)

for(i in 1:length(swarm_data[1,])){
  term <- (1:length(swarm_data[,i]))[swarm_data[,i]>0]
  print( max(term) - min(term) + 1 )
}


# I?w??
# ?????????A?e?R???j?[?őS?Ă̌̂????яI?????܂ł̊??ԂŁA???P?ʂ?I?w?????v?Z
L <- length(swarm_data[,1])
flight_term <- swarm_data[1:L,]
for(i in 1:length(swarm_data[1,])){
  fd <- flight_term[,i]
  end_flight <- max( (1:length(fd))[fd>0] )
  fd <- fd[1:end_flight]
  print(Idelta(fd))
}

# ???L?̊?ł̌Q???\???ɂ??????A???Č̐??̕??z??fitting
# Exponential: ???ۂ̐??N?m????????
# Weibull: shape parameter > 1 ?̎??A???ۂ̐??????m???????ԂƋ??ɑ???
#   ???܂??̊Ԃœ????悤?ȃ^?C?~???O??flight???N?????B
par(mfrow=c(4,2), pin=c(3,1))
for(i in 1:length(swarm_data[1,])){
  fd <- flight_term[,i]
  end_flight <- max( (1:length(fd))[fd>0] )
  fd <- fd[1:end_flight]
  L2 <- length(fd)
  fit_data <- rep(1:L2, flight_term[1:L2,i])
  id <- colnames(flight_term)[i]
  matplot(cuml_swarm_data[1:L,i]/max(cuml_swarm_data[,i]), type="o", pch=19,
          xlab="day", ylab="Cumulative proportion", main=id, axes=F)
  axis(1)
  axis(2, at=c(0,0.5,1), las=1)
  box()
  
  if(i==1||i==4){print(NA);next;}
  
  #weibull_fit <- fitdistr(fit_data, densfun="weibull")
  weibull_fit <- fitdist(fit_data, distr="weibull", method="mle")
  r <- summary(weibull_fit)
  #print(r$aic)
  #print(r$estimate[1])
  #print(r$estimate[2])
  points(1:L2, pweibull(1:L2, weibull_fit$estimate[1], weibull_fit$estimate[2]), type="l", col=2 )
  
  #exp_fit <- fitdistr(fit_data, densfun="exponential")
  exp_fit <- fitdist(fit_data, distr="exp", method="mle")
  r <- summary(exp_fit)
  print(r$aic)
  points(1:L2, pexp(1:L2, exp_fit$estimate[1]), type="l", col=4 )
  
}





###############
## try to CCM analysis
ccm_data <- data.frame(num_col = num_swarm_col, sum_swarm = sum_swarm_data)
ccm_data <- cbind(ccm_data, dt[,3:5])[1:36,]

## make them only analyzing data
varlst <- colnames(ccm_data)

# Embedding Dimension by simplex projection:
par(mfrow = c(3,3), pin=c(1,1))
simplex_output_list <- NULL

lib = pred <- c(1, length(ccm_data[,1]))
for (i in 1:length(varlst)) {
  simplex_output_list[[i]] <- simplex(ccm_data[,i], lib, pred, tau = -1, E=1:6)
  plot(simplex_output_list[[i]]$E, simplex_output_list[[i]]$rho, type = "l", 
       xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)", main = varlst[i])
}
names(simplex_output_list) <- c(varlst)
bestE <- sapply(simplex_output_list, function(simplex_output) {
  simplex_output$E[which.max(simplex_output$rho)]
})
bestE

# Using these embedding dimensions, we can now apply S-maps to identify nonlinearity:
smap_output_list <- NULL
par(mfrow = c(3, 3))
for (i in 1:length(varlst)) {
  if(sum(ccm_data[,i])==0){ plot(0,main = paste(varlst[i],"is not validate")); next; }
  smap_output_list[[i]] <- s_map(ccm_data[,i], lib, pred, E = as.numeric(bestE[i]))
  plot(smap_output_list[[i]]$theta, smap_output_list[[i]]$rho, type = "l", 
       xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)", main = varlst[i])
}

names(smap_output_list) <- c(varlst)
nonlinear <- sapply(smap_output_list, function(smap_output) {
  smap_output$theta[which.max(smap_output$rho)]
})
nonlinear


## Cross-map matrix
ncol <- dim(ccm_data)[2]
M_rho <- array(NA, dim = c(ncol, ncol), 
               dimnames = list(varlst, varlst))

for (i in 1:ncol) {
  for (j in 1:ncol) {
    if (i != j) {
      if(nonlinear[j]>0){
        out <- ccm(ccm_data, E = as.numeric(bestE[j]), lib_column = varlst[i], target_column = varlst[j],
                   lib_sizes = dim(ccm_data)[1], replace = FALSE, silent = TRUE)
        M_rho[i, j] <- round(out$rho,7)
      } else {
        M_rho[i, j] <- NA
      }
    }
  }
}

M_rho


M_corr <- array(NA, dim = c(ncol, ncol), dimnames = list(varlst, varlst))
for (i in 1:ncol) {
  for (j in 1:ncol) {
    if (i != j) {
      cf_temp <- ccf(x = ccm_data[, i], y = ccm_data[, j]
                     , type = "correlation", lag.max = 6, plot = FALSE)$acf
      M_corr[i, j] <- max(abs(cf_temp))
    }
  }
}

M_corr


## Convergent Cross-Mapping
#?S?Ăɂ???

par(mfrow=c(5,5))
count <- 1
for(j in 1:length(varlst))
  for(i in 1:length(varlst)){
    if(i>j){
      lib_list <- varlst[j]
      tar_list <- varlst[i]
      
      if( nonlinear[which(varlst == tar_list)] == 0 ){ nonL1 <- "(linear)"} else {nonL1 <- ""}
      if( nonlinear[which(varlst == lib_list)] == 0 ){ nonL2 <- "(linear)"} else {nonL2 <- ""}
      lib_xmap_tar <- ccm(ccm_data, E = as.numeric(bestE[which(varlst == tar_list)]), random_libs = F, lib_column = lib_list, 
                          target_column = tar_list, lib_sizes = seq(5, 35, by = 5),
                          num_samples = 35)
      tar_xmap_lib <- ccm(ccm_data, E = as.numeric(bestE[which(varlst == lib_list)]), random_libs = F, lib_column = tar_list, 
                          target_column = lib_list, lib_sizes = seq(5, 35, by = 5),
                          num_samples = 35)
      lib_xmap_tar_means <- ccm_means(lib_xmap_tar)
      tar_xmap_lib_means <- ccm_means(tar_xmap_lib)
      plot(rho ~ lib_size, type="l", col=2, xlab="Library size", ylab="Cross Map Skill (rho)",
           xlim=c(5,35), ylim=c(0,1), data = lib_xmap_tar_means)
      lines(tar_xmap_lib_means$lib_size, tar_xmap_lib_means$rho, col = "blue")
      lib_name <- gsub("cells", "", lib_list)
      tar_name <- gsub("cells", "", tar_list)
      
      if(1-max(c(lib_xmap_tar_means$rho,tar_xmap_lib_means$rho), na.rm = TRUE) < min(c(lib_xmap_tar_means$rho,tar_xmap_lib_means$rho), na.rm = TRUE)){
        legend(x = "bottomright", legend = c(paste(lib_name,"-",tar_name,nonL1), paste(tar_name,"-",lib_name,nonL2)),
               col = c("red","blue"), lwd = 1, inset = 0.02, cex = 0.8)
      } else {
        legend(x = "topleft", legend = c(paste(lib_name,"-",tar_name,nonL1), paste(tar_name,"-",lib_name,nonL2)),
               col = c("red","blue"), lwd = 1, inset = 0.02, cex = 0.8)
      }
      
      #abline(h=M_corr[lib_list, tar_list], lty=2)
      
      
      if(count%%9 == 0){
        # dev.off()
      } 
      count <- count + 1
      
    }
  }

## surrogate test

num_surr <- 100
par(mfrow=c(3,3))
surrogate_list <- NULL
for (i in 1:length(varlst)) {
  surrogate_list[[i]] <- make_surrogate_data(ccm_data[,varlst[i]], method = "random_shuffle", num_surr = num_surr)
}
names(surrogate_list) <- c(varlst)

rho_surr <- NULL

for(j in 1:length(varlst)){
  for(i in 1:length(varlst)){
    if(i!=j){
      lib_list <- varlst[j]
      tar_list <- varlst[i]
      
      ## ?ʏ??ʂ?CCM
      if( nonlinear[which(varlst == tar_list)] == 0 ){ next; }
      
      lib_xmap_tar <- ccm(ccm_data, E = as.numeric(bestE[which(varlst == tar_list)]), random_libs = TRUE, lib_column = lib_list, 
                          target_column = tar_list, lib_sizes = seq(4, 36, by = 4),
                          num_samples = 300)
      lib_xmap_tar_means <- ccm_means(lib_xmap_tarna.rm=T)
      lib_name <- lib_list
      tar_name <- tar_list
      
      rho_surr <- matrix(0, nrow = num_surr, ncol=length(seq(5,35,5)))
      ## surrogate test
      for(h in 1:num_surr){
        surr_ccm <- ccm(cbind(ccm_data[,tar_list], surrogate_list[[lib_list]][,h]),
                        E = as.numeric(bestE[which(varlst == tar_list)]), lib_column = 2, target_column = 1, 
                        lib_sizes = seq(5, 35, by = 5), random_libs = TRUE,  num_samples = 300,
                        replace = T)
        rho_surr[h,] <- ccm_means(surr_ccm)$rho
      }
      surr_res <- apply(rho_surr, 2, function(vec){ quantile(vec,p=0.95,na.rm = T) } )
      surr_mean <- apply(rho_surr, 2, mean )
      
      xs <- ccm_means(surr_ccm)$lib_size
      
      plot(0, type="n", xlim=c(5,35), ylim=c(-0.1,1), main = paste("from", lib_name,"to",tar_name),
           xlab="Library size", ylab="Cross Map Skill (rho)")
      polygon( c(xs,rev(xs)), c(rep(-0.1,length(xs)),rev(surr_res)), col="#ffd1d1", border = NA)
      lines(xs, surr_mean, col = 1, lty=2, lwd=3)
      lines(rho ~ lib_size, col=2, data = lib_xmap_tar_means)
      
      
      if(count%%9 == 0){
        # dev.off()
      } 
      count <- count + 1
    }
    
  }
}


## Analysis for foundation (Rk)
## 170905 N.Mizumoto 

## packages
library(car)
library(lme4)
library(multcomp)

## functions
se  <-  function(x){
  y  <-  x[!is.na(x)]  #  remove  the  missing  values
  sqrt(var(as.vector(y))/length(y))
}


## data
d <- read.delim("clipboard",header=T)
d_dish <- d[d$case=="dish",] # glass cell ?̃f?[?^?����O


## main
# 1. analysis for foundation success
found <- tapply(d_dish$foundation, d_dish$treat, mean)
par(pin=c(3,3))
Fig <- barplot(found[c(3,1,2)], ylim=c(0,1), col="#545454", las=1,
               ylab="foundation ratio", xlab="treatment")

r <- glmer(foundation ~ treat + (1|colony), family=binomial(link="logit"), data=d_dish)
Anova(r)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: foundation
# Chisq Df Pr(>Chisq)  
# treat 6.3168  2    0.04249 *

multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
summary(multicomparison)
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = foundation ~ treat + (1 | colony), data = d_dish, 
#            family = binomial(link = "logit"))
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# nat_flight - lab_flight == 0  -2.1525     0.8584  -2.508   0.0312 *
# wood - lab_flight == 0        -0.4037     1.2199  -0.331   0.9399  
# wood - nat_flight == 0         1.7488     1.3453   1.300   0.3872  


# 2. brood (offspring, total:include eggs)
# ?n?ݐ????????y?A?̃f?[?^?̂ݎg?p
d_suc <- d_dish[d_dish$foundation == 1,]
head(d_suc)
offspring_mean <- tapply(d_suc$offspring, d_suc$treat, mean)[c(3,1,2)]
offspring_sd <- tapply(d_suc$offspring, d_suc$treat, sd)[c(3,1,2)]

par(pin=c(3,3))
Fig <- barplot(offspring_mean, ylim=c(0,15), col="#545454", las=1,
               ylab="number of offspring (worker, larva)", xlab="treatment")
arrows(Fig, offspring_mean, Fig, offspring_mean+offspring_sd, angle=90,
       length=0.1)
r <- glmer(offspring ~ treat + (1|colony), family=poisson, data=d_suc)
Anova(r)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: offspring
# Chisq Df Pr(>Chisq)
# treat 2.3857  2     0.3033

total_mean <- tapply(d_suc$total, d_suc$treat, mean)[c(3,1,2)]
total_sd <- tapply(d_suc$total, d_suc$treat, sd)[c(3,1,2)]

par(pin=c(3,3))
Fig <- barplot(total_mean, ylim=c(0,22), col="#545454", las=1,
               ylab="number of offspring (worker, larva, egg)", xlab="treatment")
arrows(Fig, total_mean, Fig, total_mean+total_sd, angle=90,
       length=0.1)
r <- glmer(total ~ treat + (1|colony), family=poisson, data=d_suc)
Anova(r)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: total
# Chisq Df Pr(>Chisq)
# treat 3.8025  2     0.1494

# 3. Hybrid
d <- read.delim("clipboard",header=T)

found <- tapply(d$foundation, d$unit, mean)[c(2,3,1,4)]
par(pin=c(3,3))
Fig <- barplot(found, ylim=c(0,1), col="#545454", las=1,
               ylab="foundation ratio", xlab="treatment")

d2 <- d[d$foundation==1,]
total_mean <- tapply(d2$total, d2$unit, mean)[c(2,3,1,4)]
total_sd <- tapply(d2$total, d2$unit, sd)[c(2,3,1,4)]
par(pin=c(3,3))
Fig <- barplot(total_mean, ylim=c(0,25), col="#545454", las=1,
               ylab="foundation ratio", xlab="treatment")
arrows(Fig, total_mean, Fig, total_mean+total_sd, angle=90,
       length=0.1)








