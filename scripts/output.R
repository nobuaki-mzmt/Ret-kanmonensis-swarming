## Analysis for R.kanmonensis swarming
## Nobuaki Mizumoto

#------------------------------------------------------------------------------#
{
  library(MASS)
  
  library(car)
  library(lme4)
  library(multcomp)
  
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
d_colony_found <- read.csv("data/raw/colony_foundation.csv",header=T)
d_hybrid_found <- read.csv("data/raw/hybrid_foundation.csv",header=T)

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



r <- glmer(cbind(total,remain)~max+(1|Colony_num), 
           family=binomial(link="logit"), data=d_field[d_field$end ==0,])
summary(r)
Anova(r)


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
Anova(r)



apply(swarm_data>0, 2, sum)

for(i in 1:length(swarm_data[1,])){
  term <- (1:length(swarm_data[,i]))[swarm_data[,i]>0]
  print( max(term) - min(term) + 1 )
}


L <- length(swarm_data[,1])
flight_term <- swarm_data[1:L,]
for(i in 1:length(swarm_data[1,])){
  fd <- flight_term[,i]
  end_flight <- max( (1:length(fd))[fd>0] )
  fd <- fd[1:end_flight]
  print(Idelta(fd))
}

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





d_colony_found
d_dish <- d_colony_found[d_colony_found$case=="dish",]

## main
# 1. analysis for foundation success
found <- tapply(d_dish$foundation, d_dish$treat, mean)
Fig <- barplot(found[c(3,1,2)], ylim=c(0,1), col="#545454", las=1,
               ylab="foundation ratio", xlab="treatment")
r <- glmer(foundation ~ treat + (1|colony), family=binomial(link="logit"), data=d_dish)
Anova(r)
multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
summary(multicomparison)


# 2. brood (offspring, total:include eggs)
d_suc <- d_dish[d_dish$foundation == 1,]
offspring_mean <- tapply(d_suc$offspring, d_suc$treat, mean)[c(3,1,2)]
offspring_sd <- tapply(d_suc$offspring, d_suc$treat, sd)[c(3,1,2)]

Fig <- barplot(offspring_mean, ylim=c(0,15), col="#545454", las=1,
               ylab="number of offspring (worker, larva)", xlab="treatment")
arrows(Fig, offspring_mean, Fig, offspring_mean+offspring_sd, angle=90,
       length=0.1)
r <- glmer(offspring ~ treat + (1|colony), family=poisson, data=d_suc)
Anova(r)

total_mean <- tapply(d_suc$total, d_suc$treat, mean)[c(3,1,2)]
total_sd <- tapply(d_suc$total, d_suc$treat, sd)[c(3,1,2)]

par(pin=c(3,3))
Fig <- barplot(total_mean, ylim=c(0,22), col="#545454", las=1,
               ylab="number of offspring (worker, larva, egg)", xlab="treatment")
arrows(Fig, total_mean, Fig, total_mean+total_sd, angle=90,
       length=0.1)
r <- glmer(total ~ treat + (1|colony), family=poisson, data=d_suc)
Anova(r)

# 3. Hybrid
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








