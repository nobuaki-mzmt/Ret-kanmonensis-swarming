## Analysis for R.kanmonensis swarming
## Nobuaki Mizumoto

#------------------------------------------------------------------------------#
{
  library(MASS)
  
  library(car)
  library(lme4)
  library(multcomp)
  
  library(fitdistrplus)
  
  library(ggplot2)
  library(viridis)
  
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

d_lab_5 <- d_lab[d_lab$temp_treat==5,]
d_lab_20 <- d_lab[d_lab$temp_treat==20,]

## time development
ggplot(d_field, aes(x = day, y=alates))  +
  geom_path() +
  geom_point() +
  facet_grid(colony_id~., scales = "free") +
  theme_classic() +
  theme(strip.background = element_blank(),
    strip.text.y = element_blank()  ) +
  xlab("Day")+
  ylab("Alates")

ggplot(d_field[d_field$colony_id == "col-1",])  +
  geom_path(aes(x = day, y=ave_temp)) +
  geom_path(aes(x = day, y=min_temp)) +
  geom_path(aes(x = day, y=max_temp)) +
  ylim(-5,45) +
  theme_classic()+
  xlab("Day")+
  ylab("Temperature (°C)")

## Logistic fitting max temp vs alate fraction
total_nest_alate = tapply(d_field$cumulative_alate, d_field$colony_id, max)
d_field$alate_remaining_nest = rep(total_nest_alate, each = 40) - d_field$cumulative_alate

r <- glmer(cbind(alates, alate_remaining_nest)~max_temp+(1|colony_id), 
           family=binomial(link="logit"), data=d_field[d_field$end ==0,])
summary(r)
Anova(r)

logistic <- function(x){ 1/(1+exp(-x)) }
Est <- summary(r)$coefficient[,1]
fit_x <- seq(min(d_field[d_field$end==0,]$max_temp),
             max(d_field[d_field$end==0,]$max_temp),length.out=120) 
fit_y <- logistic(Est[1]+Est[2]*xdata)
logistic_fitting <- data.frame(fit_x, fit_y)

ggplot(d_field[d_field$end==0,])+
  geom_point( aes(x=max_temp, y=alate_fraction),
              alpha=1) +
  theme_classic()+
  #scale_color_viridis(discrete=F)+
  xlab("Maximum temperature (°C)")+
  ylab("Fraction of alates swarmed")+
  geom_line(aes(fit_x, fit_y), data=logistic_fitting)

r <- glmer(cbind(alates, alate_remaining_nest)~ave_temp+(1|colony_id), 
           family=binomial(link="logit"), data=d_field[d_field$end ==0,])
summary(r)
Anova(r)

r <- glmer(cbind(alates, alate_remaining_nest)~min_temp+(1|colony_id), 
           family=binomial(link="logit"), data=d_field[d_field$end ==0,])
summary(r)
Anova(r)



## lab time development
ggplot(d_lab_20, aes(x = day, y=alates))  +
  geom_path() +
  geom_point() +
  facet_grid(ID~., scales = "free") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()  ) +
  xlab("Day")+
  ylab("Alates")

ggplot(d_lab_20[d_lab_20$ID == "A1",])  +
  geom_path(aes(x = day, y=ave_temp)) +
  geom_path(aes(x = day, y=min_temp)) +
  geom_path(aes(x = day, y=max_temp)) +
  theme_classic()+
  ylim(-5,45) +
  xlab("Day")+
  ylab("Temperature (°C)")


## largest swarming
df_largest_swarm = data.frame(
  treat = c(rep("field", 7), rep("lab", 6)),
  colony = c( unique(d_field$colony_id), unique(d_lab_20$ID)),
  largest_swarm = c( tapply(d_field$alates, d_field$colony_id, max),
                     tapply(d_lab_20$alates, d_lab_20$ID, max)),
  total_swarm = c( tapply(d_field$cumulative_alate, d_field$colony_id, max),
                   tapply(d_lab_20$cumlative, d_lab_20$ID, max))
)

r <- glmer(cbind(largest_swarm, total_swarm-largest_swarm) ~ 
             treat+(1|colony), 
           family=binomial(link="logit"), data=df_largest_swarm)
summary(r)
Anova(r)

ggplot(df_largest_swarm, aes(x=treat, y = largest_swarm/total_swarm)) +
  geom_boxplot() + 
  geom_point() +
  ylim(c(0,1))









