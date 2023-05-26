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
# output of swarming observations
#------------------------------------------------------------------------------#
swarming_output <- function(){
  
  # read data
  {
    d_lab <- read.csv("data/raw/lab_swarm.csv",header=T)
    d_field <- read.csv("data/raw/semifield_swarm.csv",header=T)
    
    d_lab_5 <- d_lab[d_lab$temp_treat==5,]
    d_lab_20 <- d_lab[d_lab$temp_treat==20,]
  }

  # define swarming events, by obtaining outlines using z-score method
  {
    outlier = NULL
    for (i in unique(d_field$colony_id) ){
      x = subset(d_field, colony_id == i)$alates
      outlier = c(outlier, abs(scale(x)) > 3)
    }
    d_field$outlier = outlier
    
    outlier = NULL
    for (i in unique(d_lab_20$ID) ){
      x = subset(d_lab_20, ID == i)$alates
      outlier = c(outlier, abs(scale(x)) > 3)
    }
    d_lab_20$outlier = outlier
  }
  
  # plot time development
  {
    ## field
    ggplot(d_field, aes(x = day, y=alates))  +
      geom_path() +
      geom_point(data=d_field[d_field$outlier,], col= 2) +
      #geom_point(aes(col= outlier)) +
      facet_grid(colony_id~., scales = "free") +
      theme_classic() +
      theme(strip.background = element_blank(),
        strip.text.y = element_blank()  ) +
      xlab("Day")+
      ylab("Alates")
    ggsave("output/field_swarm_timedevelopment.pdf",
           width=3, height=5)  
    
    temp_cols = viridis(4, option = "C")
    ggplot(d_field[d_field$colony_id == "col-1",])  +
      geom_path(aes(x = day, y=ave_temp), col=temp_cols[3]) +
      geom_path(aes(x = day, y=min_temp), col=temp_cols[2]) +
      geom_path(aes(x = day, y=max_temp), col=temp_cols[1]) +
      ylim(-5,45) +
      theme_classic()+
      xlab("Day")+
      ylab("Temperature (°C)")
    ggsave("output/field_temp.pdf",
           width=3, height=1.5)  
    
    ## lab
    ggplot(d_lab_20, aes(x = day, y=alates))  +
      geom_path() +
      geom_point(data=d_lab_20[d_lab_20$outlier,], col= 2) +
      facet_grid(ID~., scales = "free") +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.text.y = element_blank()  ) +
      xlab("Day")+
      ylab("Alates")
    ggsave("output/lab_swarm_timedevelopment.pdf",
           width=3, height=5)  
    
    temp_cols = viridis(4, option = "C")
    ggplot(d_lab_20[d_lab_20$ID == "A1",])  +
      geom_path(aes(x = day, y=ave_temp), col=temp_cols[3]) +
      geom_path(aes(x = day, y=min_temp), col=temp_cols[2]) +
      geom_path(aes(x = day, y=max_temp), col=temp_cols[1]) +
      theme_classic()+
      ylim(-5,45) +
      xlab("Day")+
      ylab("Temperature (°C)")
    ggsave("output/lab_temp.pdf",
           width=3, height=1.5)  
  }
  
  # Logistic fitting max temp vs alate fraction in the field
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
  fit_y <- logistic(Est[1]+Est[2]*fit_x)
  logistic_fitting <- data.frame(fit_x, fit_y)
  
  ggplot(d_field[d_field$end==0,])+
    geom_point( aes(x=max_temp, y=alate_fraction, col=alates),
                alpha=1) +
    scale_color_viridis(option = "A", end=0.8) +
    theme_classic()+
    #scale_color_viridis(discrete=F)+
    xlab("Maximum temperature (°C)")+
    ylab("Fraction of alates swarmed")+
    geom_line(aes(fit_x, fit_y), data=logistic_fitting)+
    theme(legend.position = "top")
  ggsave("output/field_swarm_fitting_temp.pdf",
         width=3, height=5)  
  
  
  r <- glmer(cbind(alates, alate_remaining_nest)~ave_temp+(1|colony_id), 
             family=binomial(link="logit"), data=d_field[d_field$end ==0,])
  summary(r)
  Anova(r)
  
  r <- glmer(cbind(alates, alate_remaining_nest)~min_temp+(1|colony_id), 
             family=binomial(link="logit"), data=d_field[d_field$end ==0,])
  summary(r)
  Anova(r)
  
  hist(d_field[d_field$flight>0,"max_temp",])
  hist(d_field[d_field$flight>0,"ave_temp",])
  hist(d_field[d_field$flight>0,"min_temp",])
  
  ## lab time development
  
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
}


##
d_colony_found <- read.csv("data/raw/colony_foundation.csv",header=T)
d_hybrid_found <- read.csv("data/raw/hybrid_foundation.csv",header=T)

d_dish <- d_colony_found[d_colony_found$case=="dish",] # removal the results of glass cell

#
ggplot(d_dish, aes(x = treat, fill = as.factor(foundation))) + 
  geom_bar(position = "fill", alpha = 0.8)+
  scale_fill_viridis(discrete = T, direction = -1) +
  theme_classic()+
  theme(legend.position = "none")+
  xlab("")+
  ylab("Colony foundation success")

r <- glmer(foundation ~ treat + (1|colony), family=binomial(link="logit"),
           data=d_dish)
Anova(r)
multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
summary(multicomparison)

d_suc <- d_dish[d_dish$foundation == 1,]
ggplot(data=d_suc, aes(x=treat, y=offspring))+
  geom_boxplot()
r <- glmer(total ~ treat + (1|colony), family=poisson, data=d_suc)
Anova(r)


## 

ggplot(d_hybrid_found, aes(x = unit, fill = as.factor(foundation))) + 
  geom_bar(position = "fill", alpha = 0.8)+
  scale_fill_viridis(discrete = T, direction = -1) +
  theme_classic()+
  theme(legend.position = "none")+
  xlab("")+
  ylab("Colony foundation success")

d_hybrid_found$unit <- as.factor(d_hybrid_found$unit)
r <- glm(foundation ~ unit, family=binomial(link="logit"),
           data=d_hybrid_found)
Anova(r)

d_suc <- d_hybrid_found[d_hybrid_found$foundation == 1,]
ggplot(data=d_suc, aes(x=unit, y=offspring))+
  geom_boxplot()
r <- glm(total ~ unit, family=poisson, data=d_suc)
Anova(r)
multicomparison<-glht(r,linfct=mcp(unit="Tukey"))
summary(multicomparison)

