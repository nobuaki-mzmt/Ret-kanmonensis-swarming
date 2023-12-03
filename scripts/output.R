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
{
  swarming_output()
  foundation_output()
  hybrid_output()
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
           width=3, height=4)  
    
    temp_cols = viridis(4, option = "C")
    ggplot(d_field[d_field$colony_id == "col-1",])  +
      geom_path(aes(x = day, y=ave_temp), col=temp_cols[3]) +
      geom_path(aes(x = day, y=min_temp), col=temp_cols[2]) +
      geom_path(aes(x = day, y=max_temp), col=temp_cols[1]) +
      ylim(-5,45) +
      theme_classic()+
      xlab("Day")+
      ylab("Temperature (°C)")+
      geom_vline(xintercept  = (d_field$day[d_field$outlier]))
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
           width=3, height=4)  
    
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
  sink("output/swarming_stat.txt")
  {
    cat("----------------------------------------------\n")
    cat("Logistic fitting max temp vs alate fraction in semi-field\n")
    cat("----------------------------------------------\n")
    cat("##### max temp #####\n")
    total_nest_alate = tapply(d_field$cumulative_alate, d_field$colony_id, max)
    d_field$alate_remaining_nest = rep(total_nest_alate, each = 40) - d_field$cumulative_alate
    
    r <- glmer(cbind(alates, alate_remaining_nest)~max_temp+(1|colony_id), 
               family=binomial(link="logit"), data=d_field[d_field$end ==0,])
    cat('glmer(cbind(alates, alate_remaining_nest)~max_temp+(1|colony_id), 
        family=binomial(link="logit")\n')
    cat('summary(r)\n')
    print(summary(r))
    cat('Anova(r)\n')
    print(Anova(r))
    
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
    
    cat("##### ave temp #####\n")
    r <- glmer(cbind(alates, alate_remaining_nest)~ave_temp+(1|colony_id), 
               family=binomial(link="logit"), data=d_field[d_field$end ==0,])
    print(summary(r))
    print(Anova(r))
    
    cat("##### min temp #####\n")
    r <- glmer(cbind(alates, alate_remaining_nest)~min_temp+(1|colony_id), 
               family=binomial(link="logit"), data=d_field[d_field$end ==0,])
    print(summary(r))
    print(Anova(r))
  }

  # Comparison of the largest swarming
  {
    df_largest_swarm = data.frame(
      treat = c(rep("field", 7), rep("lab", 6)),
      colony = c( unique(d_field$colony_id), unique(d_lab_20$ID)),
      largest_swarm = c( tapply(d_field$alates, d_field$colony_id, max),
                         tapply(d_lab_20$alates, d_lab_20$ID, max)),
      total_swarm = c(tapply(d_field$alates * d_field$outlier, d_field$colony_id, sum),
                      tapply(d_lab_20$alates * d_lab_20$outlier, d_lab_20$ID, sum)),
      all_alates = c( tapply(d_field$cumulative_alate, d_field$colony_id, max),
                       tapply(d_lab_20$cumlative, d_lab_20$ID, max))
      )
    
    # large swarm is prediced by treatment and/or colony size?
    cat("----------------------------------------------\n")
    cat("large swarm is prediced by treatment and/or colony size?\n")
    cat("----------------------------------------------\n")
    cat("glmer(cbind(total_swarm, all_alates-total_swarm) ~ 
    log10(total_swarm)*treat+(1|colony),family=binomial(link='logit')\n")
    r <- glmer(cbind(total_swarm, all_alates-total_swarm) ~ 
                 log10(total_swarm)*treat+(1|colony), 
               family=binomial(link="logit"), data=df_largest_swarm)
    print(summary(r))
    print(Anova(r))
    
    logistic <- function(x){ 1/(1+exp(-x)) }
    Est <- summary(r)$coefficient[,1]
    fit_x <- seq(log10(min(df_largest_swarm$all_alates)),
                 log10(max(df_largest_swarm$all_alates)),length.out=120) 
    fit_y1 <- logistic(Est[1]+(Est[2])*fit_x)
    fit_y2 <- logistic(Est[1]+Est[3]+(Est[2]+Est[4])*fit_x)
    logistic_fitting <- data.frame(fit_x, fit_y1, fit_y2)
    
    ggplot(df_largest_swarm, aes(x=log10(all_alates), 
                                 y =(total_swarm)/(all_alates),
                                 color=treat))+
      geom_point()+
      scale_color_viridis(discrete = T, end=0.5)+
      ylim(c(0,1)) +
      xlim(c(1,3.75)) +
      theme_classic()+
      theme(legend.position = "none")+
      xlab("Log10 (All alates)")+
      ylab("Fraction of alates swarmed")+
      geom_line(aes(fit_x, fit_y1), data=logistic_fitting, color=viridis(3)[1])+
      geom_line(aes(fit_x, fit_y2), data=logistic_fitting, color=viridis(3)[2])
    ggsave("output/colonysize_synchronization.pdf",
           width=3, height=3)
      
    
    
    ggplot(df_largest_swarm, aes(x=treat, y = total_swarm/all_alates)) +
      geom_point() +
      ylim(c(0,1)) +
      theme_classic()+
      xlab("")+
      ylab("Fraction of alates swarmed")
    ggsave("output/largest_swarm.pdf",
           width=3, height=3)  
  }
  sink()
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# output of colony foundation experiments
#------------------------------------------------------------------------------#
foundation_output <- function(){
  # data
  d_colony_found <- read.csv("data/raw/colony_foundation.csv",header=T)
  
  d_dish <- d_colony_found[d_colony_found$case=="dish",] # removal the results of glass cell
  
  # colony foundation success
  ggplot(d_dish, aes(x = treat, fill = as.factor(foundation))) + 
    geom_bar(position = "fill", alpha = 0.8)+
    scale_fill_viridis(discrete = T, direction = -1) +
    theme_classic()+
    theme(legend.position = "none")+
    xlab("")+
    ylab("Colony foundation success")
  ggsave("output/foundation_success.pdf")
  
  r <- glmer(foundation ~ treat + (1|colony), family=binomial(link="logit"),
             data=d_dish)
  Anova(r)
  multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
  summary(multicomparison)
  
  # offspring production
  d_suc <- d_dish[d_dish$foundation == 1,]
  
  ggplot(data=d_suc, aes(x=treat, y=offspring, fill=treat))+
    geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 1, dotsize = .5, alpha = .5) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red") +
    scale_fill_viridis(discrete = T, direction = -1) +
    theme_classic()+
    theme(legend.position = "none")+
    xlab("")+
    ylab("Number of offspring")
  ggsave("output/foundation_offspring.pdf")
  
  r <- glmer(total ~ treat + (1|colony), family=poisson, data=d_suc)
  Anova(r)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# output of hybrid experiments
#------------------------------------------------------------------------------#
hybrid_output <- function(){
  d_hybrid_found <- read.csv("data/raw/hybrid_foundation.csv",header=T)

  ggplot(d_hybrid_found, aes(x = unit, fill = as.factor(foundation))) + 
    geom_bar(position = "fill", alpha = 0.8)+
    scale_fill_viridis(discrete = T, direction = -1) +
    theme_classic()+
    theme(legend.position = "none")+
    xlab("")+
    ylab("Colony foundation success")
  ggsave("output/hybrid_success.pdf")
  
  d_hybrid_found$unit <- as.factor(d_hybrid_found$unit)
  r <- glm(foundation ~ unit, family=binomial(link="logit"),
             data=d_hybrid_found)
  Anova(r)
  
  d_suc <- d_hybrid_found[d_hybrid_found$foundation == 1,]
  ggplot(data=d_suc, aes(x=unit, y=offspring, fill=unit))+
    geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 1, dotsize = .5, alpha = .5) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red") +
    scale_fill_viridis(discrete = T, direction = -1) +
    theme_classic()+
    theme(legend.position = "none")+
    xlab("")+
    ylab("Number of offspring")
  ggsave("output/hybgrid_offspring.pdf")
  
  r <- glm(offspring ~ unit, family=poisson, data=d_suc)
  Anova(r)
  multicomparison<-glht(r,linfct=mcp(unit="Tukey"))
  summary(multicomparison)

}
#------------------------------------------------------------------------------#