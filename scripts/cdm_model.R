# Simple simulation model for collective decision making
# based on Jeanson et al 2012
# N. Mizumoto

rm(list = ls())

#------------------------------------------------------------------------------#
{
  library(plyr)
  library(ggplot2)
  library(viridis)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  param <- list(a = 0.01, k  = c(2, 1), r  = 2)
  # a: state change speed;
  # k: propensity of the state (non-active vs active)
  # r: strength of positive feedback
  
  sim_time     <- 1000
  sim_rep      <- 100
  num_ind_vec  <- c(10,100,1000)
  prop_vec     <- (0:10)/10
  RunSimulations(sim_rep, param)
  Sim_output()
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Simulation function
#------------------------------------------------------------------------------#
SwarmingActivationSim <- function(num_ind, active_prop, param=param){
  param_a = param$a
  param_k = param$k
  param_r = param$r
  ind_X <- rep(c(1, 2), c((num_ind-round(num_ind*active_prop)), round(num_ind*active_prop)))
  X     <- matrix(num_ind/2, sim_time, 2)
  for(i_time in 2:sim_time){
    for(i_ind in 1:num_ind){
      prob_change = param_a / (param_k[ind_X[i_ind]] + 
                               (X[i_time - 1, ind_X[i_ind]])^param_r/
                                 (X[i_time - 1, 3-ind_X[i_ind]])^param_r)
      rnd = runif(1,0,1)
      if(rnd < prob_change){
        ind_X[i_ind] = 3 - ind_X[i_ind]
      }
    }
    X[i_time, 1] = sum(ind_X == 1)
    X[i_time, 2] = sum(ind_X == 2)
  }
  return(X[, 1] / num_ind)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Run simualtions
#------------------------------------------------------------------------------#
RunSimulations <- function(sim_rep, param){
  res_sum = NULL
  for(iter_num_ind in num_ind_vec){
    for(iter_prop in prop_vec){
      for(iter in 1:sim_rep){
        print(paste(iter_num_ind, iter_prop, iter))
        res = SwarmingActivationSim(iter_num_ind, iter_prop, param = param)
        res_sum = c(res_sum, res[sim_time])
      }
    }
  }
  
  df_sim_res = data.frame(
    num_ind = rep(num_ind_vec, each = 11*sim_rep),
    active_prop = rep(rep(prop_vec, each=sim_rep), length(num_ind_vec)),
    inactive = res_sum
  )
  save(df_sim_res, file = "data/df_sim_res.rda")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Output
#------------------------------------------------------------------------------#
Sim_output <- function(){
  load("data/df_sim_res.rda")
  
  df_sim_sum <- ddply(df_sim_res, .(active_prop, num_ind), summarize, 
                      inactive=mean(inactive))
  df_sim_sum <- cbind(df_sim_sum,  
                      inactive_error=as.vector(tapply(df_sim_res$inactive, df_sim_res[,1:2], sd)))
  
  ggplot(df_sim_res, aes(x=active_prop, y=1-inactive, col=log10(num_ind)))+
    #scale_color_viridis(discrete=T, end=0, begin = 1) +
    geom_jitter(position = position_jitterdodge(dodge.width=.05, jitter.width = .02),
                alpha=0.1)+
    theme_classic() +
    theme(legend.position = "none")+
    geom_point(data = df_sim_sum, aes(x=active_prop, y=1-inactive, 
                                      group=as.factor(num_ind), col=log10(num_ind)),
               position = position_dodge(.05))+
    geom_path(data = df_sim_sum, aes(x=active_prop, y=1-inactive,
                                     group=as.factor(num_ind), col=log10(num_ind)),
              position = position_dodge(.05), )+
    geom_abline(linetype = "dashed", col="grey", )+
    xlab("Propotion of active alates at the initial state (~temperature activation)")+
    ylab("Propotion of active alates at the end (~dispersed individuals)")
  ggsave("output/model_resutls.pdf", width=5, height=3)
}
#------------------------------------------------------------------------------#