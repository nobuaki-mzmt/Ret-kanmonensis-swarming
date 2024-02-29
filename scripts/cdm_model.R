# Simple simulation model for collective decision making
# based on Jeanson et al., 2012 10.3389/fnins.2012.00121

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
  sim_time     <- 100
  sim_rep      <- 1000
  num_ind_vec  <- c(10,100,1000)
  prop_vec     <- (0:10)/10
  
  run_all <- function(a, k, r){
    param <- list(a=a, k=k, r=r)
    print(param)
    # a: state change speed;
    # k: propensity of the state (non-active vs active)
    # r: strength of positive feedback
    
    RunSimulations(sim_rep, param)
    Sim_output(param)
  }

  # default
  a = 0.1
  k = c(2, 1)
  r = 2
  
  for(a in c(0.05, 0.1, 0.5)){
    run_all(a, k, r)
  }
  
  a = 0.1
  k = c(2, 1)
  r = 2
  for(ki in c(0.5, 1, 2)){
    k[1] = ki
    run_all(a, k, r)
  }
  
  a = 0.1
  k = c(2, 1)
  r = 2
  for(r in seq(1,3)){
    run_all(a, k, r)
  }
  SA_output()
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
  param_a = param$a
  param_k = param$k
  param_r = param$r
  
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
  file_name = paste0("data/df_sim_res_a-", param_a, 
         "_k1-", param_k[1], "_k2-", param_k[2], 
         "_r-", param_r, ".rda")
  save(df_sim_res, param, file = file_name)
  print(paste(file_name , " saved"))
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Output
#------------------------------------------------------------------------------#
Sim_output <- function(param){
  param_a = param$a
  param_k = param$k
  param_r = param$r
  load(paste0("data/df_sim_res_a-", param_a, 
              "_k1-", param_k[1], "_k2-", param_k[2], 
              "_r-", param_r, ".rda"))
  
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
  ggsave(paste0("output/model_resutls_a-", param_a, 
         "_k1-", param_k[1], "_k2-", param_k[2], 
         "_r-", param_r, 
         ".pdf"), width=5, height=3)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Output sensitivity analysis
#------------------------------------------------------------------------------#
SA_output <- function(){
  
  ## SA for α
  {
    df_plot <- NULL
    param_a = 0.1
    param_k = c(2, 1)
    param_r = 2
    for(param_a in c(0.05, 0.1, 0.5)){
      param_k[1] = ki
      load(paste0("data/df_sim_res_a-", param_a, 
                  "_k1-", param_k[1], "_k2-", param_k[2], 
                  "_r-", param_r, ".rda"))
      
      df_sim_sum <- ddply(df_sim_res, .(active_prop, num_ind), summarize, 
                          inactive=mean(inactive))
      df_sim_sum <- cbind(df_sim_sum,  
                          inactive_error = as.vector(tapply(df_sim_res$inactive,
                                                            df_sim_res[,1:2], sd)),
                          a = param_a,
                          k1 = param_k[1],
                          k2 = param_k[2],
                          r = param_r)
      df_plot <- rbind(df_plot, df_sim_sum)
    }
    
    ggplot(df_plot, aes(x=active_prop, y=1-inactive, col=log10(num_ind)))+
      geom_jitter(position = position_jitterdodge(dodge.width=.05, 
                                                  jitter.width = .02),
                  alpha=0.1)+
      geom_point(data = df_plot, aes(x=active_prop, y=1-inactive, 
                                     group=as.factor(num_ind), col=log10(num_ind)),
                 position = position_dodge(.05))+
      geom_path(data = df_plot, aes(x=active_prop, y=1-inactive,
                                    group=as.factor(num_ind), col=log10(num_ind)),
                position = position_dodge(.05), )+
      geom_abline(linetype = "dashed", col="grey", )+
      scale_x_continuous(breaks = c(0,0.5,1)) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_color_continuous(breaks = c(1,2,3), labels = c("10", "100", "1,000")) +
      guides(col = guide_colourbar(title = "Group size")) +
      xlab("Prop. active alates at the initial state (~temperature activation)")+
      ylab("Prop. active alates at the end \n (~dispersed individuals)")+
      facet_wrap(~ a,
                 nrow = 1, ncol = 3,
                 strip.position = "top",
                 labeller = labeller(
                   a = c("0.05" = "α = 0.05",
                         "0.1" = "α = 0.1",
                         "0.5" = "α = 0.5")))+
      theme_classic() +
      theme(strip.placement = "outside",
            strip.background = element_blank(),
            legend.position = c(0.06,0.8),
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            legend.key.height= unit(0.1, 'inch'),
            legend.key.width= unit(0.1, 'inch'),
            panel.grid = element_blank(),
            text = element_text(size = 10))
    ggsave("output/sim_sa_α.png", height = 2.5, width = 7)
    ggsave("output/sim_sa_α.pdf", height = 2.5, width = 7)
  }
  
  ## SA for k1
  {
    df_plot <- NULL
    param_a = 0.1
    param_k = c(2, 1)
    param_r = 2
    for(ki in c(0.5, 1, 2)){
      param_k[1] = ki
      load(paste0("data/df_sim_res_a-", param_a, 
                  "_k1-", param_k[1], "_k2-", param_k[2], 
                  "_r-", param_r, ".rda"))
      
      df_sim_sum <- ddply(df_sim_res, .(active_prop, num_ind), summarize, 
                          inactive=mean(inactive))
      df_sim_sum <- cbind(df_sim_sum,  
                          inactive_error = as.vector(tapply(df_sim_res$inactive,
                                                            df_sim_res[,1:2], sd)),
                          a = param_a,
                          k1 = param_k[1],
                          k2 = param_k[2],
                          r = param_r)
      df_plot <- rbind(df_plot, df_sim_sum)
    }
    
    ggplot(df_plot, aes(x=active_prop, y=1-inactive, col=log10(num_ind)))+
      geom_jitter(position = position_jitterdodge(dodge.width=.05, 
                                                  jitter.width = .02),
                  alpha=0.1)+
      geom_point(data = df_plot, aes(x=active_prop, y=1-inactive, 
                                     group=as.factor(num_ind), col=log10(num_ind)),
                 position = position_dodge(.05))+
      geom_path(data = df_plot, aes(x=active_prop, y=1-inactive,
                                    group=as.factor(num_ind), col=log10(num_ind)),
                position = position_dodge(.05), )+
      geom_abline(linetype = "dashed", col="grey", )+
      scale_x_continuous(breaks = c(0,0.5,1)) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_color_continuous(breaks = c(1,2,3), labels = c("10", "100", "1,000")) +
      guides(col = guide_colourbar(title = "Group size")) +
      xlab("Prop. active alates at the initial state (~temperature activation)")+
      ylab("Prop. active alates at the end \n (~dispersed individuals)")+
      facet_wrap(~ k1,
                 nrow = 1, ncol = 3,
                 strip.position = "top",
                 labeller = labeller(
                   k1 = c("0.5" = "k0 = 0.5",
                         "1" = "k0 = 1",
                         "2" = "k0 = 2")))+
      theme_classic() +
      theme(strip.placement = "outside",
            strip.background = element_blank(),
            legend.position = c(0.06,0.8),
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            legend.key.height= unit(0.1, 'inch'),
            legend.key.width= unit(0.1, 'inch'),
            panel.grid = element_blank(),
            text = element_text(size = 10))
    ggsave("output/sim_sa_k1.png", height = 2.5, width = 7)
    ggsave("output/sim_sa_k1.pdf", height = 2.5, width = 7)
  }
  
  ## SA for r
  {
    df_plot <- NULL
    param_a = 0.1
    param_k = c(2, 1)
    param_r = 2
    for(param_r in seq(1,3)){
      load(paste0("data/df_sim_res_a-", param_a, 
                  "_k1-", param_k[1], "_k2-", param_k[2], 
                  "_r-", param_r, ".rda"))
      
      df_sim_sum <- ddply(df_sim_res, .(active_prop, num_ind), summarize, 
                          inactive=mean(inactive))
      df_sim_sum <- cbind(df_sim_sum,  
                          inactive_error = as.vector(tapply(df_sim_res$inactive,
                                                            df_sim_res[,1:2], sd)),
                          a = param_a,
                          k1 = param_k[1],
                          k2 = param_k[2],
                          r = param_r)
      df_plot <- rbind(df_plot, df_sim_sum)
    }
    
    ggplot(df_plot, aes(x=active_prop, y=1-inactive, col=log10(num_ind)))+
      geom_jitter(position = position_jitterdodge(dodge.width=.05, 
                                                  jitter.width = .02),
                  alpha=0.1)+
      geom_point(data = df_plot, aes(x=active_prop, y=1-inactive, 
                                     group=as.factor(num_ind), col=log10(num_ind)),
                 position = position_dodge(.05))+
      geom_path(data = df_plot, aes(x=active_prop, y=1-inactive,
                                    group=as.factor(num_ind), col=log10(num_ind)),
                position = position_dodge(.05), )+
      geom_abline(linetype = "dashed", col="grey", )+
      scale_x_continuous(breaks = c(0,0.5,1)) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      guides(col = guide_colourbar(title = "Group size")) +
      xlab("Prop. active alates at the initial state (~temperature activation)")+
      ylab("Prop. active alates at the end \n (~dispersed individuals)")+
      scale_color_continuous(breaks = c(1,2,3), labels = c("10", "100", "1,000")) +
      facet_wrap(~ r,
                 nrow = 1, ncol = 3,
                 strip.position = "top",
                 labeller = labeller(
                   r = c("1" = "η = 1",
                          "2" = "η = 2",
                          "3" = "η = 3")))+
      theme_classic() +
      theme(strip.placement = "outside",
            strip.background = element_blank(),
            legend.position = c(0.06,0.8),
            legend.title = element_text(size=7),
            legend.text = element_text(size=6),
            legend.key.height= unit(0.1, 'inch'),
            legend.key.width= unit(0.1, 'inch'),
            panel.grid = element_blank(),
            text = element_text(size = 10))
    ggsave("output/sim_sa_η.png", height = 2.5, width = 7)
    ggsave("output/sim_sa_η.pdf", height = 2.5, width = 7)
  }
}
#------------------------------------------------------------------------------#