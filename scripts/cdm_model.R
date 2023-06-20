# Simple model for collective decision making
# based on Jeanson et al 2012
# N. Mizumoto

param_a  <- 0.01    # state change speed
param_k  <- c(2, 1) # propensity of the state (non-active vs active)
param_r  <- 2       # strength of positive feedback
sim_time <- 1000
sim_rep  <- 1

num_ind_vec  <- c(10,20,50,100,200,500,1000,2000,5000,10000)

# Simulation
SwarmingActivationSim <- function(num_ind, active_prop){
  ind_X <- rep(c(1, 2), c((num_ind-num_ind*active_prop), (num_ind*active_prop)))
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

# Show results
res_sum = NULL
for(iter_num_ind in num_ind_vec){
  for(iter_prop in 0:10){
    for(iter in 1:sim_rep){
      print(paste(iter_num_ind, iter_prop, iter))
      res = SwarmingActivationSim(iter_num_ind, iter_prop/10)
      res_sum = c(res_sum, res[sim_time])
    }
  }
}

df = data.frame(
  num_ind = rep(num_ind_vec, each = 11*sim_rep),
  inactive_prop = rep(0:10/10, sim_rep*length(num_ind_vec)),
  inactive = res_sum
)
library(MASS)
image.plot(matrix(res_sum, 11, 10))
ggplot(df, aes(as.factor(inactive_prop), as.factor(num_ind), fill= inactive)) + 
  geom_tile()+
  scale_fill_viridis()
res_sum
library(viridis)
SwarmingActivationSim(10, 0.6)
