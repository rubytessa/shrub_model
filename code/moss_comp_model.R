# description: Physiological Leaf type simulations based on Weng et al. equations "big leaf" model 
# author: Ruby An (rubya@princeton.edu)
# date: 2023-04-06

# Load Packages
library(tidyverse)
library(deSolve)

# Functions 

decompose_soil <- function(temp, moisture) {
  # not explored yet, so set to temp = 1, moisture = 1 below for a linear relationship
  rate <- temp*moisture
  return(rate)
}

link_traits <- function(sigma, parameter_list) {
  ## INPUTS:
  # sigma = LMA
  # parameters = parameter list 
  
  # assign parameters within local function environment
  list2env(parameter_list, envir = environment())
  
  # soil conditions -- not explored, so set to 1
  temp <- 1 # deg C
  moisture <- 1 # random soil moisture value
  
  ## Calculate traits derived from LMA
  # leaf lifespan
  lambda = c*sigma
  
  # leaf N content
  n = A + B*sigma
  
  # leaf respiration
  R = n*r
  
  # SOM residence time 
  tau_s = decompose_soil(temp, moisture)*s*sigma
  
  ## Create list of leaf traits
  leaf_traits <- tibble::lst(sigma, lambda, n, R, tau_s)
  
  return(leaf_traits)
}

grow_moss <- function(t,y, p) {
  N1 = y[1] # leaf N
  N2 = y[2] # litter N 
  N3 = y[3] # soil N
  L1 = y[4] # leaf area
  L2 = y[5]
  
  with(as.list(p),{
    L1_new = min(1/(sigma + gamma)*V/k*(1-exp(-k*L) - R*L), # carbon limited growth 
                u/(n*sigma)*N3*L) # n limited growth
    
    dN1.dt = L_new*sigma*n - N1 / lambda
    dN2.dt = N1/lambda - N2/tau_s
    dN3.dt = N2/tau_s - L_new*sigma*n
    
    dL.dt = L_new - L/lambda
    
    return(list(c(dN1.dt, dN2.dt, dN3.dt, dL1.dt, dL2.dt)));
  })
}

# Parameters 
param_file <- "parameters/moss.p"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) # list
p_derived <- c(p_list, link_traits(0.1, p_list)) # parameters and derived parameters 

# Initial Conditions
L.0 = 1
N1.0 = L.0*p_derived$n
N2.0 = 1
N3.0 = 1
y.0 = c(L.0, N1.0, N2.0, N3.0)

# Integrate
t = seq(from=0,to=50,by=0.01)
out = ode(y=y.0,times=t,func=grow_moss,parms=p_derived) %>% as_tibble() 
colnames(out) <- c("time", "L", "N_leaf", "N_litter", "N_soil")

# Viz
out %>% pivot_longer(L:N_soil, names_to = "variable", values_to = "value") %>% 
  ggplot(aes(x=time, y=value, color=variable)) + 
  geom_line()

# Limitation Diagnostics (C or N?)

new_leaf <- with(as.list(p_derived),{
  new_leaf <- out %>% mutate(c_lim_leaf = 1/(sigma + gamma)*V/k*(1-exp(-k*L) - R*L),
                             n_lim_leaf = u/(n*sigma)*N_soil*L,
                             leaf_death = L/lambda)
  return(new_leaf)
})

new_leaf %>% pivot_longer(-time, names_to = "variable", values_to = "value") %>% 
  mutate(class = if_else(variable %in% c("c_lim_leaf", "n_lim_leaf", "leaf_death"), "rate", "state")) %>% 
  ggplot(aes(x=time, y=value, color=variable)) + 
  geom_line() + #ylim(c(0,20)) + 
  facet_wrap(vars(class))

