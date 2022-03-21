# Description: Reproducing figures from Weng et al. 2016 GCB paper on predicting vegetation type
# Author: Ruby
# Date: 2022-03-21

# Packages
library(tidyverse)


## Parameters 

# Read parameter file
param_file <- "parameters/default_species.p"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) # list


## Trait Linkages 

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

cycle_carbon <- function(sigma, N_min, parameter_list = p_list) {
  ##INPUT
  ##OUTPUT: list with the following
  # L : leaf area index
  # gain : carbon gain
  # cost :carbon cost
  # net : net carbon
  
  leaf_traits <- link_traits(sigma, parameter_list)
  list2env(parameter_list, envir = environment())
  list2env(leaf_traits, envir = environment())
  
  # leaf area
  L = N_min*lambda/n
  
  # carbon balance
  gain = V/k * (1 - exp(-k*L))
  cost = (R + G * sigma/lambda)*L
  net = gain - cost
  
  carbon_balance <- tibble(sigma, N_min, L, gain, cost, net)
  
  return(carbon_balance)
}

## testing 
#traits <- link_traits(0.05, p_list)
#cycle_carbon(0.02, p_list, 0.05)

## Figure 1

# generate data
lma <- rep(seq(0,0.25, 0.001), 5)
N_min <- rep(c(1,5,10,15,20), each = length(lma)/5)
carbon_balance <- pmap_dfr(list(lma, N_min), cycle_carbon) %>% 
  pivot_longer(cols = gain:net, names_to = "carbon", values_to = "flux")

# viz
ggplot(carbon_balance) + 
  geom_line(aes(x = sigma, y = flux, color = carbon)) + 
  facet_wrap(vars(N_min))

ggplot(carbon_balance) + 
  geom_line(aes(x = sigma, y = L))  + 
  facet_wrap(vars(N_min))

## Scrap code

# Add parameters to global environment
list2env(setNames(as.list(params$value), params$parameter), envir = .GlobalEnv)
