## Shrub Model V0
## Date: 2023-10-09
## Author: Ruby An - rubya@princeton.edu

## packages
library(tidyverse)
library(deSolve)

## Functions ----
source("shrub_model_code.R")

## Parameters -----
param_file <- "parameters/weng_units.p"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) # list
p_derived <- c(p_list, link_traits(0.002, p_list)) # parameters and derived parameters 

## Equilibria ----

# N Equilibria
list2env(p_derived, envir = environment())
NS.eq = 1/(lambda*u*a_r)
NL.eq = (N_tot*lambda - 1/(u*a_r))/(tau + lambda)
N_lit.eq = N_tot - NS.eq - NL.eq

## Timseries Simulation -----

# Initial Conditions 

## N_tot 
NL.0 <- N_tot-N_tot*0.5
NS.0 <- p_derived$N_tot - NL.0
y.0 = c(NL.0, NS.0)
y.star = c(NL.eq, NS.eq)

# Integrate
t = seq(from=0,to=50,by=1/365)
out = ode(y=y.0,times=t,func=grow_plant_n_limited,parms=p_derived) %>% as_tibble() 
colnames(out) <- c("time", "NL", "NS")

# Viz 
out %>% 
  ## calculate litter
  mutate(N_lit = p_derived$N_tot - NL - NS) %>% 
  ## calculate flux rates 
  mutate(photo = g*n*(V/k)*(1-exp(-k*NL/(n*sigma))),
         resp = g*n*(r*n*NL/(n*sigma)),
         net_gain = photo-resp) %>% 
  ## calculate N-limitation
  mutate(n_uptake = u*a_r*NS*NL) %>% 
  
  ## plotting 
  pivot_longer(NL:net_gain, names_to = "variable", values_to = "value") %>% 
  filter(variable %in% c("N_lit", "NL", "NS")) %>% 
  
  ## VIZ
  ggplot(aes(x=time, y=value, color=variable)) + 
  geom_line()  +
  
  ## simple equilibria 
  geom_hline(yintercept = NL.eq, color = "darkgreen") + 
  geom_hline(yintercept = N_lit.eq, color = "darkblue") + 
  geom_hline(yintercept = NL.eq, color = "darkblue") + 
  
  labs(y = "gN in each state (gN)", x = "years")


## plot uptake rates 
potential_growth <- out %>% 
  ## calculate litter
  mutate(N_lit = p_derived$N_tot - NL - NS) %>% 
  ## calculate flux rates 
  mutate(photo = g*n*(V/k)*(1-exp(-k*NL/(n*sigma))),
         resp = g*n*(r*n*NL/(n*sigma)),
         net_gain = photo-resp) %>% 
  ## calculate N-limitation
  mutate(n_uptake = u*a_r*NS*NL)

potential_growth %>% 
  
  ## plotting 
  pivot_longer(NL:n_uptake, names_to = "variable", values_to = "value") %>% 
  filter(variable %in% c("photo", "resp", "net_gain", "n_uptake")) %>% 
  
  ## VIZ
  ggplot(aes(x=time, y=value, color=variable)) + 
  geom_line()  +
  lims(y = c(0,250)) + 
  
  labs(y = "gN flux", x = "year")




## DIAGNOSING NO GROWTH -----

L <- seq(0,2, 0.001)
NL <- L/(n*sigma)

calculate_leaf_flux <- function(NL) {
  NS = N_tot
  photo <- (V/k)*(1-exp(-k*NL/(n*sigma)))
  resp <- r*n*NL/(n*sigma)
  
  leaf_flux <- tibble(NL, photo, resp) 
  
  return(leaf_flux)
}

leaf_flux <- calculate_leaf_flux(NL) %>% 
  mutate(net_gain = photo-resp,
         L = NL*n*sigma)

leaf_flux %>% 
  pivot_longer(cols = photo:net_gain, names_to = "variable", values_to = "value") %>% 
  ggplot(aes(x=L, y = value, color = variable)) +
  geom_line() + 
  geom_hline(yintercept = V/k, aes(color = "black")) +
  labs(x = "Leaf Area per m2 ground (m2)", y = "Carbon Flux (kgC m-2 yr-1)")
