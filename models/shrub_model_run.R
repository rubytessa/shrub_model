## Moss Model V2 
## Date: 2023-06-10
## Ruby & Helen

## packages
library(tidyverse)
library(deSolve)

## Functions ----

grow_moss_v2 <- function(t,y,p) {
  NL = y[1] # nitrogen in leaves
  NS = y[2] # nitrogen in soil 
  
  with(as.list(p),{
    U = min(u*a_r*NS*NL, # n limited uptake 
              g*n*((V/k)*(1-exp(-k*NL/(n*sigma))) - r*n*NL/(n*sigma))) # c-limited uptake
    
    # U = min(u*a_r*NS*NL, # n limited uptake 
    #        g*n*(V/k*NL/(n*sigma) - r*n*NL/(n*sigma))) # c-limited uptake
    # 
    # simple_photo = V/k*NL/(n*sigma)
    # #photo = g*n*((V/k)*(1-exp(-k*NL/(n*sigma))))
    #U = u*a_r*NS*NL
    #print(simple_photo)
    
    dNL.dt = U - NL/(lambda)
    dNS.dt = (N_tot - NL - NS)/tau - U
    
    return(list(c(dNL.dt, dNS.dt)));
  })
}

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
  tau = decompose_soil(temp, moisture)*s*sigma
  
  ## Create list of leaf traits
  leaf_traits <- tibble::lst(sigma, lambda, n, R, tau)
  
  return(leaf_traits)
}

##Parameters -----
#param_file <- "parameters/moss_v2.p"
param_file <- "parameters/weng_units.p"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) # list
p_derived <- c(p_list, link_traits(0.002, p_list)) # parameters and derived parameters 

## Equilibria ----

# N Equilibria
list2env(p_derived, envir = environment())
NS.eq = 1/(lambda*u*a_r)
NL.eq = (N_tot*lambda - 1/(u*a_r))/(tau + lambda)


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

## Timseries Simulation -----

# Initial Conditions 
NL.0 <- N_tot-N_tot*0.5
NS.0 <- p_derived$N_tot - NL.0
y.0 = c(NL.0, NS.0)
y.star = c(NL.eq, NS.eq)

# Integrate
t = seq(from=0,to=50,by=0.01)
out = ode(y=y.0,times=t,func=grow_moss_v2,parms=p_derived) %>% as_tibble() 
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
  geom_hline(yintercept = V/k, color = "gray") + 
  geom_hline(yintercept = NL.eq, color = "red") + 
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
