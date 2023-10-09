# description: Physiological Leaf type simulations based on Weng et al. equations
# author: Ruby An (rubya@princeton.edu)
# date: 2022-11-10

# Load Packages
library(tidyverse)
library(deSolve)

# Plotting 
theme_set(theme_cowplot())

# Functions --------
decompose_soil <- function(temp, moisture) {
  # not explored yet, so set to temp = 1, moisture = 1 below for a linear relationship
  rate <- temp*moisture
  return(rate)
}

cycle_carbon <- function(sigma, N_min, parameter_list = p_list) {
  ##INPUT
  ##OUTPUT: list with the following
  # L : leaf area index
  # gain : carbon gain
  # cost :carbon cost
  # net : net carbon

  # soil conditions -- not explored, so set to 1
  temp <- 1 # deg C
  moisture <- 1 # random soil moisture value
  
  with(as.list(p),{
    ## Calculate traits derived from LMA
    # leaf lifespan
    lambda = c*sigma
    
    # leaf N content
    n = A + B*sigma
    
    # leaf respiration
    R = n*r
    
    # SOM residence time 
    tau_s = decompose_soil(temp, moisture)*s*sigma
    
    # N mineralization rate 
    #N_min = N_tot / (lambda + tau_s)
    
    # leaf area
    L = N_min*c*sigma/(A + B*sigma)
    
    # carbon balance
    gain = V/k * (1 - exp(-k*L))
    cost = (R + G * sigma/lambda)*L
    net = gain - cost
    
    return(net)
    
  })
}

mineralize_N <- function(sigma, total_N, parameter_list) {
  
  # link traits 
  leaf_traits <- link_traits(sigma, parameter_list)
  list2env(parameter_list, envir = environment())
  list2env(leaf_traits, envir = environment())
  
  N_min = total_N/(lambda + tau_s)
  
  return(N_min)
}

list2env(p, envir = environment())


grow_two_plants <- function(t, y, p) { 
  B1 = y[1]
  B2 = y[2]
  
  with(as.list(p),{
    sigma_avg = B1/(B1 + B2)*sigma_1 + B2/(B1 + B2)*sigma_2
    N_min = N_tot / ((c + s)*sigma_avg)
    L = N_min*c*sigma_1/(A + B*sigma_1)

    #dB1.dt = cycle_carbon(sigma_1, N_min, p);
    #dB2.dt = cycle_carbon(sigma_2, N_min, p);
    dB1.dt = V/k * (1 - exp(-k*N_min*c*sigma_1/(A + B*sigma_1))) - 
      ((A + B*sigma_1)*r + G * sigma_1/(c*sigma_1))*N_min*c*sigma_1/(A + B*sigma_1);
    dB2.dt = V/k * (1 - exp(-k*N_min*c*sigma_2/(A + B*sigma_2))) - 
      ((A + B*sigma_2)*r + G * sigma_2/(c*sigma_2))*N_min*c*sigma_2/(A + B*sigma_2)
    return(list(c(dB1.dt, dB2.dt)));
  })
}

grow_two_plants_with_N <- function(t, y, p) { 
  B1 = y[1]
  B2 = y[2]
  N_tot = y[3]
  
  with(as.list(p),{
    sigma_avg = B1/(B1 + B2)*sigma_1 + B2/(B1 + B2)*sigma_2
    N_min = N_tot / ((c + s)*sigma_avg)
    L = N_min*c*sigma_1/(A + B*sigma_1)
    
    #dB1.dt = cycle_carbon(sigma_1, N_min, p);
    #dB2.dt = cycle_carbon(sigma_2, N_min, p);
    dB1.dt = V/k * (1 - exp(-k*N_min*c*sigma_1/(A + B*sigma_1))) - 
      ((A + B*sigma_1)*r + G * sigma_1/(c*sigma_1))*N_min*c*sigma_1/(A + B*sigma_1);
    dB2.dt = V/k * (1 - exp(-k*N_min*c*sigma_2/(A + B*sigma_2))) - 
      ((A + B*sigma_2)*r + G * sigma_2/(c*sigma_2))*N_min*c*sigma_2/(A + B*sigma_2)
    return(list(c(dB1.dt, dB2.dt)));
  })
}


# Parameters --------

# Read parameter file
param_file <- "parameters/default_species.p"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) # list
p_add <- setNames(as.list(c(0.02, 2, 40)), 
                  c("sigma_1", "sigma_2", "N_tot"))

# Run the model -------- 
p = c(p_list, p_add)
t = seq(from=0,to=250,by=0.01)

# initial conditions
B1_0 = 1
B2_0 = 1
y0 = c(B1_0, B2_0)

# Integrate
out = ode(y=y0,times=t,func=grow_two_plants,parms=p);

plot(out) 

# Viz 
list2env(p, envir = environment())

data <- as.data.frame(out) 
names(data) <- c("time", "B1", "B2")
plant_data <- data %>% 
  mutate(N_min = p$N_tot/(B1/(B1 + B2)*sigma_1 + B2/(B1 + B2)*sigma_2)) %>% 
  mutate(L1 = N_min*c*sigma_1/(A + B*sigma_1), 
         L2 = N_min*c*sigma_2/(A + B*sigma_2), 
         gain1 = V/k * (1 - exp(-k*L1)),
         gain2 = V/k * (1 - exp(-k*L2)),
         cost1 = (r*(A+B*sigma_1) + G * sigma_1/(c*sigma_1))*L1,
         net1 = gain1-cost1,
         cost2 = (r*(A+B*sigma_2) + G * sigma_2/(c*sigma_2))*L2,
         net2 = gain2-cost2, 
         dB1.dt = cycle_carbon(sigma_1, N_min, p),
         dB2.dt = cycle_carbon(sigma_2, N_min, p))
  
plot_data <- plant_data %>% 
  pivot_longer(B1:dB2.dt,names_to = "variable", values_to = "biomass")

ggplot(plot_data) + 
  geom_line(aes(x=time, y = biomass)) + 
  facet_wrap(vars(variable), scales = "free")
