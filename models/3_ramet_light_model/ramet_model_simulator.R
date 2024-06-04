# Description: R script for running tundra competition model
# Date: 2024-03-15
# Author: Ruby An

# Load Packages
library(tidyverse)
library(deSolve)


# discrete time makes more sense. choose one ramet, accumulates in full sun. next year, makes more ramets. can be a density 
# can do it as units of leaf area 
#A_i (leaf area in one year) => G_i (carbon gain, t). A_t+1 = G(t) 

# choose the ones that 

# Model functions ------

simulate_ramets <- function(t,y,p) {
  
  # number of species
  S <- length(y) 
  
  # p is a list of species-specific parameters: 
  # f : fecundity
  # m : mortality
  # k : light capture per ramet
  
  dy.dt <- rep(NA, S)
  with(as.list(p), {
    #first species 
    dy.dt[1] = L_above*f[1]*(1-exp(-k[1]*y[1]))-m[1]*y[1] 
    
    if (S > 1) {
      for (i in 2:S) {
        dy.dt[i] = L_above*f[i]*exp(-k[i]*sum(y[1:(i-1)]))*(1-exp(-k[i]*y[i])) - m[i]*y[i]
      }
      return(list(c(dy.dt)))
    } else { # only one species
      return(list(c(dy.dt)))
    }
  })
}

calculate_ramet_eq <- function(u,p) {
  
  # u is the light requirement per species that satisfies the invasion condition
  
  with(as.list(p), {
    # set vectors
    S = length(u) # f, m, and k should be the same dimension, number of species
    y = rep(NA, S)
    L = rep(NA, S) # light level
    
    # First species 
    L_total = p$L_above # light level at top of canopy
    
    # Search interval, with total light
    lower_bound <- (1/k[1])*log(L_total/u[1]) # first derivative of f_x
    upper_bound <- L_total/(k[1]*u[1]) 
    
    # numerical root find implicit solution for eq 
    num_solve <- uniroot(function(x) L_total*(1-exp(-k[1]*x)) - k[1]*u[1]*x, 
                         lower = lower_bound,
                         upper = upper_bound)
    
    # First Species EQ density
    y[1] <- max(0,num_solve$root)
    
    # Light environment for shorter species
    L[1] <- L_total*exp(-k[1]*y[1])
    
    # Species 2:S
    
    if(S > 1) {
      for (i in 2:S) {
        
        # define search interval
        lower_bound <- (1/k[i])*log(L[i-1]/u[i]) # first derivative of f_x
        upper_bound <- L[i-1]/(k[i]*u[i])
        
        # numerical root find
        num_solve <- uniroot(function(x) L[i-1]*(1-exp(-k[i]*x)) - k[i]*u[i]*x, 
                             lower = lower_bound,
                             upper = upper_bound)
        y[i] <- max(0, num_solve$root)
        
        # set light below canopy i
        L[i] <- L[i-1]*exp(-k[i]*y[i])
      }
    }  else { # do nothing, only one specie
    }
    
    L_above <- c(L_total, L[1:(i-1)])
    
    env_eq <- tibble(L_above = L_above, u = u, y = y, L_below = L) %>% 
      mutate(dL = -(L - lag(L_above,default = L_total)))
    
    return(env_eq)
  })
}

# Parameters -------

# read species-specific and species-agnostic parameters
param_file <- "3_ramet_light_model/ramet_parameters.csv"
params <- read_csv(param_file)

all_sp <- params %>% filter(species == "all")
all_sp_list <- setNames(as.list(all_sp$value), all_sp$parameter) 
sp_heights <- params %>% filter(parameter == "h") %>% arrange(desc(value))
num_species <- nrow(sp_heights)

all_env <- params %>% filter(species == "env") 

# calculate f, m, kc, and environmental parameters

# p is a list of species-specific parameters: 
# f : fecundity
# m : mortality
# kc : light capture per ramet, k is the extinction coefficient and c is canopy leaf area

p_sp <- with(as.list(all_sp_list), {
  k <- rep(k*c, num_species) # light capture, for future crown size differences
  fi <- a/(b*sp_heights$value^beta) # fecundity
  mi <- r/(b*sp_heights$value^beta) + m
  
  p <- setNames(list(fi,mi,k), c("f", "m", "k"))
  return(p)
})

# calculate light requirement (u) per species
ui <- with(as.list(p), {
  ui <- m/(f*k)
  return(ui)
})

# p_env is a list of environmental parameters
# L_above : light available above the canopy 

p_env <- setNames(as.list(all_env$value), all_env$parameter) 

# full parameter list: species + environment parameters 
p <- c(p_sp, p_env)

# Solve EQ analytically

eq_sols <- calculate_ramet_eq(ui,p) %>% mutate(species = paste0("y",1:num_species))

eq_heights <- sp_heights %>% mutate(species = paste0("y",1:num_species)) %>% 
  rename(height = value)

eq_sols <- eq_sols %>% left_join(eq_heights) %>% 
  rename(y_eq=y)

# Simulations -------

# Set timespan
t = seq(from=0,to=1/p$m[1]*200,by=0.01)

# Initial conditions
y0 = rep(0.01,num_species)

# Integrate
out = ode(y=y0,times=t,func=simulate_ramets,parms=p);

# Arrange Data
data <- as_tibble(as.data.frame(out))
colnames(data) <- c("time", paste0("y",1:num_species))
viz_data <- data %>% pivot_longer(cols = contains("y"), 
                                  names_to = "species", values_to = "value") %>% 
  left_join(eq_sols)

# Plot 
ggplot(viz_data, aes(x = time, y = value)) + 
  geom_line(aes(color = factor(height), group = factor(height)))+ 
  geom_hline(aes(yintercept = y_eq, color = factor(height)), linetype = 3)+
  theme_bw() 

