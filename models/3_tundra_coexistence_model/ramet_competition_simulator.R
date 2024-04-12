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
one_species_light_comp <- function(t,y,p) {
  # This model has one species competing for light
  
  with(as.list(p),{
    
    B = b*H^beta # allometry 
    g = (a*L0*(1-exp(-k*y))/y-r_l) # carbon economy for all ramets of species
    
    dy.dt = (g/B-m)*y; 
    
    return(list(c(dy.dt)));
  })
}

n_species_light_comp <- function(t, y, p) {
  # This model has S species competing for light 
  
  S <- length(y) # of species 
  g <- rep(NA, S)
  dx.dt <- rep(NA, S)
  
  with(as.list(p),{
    
    # light levels
    y.above <- cumsum(y) # cumulative leaves in each layer
    L.beneath = exp(-k*y.above)
    L.r = -diff(c(1, L.beneath))/y
    
    # vectors per species
    B = b*H^beta # allometry
    
    # carbon balance
    g[1] = a*L0/y[1]-r_l # for first species
    for (i in 2:S) { # for rest of species
      g[i] = a*L.r[i-1] - r_l
    }
    
    # fecundity
    f = g/B 
    
    # derivatives (vector of length S)
    for (i in 1:S) {
      dx.dt[i] <- (f[i]-m)*y[i]
    }
    
    return(list(c(dx.dt)))
  })
}

light_extinction_given_u <- function(u,k) {
 
  # set vectors 
  S = length(u) 
  x = rep(NA,times=S) # density
  L = rep(NA,times=S)
  
  # Solve first species, to get first light level
  L_above = 1 # light level at top of canopy
  
  # Search interval, with total light
  lower_bound <- (1/k)*log(L_above/u[1]) # first derivative of f_x
  upper_bound <- L_above/(k*u[1]) 
  
  # numerical root find
  num_solve <- uniroot(function(x) L_above*(1-exp(-k*x)) - k*u[1]*x, 
                       lower = lower_bound,
                       upper = upper_bound)
  x[1] <- num_solve$root
  
  # set light
  L[1] <- exp(-k*x[1])
  
  if(S > 1) {
    for (i in 2:S) {
      
      # check trait is outside shadow
      
      if(u[i]<L[i-1]) {
        # define search interval
        lower_bound <- (1/k)*log(L[i-1]/u[i]) # first derivative of f_x
        upper_bound <- L[i-1]/(k*u[i])
        
        # numerical root find
        num_solve <- uniroot(function(x) L[i-1]*(1-exp(-k*x)) - k*u[i]*x, 
                             lower = lower_bound,
                             upper = upper_bound)
        x[i] <- num_solve$root
        
        # set light
        L[i] <- L[i-1]*exp(-k*x[i])
      } 
      else{ 
        x[i] = 0
        # set light
        L[i] <- L[i-1]
      }
    }
  } else {
    # do nothing
  }
  
  community_eq  <- tibble(L = L, u = u, x = x) %>% 
    mutate(dL = -(L - lag(L,default = 1)))
  
  return(community_eq)
}

# Run model for one species -------- 


# Read parameter file
param_file <- "3_tundra_coexistence_model/ramet_parameters.csv"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) 

# Set species specific parameters
S <- 1 # number of species 
H <- sort(0.4, decreasing = T) # species heights
p_add <- setNames(as.list(c(1, list(H))), 
                  c("L0", "H"))
p = c(p_list, p_add)

# Calculate light requirement
u <- with(as.list(p), {
  u = (m*b*H^beta + r_l)/(a*k)
  return(u)
})

# Solve EQ analytically
light_extinction_given_u(u,p$k)

# Set timespan
t = seq(from=0,to=1/p$m*100,by=0.01)

# Initial conditions
y0 = rep(0.1,S)

# Integrate
out = ode(y=y0,times=t,func=one_species_light_comp,parms=p);

data <- as_tibble(as.data.frame(out))
colnames(data) <- c("time", paste0("y",1:S))

# Viz data
(viz_data <- data %>% pivot_longer(cols = contains("y"), 
                                  names_to = "variable", values_to = "value") %>% 
ggplot(aes(x = time, y = value)) + 
  geom_line(aes(color = variable)) + 
  # geom_hline(yintercept = n1, color = "red") + 
  # geom_hline(yintercept = n2, color = "green") + 
  theme_bw() + 
  lims(y = c(0,max(data$y1))))


# Run model for N species -------- 

## put S in the parameters 
s_species_light_comp <- function(t, y, p) {
  # This model has S species competing for light 
  
  g <- rep(NA, S)
  dx.dt <- rep(NA, S)
  
  with(as.list(p),{
    
    # light levels
    y.above <- cumsum(y) # cumulative leaves in each layer
    L.beneath = exp(-k*y.above)
    L.r = -diff(c(1, L.beneath))/y
    
    # vectors per species
    B = b*H^beta # allometry
    
    # carbon balance
    g[1] = a*L0/y[1]-r_l # for first species
    for (i in 2:S) { # for rest of species
      g[i] = a*L.r[i-1] - r_l
    }
    
    # fecundity
    f = g/B 
    print(g)
    print(B)
    
    # derivatives (vector of length S)
    for (i in 1:S) {
      dx.dt[i] <- (f[i]-m)*y[i]
    }
    
    return(list(c(dx.dt)))
  })
}

## Integrate ----

# Read parameter file
param_file <- "3_tundra_coexistence_model/ramet_parameters.csv"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) 

# Choose number of species
number_sp <- 8

# Set species specific parameters
set.seed(1)
H <- sort(runif(number_sp)+1, decreasing = T)
p_add <- setNames(as.list(c(1, list(H))), 
                  c("L0", "H"))
p = c(p_list, p_add)

# Set timespan
t = seq(from=0,to=100,by=0.01)

# Initial conditions
y0 = rep(1,number_sp)

n_species_light_comp(t = t, y = y0, p = p)

# Integrate
out = ode(y=y0,times=t,func=n_species_light_comp,parms=p);

data <- as_tibble(as.data.frame(out))
colnames(data) <- c("time", paste0("y",1:number_sp))

viz_data <- data %>% pivot_longer(cols = contains("y"), names_to = "variable", values_to = "value")

viz_data <- left_join(viz_data, sp_height)

# Viz 

ggplot(viz_data, aes(x = time, y = value)) + 
  geom_line(aes(color = factor(height))) + 
  # geom_hline(yintercept = n1, color = "red") + 
  # geom_hline(yintercept = n2, color = "green") + 
  theme_bw() 
  #lims(x= c(0,25),y = c(0, 0.001))
  

# Does it match analytical expectations?

# Calculate light requirements per species
u <- with(as.list(p), { # light traits 
  u = (m*b*H^beta + r_l)/(a*k)
  sp_traits <- tibble(height = p$H, u = u, variable = viz_data$variable %>% unique())
  return(sp_traits)
})

light_extinction_given_u(u$u,p$k)

list2env(p, envir = environment())

a/(m*b*H^beta)

H = 0.8
n <- seq(0,10,by = 0.01)
g <- (a*(1-exp(-k*n)))/(b*H^beta)
        
curves = tibble(n=n, g=g, one = n) %>% 
  pivot_longer(g:n, values_to = "value", names_to = "variable")

ggplot(curves, aes(x=n, y = values)) + 
  geom_line()

H1 = p$H[1]
B1 = b*H1^beta # allometry - species biomass 
g1 = a*L0/tail(data$y1, n=1)-r_l # carbon economy - growth rate
T1 = B1/g1 # time to reproduction
m1 = m # demography - mortality
f1 = 1/T1 # fecundity - rate of reproduction

