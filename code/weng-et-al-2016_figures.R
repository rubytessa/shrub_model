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

## Figure 1 Variations ----------

## Trait Linkages 

decompose_soil <- function(temp, moisture) {
  # not explored yet, so set to temp = 1, moisture = 1 below for a linear relationship
  rate <- temp*moisture
  return(rate)
}

mineralize_N <- function(sigma, total_N, parameter_list) {
  
  # assign parameters within local function environment
  list2env(parameter_list, envir = environment())
  tau_s = decompose_soil(temp, moisture)*s*sigma
  
  N_min = total_N/(lambda + tau_s)
  
  return(N_min)
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

# generate plants and environment
lma <- seq(0,1.2, length = 20)
N_min <- round(seq(0.05, 40, length = 20), digits = 2)
grid <- expand.grid(lma = lma, N_min = N_min) %>% tibble()

# grow plants
carbon_balance <- pmap_dfr(list(grid[,1], grid[,2]), cycle_carbon) %>% 
  pivot_longer(cols = gain:net, names_to = "carbon", values_to = "flux") %>% 
  ## for plotting aesthetics 
  mutate(N_min_short = round(N_min, digits = 2 ))

# viz carbon balance
ggplot(carbon_balance) + 
  geom_line(aes(x = sigma, y = flux, color = carbon)) + 
  facet_wrap(vars(N_min))

# viz LAI
ggplot(carbon_balance) + 
  geom_line(aes(x = L, y = flux, color = carbon))  + 
  facet_wrap(vars(N_min))

## Figure 3 Variations ---------

list2env(p_list, envir = .GlobalEnv)

calculate_n_ref_optimal <- function(y) {
  # y = optimal lma
  (A + B*y)/(k*c*y)*log(V/((A+B*y)^2*r/A + G/c))
}

calculate_n_ref_ess <- function(y){ ## S36, #EQ13
  # y = ess lma
  (A + B*y)/(k*c*sigma_min)*log(V/((A+B*y)^2*r/A + G/c))
}

# ESS Thresholds
sigma_max <- (sqrt((V-G/c)/(r*A)) - 1)*A/B
sigma_ess_min <- (sqrt(V*A/((exp(1)^2*r)))- A)/B
## equivalent expression: (1/exp(1)*sqrt(V/(r*A))-1)*A/B

n_1 <- (A + B*sigma_min)/(k*c*sigma_min)*log(V/((A+B*sigma_min)^2*r/A+ G/c))
n_2 <- A/(exp(1)*k*c*sigma_min)*sqrt(V/(r*A))*log(V/(V/exp(1)^2 + G/c))

## VIZ Plot 
y <- seq(-sigma_max,sigma_max,length = 20000) # sequence of LMA values 

fig_3_transpose <- tibble(lma = y, 
                          n_ref_ess = calculate_n_ref_ess(y), 
                          n_ref_opt = calculate_n_ref_optimal(y)) %>% 
  pivot_longer(cols = starts_with("n"), names_to = "strategy", values_to = "n_min")

fig_3_transpose$strategy = factor(fig_3_transpose$strategy,levels = c("n_ref_opt", "n_ref_ess")
)

## ESS LMA
ggplot(fig_3_transpose) + 
  # Thresholds
  geom_hline(yintercept = n_1, color = "blue2", linetype =2) + 
  geom_hline(yintercept = n_2, color = "red", linetype = 2) + 
  geom_vline(xintercept = 0.02, color = "green4") + 
  geom_vline(xintercept = sigma_max, color = "green4") + 
  geom_vline(xintercept = sigma_ess_min, color = "green4", linetype = 2) + 
  # LMA
  geom_line(aes(x=lma, y= n_min)) +
  #geom_segment(aes(x = sigma_min, xend = sigma_min, y = n_1, yend = 40)) + 
  # Axes
  ylim(0,40) +
  xlim(0, sigma_max) + 
  coord_flip() + 
  facet_grid(rows = vars(strategy), 
             labeller = labeller(strategy = c(n_ref_ess = "ESS", n_ref_opt = "optimal"))) + 
  # Formatting 
  theme_bw() + 
  ggtitle("Figure 3") + 
  labs(y = "Reference N mineralization rate", x = "LMA")


### Fig 3a and 3b alone -----------
y = seq(sigma_min,2,length = 2000)
fig_3_transpose <- tibble(lma = y, 
                          n_ref_ess = calculate_n_ref_ess(y), 
                          n_ref_opt = calculate_n_ref_optimal(y)) 
## ESS LMA
(fig_3 <- ggplot(fig_3_transpose) + 
  #thresholds
  geom_hline(yintercept = n_1, color = "blue2", linetype = 2) + 
  geom_hline(yintercept = n_2, color = "red", linetype = 2) + 
  geom_vline(xintercept = sigma_min, color = "green4") + 
  geom_vline(xintercept = sigma_max, color = "green4") + 
  geom_line(aes(x=lma, y= n_ref_ess), size = 1) +
  geom_segment(aes(x = sigma_min, xend = sigma_min, y = n_1, yend = 40), size = 1) + 
  ylim(0,40) +
  coord_flip() + 
  theme_bw())

## Optimal LMA
ggplot(fig_3_transpose) + 
  geom_line(aes(x =y, y= n_ref_opt)) + 
  geom_hline(yintercept = n_1, color = "blue2") + 
  geom_hline(yintercept = n_2, color = "red") + 
  geom_vline(xintercept = 0.02, color = "green4") + 
  geom_vline(xintercept = sigma_max, color = "green4") + 
  ylim(0,40) + 
  coord_flip() + 
  theme_bw()

## Simulation Model -----

p = p_list # parameters 
t = seq(from=0,to=1000,by=0.1) # time

Bd_0 = 1;
Be_0 = 1;
Nm = Ntot/
y0 = c(Bd_0, Be_0, Nm_0) 

