# Description: R script for running simplified light competition model
# Date: 2024-03-15
# Author: Ruby An

# Load Packages
library(tidyverse)
library(deSolve)

# Functions ----
eq_light_extinction <- function(u,p) {
  
  with(as.list(p), {
    # set vectors
    S = length(u) # f, m, and k should be the same dimension, number of species
    y = rep(NA, S)
    L = rep(NA, S)
    
    # First species 
    L_above = p$L_above # light above the canopy
    
    # Search interval, with total light
    lower_bound <- (1/k[1])*log(L_above/u[1]) # first derivative of f_x
    upper_bound <- L_above/(k[1]*u[1]) 
    
    # numerical root find implicit solution for eq 
    num_solve <- uniroot(function(x) L_above*(1-exp(-k[1]*x)) - k[1]*u[1]*x, 
                         lower = lower_bound,
                         upper = upper_bound)
    
    # First Species EQ density
    y[1] <- max(0,num_solve$root)
    
    # Light environment for shorter species
    L[1] <- exp(-k[1]*y[1])
    
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
    
    L_above <- c(1, L[1:(i-1)])
    
    env_eq <- tibble(L_above = L_above, u = u, y = y, L_below = L) %>% 
      mutate(dL = -(L - lag(L_above,default = 1)))
    
    return(env_eq)
  })
}

simulate_light_extinction <- function(t,y,p) {
  
  # number of species
  S <- length(y) 
  
  # p is a list of parameters: 
  # f
  # m
  # k
  
  dy.dt <- rep(NA, S)
  
  #first species 
  dy.dt[1] = f[1]*(1-exp(-k[1]*y[1]))-m[1]*y[1] 
  
  if (S > 1) {
    for (i in 2:S) {
      dy.dt[i] = f[i]*exp(-k[i]*sum(y[1:(i-1)]))*(1-exp(-k[i]*y[i])) - m[i]*y[i]
    }
    return(list(c(dy.dt)))
  } else { # only one species
    return(list(c(dy.dt)))
  }
}

eigens_of_eq <- function(y,p) {
  S <- length(y)
  lambda <- rep(NA,S)
  lambda <- with(as.list(p), {
    for (i in 1:S) {
      lambda[i] = f[i]*k[i]*exp(-k[i]*sum(y[1:i])) - m[i]
    }
    return(lambda)
  })
}

eigens_of_S <- function(S) {
  
  # Set parameters
  m <- rep(1,S)
  k <- rep(1,S)
  pre_u <- 0.01 + sort(runif(S), decreasing =T) # u_min = 0.01
  f <- m/(pre_u*k)
  p <- setNames(list(f,m,k), c("f","m","k"))
  
  # Light requirement
  u <- with(as.list(p), { 
    u <- m/(f*k)
    return(u)
  })
  
  # EQ solutions
  eq_density <- eq_light_extinction(u,p)
  y <- eq_density$y
  #eq_sols <- tibble(species = paste0("y",1:S), u = u) %>% left_join(eq_density) %>% rename(y_eq = y) 
  
  # Eigenvalues of EQ
  lambdas <- with(as.list(p), {
    lambdas <- rep(NA,S)
    
    for (i in 1:S) {
      lambdas[i] = f[i]*k[i]*exp(-k[i]*sum(y[1:i])) - m[i]
    }
    return(lambdas)
  })
  
  # Output
  return(c(max(lambdas), mean(lambdas), var(lambdas), sum(y)))
}

# Parameters ---- 
# Choose number of species + initial conditions
S <- 10
y0 <- rep(0.01, S)

# Set parameters
m <- rep(1,S)
k <- rep(1,S)
pre_u <- 0.01 + sort(runif(S), decreasing =T) # u_min = 0.01
f <- m/(pre_u*k)
p <- setNames(list(f,m,k), c("f","m","k"))

# EQ Solutions ---- 

# Analytical solutions
u <- with(as.list(p), { # redundant, but if you choose parameters differently
  u <- m/(f*k)
  return(u)
})

(eq_density <- eq_light_extinction(u,p))
eq_sols <- tibble(species = paste0("y",1:S), u = u) %>% left_join(eq_density) %>% 
  rename(y_eq = y)  %>% 
  mutate(light_surplus = L_above - u)

# plot numerical solution
(eq_plot <- eq_sols %>%
  ggplot(aes(x=u, y=y_eq, color = (L_above-u)/u))+ 
    geom_point() + 
  #geom_vline(xintercept = eq_state$L, color = "grey") + 
 # scale_color_manual(values = c("red", "black")) +
  geom_hline(yintercept=0, color = "red", linetype = 2))
  #geom_rect(aes(xmin = L, xmax = u,ymin = -0.02,ymax = 0.02, alpha = if_else(x==0, Inf, 0.4))) + 

# Simulations ----

# Set timespan
t = seq(from=0,to=300,by=0.01)

# Integrate
out = ode(y=y0,times=t,func=simulate_light_extinction,parms=p);

# Arrange Data
data <- as_tibble(as.data.frame(out))
colnames(data) <- c("time", paste0("y",1:S))
viz_data <- data %>% pivot_longer(cols = contains("y"), 
                                  names_to = "species", values_to = "value") %>% 
  left_join(eq_sols)

# Plot 
ggplot(viz_data, aes(x = time, y = value)) + 
  geom_line(aes(color = u, group = factor(u)))+ 
  geom_hline(aes(yintercept = y_eq, color = u), linetype = 3)+
  theme_bw() 

# Eigenvalues -----
# Jacobian of this model is a lower triangular matrix (with zero's in the upper bit)
# This means the eigenvalues are the elements along the diagonal

lambda <- rep(NA,S)
y <- eq_sols$y_eq

lambda <- with(as.list(p), {
  for (i in 1:S) {
    lambda[i] = f[i]*k[i]*exp(-k[i]*sum(y[1:i])) - m[i]
  }
  return(lambda)
})

# Number Species & Stability ----

iter = 100
S <- rep(seq(10,10000, by = 10), each=iter)
max_lambda <- rep(NA, length(S))
mean_lambda <- rep(NA, length(S))
var_lambda <- rep(NA, length(S))
sum_x <- rep(NA, length(S))
for (i in 1:length(S)) {
  ls <- eigens_of_S(S[i])
  max_lambda[i] <- ls[1]
  mean_lambda[i] <- ls[2]
  var_lambda[i] <- ls[3]
  sum_x[i] <- ls[4]
}

# plot eigenvalues
k<- 1
tibble(S, max_lambda, mean_lambda, var_lambda, sum_x) %>% 
  ggplot(aes(x=S, y= max_lambda)) +
  geom_point(alpha=0.2)

#histograms
hist(mean_lambda)
hist(exp(-k*sum_x))
  
tibble(S,lambda) %>% ggplot(aes(x=S, y = lambda)) +
  geom_point(alpha = 0.2)

df_lambda <- tibble(S = as.numeric(S),lambda)
ggplot(df_lambda) + 
  geom_point(aes(x=lambda, y=S, color = S), alpha = 0.2)

# three community sizes
N <- rep(c(10,100,1000), each=iter)
lambda <- rep(NA, length(N))
for (i in 1:length(N)) {
  lambda[i] <- eigens_of_S(N[i])
}

tibble(N,lambda) %>% ggplot(aes(x=N, y = lambda)) +
  geom_point(alpha = 0.2)

# Per capita growth rates ----- 

eigens_of_S(100)
