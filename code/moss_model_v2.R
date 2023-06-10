## Moss Model V2 
## Date: 2023-06-10
## Ruby & Helen


## packages
library(tidyverse)
library(deSolve)

# Functions 

grow_moss_v2 <- function(t,y,p) {
  NL = y[1] # nitrogen in leaves
  NS = y[2] # nitrogen in soil 
  
  with(as.list(p),{
    dNL.dt = k*a_r*NS*NL - NL/lambda
    dNS.dt = (N_tot - NL - NS)/tau - k*a_r*NS*NL
    
    return(list(c(dNL.dt, dNS.dt)));
  })
}

# Parameters 
param_file <- "parameters/moss_v2.p"
params <- read_csv(param_file)
p_list <- setNames(as.list(params$value), params$parameter) # list
p_derived <- p_list # c(p_list, link_traits(0.1, p_list)) # parameters and derived parameters 

## Equilibria

list2env(p_derived, envir = environment())
NS.eq = 1/(lambda*k*a_r)
NL.eq = (N_tot*lambda - 1/(k*a_r))/(tau + lambda)

# Initial Conditions
NL.0 <- 10
NS.0 <- p_derived$N_tot - NL.0
y.0 = c(NL.0, NS.0)
y.star = c(NL.eq, NS.eq)

# Integrate
t = seq(from=0,to=500,by=0.1)
out = ode(y=y.star,times=t,func=grow_moss_v2,parms=p_derived) %>% as_tibble() 
colnames(out) <- c("time", "NL", "NS")

# Viz 
out %>% 
  mutate(N_lit = p_derived$N_tot - NL - NS) %>% 
  pivot_longer(NL:N_lit, names_to = "variable", values_to = "value") %>% 
  ggplot(aes(x=time, y=value, color=variable)) + 
  geom_line() 


