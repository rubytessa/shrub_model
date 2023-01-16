# description: Physiological Leaf type simulations based on Weng et al. equations
# author: Ruby An (rubya@princeton.edu)
# date: 2022-11-10

# Load Packages
library(tidyverse)
library(deSolve)

one_consumer_one_resource_decomp <- function(t,y,p) {
  # This system tends to go to oscillations of one consumer. If you add carrying 
  # capacity to the resource, it will stabilize. 
  
  n1 = y[1]
  R = y[2]

  
  with(as.list(p),{
    dR.dt = r*d1*n1-a1*c1*n1*R;
    dn1.dt = (e1*c1*a1*R - d1)*n1; 
    return(list(c(dn1.dt,  dR.dt)));
  })
}

two_consumer_one_resource <- function(t, y, p) {
  # This system tends to go to oscillations of one consumer. If you add carrying 
  # capacity to the resource, it will stabilize. 
  
  n1 = y[1]
  n2 = y[2]
  R = y[3]
  
  with(as.list(p),{
    dR.dt = (r-a1*c1*n1 - a2*c2*n2)*R;
    dn1.dt = (e1*c1*a1*R - d1)*n1; 
    dn2.dt = (e2*c2*a2*R - d2)*n2; 
    return(list(c(dn1.dt, dn2.dt, dR.dt)));
  })
}

two_consumer_one_resource_capped <- function(t, y, p) {
  # This system tends to go to oscillations of one consumer. If you add carrying 
  # capacity to the resource, it will stabilize. 
  
  n1 = y[1]
  n2 = y[2]
  R = y[3]
  
  with(as.list(p),{
    dR.dt = (r/k*(k-R)-a1*c1*n1 - a2*c2*n2)*R;
    dn1.dt = (e1*c1*a1*R - d1)*n1; 
    dn2.dt = (e2*c2*a2*R - d2)*n2; 
    return(list(c(dn1.dt, dn2.dt, dR.dt)));
  })
}

two_consumer_one_resource_decomp <- function(t, y, p) {
  # This system tends to go to oscillations of one consumer. If you add carrying 
  # capacity to the resource, it will stabilize. 
  
  n1 = y[1]
  n2 = y[2]
  R = y[3]
  
  with(as.list(p),{
    dR.dt = r + d1*r*n1 + d2*r*n2 + (-a1*c1*n1 - a2*c2*n2)*R;
    dn1.dt = (e1*c1*a1*R - d1)*n1; 
    dn2.dt = (e2*c2*a2*R - d2)*n2; 
    return(list(c(dn1.dt, dn2.dt, dR.dt)));
  })
}


# Read parameter file
param_file <- "parameters/two_consumers_one_resource_parameters.csv"
params <- read_csv(param_file)
p <- setNames(as.list(params$value), params$parameter) 

# Run the model -------- 
t = seq(from=0,to=100,by=0.01)

# initial conditions
n1_0 = 5
n2_0 = 5
R_0 = 10
y0 = c(n1_0, n2_0, R_0)

# Integrate
out = ode(y=y0,times=t,func=two_consumer_one_resource_decomp,parms=p);

plot(out) 

# Viz 
list2env(p, envir = environment())

R1 = d1/(e1*c1*a1)
R2 = d2/(e2*c2*a2)
n1 = r/(a1*c1)
n2 = r/(a2*c2)

R1
R2
n1
n2

data <- as_tibble(as.data.frame(out))
colnames(data) <- c("time", "n1", "n2", "R")

viz_data <- data %>% pivot_longer(cols = n1:R, names_to = "variable", values_to = "value")

ggplot(viz_data, aes(x = time, y = value)) + 
  geom_line(aes(color = variable)) + 
  # geom_hline(yintercept = n1, color = "red") + 
  # geom_hline(yintercept = n2, color = "green") + 
  theme_minimal()
