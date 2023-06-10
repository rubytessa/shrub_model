# Resource Competition Isoclines 
# author: Ruby An
# date: 2023-01-26

# Packages
library(tidyverse)
library(deSolve)

# Functions ---- 

## Define Model with One consumer, One Resource 

one_sp_one_resource <- function(t, y, p) {
  ## inputs 
  # t: time
  # y: initial conditions [B, R] 
  # p: parameters : r (resource supply rate), 
  #                 g (plant growth rate), 
  #                 f (resource conversion efficiency)
  #                 L (plant loss rate) 
  
  
  B = y[1]
  R = y[2]
  
  ## outputs: B (biomass), R (resource level) 
  with(as.list(p),{
    dR.dt = r - g*B*R;
    dB.dt = f*g*R*B/(1 + h*g*R) - L*B; 
    return(list(c(dB.dt,  dR.dt)));
  })
}

# Run Models  ---- 

## Integrate 1 instance ----

# Assign Parameter Values 
p_values <- list(10, 5, 1, 1, 0)
p_names <- c("r", "g", "f", "L", "h") 
p <- setNames(p_values, p_names) 

# Initial conditions 
t = seq(from=0,to=100,by=0.01)
B_0 = 5
R_0 = 10
y0 = c(B_0, R_0)

# Integrate
out = ode(y=y0,times=t,func=one_sp_one_resource,parms=p);
out <- as_tibble(out) 
colnames(out) <- c("time", "B", "R")

## Visualize and analyze results 

# analytical equilibrium resource levels 
R_star = p$L/(p$f*p$g) 
B_star = p$r/(p$g*R_star)
# Viz 
out_viz <- out %>% pivot_longer(cols = B:R, names_to = "variable", values_to = "value")
ggplot(out_viz) + geom_line(aes(x = time, y = value, color = variable)) + 
  theme_bw() + 
  geom_hline(aes(yintercept = R_star, color = "R"), linetype = 2) + 
  geom_hline(aes(yintercept = B_star, color = "B"), linetype = 2)

## Scan over 2 parameters ---- 

# Set-up storage matrices 
outlist <- list()
n = 20 # number of steps to break down the variables
out_end <- tibble(B_end = rep(NA, n), R_end = rep(NA,n)) # data frame of outputs

## Assign parameters in grid - choose two variables to scan over (r, L)
r = seq(0,5, length.out = n) #resource supply rate 
L = seq(0,5, length.out = n) #loss rate
plist <- expand.grid(r=r,L=L) %>% bind_cols(g = rep(1, n = n^2), f = rep(1, n = n^2))

# Initial conditions 
t = seq(from=0,to=100,by=0.02)
B_0 = 5
R_0 = 10
y0 = c(B_0, R_0)

## Integrate & store output 
for (i in 1:nrow(plist)) {
  outlist[[i]] <- ode(y = y0, parms = plist[i,], func = one_sp_one_resource, times = t)
  out_end[i, 1] <- outlist[[i]][,2] %>% tail(n=1)
  out_end[i, 2] <- outlist[[i]][,3] %>% tail(n=1)
}

out_full <- bind_cols(plist, out_end) %>% 
  # add analytical solutions
  mutate(R_star = L/(g*f)) %>% 
  mutate(B_star = r/(g*R_star))

## VIZ
# Heatmap of bimoass for r vs L 
ggplot(out_full, aes(x=L, y=r)) + 
  geom_raster(aes(fill = B_end)) + 
  theme_minimal()

# Line plot of B_end vs L 
ggplot(out_full, aes(x=L, y = B_end)) + 
  geom_line(aes(color = factor(r))) + 
  lims(y=c(0,20))

# Line plot of B_end vs r 
ggplot(out_full, aes(x=r, y = B_end)) + 
  geom_line(aes(color = factor(L))) 


## Plot timeseries for a specific value
m <- 50
out_plot_p <- out_full[m,] # parameters
out_plot <- as_tibble(outlist[[m]]) # timeseries
colnames(out_plot) <- c("time", "B", "R")

# Viz 
out_viz <- out_plot %>% pivot_longer(cols = B:R, names_to = "variable", values_to = "value")
ggplot(out_viz) + geom_line(aes(x = time, y = value, color = variable)) + 
  theme_bw() + 
  geom_hline(aes(yintercept = out_plot_p$R_star, color = "R"), linetype = 2) + 
  geom_hline(aes(yintercept = out_plot_p$B_star, color = "B"), linetype = 2)


## ANALYTICAL ----- 

## pmap version ------- should be possible but still #TODO
# p_list 
# args0 <- list(times = list(t), yini = list(y0), params = list(p))
# args <- list(times = rep(list(t), n), 
#              yini = rep(list(y0), n), 
#              params = rep(list(p), n))
# 
# one_sp_one_resource(t=args0$times, y=args0$yini, )

