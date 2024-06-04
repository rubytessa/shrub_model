
### other stuff
```{r}
## plot uptake rates 
output %>% 
  
  ## plotting 
  pivot_longer(NL:n_lim_uptake, names_to = "variable", values_to = "value") %>% 
  filter(variable %in% c("photo", "resp", "net_gain")) %>% 
  
  ## VIZ
  ggplot(aes(x=time, y=value, color=variable)) + 
  geom_line()  +
  #lims(y = c(0,250)) + 
  
  labs(y = "gN flux", x = "year")


## Plot Leaf Area as a function of NL

tibble(NL.seq = seq(0, N_tot, by = 1)) %>% 
  mutate(LeafArea = NL.seq/n) %>% 
  
  ggplot() + 
  geom_line(aes(x=NL.seq, y = LeafArea)) + 
  labs(y = "Leaf Area (m2 m-2)", x = "Leaf N (gN m-2)")

## Plot of photosynthesis and respiration as a function of L

tibble(L = seq(0, 5, by = 0.01)) %>% 
  mutate(photo = (V/k)*(1-exp(-k*L)),
         resp = (R + G*sigma/lambda)*L) %>% 
  mutate(net_gain = photo - resp) %>% 
  pivot_longer(photo:net_gain, names_to = "variable", values_to = "value") %>% 
  filter(variable %in% c("photo", "resp", "net_gain")) %>% 
  
  ggplot() + 
  geom_line(aes(x=L, y = value, color = variable)) + 
  labs(y = "carbon flux", x = "leaf area")


calculate_leaf_flux <- function(NL) {
  NS = N_tot
  photo <- (V/k)*(1-exp(-k*NL/(n*sigma)))
  resp <- r*n*NL/(n*sigma)
  
  leaf_flux <- tibble(NL, photo, resp) 
  
  return(leaf_flux)
}
```
## simulating 5 species 2024-05-30
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
param_file <- "3_ramet_light_model/ramet_parameters.csv"
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

