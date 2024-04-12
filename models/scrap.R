
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
