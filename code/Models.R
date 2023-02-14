## Model scripts for "nested_cr_models.Rmd" file 

one_sp_one_resource <- function(t, y, p) { # one sp/one resource consumer resource model w/handling time
  ## inputs 
  # t: time
  # y: initial conditions [B, R] 
  # p: parameters : r (resource supply rate), 
  #                 g (plant growth rate), 
  #                 f (resource conversion efficiency)
  #                 L (plant loss rate) 
  #                 h (handling time, or saturating, set h = 0 to get linear growth function)
  
  B = y[1]
  R = y[2]
  
  ## outputs: B (biomass), R (resource level) 
  with(as.list(p),{
    dR.dt = r - g*B*R;
    dB.dt = f*g*R*B/(1 + h*g*R) - L*B; 
    return(list(c(dB.dt,  dR.dt)));
  })
}
