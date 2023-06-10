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
  
  R = y[1]
  B = y[2]
  
  ## outputs: B (biomass), R (resource level) 
  with(as.list(p),{
    dR.dt = r - g*R/(1+g*h*R)*B;
    dB.dt = f*g*R*B/(1 + g*h*R) - L*B; 
    return(list(c(dR.dt,  dB.dt)));
  })
}

two_sp_one_resource <- function(t, y, p) { # one sp/one resource consumer resource model w/handling time
  ## inputs 
  # t: time
  # y: initial conditions [R, B1, B2] 
  # p: parameters : r (resource supply rate), 
  #                 g_n (plant growth rate), 
  #                 f_n (resource conversion efficiency)
  #                 L_n (plant loss rate) 
  #                 h_n (handling time, or saturating, set h = 0 to get linear growth function)
  
  R = y[1]
  B1 = y[2]
  B2 = y[3]
  
  ## outputs: B (biomass), R (resource level) 
  with(as.list(p),{
    dR.dt = r - g1*R*B1/(1 + h1*g1*R) - g2*R*B2/(1 + h2*g2*R)
    dB1.dt = f1*g1*R*B1/(1 + h1*g1*R) - L1*B1; 
    dB2.dt = f2*g2*R*B2/(1 + h2*g2*R) - L2*B2; 
    
    return(list(c(dR.dt, dB1.dt, dB2.dt )));
  })
}