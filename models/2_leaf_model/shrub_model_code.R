## Model scripts for "shrub_model_run.R" file 


grow_plant_n_limited <- function(t,y,p) {
  NL = y[1] # nitrogen in leaves
  NS = y[2] # nitrogen in soil 
  
  with(as.list(p),{
    
    U = u*a_r*NS*NL # n-limited uptake
    
    dNL.dt = U - NL/(lambda)
    dNS.dt = (N_tot - NL - NS)/tau - U
    
    return(list(c(dNL.dt, dNS.dt)));
  })
}
 
grow_plant_c_limited <- function(t,y,p) {
  NL = y[1] # nitrogen in leaves
  NS = y[2] # nitrogen in soil 
  L = NL/n  # leaf area in m2 per m2 ground area
  
  with(as.list(p),{

    U = (n/sigma)*((V/k)*(1-exp(-k*L)) - (R + G*sigma/lambda)*L) # c-limited uptake
    
    dNL.dt = U - NL/(lambda)
    dNS.dt = (N_tot - NL - NS)/tau - U
    
    return(list(c(dNL.dt, dNS.dt)));
  })
}

grow_plant_co_limited <- function(t,y,p) {
  NL = y[1] # nitrogen in leaves
  NS = y[2] # nitrogen in soil 
  L = NL/n  # leaf area in m2 per m2 ground area
  
  with(as.list(p),{
    U = min(max(0, u*a_r*NS*NL), #n limited uptake
            max(0,(n/sigma)*((V/k)*(1-exp(-k*L)) - (R + G*sigma/lambda)*L))) # c-limited uptake
        
    dNL.dt = U - NL/(lambda)
    dNS.dt = (N_tot - NL - NS)/tau - U
    
    return(list(c(dNL.dt, dNS.dt)));
  })
}

decompose_soil <- function(temp, moisture) {
  # not explored yet, so set to temp = 1, moisture = 1 below for a linear relationship
  rate <- temp*moisture
  return(rate)
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
  tau = decompose_soil(temp, moisture)*s*sigma
  
  ## Create list of leaf traits
  leaf_traits <- tibble::lst(sigma, lambda, n, R, tau)
  
  return(leaf_traits)
}
