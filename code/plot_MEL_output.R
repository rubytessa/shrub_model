## Description: Plot MELVI output
## Creator: Ruby An (rubya@princeton.edu) 
## Date: 2022-11-01


## Required Packages
library(tidyverse)


## Load Data
dir <- "/home/ruby/MELVI/"
setwd(dir)
output <- read_csv("output.csv", skip = 2)


## Explore data
names(output)
str(output)

## Plot data

output$day...1

ggplot(output) + 
  geom_line(aes(x="Time", y = "biomass C"))
