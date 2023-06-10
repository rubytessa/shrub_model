## Description: Plot ITEX trait data from Betway Tundra Traits papers 
## Author: Ruby An
## Date: 2022-02-14


## Required Packages
library(tidyverse)


itex_traits <- read_csv("data/itex-aon_Betway_PlantTraits_2018.csv",skip =1)

names(itex_traits)

ggplot(itex_traits, aes(x = SLA, y = PlantHeight_Veg)) + 
  geom_point(aes(color = Taxon)) + 
  scale_x_continuous(limits = c(0,500))


itex_traits %>% filter(Taxon %in% c("BETNAN", "SALPUL", "LEDPAL", "VACVIT")) %>% 
  group_by(Taxon, SiteName) %>%
  ## units cm2/g to m2/kg
  mutate(SLA = .1*SLA) %>% 
  summarize(mean_LMA = 1/mean(SLA), mean_SLA = mean(SLA), sd_SLA = sd(SLA))
