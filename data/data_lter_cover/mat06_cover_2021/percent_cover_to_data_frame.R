rm(list=ls())

library(fs)
library(data.table)
library(readxl)
library(tidyverse)
library(lubridate)

# This is the directory/folder where you saved the percent cover excel files
data.dir <- ("/home/ruby/tundra_coexistence/data/mat06_cover_2021/data/")

# This creates a list of paths to the data sheets.
file.paths <- data.dir %>% 
  list.files(pattern = "*.xlsx", full.names = T)

# Create a 3-step function
join_data_sheets <- function(file.name){
  # Step 1: create unique ID
  # pull out site, block, treatment, date information from the data sheet

  date <- as.character(read_excel(file.name, col_names = F, range = "A2"))

  id.data <- as.character(read_excel(file.name, col_names = F, n_max = 1))
  # identify strings that contain the words site ^[S], block ^[B], treatment ^[T]
  id <- tibble(DATE = mdy(unlist(strsplit(date[grep("^[Dd]", date)], "\\: "))[2]),
               SITE = unlist(strsplit(id.data[grep("^[Ss]", id.data)], "\\: "))[2],
               BLOCK = unlist(strsplit(id.data[grep("^[Bb]", id.data)], "\\: "))[2],
               TREATMENT = unlist(strsplit(id.data[grep("^[Tt]", id.data)], "\\: "))[2])
  
  # Step 2: Cover Data
  # insert the percent cover data
  cover.data <- read_excel(file.name,
                           skip = 3, n_max = 31, 
                           col_names = T,
                           col_types = "text") %>% 
    select(-...10) %>% 
    pivot_longer('1':'8',
                 names_to = 'QUADRAT',
                 values_to = 'percent_cover') %>% 
    
    pivot_wider('QUADRAT',
                names_from = '...1', values_from = 'percent_cover')
  
  # Step 3: Merge the data and remove NAs
  data <- as_tibble(merge(id, cover.data))
}

# apply this new function to the list of file paths to create one long data frame
final <- map_dfr(file.paths, join_data_sheets) %>%
  
# Remove columns that just helped the data sheet look nice
  select(-'NA', -'VOLES')

colnames(final) <- tolower(colnames(final))

# These three lines check to see if there are any naming inconsistencies
colnames(final)
unique(final$site)
unique(final$treatment)

# Write the final file to a .csv in the location of the data sheets
write_csv(x = final, file = paste0(data.dir, 'percentCover', year(today()), '.csv'))

