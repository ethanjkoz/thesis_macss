# Load packages
library(tidyverse)
library(janitor)
library(haven)
library(arrow)

# The file path for the CPS CEV columns available
cev_col_file <- "D:/Thesis/data/2023_CEV_Data__Current_Population_Survey_Civic_Engagement_and_Volunteering_Supplement_20250204.csv"

# Read the data
cev_cols <- read_csv(cev_col_file, show_col_types = FALSE)

# Filter columns of interest to ensure that variables are measured for all years
filtered_data <- cev_cols %>% janitor::clean_names() %>% 
  filter(in_2017_dataset == TRUE & 
         in_2019_dataset == TRUE & 
         in_2021_dataset == TRUE &
         in_2023_dataset == TRUE) %>% 
  select(
    c(
      analysis_variable_name,
      original_variable_name,
      variable_type,
      label
    )
  )

####

cps_folder <- "D:/Thesis/data/cps_cev_data"

# Get all .dta files in the directory
cps_files <- list.files(cps_folder, pattern = "\\.dta", full.names = TRUE)

# Function to read and filter .dta files
read_cps_dta <- function(file, selected_cols) {
  read_dta(file) %>%
    # change to chr for now
    mutate(across(c("penatvty", "pemntvty", "pefntvty"), as.character)) %>%
    select(any_of(c(selected_cols,
                    # Variables to merge with ACS
                    "hryear4", "state", "gtcbsa",
                    # Immigration and citizenship questions
                    "prcitshp", "penatvty", "pemntvty", "pefntvty"
    )))
}

# Load all .dta files and filter columns
cps_data_list <- map(cps_files, ~ read_cps_dta(.x, filtered_data$analysis_variable_name))

# Combine all datasets into a single tibble
cps_data <- bind_rows(cps_data_list)

cps_data <- cps_data %>%
  mutate(
    citizenship_status = case_when(
      prcitshp %in% c(1, 2, 3) ~ "Native",
      prcitshp == 4 ~ "Naturalized Citizen",
      prcitshp == 5 ~ "Non-Citizen",
      TRUE ~ NA_character_),
    foreign_born = case_when(
      penatvty %in% 100:554 ~ 1,
      # U.S. territories also counted
      penatvty %in% c(57, 66, 73, 78, 96) ~ 0,  
      TRUE ~ NA_real_),
    mother_foreign_born = case_when(
      pemntvty %in% 100:554 ~ 1,
      pemntvty %in% c(57, 66, 73, 78, 96) ~ 0,
      TRUE ~ NA_real_),
    father_foreign_born = case_when(
      pefntvty %in% 100:554 ~ 1,
      pefntvty %in% c(57, 66, 73, 78, 96) ~ 0,
      TRUE ~ NA_real_),
    birth_region = case_when(
      penatvty %in% c(1:99) ~ "US Territories",
      penatvty %in% c(100:157, 160, 162:199) ~ "Europe",
      penatvty %in% c(158:159, 161, 200:299) ~ "Asia",
      penatvty %in% c(300:399) ~ "Americas",
      penatvty %in% c(60, 500:553) ~ "Oceania",
      penatvty == 555 ~ "Elsewhere",
      TRUE ~ NA_character_),
    father_birth_region = case_when(
      penatvty %in% c(1:99) ~ "US Territories",
      penatvty %in% c(100:157, 160, 162:199) ~ "Europe",
      penatvty %in% c(158:159, 161, 200:299) ~ "Asia",
      penatvty %in% c(300:399) ~ "Americas",
      penatvty %in% c(60, 500:553) ~ "Oceania",
      penatvty == 555 ~ "Elsewhere",
      TRUE ~ NA_character_),
    mother_birth_region = case_when(
      penatvty %in% c(1:99) ~ "US Territories",
      penatvty %in% c(100:157, 160, 162:199) ~ "Europe",
      penatvty %in% c(158:159, 161, 200:299) ~ "Asia",
      penatvty %in% c(300:399) ~ "Americas",
      penatvty %in% c(60, 500:553) ~ "Oceania",
      penatvty == 555 ~ "Elsewhere",
      TRUE ~ NA_character_))

write_feather(cps_data, "D:\\Thesis\\data\\cps_cev_data\\cps_data.feather")
