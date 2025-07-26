# Much of this code is aided by AI for ease of gathering important data attributes

library(tidyverse)
library(arrow)
library(haven)
library(flextable)
library(apc)
library(cobalt)
library(gridExtra)

# from github
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linmed.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwcde.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R")

acs = read_feather("D:\\Thesis\\data\\acs_data\\acs.feather")
cev = read_feather("D:\\Thesis\\data\\cps_cev_data\\cps_data.feather")

years = c(2017, 2019, 2021, 2023)
boot = 1000
n_sim = 1000
random_seed = 0

topmsa_gtcbsas = cev %>% drop_na(topmsa) %>% dplyr::select(topmsa, gtcbsa) %>% 
  zap_labels() %>% pull(gtcbsa) %>% unique()

### ACS ###

acs_clean = acs %>%
  # if I want to make claims about the MSAs, I need to use top 12 only
  filter(gtcbsa %in% topmsa_gtcbsas) %>%
  mutate(
    # coalescing required because of name changes in 2015–17 vs 19–23
    median_age = coalesce(sex_and_age_median_age_yrs,
                          sex_and_age_tot_pop_median_age_yrs),
    moe_median_age = coalesce(moe_sex_and_age_median_age_yrs,
                              moe_sex_and_age_tot_pop_median_age_yrs),
    pct_white = coalesce(pct_race_1_race_white,
                         pct_race_tot_pop_1_race_white),
    moe_pct_white = coalesce(pct_moe_race_1_race_white,
                             pct_moe_race_tot_pop_1_race_white),
    pct_col_ed = coalesce(pct_edu_attain_pct_bachelors_degree_up,
                          pct_edu_attain_pop_25_yrs_up_bachelors_degree_up),
    moe_pct_col_ed = coalesce(pct_moe_edu_attain_pct_bachelors_degree_up,
                              pct_moe_edu_attain_pop_25_yrs_up_bachelors_degree_up)
  ) %>%
  rename(
    med_income = income_and_benefits_infl_adj_usd_fams_median_fam_income_usd,
    moe_med_income = moe_income_and_benefits_infl_adj_usd_fams_median_fam_income_usd,
    
    unemploy_rate = pct_employment_status_civ_labor_force_unemployment_rate,
    moe_unemploy_rate = pct_moe_employment_status_civ_labor_force_unemployment_rate,
    
    # pct_hisp = pct_hisp_lat_and_race_tot_pop_hisp_lat,
    # moe_pct_hisp = pct_moe_hisp_lat_and_race_tot_pop_hisp_lat, # has many missing
    
    per_k_birth_12mo = fertility_num_women_15_to_50_yo_had_birth_in_past_12_months_per_1k_women_15_to_50_yo,
    moe_per_k_birth_12mo = moe_fertility_num_women_15_to_50_yo_had_birth_in_past_12_months_per_1k_women_15_to_50_yo,
    
    pct_fb = pct_place_birth_tot_pop_foreign_born,
    moe_pct_fb = pct_moe_place_birth_tot_pop_foreign_born,
    
    pct_renter = pct_hous_tenure_occupied_hous_units_renter_occupied,
    moe_pct_renter = pct_moe_hous_tenure_occupied_hous_units_renter_occupied
    
  ) %>%
  dplyr::select(
    year, name, gtcbsa,
    median_age, moe_median_age,
    pct_white, moe_pct_white,
    # pct_hisp, moe_pct_hisp,
    per_k_birth_12mo, moe_per_k_birth_12mo,
    pct_col_ed, moe_pct_col_ed,
    unemploy_rate, moe_unemploy_rate,
    med_income, moe_med_income,
    pct_fb, moe_pct_fb, 
    pct_renter, moe_pct_renter
  )

vars = c(
  "median_age",
  "pct_white",
  # "pct_hisp",
  "per_k_birth_12mo",
  "pct_col_ed",
  "unemploy_rate",
  "med_income",
  "pct_fb",
  "pct_renter"
)

acs_new = acs_clean %>%
  arrange(gtcbsa, year) %>%
  group_by(gtcbsa) %>%
  
  # lag vars
  mutate(across(all_of(c(vars, paste0("moe_", vars))), lag, .names = "lag_{.col}")) %>%
  
  # find diffs
  mutate(across(all_of(vars), ~.x - get(paste0("lag_", cur_column())), .names = "diff_{.col}")) %>%
  
  # find SE for diffs
  mutate(across(all_of(vars), ~{
    moe_now = get(paste0("moe_", cur_column()))
    moe_lag = get(paste0("lag_moe_", cur_column()))
    sqrt(moe_now^2 + moe_lag^2) / 1.645
  }, .names = "se_diff_{.col}")) %>%
  
  # determine significance
  mutate(across(all_of(vars), ~as.factor(ifelse(
    abs(get(paste0("diff_", cur_column()))) > 1.96 * get(paste0("se_diff_", cur_column())), 
    1, 0)), .names = "sig_diff_{.col}")) %>%
  
  ungroup() %>%
  filter(year %in% years)

rm(acs_clean, acs)

### CEV ###

cev_new = cev %>% 
  mutate(
    employed = as.numeric(employed),
    social_capital = as.numeric(rowSums(
      across(c(yesfftalk, yesffissues, yesntalk, yesnissues, yesfavors, action), 
             ~ifelse(. == 1, 1, 0)),
      na.rm = FALSE)),
    high_social_capital = as.factor(ifelse(social_capital >= 3, 1, 0)),
    age = prtage,
    is_female = as.factor(female),
    race = white,
    college = as.factor(
      ifelse(simedu >= 4, 1, 0)),
    
    foreign_born = as.factor(foreign_born),
    # does not distinguish between native and naturalized citizens
    citizen = as.factor(
      ifelse(citizenship_status %in% c("Native", "Naturalized Citizen"), 1, 0)),
    
    year = hryear4) %>%
  # Remove missing values on key variables
  drop_na(topmsa, # drop rows not in top msa regions
          employed, high_social_capital, is_female, age, hispanic, race, 
          college, foreign_born, citizen) %>%
  dplyr::select(
    gtcbsa,
    
    employed, high_social_capital, social_capital,
    is_female, race, hispanic, college, 
    foreign_born, citizen, mother_foreign_born, father_foreign_born,
    
    yesfftalk, yesffissues, yesntalk, yesnissues, yesfavors, action,
    fftalk, ffissues, ntalk, nissues, favors,
    
    age, year, generation, # age, period, and cohort
    
    supwgt # weighting
  )
rm(cev)

### JOINING CEV AND ACS###
joined_df = cev_new %>% 
  left_join(acs_new, by = c("year", "gtcbsa")) %>%
  dplyr::select(-starts_with("lag"))




#### stuff #####

# Summary Statistics for Data Sources Section
library(tidyverse)
library(flextable)

# 1. Sample Size and Geographic Coverage
sample_summary <- joined_df %>%
  summarise(
    total_observations = n(),
    unique_metros = n_distinct(gtcbsa),
    years_covered = paste(sort(unique(year)), collapse = ", "),
    .groups = 'drop'
  )

cat("Sample Overview:\n")
cat("Total observations:", sample_summary$total_observations, "\n")
cat("Number of metropolitan areas:", sample_summary$unique_metros, "\n") 
cat("Years covered:", sample_summary$years_covered, "\n\n")

# 2. Sample by Year
year_summary <- joined_df %>%
  group_by(year) %>%
  summarise(
    n_respondents = n(),
    n_metros = n_distinct(gtcbsa),
    .groups = 'drop'
  )

print("Sample Distribution by Year:")
print(year_summary)

# 3. Demographic Characteristics
demo_vars <- joined_df %>%
  summarise(
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    pct_female = mean(as.numeric(is_female) - 1, na.rm = TRUE) * 100,
    pct_white = mean(race, na.rm = TRUE) * 100,
    pct_hispanic = mean(hispanic, na.rm = TRUE) * 100,
    pct_college = mean(as.numeric(college) - 1, na.rm = TRUE) * 100,
    pct_foreign_born = mean(as.numeric(foreign_born) - 1, na.rm = TRUE) * 100,
    pct_citizen = mean(as.numeric(citizen) - 1, na.rm = TRUE) * 100,
    pct_employed = mean(employed, na.rm = TRUE) * 100,
    .groups = 'drop'
  )

cat("\nIndividual-Level Demographic Characteristics:\n")
cat(sprintf("Age: Mean = %.1f years (SD = %.1f)\n", demo_vars$age_mean, demo_vars$age_sd))
cat(sprintf("Female: %.1f%%\n", demo_vars$pct_female))
cat(sprintf("White: %.1f%%\n", demo_vars$pct_white))
cat(sprintf("Hispanic: %.1f%%\n", demo_vars$pct_hispanic))
cat(sprintf("College educated: %.1f%%\n", demo_vars$pct_college))
cat(sprintf("Foreign-born: %.1f%%\n", demo_vars$pct_foreign_born))
cat(sprintf("U.S. citizen: %.1f%%\n", demo_vars$pct_citizen))
cat(sprintf("Employed: %.1f%%\n", demo_vars$pct_employed))

# 4. Social Capital Measures
social_capital_summary <- joined_df %>%
  summarise(
    social_capital_mean = mean(social_capital, na.rm = TRUE),
    social_capital_sd = sd(social_capital, na.rm = TRUE),
    high_social_capital_pct = mean(as.numeric(high_social_capital) - 1, na.rm = TRUE) * 100,
    
    # Individual components
    talk_family_friends_pct = mean(yesfftalk, na.rm = TRUE) * 100,
    discuss_issues_family_friends_pct = mean(yesffissues, na.rm = TRUE) * 100,
    talk_neighbors_pct = mean(yesntalk, na.rm = TRUE) * 100,
    discuss_issues_neighbors_pct = mean(yesnissues, na.rm = TRUE) * 100,
    do_favors_pct = mean(yesfavors, na.rm = TRUE) * 100,
    community_action_pct = mean(action, na.rm = TRUE) * 100,
    .groups = 'drop'
  )

cat("\nSocial Capital Measures:\n")
cat(sprintf("Social capital index: Mean = %.2f (SD = %.2f, Range: 0-6)\n", 
            social_capital_summary$social_capital_mean, social_capital_summary$social_capital_sd))
cat(sprintf("High social capital (≥3): %.1f%%\n", social_capital_summary$high_social_capital_pct))
cat("\nComponent behaviors:\n")
cat(sprintf("  Talk with family/friends: %.1f%%\n", social_capital_summary$talk_family_friends_pct))
cat(sprintf("  Discuss issues with family/friends: %.1f%%\n", social_capital_summary$discuss_issues_family_friends_pct))
cat(sprintf("  Talk with neighbors: %.1f%%\n", social_capital_summary$talk_neighbors_pct))
cat(sprintf("  Discuss issues with neighbors: %.1f%%\n", social_capital_summary$discuss_issues_neighbors_pct))
cat(sprintf("  Do favors for neighbors: %.1f%%\n", social_capital_summary$do_favors_pct))
cat(sprintf("  Participate in community action: %.1f%%\n", social_capital_summary$community_action_pct))

# 5. Metropolitan-Level Characteristics (ACS data)
metro_summary <- joined_df %>%
  group_by(gtcbsa, year) %>%
  slice(1) %>% # One observation per metro-year
  ungroup() %>%
  summarise(
    median_age_mean = mean(median_age, na.rm = TRUE),
    median_age_sd = sd(median_age, na.rm = TRUE),
    pct_white_mean = mean(pct_white, na.rm = TRUE),
    pct_white_sd = sd(pct_white, na.rm = TRUE),
    pct_foreign_born_mean = mean(pct_fb, na.rm = TRUE),
    pct_foreign_born_sd = sd(pct_fb, na.rm = TRUE),
    unemployment_rate_mean = mean(unemploy_rate, na.rm = TRUE),
    unemployment_rate_sd = sd(unemploy_rate, na.rm = TRUE),
    median_income_mean = mean(med_income, na.rm = TRUE) / 1000, # Convert to thousands
    median_income_sd = sd(med_income, na.rm = TRUE) / 1000,
    pct_college_mean = mean(pct_col_ed, na.rm = TRUE),
    pct_college_sd = sd(pct_col_ed, na.rm = TRUE),
    pct_renter_mean = mean(pct_renter, na.rm = TRUE),
    pct_renter_sd = sd(pct_renter, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nMetropolitan-Level Characteristics (ACS):\n")
cat(sprintf("Median age: Mean = %.1f years (SD = %.1f)\n", 
            metro_summary$median_age_mean, metro_summary$median_age_sd))
cat(sprintf("Percent white: Mean = %.1f%% (SD = %.1f)\n", 
            metro_summary$pct_white_mean, metro_summary$pct_white_sd))
cat(sprintf("Percent foreign-born: Mean = %.1f%% (SD = %.1f)\n", 
            metro_summary$pct_foreign_born_mean, metro_summary$pct_foreign_born_sd))
cat(sprintf("Unemployment rate: Mean = %.1f%% (SD = %.1f)\n", 
            metro_summary$unemployment_rate_mean, metro_summary$unemployment_rate_sd))
cat(sprintf("Median family income: Mean = $%.0fk (SD = $%.0fk)\n", 
            metro_summary$median_income_mean, metro_summary$median_income_sd))
cat(sprintf("Percent college educated: Mean = %.1f%% (SD = %.1f)\n", 
            metro_summary$pct_college_mean, metro_summary$pct_college_sd))
cat(sprintf("Percent renter-occupied: Mean = %.1f%% (SD = %.1f)\n", 
            metro_summary$pct_renter_mean, metro_summary$pct_renter_sd))

# 6. Missing Data Assessment
missing_summary <- joined_df %>%
  summarise(
    across(c(age, employed, high_social_capital, is_female, race, hispanic, 
             college, foreign_born, citizen, median_age, pct_white, pct_fb, 
             unemploy_rate, med_income, pct_col_ed, pct_renter),
           ~sum(is.na(.)) / n() * 100, .names = "missing_{.col}")
  )

cat("\nMissing Data Summary (% missing):\n")
missing_long <- missing_summary %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  mutate(variable = str_remove(variable, "missing_")) %>%
  arrange(desc(pct_missing))

for(i in 1:nrow(missing_long)) {
  if(missing_long$pct_missing[i] > 0) {
    cat(sprintf("%s: %.1f%%\n", missing_long$variable[i], missing_long$pct_missing[i]))
  }
}

# 7. Create a formal table for the paper
create_summary_table <- function() {
  # Individual characteristics
  indiv_stats <- joined_df %>%
    summarise(
      Variable = c("Age (years)", "Female (%)", "White (%)", "Hispanic (%)", 
                   "College educated (%)", "Foreign-born (%)", "U.S. citizen (%)", 
                   "Employed (%)", "Social capital index (0-6)", "High social capital (%)"),
      `Mean/Percent` = c(
        sprintf("%.1f", mean(age, na.rm = TRUE)),
        sprintf("%.1f", mean(as.numeric(is_female) - 1, na.rm = TRUE) * 100),
        sprintf("%.1f", mean(race, na.rm = TRUE) * 100),
        sprintf("%.1f", mean(hispanic, na.rm = TRUE) * 100),
        sprintf("%.1f", mean(as.numeric(college) - 1, na.rm = TRUE) * 100),
        sprintf("%.1f", mean(as.numeric(foreign_born) - 1, na.rm = TRUE) * 100),
        sprintf("%.1f", mean(as.numeric(citizen) - 1, na.rm = TRUE) * 100),
        sprintf("%.1f", mean(employed, na.rm = TRUE) * 100),
        sprintf("%.2f", mean(social_capital, na.rm = TRUE)),
        sprintf("%.1f", mean(as.numeric(high_social_capital) - 1, na.rm = TRUE) * 100)
      ),
      `SD/Range` = c(
        sprintf("%.1f", sd(age, na.rm = TRUE)),
        sprintf("%d-%d", min(as.numeric(is_female) - 1), max(as.numeric(is_female) - 1)),
        sprintf("%d-%d", 0, 1),
        sprintf("%d-%d", 0, 1),
        sprintf("%d-%d", 0, 1),
        sprintf("%d-%d", 0, 1),
        sprintf("%d-%d", 0, 1),
        sprintf("%d-%d", 0, 1),
        sprintf("%.2f", sd(social_capital, na.rm = TRUE)),
        sprintf("%d-%d", 0, 1)
      )
    )
  
  # Metropolitan characteristics
  metro_stats <- joined_df %>%
    group_by(gtcbsa, year) %>%
    slice(1) %>%
    ungroup() %>%
    summarise(
      Variable = c("Median age (years)", "White (%)", "Foreign-born (%)", 
                   "Unemployment rate (%)", "Median income ($1000s)", 
                   "College educated (%)", "Renter-occupied (%)"),
      `Mean/Percent` = c(
        sprintf("%.1f", mean(median_age, na.rm = TRUE)),
        sprintf("%.1f", mean(pct_white, na.rm = TRUE)),
        sprintf("%.1f", mean(pct_fb, na.rm = TRUE)),
        sprintf("%.1f", mean(unemploy_rate, na.rm = TRUE)),
        sprintf("%.0f", mean(med_income, na.rm = TRUE) / 1000),
        sprintf("%.1f", mean(pct_col_ed, na.rm = TRUE)),
        sprintf("%.1f", mean(pct_renter, na.rm = TRUE))
      ),
      `SD/Range` = c(
        sprintf("%.1f", sd(median_age, na.rm = TRUE)),
        sprintf("%.1f", sd(pct_white, na.rm = TRUE)),
        sprintf("%.1f", sd(pct_fb, na.rm = TRUE)),
        sprintf("%.1f", sd(unemploy_rate, na.rm = TRUE)),
        sprintf("%.0f", sd(med_income, na.rm = TRUE) / 1000),
        sprintf("%.1f", sd(pct_col_ed, na.rm = TRUE)),
        sprintf("%.1f", sd(pct_renter, na.rm = TRUE))
      )
    )
  
  # Combine and create flextable
  combined_stats <- bind_rows(
    indiv_stats %>% mutate(Category = "Individual Characteristics"),
    metro_stats %>% mutate(Category = "Metropolitan Characteristics")
  ) %>%
    dplyr::select(Category, Variable, `Mean/Percent`, `SD/Range`)
  
  return(combined_stats)
}

# Generate the table
summary_table <- create_summary_table()
print(summary_table)

ft <- flextable(summary_table) %>%
  merge_v(j = "Category") %>%
  theme_vanilla() %>%
  autofit()
print(ft)



# Panel Structure Statistics Calculator
# Calculate the key numbers for your thesis sentence

# 1. Total metropolitan area-year observations (X)
total_metro_year_obs <- acs_new %>% 
  filter(year %in% years) %>%
  nrow()

# 2. Metro-year observations with significant increases in foreign-born (Y)
sig_increase_obs <- acs_new %>% 
  filter(year %in% years, 
         sig_diff_pct_fb == 1) %>%
  nrow()

# 3. Percentage of observations with significant increases (Z%)
pct_obs_with_increase <- round((sig_increase_obs / total_metro_year_obs) * 100, 1)

# 4. Individual respondents in areas with significant increases (A%)
# First, get total individual respondents
total_respondents <- joined_df %>% nrow()

# Respondents in areas with significant increases
respondents_in_increase_areas <- joined_df %>%
  filter(sig_diff_pct_fb == 1) %>%
  nrow()

# Percentage of respondents in areas with increases
pct_respondents_in_increase_areas <- round((respondents_in_increase_areas / total_respondents) * 100, 1)

# Create summary table
panel_summary <- data.frame(
  Statistic = c("Total Metro-Year Observations (X)",
                "Metro-Years with Significant FB Increases (Y)", 
                "Percentage of Observations with Increases (Z%)",
                "Percentage of Respondents in Increase Areas (A%)"),
  Value = c(total_metro_year_obs,
            sig_increase_obs,
            paste0(pct_obs_with_increase, "%"),
            paste0(pct_respondents_in_increase_areas, "%"))
)

# Display results
cat("PANEL STRUCTURE STATISTICS\n")


for(i in 1:nrow(panel_summary)) {
  cat(paste0(panel_summary$Statistic[i], ": ", panel_summary$Value[i], "\n"))
}

cat("FILLED SENTENCE FOR YOUR THESIS:\n\n")

filled_sentence <- paste0(
  "The panel structure creates ", total_metro_year_obs,
  " total metropolitan area-year observations, of which ", 
  sig_increase_obs, " experienced statistically significant increases in foreign-born population ",
  "compared to the previous survey year, representing ", pct_obs_with_increase,
  "% of observations and encompassing ", pct_respondents_in_increase_areas, 
  "% of individual respondents."
)

cat(filled_sentence)

# Additional descriptive statistics for context

cat("ADDITIONAL CONTEXT:\n\n")

cat("Years included:", paste(years, collapse = ", "), "\n")
cat("Metropolitan areas included:", length(unique(acs_new$gtcbsa)), "\n")
cat("Total individual respondents:", total_respondents, "\n")

# Show distribution by year
year_breakdown <- acs_new %>%
  filter(year %in% years) %>%
  group_by(year) %>%
  summarise(
    total_metros = n(),
    sig_increases = sum(sig_diff_pct_fb == 1, na.rm = TRUE),
    pct_with_increases = round((sig_increases/total_metros)*100, 1)
  )

cat("\nBreakdown by year:\n")
print(year_breakdown)

# Show which metros had increases and when
metros_with_increases <- acs_new %>%
  filter(year %in% years, sig_diff_pct_fb == 1) %>%
  dplyr::select(year, name, gtcbsa, diff_pct_fb, se_diff_pct_fb) %>%
  arrange(year, name)

cat("\nMetropolitan areas with significant increases:\n")
print(metros_with_increases)

# Return values for use in other parts of analysis
return_values <- list(
  total_metro_year_obs = total_metro_year_obs,
  sig_increase_obs = sig_increase_obs,
  pct_obs_with_increase = pct_obs_with_increase,
  pct_respondents_in_increase_areas = pct_respondents_in_increase_areas,
  filled_sentence = filled_sentence,
  panel_summary = panel_summary,
  year_breakdown = year_breakdown,
  metros_with_increases = metros_with_increases
)

# Save to file for later reference
write.csv(panel_summary, "panel_structure_summary.csv", row.names = FALSE)
write.csv(year_breakdown, "panel_year_breakdown.csv", row.names = FALSE)
write.csv(metros_with_increases, "metros_with_increases.csv", row.names = FALSE)

cat("\n\nFiles saved: panel_structure_summary.csv, panel_year_breakdown.csv, metros_with_increases.csv\n")



### alpha


library(psych)

# Select your social capital items
social_capital_items <- cev_new %>%
  dplyr::select(yesfftalk, yesffissues, yesntalk, yesnissues, yesfavors, action) %>%
  # Convert to numeric if needed
  mutate(across(everything(), ~as.numeric(.))) %>%
  # Remove rows with any missing values
  drop_na()

# Calculate Cronbach's alpha
cronbach_result <- psych::alpha(social_capital_items)

# Print the result
print(cronbach_result)

# Just get the alpha coefficient
cronbach_alpha <- cronbach_result$total$raw_alpha
cat("Cronbach's alpha:", round(cronbach_alpha, 3))


### table ###

# 1) Load packages
library(dplyr)
library(knitr)
library(kableExtra)

# 2) Define your table rows: variable name, label, grouping, and scale factor
table_rows <- tibble::tribble(
  ~var,                  ~label,~group, ~scale,
  "social_capital",      "Social Capital (Continuous)","Outcome Variables", 1,
  "high_social_capital", "High Social Capital (Binary)","Outcome Variables", 1,
  
  "employed",            "Employment Status", "Individual Characteristics",  1,
  "age",                 "Age (years)", "Individual Characteristics",    1,
  "is_female",           "Female", "Individual Characteristics",  1,
  "race",                "White", "Individual Characteristics",  1,
  "hispanic",            "Hispanic", "Individual Characteristics",  1,
  "college",             "College Education", "Individual Characteristics",  1,
  "foreign_born",        "Foreign‑Born", "Individual Characteristics",  1,
  "citizen",             "U.S. Citizen", "Individual Characteristics",  1,
  "mother_foreign_born", "Mother Foreign‑Born", "Individual Characteristics",  1,
  "father_foreign_born", "Father Foreign‑Born", "Individual Characteristics",  1,
  
  "pct_white",           "Percent White (%)",               "Metropolitan Characteristics", 1,
  "median_age",          "Median Age (years)",              "Metropolitan Characteristics", 1,
  "med_income",          "Median Income ($1{,}000s)",       "Metropolitan Characteristics", 1/1000,
  "pct_renter",          "Percent Renters (%)",             "Metropolitan Characteristics", 1,
  "unemploy_rate",       "Unemployment Rate (%)",           "Metropolitan Characteristics", 1,
  "per_k_birth_12mo",    "Births per 1,000 Women",          "Metropolitan Characteristics", 1
)

# 3) Helper to compute mean and SD on the scaled variable
compute_stats <- function(df, var, scale) {
  x <- as.numeric(df[[var]]) * scale
  c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
}

# 4) Helper to compute SMD for a scaled variable
compute_smd <- function(df, var, treat, scale) {
  x0 <- as.numeric(df[[var]][df[[treat]]==0]) * scale
  x1 <- as.numeric(df[[var]][df[[treat]]==1]) * scale
  m0 <- mean(x0, na.rm=TRUE); m1 <- mean(x1, na.rm=TRUE)
  s0 <- sd(x0, na.rm=TRUE);  s1 <- sd(x1, na.rm=TRUE)
  n0 <- sum(!is.na(x0));     n1 <- sum(!is.na(x1))
  pooled_sd <- sqrt(((n0-1)*s0^2 + (n1-1)*s1^2)/(n0+n1-2))
  (m1 - m0)/pooled_sd
}

# 5) Build the table
table_data <- table_rows %>% rowwise() %>%
  mutate(
    # stats for control (treat=0) and treated (treat=1)
    stats0 = list(compute_stats(joined_df %>% filter(sig_diff_pct_fb==0), var, scale)),
    stats1 = list(compute_stats(joined_df %>% filter(sig_diff_pct_fb==1), var, scale)),
    mean_control = stats0["mean"],
    sd_control   = stats0["sd"],
    mean_treated = stats1["mean"],
    sd_treated   = stats1["sd"],
    smd          = compute_smd(joined_df, var, "sig_diff_pct_fb", scale),
    Imbalance    = ifelse(abs(smd) > 0.1, "Yes", "No")
  ) %>% 
  ungroup() %>%
  # Round nicely
  mutate(across(c(mean_control:smd), ~round(.x, ifelse(group=="Outcome Variables" & var=="social_capital",2, ifelse(grepl("%", label),1,1)) ))) %>%
  # Format SD columns: if var is the binary outcomes, leave SD blank
  mutate(
    sd_control = ifelse(var %in% c("high_social_capital", "is_female","race","hispanic","college","foreign_born","citizen","mother_foreign_born","father_foreign_born"),
                        NA, sd_control),
    sd_treated = ifelse(var %in% c("high_social_capital", "is_female","race","hispanic","college","foreign_born","citizen","mother_foreign_born","father_foreign_born"),
                        NA, sd_treated)
  ) %>%
  dplyr::select(group, label, mean_control, sd_control, mean_treated, sd_treated, smd, Imbalance)

# 6) Print with kableExtra
table_data %>%
  arrange(factor(group, levels=c(
    "Outcome Variables",
    "Individual Characteristics",
    "Metropolitan Characteristics"
  ))) %>%
  dplyr::select(-group) %>%
  kbl(
    caption = "Table 2: Descriptive Statistics and Balance Assessment",
    booktabs = TRUE,
    align = c("l","r","r","r","r","r","c")
  ) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, "Stable Diversity Change" = 2, "High Diversity Change" = 2, " " = 2)) %>%
  pack_rows("Outcome Variables",              1,  2) %>%
  pack_rows("Individual Characteristics",     3, 12) %>%
  pack_rows("Metropolitan Characteristics",  13, 18)

write.csv(table_data, "table2_desc_balance.csv", row.names = FALSE)
