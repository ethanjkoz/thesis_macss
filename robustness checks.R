#### INITIALIZATION ####
library(tidyverse)
library(arrow)
library(MatchIt)
library(cem)
library(ggplot2)
library(mediation)
library(cobalt)

set.seed(0)
acs <- read_feather("D:/Thesis/data/acs_data/acs.feather")
cev <- read_feather("D:/Thesis/data/cps_cev_data/cps_data.feather")
years <- c(2017, 2019, 2021, 2023)

topmsa_gtcbsas <- cev %>% drop_na(topmsa) %>% pull(gtcbsa) %>% unique()


library(tidyverse)
library(arrow)
library(haven)
library(flextable)
library(apc)
library(cobalt)
library(gridExtra)
library(mediation)


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

# treatment
D = "sig_diff_pct_fb"
# mediator
M = "employed"

# covars
C = c(
  # age, period, cohort effects
  "year",
  # "generation"
  "age",
  
  # individual lvl covars
  "is_female",
  "foreign_born",
  "race",
  "hispanic",
  "college",
  "citizen",
  "mother_foreign_born",
  "father_foreign_born",
  
  # county lvl covars
  "pct_white",
  "median_age",
  "med_income",
  "pct_renter",
  "unemploy_rate",
  "per_k_birth_12mo"
)

# outcome
Y_bin= "high_social_capital"
Y_cont = "social_capital"

analysis_df = joined_df %>%
  mutate(
    across(all_of(c(D, M, C, Y_bin, Y_cont)),
           ~as.numeric(as.character(.)))
  ) %>%
  drop_na(all_of(c(D, M, C, Y_bin, Y_cont)))



#### 3. PROPENSITY SCORE MATCHING ####
ps_formula <- as.formula(paste(D, "~", paste(C, collapse = " + ")))
ps_model <- glm(ps_formula, data = analysis_df, family = binomial)
analysis_df$pscore <- predict(ps_model, type = "response")

ps_match <- matchit(ps_formula, data = analysis_df, method = "nearest", distance = "logit", ratio = 1)
ps_matched <- match.data(ps_match)
print(summary(ps_match))

ps_plot <- ggplot(analysis_df, aes(x = pscore, fill = factor(get(D)))) +
  geom_density(alpha = 0.7) +
  labs(title = "Propensity Score Distribution by Treatment Status",
       x = "Propensity Score", 
       fill = "Treatment") +
  theme_minimal()

print(ps_plot)

#### 4. COARSENED EXACT MATCHING ####
run_cem <- function(data, treat, vars) {
  df <- data %>%
    dplyr::select(all_of(c(treat, vars))) %>%
    mutate(across(everything(), as.numeric)) %>%
    drop_na()
  names(df)[1] <- "treatment"
  df[paste0(vars,"_c")] <- lapply(df[vars], function(x) as.numeric(cut(x,4,labels=FALSE)))
  cem_res <- cem("treatment", data = df[c("treatment",paste0(vars,"_c"))], drop = "treatment")
  matched_ids <- as.numeric(rownames(df))[cem_res$w > 0]
  list(summary = cem_res, data = analysis_df[matched_ids, ])
}
cem_out <- run_cem(analysis_df, D, c("med_income","unemploy_rate","median_age","pct_white"))
print(cem_out$summary)


#### 5. ALTERNATIVE MEDIATOR ANALYSIS ####

alternative_mediators <- c("college", "foreign_born", "med_income",
                           "pct_white", "pct_renter", "unemploy_rate")

#### 0. Prep: make sure everything is numeric ####
# (linmed() requires D, M, Y, and C all to be numeric)
all_vars <- unique(c(D, Y_bin, Y_cont, alternative_mediators, C))
analysis_df <- analysis_df %>%
  mutate(across(all_of(all_vars),
                ~ as.numeric(as.character(.))))

#### 1. A small wrapper around linmed() ####
run_linmed_summary <- function(data, D, M, Y, C,
                               d = 1, dstar = 0,
                               boot = TRUE, boot_reps = 1000,
                               boot_seed = 1234, boot_parallel = TRUE) {
  res <- linmed(
    data          = data,
    D             = D,
    M             = M,
    Y             = Y,
    C             = C,
    d             = d,
    dstar         = dstar,
    boot          = boot,
    boot_reps     = boot_reps,
    boot_seed     = boot_seed,
    boot_parallel = boot_parallel
  )
  
  # Extract the point estimates + bootstrap CIs from the returned list
  tibble::tibble(
    Mediator      = M,
    Outcome       = Y,
    ATE           = res$ATE,
    ATE_CI        = paste0("[", round(res$ci_ATE[1],3), ", ", round(res$ci_ATE[2],3), "]"),
    NDE           = res$NDE,
    NDE_CI        = paste0("[", round(res$ci_NDE[1],3), ", ", round(res$ci_NDE[2],3), "]"),
    NIE           = res$NIE,
    NIE_CI        = paste0("[", round(res$ci_NIE[1],3), ", ", round(res$ci_NIE[2],3), "]"),
    Prop_Mediated = res$NIE / res$ATE
  )
}

#### 2. Loop over mediators × outcomes ####

robustness_results <- bind_rows(
  lapply(c(Y_bin, Y_cont), function(Yv) {
    bind_rows(lapply(alternative_mediators, function(Mv) {
      if (!Mv %in% names(analysis_df)) return(NULL)
      message("Running linmed for mediator = ", Mv, " | outcome = ", Yv)
      run_linmed_summary(
        data      = analysis_df,
        D         = D,
        M         = Mv,
        Y         = Yv,
        C         = C,
        boot_reps = boot,
        boot_seed = random_seed,
        boot_parallel = TRUE
      )
    }))
  })
)

#### 3. View your final table ####
print(robustness_results)
write.csv(robustness_results, "robustness_results.csv", row.names = FALSE)


#### 6. ALTERNATIVE OUTCOME ANALYSIS ####
social_capital_components <- c("yesfftalk","yesffissues","yesntalk","yesnissues","yesfavors","action")
run_outcome <- function(data, treat, outcome, mediator, covars) {
  lm(as.formula(paste(outcome, "~", treat, "+", mediator, "+", paste(covars, collapse = " + "))), data)
}
component_results <- lapply(social_capital_components, function(out) {
  cat("Testing", out, "as outcome\n")
  summary(run_outcome(analysis_df, D, out, M, C))
})
names(component_results) <- social_capital_components

#### 7. TEMPORAL ROBUSTNESS ####
pre_covid <- analysis_df %>% filter(year <= 2019)
covid_after <- analysis_df %>% filter(year >= 2021)

run_temporal <- function(data, treat, outcome, mediator, covars, period) {
  if (nrow(data) < 100) return(NULL)
  model <- lm(as.formula(paste(outcome, "~", treat, "+", mediator, "+", paste(covars, collapse = " + "))), data)
  list(period = period, coef = coef(model)[treat], se = summary(model)$coefficients[treat, "Std. Error"], n = nrow(data))
}
temporal_results <- list(
  pre_covid = run_temporal(pre_covid, D, Y_cont, M, C, "Pre-COVID"),
  covid_after = run_temporal(covid_after, D, Y_cont, M, C, "COVID and After")
)
for (res in temporal_results) if (!is.null(res)) cat(res$period, ": Coef =", round(res$coef,4), "SE =", round(res$se,4), "N =", res$n, "\n")

#### 8. BALANCE ASSESSMENT ####
calculate_balance <- function(data, treat, covars) {
  bind_rows(lapply(covars, function(v) {
    if (!v %in% names(data)) return(NULL)
    control <- data[data[[treat]] == 0, v]; treated <- data[data[[treat]] == 1, v]
    control <- control[!is.na(control)]; treated <- treated[!is.na(treated)]
    if (length(control) == 0 | length(treated) == 0) return(NULL)
    mean_c <- mean(control); mean_t <- mean(treated)
    sd_c <- sd(control); sd_t <- sd(treated)
    pooled_sd <- sqrt(((length(control)-1)*sd_c^2 + (length(treated)-1)*sd_t^2)/(length(control)+length(treated)-2))
    smd <- ifelse(pooled_sd == 0,0,(mean_t-mean_c)/pooled_sd)
    data.frame(Variable=v, Mean_Control=mean_c, Mean_Treated=mean_t, SMD=smd, Imbalanced=abs(smd)>0.1)
  }))
}
balance_table <- calculate_balance(analysis_df, D, C)
imbalanced <- balance_table %>% filter(Imbalanced)
cat("Balance Check: ", nrow(imbalanced), "imbalanced of", nrow(balance_table), "\n")
if (!is.null(cem_out$data)) {
  cat("After Matching:\n")
  balance_matched <- calculate_balance(cem_out$data, D, C)
  cat(sum(balance_matched$Imbalanced),"imbalanced of",nrow(balance_matched),"\n")
}
