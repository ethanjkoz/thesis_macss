library(tidyverse)
library(readr)
library(stringr)
library(janitor)
library(arrow)

acs_folder = "D:/Thesis/data/acs_data"

# Read metadata files
metadata_df = list.files(acs_folder, pattern = "Column-Metadata.csv$", full.names = TRUE) %>%
  tibble(file_path = .) %>%
  mutate(
    dp_category = str_extract(basename(file_path), "DP\\d{2}"),
    year = as.integer(str_extract(basename(file_path), "\\d{4}")),
    columns = map(file_path, ~ read_csv(.x, skip = 2, show_col_types = FALSE) %>%
                    transmute(
                      name = tolower(str_trim(NAME)),
                      desc = str_trim(`Geographic Area Name`)
                    ))
  ) %>%
  select(year, dp_category, columns) %>%
  unnest(columns) %>%
  distinct(year, dp_category, name, desc) %>%
  mutate(
    slug_base = desc %>%
      str_to_lower() %>%
      str_remove("^estimate!!") %>%
      str_replace_all(c(
        "(excluding units where smocapi cannot be computed)" = "w_smocapi_units_only",
        "in 20\\d{2} inflation-adjusted" = "infl_adj",
        "selected monthly owner costs" = "smoc",
        "selected characteristics" = "sel_chars",
        "noninstitutionalized" = "non_int",
        "not in labor force" = "nilf",
        "margin of error" = "moe",
        "high school" = "hs",
        "united states" = "us",
        "and over" = "up",
        "or higher" = "up",
        "years old" = "yo",
        "without" = "w_o",
        "supplemental" = "supp",
        "naturalized" = "natzd",
        "relationship" = "relship",
        "educational" = "edu",
        "attainment" = "attain",
        "government" = "gov",
        "household" = "hshld",
        "population" = "pop",
        "different" = "diff",
        "structure" = "struct",
        "civilian" = "civ",
        "number" = "num",
        "language" = "lang",
        "english" = "eng",
        "per 1,000" = "per_1k",
        "spanish" = "span",
        "children" = "kids",
        "grand" = "g",
        "one or more" = "oneplus",
        "child" = "kid",
        "less than very well" = "less_v_well",
        "asian and pacific islander" = "aapi",
        "residence" = "res",
        "families" = "fams",
        "average" = "avg",
        "veteran" = "vet",
        "graduate" = "grad",
        "poverty" = "pov",
        "housing" = "hous",
        "dollars" = "usd",
        "percent" = "pct",
        "people" = "ppl",
        "family" = "fam",
        "hispanic or latino" = "hisp_lat",
        "(of any race)" = "",
        "(excluding Hispanic origin groups)" = "",
        "total" = "tot",
        "year" = "yr",
        "with" = "w",
        "\\b(one)\\b" = "1",
        "\\b(two)\\b" = "2",
        "\\b(three)\\b" = "3",
        "widowed, divorced, and never married" = "",
        "\\b(who|a|an|the|of)\\b" = "",
        "\\s+" = "_",
        "-" = "_",
        "!!" = "_",
        "tached" = ""
      )) %>%
      str_remove_all("[^A-Za-z0-9_\\.]") %>%
      str_replace_all("_+", "_") %>% 
      str_replace_all("^|_$", "")
  ) %>%
  group_by(year, dp_category, slug_base) %>%
  mutate(slug = if(n() > 1) paste0(slug_base, "_", name) else slug_base) %>%
  ungroup() %>%
  select(year, dp_category, name, desc, slug)

write_csv(metadata_df, file.path(acs_folder, "metadata.csv"))

# Read ACS data function
read_acs_csv = function(file_path, year, dp_category) {
  yr = as.integer(year)
  
  # Read column names from first row
  codes = read_csv(file_path, n_max = 1, col_names = FALSE, show_col_types = FALSE) %>%
    as.character() %>%
    tolower()
  
  # Read data
  df = read_csv(file_path, col_names = codes, skip = 2, show_col_types = FALSE)
  names(df)[1:2] = c("geo_id", "name")
  
  # Get rename mapping
  rename_map = metadata_df %>%
    filter(year == yr, dp_category == dp_category, name %in% names(df)) %>%
    select(name, slug) %>%
    deframe()
  
  # Apply renames
  if(length(rename_map) > 0) {
    col_idx = match(names(rename_map), names(df))
    names(df)[col_idx] = rename_map
  }
  
  df %>%
    clean_names() %>%
    mutate(
      across(-c(geo_id, name), as.numeric),
      year = yr
    ) %>%
    select(geo_id, name, year, everything())
}

# Process all files
final_acs_data = list.files(acs_folder, pattern = "Data\\.csv$", full.names = TRUE) %>%
  tibble(file_path = .) %>%
  mutate(
    year = str_extract(basename(file_path), "\\d{4}"),
    dp_category = str_extract(basename(file_path), "DP\\d{2}"),
    data = pmap(list(file_path, year, dp_category), read_acs_csv)
  ) %>%
  group_by(year) %>%
  summarize(df = list(reduce(data, full_join, by = c("geo_id", "name", "year"))), .groups = "drop") %>%
  pull(df) %>%
  bind_rows() %>%
  mutate(gtcbsa = as.numeric(str_extract(geo_id, "\\d{5}$")))

write_feather(final_acs_data, file.path(acs_folder, "acs.feather"))