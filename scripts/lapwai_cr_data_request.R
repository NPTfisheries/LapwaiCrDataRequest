# -----------------------
# Author(s): Mike Ackerman
# Purpose: Query and compile relevant data and results for a data request regarding a 
#  biological assessment for the McIntyre Bridge replacement on Lapwai Creek.
# 
# Created Date: December 30, 2025
#   Last Modified:
#
# Notes: 

# clean environment
rm(list = ls())

# load libraries
library(tidyverse)
library(readxl)
library(writexl)
library(PITcleanr)
library(PITmodelR)
library(janitor)

# iptds metadata
iptds_meta = queryInterrogationMeta() %>%
  filter(siteCode %in% c("LAP", "MIS", "SWT", "WEB")) %>%
  select(-operationsOrganizationCode,
         -lastFileName,
         -lastFileOpenedOn,
         -lastFileClosedOn,
         -lastFileProcessedOn,
         -lastTagCode,
         -lastTagTime)

# dabom site escapements
iptds_escape_df = list.files(path = "C:/Git/SnakeRiverFishStatus/output/syntheses",
                            pattern = "\\.xlsx$",
                            full.names = T) %>%
  discard(~ basename(.x) |> stringr::str_starts("~\\$")) %>%   # discard any tmp files
  map_dfr(~ read_excel(.x, sheet = "Site_Esc")) %>%
  filter(site %in% c("LAP", "MIS", "SWT", "WEB")) %>%
  select(-mean, -mode, -sd, -notes) %>%
  mutate(
    across(
      c(median, lower95ci, upper95ci),
      ~ round(.x, 0)
    ),
    cv = round(cv, 3)
  ) %>%
  mutate(
    species = recode(species, Chinook = "Sp/Sum Chinook")
  ) %>%
  mutate(origin = case_when(
    species %in% c("Sp/Sum Chinook", "Steelhead") ~ "Natural",
    species == "Coho"                             ~ "Natural + Hatchery"
  )) %>%
  select(species, origin, spawn_yr, everything())

# individual dabom detections
dabom_obs_df = list.files(path = "C:/Git/SnakeRiverFishStatus/output/life_history",
                         pattern = "\\.xlsx$",
                         full.names = T) %>%
  discard(~ basename(.x) |> stringr::str_starts("~\\$")) %>%   # discard any tmp files
  map_dfr(~ read_excel(.x, sheet = "tag_lh")) %>%
  filter(spawn_site %in% c("LAP", "MIS", "SWT", "WEB")) %>%
  mutate(
    sex = coalesce(GenSex, LGDSex)
  ) %>%
  select(species,
         spawn_year,
         tag_code,
         spawn_site,
         min_det,
         SRR,
         lgr_trap_date = CollectionDate,
         LGDFLmm,
         sex,
         path,
         fw_age,
         sw_age,
         total_age,
         brood_year)

# read in all ptagis observations at lapwai iptds, 2010-2025
lap_int_df = read_csv(file = "data/lapwai_interrogation_summary_ptagis_20251230.csv") %>%
  clean_names() %>%
  mutate(
    site_code = str_sub(site_name, 1, 3),
    srr = str_c(
      as.integer(species_code),
      as.integer(run_code),
      rear_type_code,
      sep = ""
    ),
    first_time_value = mdy_hms(first_time_value, tz = "America/Los_Angeles"),
    last_time_value  = mdy_hms(last_time_value, tz = "America/Los_Angeles"),
    year = year(first_time_value)
  ) %>%
  select(tag_code,
         site_code,
         first_time_value,
         year,
         srr,
         species_name,
         run_name,
         rear_type_name)

table(lap_int_df$year, lap_int_df$species_name)  
table(lap_int_df$year, lap_int_df$srr)

# retrieve complete tag histories
lap_cth = lap_int_df %>%
  distinct(tag_code) %>%
  pull(tag_code) %>%
  get_batch_tag_histories()



# run timing, by year
run_df = dabom_obs_df %>%
  mutate(julian = yday(min_det)) %>%
  group_by(species, spawn_site, spawn_year, julian) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species, spawn_year, spawn_site) %>%
  complete(julian = 1:366, fill = list(n = 0)) %>%
  arrange(julian, .by_group = TRUE) %>%
  mutate(
    csum = cumsum(n),
    cdf = csum / sum(n)
  ) %>%
  ungroup()

# average run timing
avg_run_df = run_df %>%
  group_by(species, spawn_site, julian) %>%
  summarise(cdf = mean(cdf), .groups = "drop")

# plot run timing
run_p = ggplot() +
  geom_line(
    data = run_df,
    aes(x = julian, y = cdf, group = spawn_year),
    linewidth = 0.4,
    linetype = "dashed",
    color = "grey60"
  ) +
  geom_line(
    data = avg_run_df,
    aes(x = julian, y = cdf),
    linewidth = 0.9,
    color = "black"
  ) +
  #scale_x_date(date_breaks = "4 weeks", date_labels = "%b-%d") +
  facet_grid(spawn_site ~ species, scales = "free_y") +
  theme_bw() +
  labs(x = "Julian Day",
       y = "CDF")
run_p

# save run timing plot
ggsave("output/figures/lapwai_iptds_run_timing.pdf")

# write some objects to excel
list("IPTDS_Metadata"   = iptds_meta,
     "IPTDS_Escapement" = iptds_escape_df,
     "IPTDS_Detections" = dabom_obs_df) %>%
  write_xlsx("output/lapwai_data_request_20251230.xlsx")

### END SCRIPT