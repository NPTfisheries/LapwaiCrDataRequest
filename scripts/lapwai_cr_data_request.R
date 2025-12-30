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
library(PITcleanr)
library(readxl)
library(writexl)
#library(PITmodelR)

# lapwai iptds metadata
iptds_meta = queryInterrogationMeta() %>%
  filter(siteCode %in% c("LAP", "MIS", "SWT", "WEB")) %>%
  select(-operationsOrganizationCode,
         -starts_with("last"),
         lastYear,
         lastDate)

# lapwai dabom adult site escapement estimates
dabom_escape_df = list.files(path = "C:/Git/SnakeRiverFishStatus/output/syntheses",
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

# summarize run timing, by year
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
ggsave("output/figures/natural_adult_run_timing_20251230.pdf")

#---------------------------
# raw lapwai interrogations

# read in all ptagis observations at lapwai iptds queried from ptagis, 2010-2025
lap_int_df = read_csv(file = "data/lapwai_interrogation_summary_ptagis_20251230.csv") %>%
  janitor::clean_names() %>%
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
         rear_type_code)

# retrieve complete tag histories
# lap_cth = lap_int_df %>%
#   distinct(tag_code) %>%
#   pull(tag_code) %>%
#   get_batch_tag_histories()

# save complete tag histories
# save(lap_cth, file = "output/lapwai_cths_20251230.rda")

# load complete tag histories
load("output/lapwai_cths_20251230.rda")

# attach mark info from complete tag histories to detections
lap_df = lap_int_df %>%
  left_join(lap_cth %>%
              filter(event_type == "Mark") %>%
              select(mark_date = event_date,
                     mark_site_code = site_code,
                     mark_site_name = site_name,
                     tag_code)) %>%
  # remove some species not of interest (all unknowns are orphan tags); single likely errant green sturgeon observation
  filter(!species_name %in% c("Northern Pikeminnow", "Unknown", "Green Sturgeon")) %>%
  # fix likely erroneous run types for steelhead
  mutate(
    run_name = if_else(species_name == "Steelhead", "Summer", run_name),
    srr = if_else(species_name == "Steelhead", str_replace(srr, "^..", "32"), srr)
  ) %>%
  # remove few records of unknown run chinook
  filter(!(species_name == "Chinook" & run_name == "Unknown")) %>%
  # column to differentiate sp/sum vs. fall run Chinook (for plotting)
  mutate(spc = case_when(
    species_name == "Steelhead"                                     ~ "Summer Steelhead",
    species_name == "Chinook" & run_name %in% c("Spring", "Summer") ~ "Sp/Sum Chinook Salmon",
    species_name == "Chinook" & run_name == "Fall"                  ~ "Fall Chinook Salmon",
    species_name == "Coho"                                          ~ "Coho Salmon",
    TRUE ~ species_name
  )) %>%
  mutate(
    julian = yday(first_time_value),
    xday = as.Date(julian - 1, origin = "2001-01-01") # dummy year for plotting
  )

table(lap_df$species_name, lap_df$run_name)
table(lap_df$year, lap_df$spc)

# jittered detections
lap_obs_p = lap_df %>%
  # just most recent 5 years
  filter(year >= 2021) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = xday, y = year, color = rear_type_code)) +
  geom_jitter(height = 0.25, alpha = 1, size = 0.8) +
  facet_grid(site_code ~ spc, scales = "free_y", space = "free_y") +
  scale_x_date(date_breaks = "3 month", date_labels = "%b") +
  theme_bw() +
  labs(title = "Timing of Detections, All Life Stages, By Species, 2021-2025",
       x = NULL,
       y = NULL, 
       color = "Rear")
lap_obs_p

# save run timing plot
ggsave("output/figures/lapwai_iptds_all_detections_timing_20251230.pdf")

# write some objects to excel
list("IPTDS_Metadata"        = iptds_meta,
     "DABOM_Escapement"      = dabom_escape_df,
     "DABOM_Detections"      = dabom_obs_df,
     "All_Lapwai_Detections" = lap_df) %>%
  write_xlsx("output/lapwai_data_request_20251230.xlsx")

# ridge plot
# library(ggridges)
# all_obs_timing_p = lap_df %>%
#   # trim to most recent 5 years, remove small # of remaining U rear_type_code
#   filter(year >= 2021) %>%
#   filter(!rear_type_code == "U") %>%
#   ggplot(aes(x = xday, y = site_code, fill = rear_type_code)) +
#   geom_density_ridges(scale = 2, alpha = 0.6, color = NA) +
#   scale_x_date(date_breaks = "1 month", date_labels = "%b") +
#   facet_wrap(~ spc, ncol = 1, scales = "free_y") +
#   theme_bw() +
#   labs(x = NULL, 
#        y = "Site", 
#        fill = "Rear Type")

### END SCRIPT