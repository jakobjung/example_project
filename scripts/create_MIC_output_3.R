##
##   Title:       Create MIC output
##
##   Author:      Jakob Jung
##
##   Date:        2023-05-15
##
##   Description: Import the table with all predictors ("../data/all_predictors.csv"). Then, go through all files in
##                "../data/raw_growth_data/k12_tables" and "../data/raw_growth_data/upec_tables" and extract the MIC
##                values from the OD values. The MIC values are then added to the table with all predictors. The
##                resulting table is saved as "../data/all_predictors_with_MIC.csv".

# Load libraries
library(tidyverse)
library(readxl)
library(stringr)

# Set working directory
setwd("~/Documents/UPEC_K12_essgenes_2023_03/scripts")

# Import table with all predictors
all_predictors <- read_csv("../data/all_predictors.csv")
View(all_predictors)

# Create a list of all files in the folder with the raw growth data
k12_files <- list.files("../data/raw_growth_data/k12_tables", full.names = TRUE)
upec_files <- list.files("../data/raw_growth_data/upec_tables", full.names = TRUE)
tiling_files <- list.files("../data/raw_growth_data/tiling_tables", full.names = TRUE)

# import one of the k12 files to try out the code
read_delim(k12_files[1])[-1,] %>%
  # make table long using pivot_longer
  pivot_longer(cols = 2:ncol(.), names_to = "concentration", values_to = "od600") %>%
  # extract the concentration using gsub:
    mutate(concentration = as.numeric(gsub("^(.*) µM.*", "\\1", concentration))) %>%
  # rename the first column to "PNA name" and cahnge - to _
    rename_with(.cols = 1, ~"PNA") %>% mutate(PNA = str_replace(PNA, "-", "_"))  %>%
  # make od600 numeric
    mutate(od600 = as.numeric(od600)) %>%
  # for rows having the same PNA and concentration, take the mean of the OD values. For NA values, keep one NA value
    group_by(PNA, concentration) %>%
    summarise(od600 = mean(od600, na.rm = TRUE)) %>%
  # make new column assessing whether or not theres growth. threshold is 0.1. If NA, classify as "growth"
    mutate(growth = ifelse(od600 > 0.1, "growth", "no_growth")) %>%
  # if NA, classify as "growth":
    mutate(growth = ifelse(is.na(growth), "growth", growth)) %>%
  # for each PNA, take the minimum concentration inhibiting growth. make sure theres no higher concentration with growth.
  # if there is only growth, add value 500
    group_by(PNA) %>%
    mutate(MIC = ifelse(all(growth == "growth"), 500, min(concentration[growth == "no_growth"]))) %>%
  # shrink table to only contain PNA and MIC, and remove duplicates
    select(PNA, MIC) %>%
    distinct()

# run the code for all k12 files
k12_MIC <- map_dfr(k12_files, ~read_delim(.x)[-1,] %>%
  # make table long using pivot_longer
  pivot_longer(cols = 2:ncol(.), names_to = "concentration", values_to = "od600") %>%
  # extract the concentration using gsub:
    mutate(concentration = as.numeric(gsub("^(.*) µM.*", "\\1", concentration))) %>%
  # rename the first column to "PNA name" and cahnge - to _
    rename_with(.cols = 1, ~"PNA") %>% mutate(PNA = str_replace(PNA, "-", "_"))  %>%
  # make od600 numeric
    mutate(od600 = as.numeric(od600)) %>%
  # for rows having the same PNA and concentration, take the mean of the OD values. For NA values, keep one NA value
    group_by(PNA, concentration) %>%
    summarise(od600 = mean(od600, na.rm = TRUE)) %>%
  # make new column assessing whether or not theres growth. threshold is 0.1. If NA, classify as "growth"
    mutate(growth = ifelse(od600 > 0.1, "growth", "no_growth")) %>%
  # if NA, classify as "growth":
    mutate(growth = ifelse(is.na(growth), "growth", growth)) %>%
  # for each PNA, take the minimum concentration inhibiting growth. make sure theres no higher concentration with growth.
  # if there is only growth, add value 500
    group_by(PNA) %>%
    mutate(MIC = ifelse(all(growth == "growth"), 500, min(concentration[growth == "no_growth"]))) %>%
  # shrink table to only contain PNA and MIC, and remove duplicates
    select(PNA, MIC) %>%
    distinct() %>%
  # add column indicating whether the MIC was extracted from a k12 or upec file
    mutate(k12_upec = "K12")) %>%
# rsummarize pnas by median MIC for each PNA
    group_by(PNA) %>%
    summarise(MIC = min(MIC, na.rm = TRUE), k12_upec = "K12")%>%
  # remove stars from PNA names
    mutate(PNA = str_replace(PNA, "\\*", ""))

# run the code for all upec files
upec_MIC <- map_dfr(upec_files, ~read_delim(.x)[-1,] %>%
  # make table long using pivot_longer
  pivot_longer(cols = 2:ncol(.), names_to = "concentration", values_to = "od600") %>%
  # extract the concentration using gsub:
    mutate(concentration = as.numeric(gsub("^(.*) µM.*", "\\1", concentration))) %>%
  # rename the first column to "PNA name" and cahnge - to _
    rename_with(.cols = 1, ~"PNA") %>% mutate(PNA = str_replace(PNA, "-", "_"))  %>%
  # make od600 numeric
    mutate(od600 = as.numeric(od600)) %>%
  # for rows having the same PNA and concentration, take the mean of the OD values. For NA values, keep one NA value
    group_by(PNA, concentration) %>%
    summarise(od600 = mean(od600, na.rm = TRUE)) %>%
  # make new column assessing whether or not theres growth. threshold is 0.1. If NA, classify as "growth"
    mutate(growth = ifelse(od600 > 0.1, "growth", "no_growth")) %>%
  # if NA, classify as "growth":
    mutate(growth = ifelse(is.na(growth), "growth", growth)) %>%
  # for each PNA, take the minimum concentration inhibiting growth. make sure theres no higher concentration with growth.
  # if there is only growth, add value 500
    group_by(PNA) %>%
    mutate(MIC = ifelse(all(growth == "growth"), 500, min(concentration[growth == "no_growth"]))) %>%
  # shrink table to only contain PNA and MIC, and remove duplicates
    select(PNA, MIC) %>%
  # add column indicating whether the MIC was extracted from a k12 or upec file
    mutate(k12_upec = "UPEC")) %>%
# rsummarize pnas by median MIC for each PNA
    group_by(PNA) %>%
    summarise(MIC = min(MIC, na.rm = TRUE), k12_upec = "UPEC") %>%
  # remove stars from PNA names
    mutate(PNA = str_replace(PNA, "\\*", ""))

# combine the k12 and upec MIC tables
MIC_table <- bind_rows(k12_MIC, upec_MIC)

# add MIC for K12, k12_MIC as column to table all_predictors
all_predictors <- all_predictors %>%
  # get MICs, and add them by pna_name column of all_predictors, and PNA column of MIC_table
  mutate(MIC_K12 = k12_MIC$MIC[match(pna_name, k12_MIC$PNA)]) %>%
  # do same for UPEC
    mutate(MIC_UPEC = upec_MIC$MIC[match(pna_name, upec_MIC$PNA)])

# Now save the updated all_predictors table.
write_csv(all_predictors, "../data/all_predictors_with_MICs.csv")


# import one of the tiling files to try out the code
read_delim(tiling_files[2])[-1,] %>%
  # make table long using pivot_longer
  pivot_longer(cols = 2:ncol(.), names_to = "concentration", values_to = "od600") %>%
  # extract the concentration using gsub:
    mutate(concentration = as.numeric(gsub("^(.*) µM.*", "\\1", concentration))) %>%
  # rename the first column to "PNA name" and cahnge - to _
    rename_with(.cols = 1, ~"PNA") %>% mutate(PNA = str_replace(PNA, "-", "_"))  %>%
  # make od600 numeric
    mutate(od600 = as.numeric(od600)) %>%
  # for rows having the same PNA and concentration, take the mean of the OD values. For NA values, keep one NA value
    group_by(PNA, concentration) %>%
    summarise(od600 = mean(od600, na.rm = TRUE)) %>%
  # make new column assessing whether or not theres growth. threshold is 0.1. If NA, classify as "growth"
    mutate(growth = ifelse(od600 > 0.1, "growth", "no_growth")) %>%
  # if NA, classify as "growth":
    mutate(growth = ifelse(is.na(growth), "growth", growth)) %>%
  # for each PNA, take the minimum concentration inhibiting growth. make sure theres no higher concentration with growth.
  # if there is only growth, add value 500
    group_by(PNA) %>%
    mutate(MIC = ifelse(all(growth == "growth"), 500, min(concentration[growth == "no_growth"]))) %>%
  # shrink table to only contain PNA and MIC, and remove duplicates
    select(PNA, MIC) %>%
    distinct()


# run the code for all tiling files
# run the code for all k12 files
tiling_MIC <- map_dfr(tiling_files, ~read_delim(.x)[-1,] %>%
  # make table long using pivot_longer
  pivot_longer(cols = 2:ncol(.), names_to = "concentration", values_to = "od600") %>%
  # extract the concentration using gsub:
    mutate(concentration = as.numeric(gsub("^(.*) µM.*", "\\1", concentration))) %>%
  # rename the first column to "PNA name" and cahnge - to _
    rename_with(.cols = 1, ~"PNA") %>% mutate(PNA = str_replace(PNA, "-", "_"))  %>%
  # make od600 numeric
    mutate(od600 = as.numeric(od600)) %>%
  # for rows having the same PNA and concentration, take the mean of the OD values. For NA values, keep one NA value
    group_by(PNA, concentration) %>%
    summarise(od600 = mean(od600, na.rm = TRUE)) %>%
  # make new column assessing whether or not theres growth. threshold is 0.1. If NA, classify as "growth"
    mutate(growth = ifelse(od600 > 0.1, "growth", "no_growth")) %>%
  # if NA, classify as "growth":
    mutate(growth = ifelse(is.na(growth), "growth", growth)) %>%
  # for each PNA, take the minimum concentration inhibiting growth. make sure theres no higher concentration with growth.
  # if there is only growth, add value 500
    group_by(PNA) %>%
    mutate(MIC = ifelse(all(growth == "growth"), 40, min(concentration[growth == "no_growth"]))) %>%
  # shrink table to only contain PNA and MIC, and remove duplicates
    select(PNA, MIC) %>%
    distinct() %>%
  # add column indicating whether the MIC was extracted from a k12 or upec file
    mutate(strain = "UPEC")) %>%
# rsummarize pnas by median MIC for each PNA
    group_by(PNA) %>%
    summarise(MIC = min(MIC, na.rm = TRUE), strain = "UPEC")%>%
  # remove stars from PNA names
    mutate(PNA = str_replace(PNA, "\\*", ""))

# write tiling mics to file
write_csv(tiling_MIC, "../data/tiling_test/tiling_mics.csv")
