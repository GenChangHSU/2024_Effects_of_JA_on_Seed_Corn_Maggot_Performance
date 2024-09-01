## -----------------------------------------------------------------------------
## Title: Data organization for the cup treatment and maggot count data
##
## Author: Gen-Chang Hsu
##
## Date: 2023-11-23
##
## Description:
## 1. Organize the cup treatment data
## 2. Organize the maggot count data 
## 3. Merge and write the two clean datasets
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(janitor)
library(magrittr)
library(readxl)


# Import files -----------------------------------------------------------------
cup_data_raw <- read_xlsx("./01_Data_Raw/Wire_Mesh_Intake_Corrected.xlsx", sheet = 1)
maggot_data_raw <- read_xlsx("./01_Data_Raw/Wire_Mesh_Data.xlsx", sheet = 1)


############################### Code starts here ###############################

# 1. Organize the cup treatment data -------------------------------------------
cup_data_clean <- cup_data_raw %>% 
  filter(Project == "M") %>% 
  separate(`Record ID`, 
           into = c("plotID_treatment", "seed_type", "seed_treatment"), 
           sep = "_") %>% 
  mutate(plot_ID = str_extract(plotID_treatment, "^[0-9]*"), 
         plot_treatment = str_remove(plotID_treatment, "^[0-9]*"),
         .after = "Project") %>% 
  mutate(seed_type = fct_recode(seed_type, `Lima bean` = "L", Corn = "C"),
         seed_treatment = fct_na_value_to_level(seed_treatment, "Water")) %>% 
  filter(seed_type == "Corn") %>% 
  mutate(plot_treatment = fct_recode(plot_treatment, 
                                     Alfafa = "A",
                                     `Alfafa&Manure` = "AM",
                                     Fertilizer = "F",
                                     Manure = "M",
                                     Control = "C"),
         plot_treatment = as.character(plot_treatment)) %>% 
  group_by(plot_ID, `Collection Date`) %>% 
  mutate(plot_treatment = ifelse(plot_treatment == "Control",
                                 paste0("Control_", str_remove(str_c(unique(plot_treatment), collapse = ""), "Control")),
                                 plot_treatment)) %>% 
  select(-plotID_treatment) %>% 
  relocate(`Collection Date`, .after = Project) %>% 
  arrange(`Collection Date`, plot_ID, plot_treatment)


# 2. Organize the maggot count data --------------------------------------------
maggot_data_clean <- maggot_data_raw %>% 
  rename(Project = `Project (M=Manure; T=Trap; R=Risk Assessment; C=Cover Crop; Intercropping = I)`) %>% 
  filter(Project == "M")


# 3. Merge and write out the two cleaned datasets ------------------------------
corn_data_clean <- cup_data_clean %>% 
  left_join(select(maggot_data_clean, -ID, -Project, -15, -16), by = join_by(`Container Number` == Cup_Nr)) %>% 
  mutate(Notes.x = fct_na_value_to_level(Notes.x, ""),
         Notes.y = fct_na_value_to_level(Notes.y, ""),
         Notes = str_c(Notes.x, Notes.y, sep = "")) %>% 
  select(-Notes.x, -Notes.y) %>% 
  clean_names() %>% 
  rename(maggot_location = "seed_type_bean_or_corn") %T>%
  write_csv(., "./03_Outputs/Data_Clean/Corn_Data_Clean.csv")
  
  
  







