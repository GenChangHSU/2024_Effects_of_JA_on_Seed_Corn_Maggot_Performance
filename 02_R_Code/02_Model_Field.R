## -----------------------------------------------------------------------------
## Title: Analysis of the relationship between the number of SCM maggots vs. 
##        plot and seed treatment
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-04
##
## Description:
## 1. Data exploration and summary 
## 2. Fit GLMMs to the maggot count data and select the best model
## 3. Perform model diagnostics and test model significance
## 4. Assess coefficient significance
## 5. Estimated marginal means for the predictors and pairwise comparisons among the treatments
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(lmtest)
library(broom)
library(broom.mixed)
library(car)
library(emmeans)
library(multcomp)


# Import files -----------------------------------------------------------------
corn_data_clean <- read_csv("./03_Outputs/Data_Clean/Corn_Data_Clean.csv")


############################### Code starts here ###############################

# 1. Data exploration ----------------------------------------------------------
### Final dataset for analysis (only look at fertilizer treatment)
corn_data_clean_for_analysis <- corn_data_clean %>% 
  filter(maggot_location == "Total") %>% 
  filter(plot_treatment %in% c("Manure", "Fertilizer", "Control_Manure", "Control_Fertilizer")) %>% 
  drop_na(total_nr_scm_maggots) %>% 
  mutate(plot_treatment = fct_relevel(plot_treatment, "Control_Fertilizer", "Fertilizer", "Control_Manure", "Manure"),
         seed_treatment = fct_relevel(seed_treatment, "Water", "JA"),
         plot_id = as.factor(plot_id),
         collection_date = as.character(collection_date))

### Number of observations
# (1) by plot treatment
corn_data_clean_for_analysis %>% 
  group_by(plot_treatment) %>% 
  summarise(n_observation = n())

# (2) by seed treatment
corn_data_clean_for_analysis %>% 
  group_by(seed_treatment) %>% 
  summarise(n_observation = n())

# (3) by collection date
corn_data_clean_for_analysis %>% 
  group_by(collection_date) %>% 
  summarise(n_observation = n())

### Proportion of zero counts
# (1) by plot treatment
corn_data_clean_for_analysis %>% 
  group_by(plot_treatment) %>% 
  summarise(Prop_zero = sum(total_nr_scm_maggots == 0)/n())

# (2) by seed treatment
corn_data_clean_for_analysis %>% 
  group_by(seed_treatment) %>% 
  summarise(Prop_zero = sum(total_nr_scm_maggots == 0)/n())

# (3) by collection date
corn_data_clean_for_analysis %>% 
  group_by(collection_date) %>% 
  summarise(Prop_zero = sum(total_nr_scm_maggots == 0)/n())

### Data summary for the data from the first collection date 2023-05-02
corn_data_clean_for_analysis_first_collection_date <- corn_data_clean_for_analysis %>% 
  filter(collection_date == "2023-05-02")

# (1) number of observations by plot treatment
corn_data_clean_for_analysis_first_collection_date %>% 
  group_by(plot_treatment) %>% 
  summarise(n_observation = n())

# (2) number of observations by seed treatment
corn_data_clean_for_analysis_first_collection_date %>% 
  group_by(seed_treatment) %>% 
  summarise(n_observation = n())

# (3) proportion of zeros by plot treatment
corn_data_clean_for_analysis_first_collection_date %>% 
  group_by(plot_treatment) %>% 
  summarise(Prop_zero = sum(total_nr_scm_maggots == 0)/n())

# (4) proportion of zeros by seed treatment
corn_data_clean_for_analysis_first_collection_date %>% 
  group_by(seed_treatment) %>% 
  summarise(Prop_zero = sum(total_nr_scm_maggots == 0)/n())

# (5) plots that did not have four observations
corn_data_clean_for_analysis_first_collection_date %>%
  group_by(plot_id) %>% 
  summarise(n = n()) %>% 
  filter(n != 4)


# 2. Fit GLMMs to the maggot count data ----------------------------------------
### (1) Poisson model
poisson_model <- glmmTMB(total_nr_scm_maggots ~ plot_treatment * seed_treatment + (1|plot_id),
                                 data = corn_data_clean_for_analysis_first_collection_date,
                                 family = poisson, 
                                 na.action = na.exclude)

summary(poisson_model)
testDispersion(poisson_model)
testZeroInflation(poisson_model)
plot(simulateResiduals(poisson_model))

### (2) Negative binomial model
nb_model <- glmmTMB(total_nr_scm_maggots ~ plot_treatment * seed_treatment + (1|plot_id),
                            data = corn_data_clean_for_analysis_first_collection_date,
                            family = nbinom2, 
                            na.action = na.exclude)

summary(nb_model)
testZeroInflation(nb_model)
plot(simulateResiduals(nb_model))

### Compare poisson and negative binomial model
lrtest(poisson_model, nb_model)  # negative binomial model is better

### (3) Truncated Poisson model
truncated_poisson_model <- glmmTMB(total_nr_scm_maggots ~ plot_treatment * seed_treatment + (1|plot_id),
                                           data = filter(corn_data_clean_for_analysis_first_collection_date, total_nr_scm_maggots > 0),
                                           family = truncated_poisson, 
                                           na.action = na.exclude)

summary(truncated_poisson_model)
testDispersion(truncated_poisson_model)
plot(simulateResiduals(truncated_poisson_model))

### (4) Truncated negative binomial model
truncated_nb_model <- glmmTMB(total_nr_scm_maggots ~ plot_treatment * seed_treatment + (1|plot_id),
                                           data = filter(corn_data_clean_for_analysis_first_collection_date, total_nr_scm_maggots > 0),
                                           family = truncated_nbinom2, 
                                           na.action = na.exclude)

summary(truncated_nb_model)
plot(simulateResiduals(truncated_nb_model))

### Test for overdispersion for the truncated part
lrtest(truncated_poisson_model, truncated_nb_model)  # the truncated part is overdispersed
AIC(truncated_poisson_model, truncated_nb_model)  # the truncated part is overdispersed

### (5) Negative binomial model without interaction term
nb_model_wo_interaction <- glmmTMB(total_nr_scm_maggots ~ plot_treatment + seed_treatment + (1|plot_id),
                            data = corn_data_clean_for_analysis_first_collection_date,
                            family = nbinom2, 
                            na.action = na.exclude)

summary(nb_model_wo_interaction)
testZeroInflation(nb_model_wo_interaction)
plot(simulateResiduals(nb_model_wo_interaction))

### Test the significance of the interaction term
lrtest(nb_model_glmmTMB, nb_model_glmmTMB_wo_interaction)  # interaction term not significant

### (6) Zero-inflated negative binomial model without interaction term
zi_nb_model_wo_interaction <- glmmTMB(total_nr_scm_maggots ~ plot_treatment + seed_treatment + (1|plot_id),
                       data = corn_data_clean_for_analysis_first_collection_date,
                       ziformula = ~ plot_treatment + seed_treatment, 
                       family = nbinom2, 
                       na.action = na.exclude)

summary(zi_nb_model_wo_interaction)
plot(simulateResiduals(zi_nb_model_wo_interaction))

### Test the significance of zero inflation
lrtest(nb_model_wo_interaction, zi_nb_model_wo_interaction)  # no zero inflation
AIC(nb_model_wo_interaction, zi_nb_model_wo_interaction)  # no zero inflation

### AIC for all models
AIC(poisson_model, 
    nb_model, 
    nb_model_wo_interaction,
    zi_nb_model_wo_interaction) %>% 
  arrange(AIC)  # best model is the negative binomial model without interaction


# 3. Model diagnostics and significance ----------------------------------------
### Model diagnostics
testZeroInflation(nb_model_wo_interaction)  # no zero inflation
plot(simulateResiduals(nb_model_wo_interaction))  # no significant heteroscedasticity

### Model global significance
lrtest(nb_model_wo_interaction)

### Deviance test for goodness of fit 
nb_model_wo_interaction_summary <- summary(nb_model_wo_interaction)
resid_deviance <- nb_model_wo_interaction_summary$AICtab["deviance"]
resid_df <- nb_model_wo_interaction_summary$AICtab["df.resid"]
pchisq(resid_deviance, resid_df, lower.tail = F)


# 4. Coefficient significance --------------------------------------------------
### Summary table of the coefficients
tidy(nb_model_wo_interaction, conf.int = TRUE) %>% view()

### Profile likelihood confidence intervals for the coefficient estimates
nb_model_wo_interaction_profile <- profile(nb_model_wo_interaction)
confint(nb_model_wo_interaction_profile) %>% view()

### LR test for the fixed effects
Anova(nb_model_wo_interaction)


# 5. Emmeans and pairwise comparisons ------------------------------------------
### EMMs for plot and seed treatment
emm_plot_treatment <- emmeans(nb_model_wo_interaction, ~ plot_treatment, type = "response")
emm_seed_treatment <- emmeans(nb_model_wo_interaction, ~ seed_treatment, type = "response")

### Pairwise comparisons between treatment levels
pairs(regrid(emm_plot_treatment))
pairs(regrid(emm_seed_treatment))

cld(emm_plot_treatment, adjust = "sidak", Letters = letters)
cld(emm_seed_treatment, adjust = "sidak", Letters = letters)



