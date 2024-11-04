## -----------------------------------------------------------------------------
## Title: Analysis of the no-choice oviposition experiment data
##
## Author: Gen-Chang Hsu
##
## Date: 2024-03-19
##
## Description:
## 1. Organize the raw data for model fitting and visualization
## 2. Fit a poisson GLMM on the number of eggs in each trial
## 3. Create a figure showing the number of eggs in the control vs. MeJA treatment
## 4. Package citations
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(car)


# ggplot theme -----------------------------------------------------------------
my_theme <- 
  theme(# Axis
        axis.text.x = element_text(size = 14, color = "black", margin = margin(t = 3)),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 8)),
        axis.ticks.length.x = unit(0.18, "cm"),
        axis.ticks.length.y = unit(0.15, "cm"),
                
        # Plot
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(colour = "transparent"),
        
        # Panel
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        
        # Legend
        legend.position = "right",
        legend.spacing.x = unit(0.2, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.box.just = "center",
        legend.justification = c(0.5, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.25, linetype = "solid", colour = "black"),
        
        # Facet strip
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13, hjust = 0.5)
        )


# Import files -----------------------------------------------------------------
nochoice_data_raw <- read_xlsx("./01_Data_Raw/Nochoice_Data.xlsx", sheet = 1)


############################### Code starts here ###############################

# 1. Data organization ---------------------------------------------------------
nochoice_data_clean <- nochoice_data_raw %>% 
  filter(Trial_ID != 0) %>% 
  select(Trial_ID, Cage_ID, Seed_treatment, N_eggs, Date) %>%
  mutate(Seed_treatment = factor(Seed_treatment, levels = c("Water", "JA"), ordered = F),
         Trial_cage_ID = str_c(Trial_ID, Cage_ID, sep = "_"), .after = Cage_ID) %>% 
  group_by(Trial_cage_ID) %>% 
  mutate(Treatment_order = order(Date),
         Treatment_order = case_match(Treatment_order, 1 ~ "First", 2 ~ "Second")) %>% 
  group_by(Trial_cage_ID) %>%  
  mutate(Total_eggs = sum(N_eggs),
         Prop_eggs = N_eggs/Total_eggs) %>% 
  mutate(Prop_eggs_modified = case_when(Prop_eggs == 0 ~ 0.001,
                                        Prop_eggs == 1 ~ 0.999,
                                        TRUE ~ Prop_eggs)) %>% 
  ungroup()

### Mean proportion of eggs in JA treatment in each trial
nochoice_data_clean %>% 
  filter(Seed_treatment == "JA") %>% 
  pull(Prop_eggs) %>% 
  mean()

### Overall proportion of eggs in JA treatment across all trials
nochoice_data_clean %>% 
  group_by(Seed_treatment) %>% 
  summarise(N_eggs = sum(N_eggs)) %>% 
  mutate(Prop_eggs = N_eggs/sum(N_eggs)) %>% 
  filter(Seed_treatment == "JA") %>% 
  pull(Prop_eggs)


# 2. A GLMM for the number of eggs in each trial ----------------------------------
### Poisson and NB GLMM using glmmTMB
glmmtmb_Poisson_out <- glmmTMB(N_eggs ~ Seed_treatment + Treatment_order + (1|Trial_cage_ID), 
                          data = nochoice_data_clean,
                          family = "poisson",
                          na.action = na.omit)

glmmtmb_nb_out <- glmmTMB(N_eggs ~ Seed_treatment + Treatment_order + (1|Trial_cage_ID), 
                               data = nochoice_data_clean,
                               family = "nbinom1",
                               na.action = na.omit)

summary(glmmtmb_Poisson_out)
Anova(glmmtmb_Poisson_out)

summary(glmmtmb_nb_out)
Anova(glmmtmb_nb_out)

testDispersion(glmmtmb_Poisson_out)
AIC(glmmtmb_nb_out, glmmtmb_Poisson_out)
library(lmtest)
lrtest(glmmtmb_nb_out, glmmtmb_Poisson_out)

### NB GLMM using glmmTMB
library(lme4)
glmer_nb_out <- glmer.nb(N_eggs ~ Seed_treatment + Treatment_order + (1|Trial_cage_ID), 
                               data = nochoice_data_clean)

summary(glmer_nb_out)
Anova(glmer_nb_out)


### Model diagnostics
sim_resid <- simulateResiduals(glmmtmb_Poisson_out, n = 1000) 
plotQQunif(sim_resid)  # QQ-plot
plotResiduals(sim_resid)  # fitted vs. residual quantile plot
plotResiduals(sim_resid, form = nochoice_data_clean$Seed_treatment)
plotResiduals(sim_resid, form = nochoice_data_clean$Treatment_order)

testUniformity(sim_resid) 
testOutliers(sim_resid)
testDispersion(sim_resid) 
testQuantiles(sim_resid)
testCategorical(sim_resid, catPred = nochoice_data_clean$Seed_treatment)
testCategorical(sim_resid, catPred = nochoice_data_clean$Treatment_order)
testZeroInflation(sim_resid)


# 3. Figures of the number and proportion of eggs in water vs. JA treatment ----
### Boxplot
ggplot(nochoice_data_clean) + 
  geom_line(aes(x = Seed_treatment, y = N_eggs, group = Trial_cage_ID),
            position = position_jitter(width = 0.01, height = 1, seed = 1), color = "grey80") + 
  geom_point(aes(x = Seed_treatment, y = N_eggs, color = Seed_treatment), 
             position = position_jitter(width = 0.01, height = 1, seed = 1), alpha = 1, size = 2) + 
  geom_boxplot(aes(x = Seed_treatment, y = N_eggs, color = Seed_treatment), 
               position = position_nudge(x = c(-0.25, 0.25)), width = 0.2, 
               outlier.color = NA, show.legend = F) +
  labs(x = "Seed treatment", y = "Number of eggs per cup", shape = "Treatment order") + 
  scale_color_manual(values = c("grey30", "grey")) + 
  scale_x_discrete(labels = c("Control", "MeJA")) +
  guides(color = "none") + 
  my_theme + 
  theme(legend.title = element_text(size = 13, margin = margin(l = -15)),
        legend.box.margin = margin(l = 5),
        legend.margin = margin(t = 8, r = 10, b = 8, l = 22.5))

ggsave("./03_Outputs/Figures/Nochoice_Boxplot.tiff", width = 6, height = 5, dpi = 600, device = "tiff")

### Barplot
nochoice_data_clean_summary <- nochoice_data_clean %>% 
  group_by(Seed_treatment) %>% 
  summarise(Mean_n_eggs = mean(N_eggs),
            n = n(),
            SD_n_eggs = sd(N_eggs),
            SE_n_eggs = SD_n_eggs/sqrt(n))

ggplot() + 
  geom_line(data = nochoice_data_clean, aes(x = Seed_treatment, y = N_eggs, group = Trial_cage_ID),
            position = position_nudge(x = c(0.25, -0.25)), color = "grey80") + 
  geom_point(data = nochoice_data_clean, aes(x = Seed_treatment, y = N_eggs, color = Seed_treatment), 
             position = position_nudge(x = c(0.25, -0.25)), alpha = 1, size = 2) + 
  geom_col(data = nochoice_data_clean_summary, 
           aes(x = Seed_treatment, y = Mean_n_eggs, fill = Seed_treatment),
           position = position_nudge(x = c(-0.1, 0.1)),
           width = 0.5, show.legend = F, color = "black") +
  geom_errorbar(data = nochoice_data_clean_summary, 
                aes(x = Seed_treatment, ymin = Mean_n_eggs - SE_n_eggs, 
                    ymax = Mean_n_eggs + SE_n_eggs), width = 0.1,
                position = position_nudge(x = c(-0.1, 0.1))) + 
  labs(x = "Seed treatment", y = "Number of eggs per cup", shape = "Treatment order") + 
  scale_color_manual(values = c("grey30", "grey")) + 
  scale_fill_manual(values = c("grey30", "grey")) +
  scale_x_discrete(labels = c("Control", "MeJA")) +
  scale_y_continuous(limits = c(-0.5, 215), expand = c(0, 0)) + 
  guides(color = "none") + 
  my_theme + 
  theme(legend.title = element_text(size = 13, margin = margin(l = -15)),
        legend.box.margin = margin(l = 5),
        legend.margin = margin(t = 8, r = 10, b = 8, l = 22.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust = c(0.8, 0)))

ggsave("./03_Outputs/Figures/Nochoice_Barplot.tiff", width = 6, height = 5, dpi = 600, device = "tiff")


# 4. Package citations ---------------------------------------------------------
capture.output(utils:::print.bibentry(citation(), style = "Bibtex"),
               utils:::print.bibentry(citation("glmmTMB"), style = "Bibtex"),
               utils:::print.bibentry(citation("car"), style = "Bibtex"),
               utils:::print.bibentry(citation("DHARMa"), style = "Bibtex"),
               file = "./06_References/R_References.bib")




        



















