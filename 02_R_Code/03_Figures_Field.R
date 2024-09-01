## -----------------------------------------------------------------------------
## Title: Visualization of the seed corn maggot data
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-04
##
## Description:
## 1. Organize the seed corn maggot data for plotting
## 2. Create barplots of the probability of zero counts in each treatment
## 3. Create boxplots of the number of maggots in each treatment (including zeros) 
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)


# ggplot theme -----------------------------------------------------------------
my_theme <- 
  theme(# Axis
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 8)),
        axis.ticks.length.x = unit(0.2, "cm"),
        
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
        legend.key.width = unit(0.5, "cm"),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.box.just = "center",
        legend.justification = c(0.5, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "black"),
        
        # Facet strip
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13, hjust = 0.5)
        )
      
  
# Import files -----------------------------------------------------------------
corn_data_clean <- read_csv("./03_Outputs/Data_Clean/Corn_Data_Clean.csv")


############################### Code starts here ###############################

# 1. Data organization ---------------------------------------------------------
corn_data_clean_for_plotting <- corn_data_clean %>% 
  filter(maggot_location == "Total") %>% 
  filter(plot_treatment %in% c("Manure", "Fertilizer", "Control_Manure", "Control_Fertilizer")) %>% 
  drop_na(total_nr_scm_maggots) %>% 
  mutate(plot_treatment = fct_relevel(plot_treatment, "Control_Fertilizer", "Fertilizer", "Control_Manure", "Manure"),
         seed_treatment = fct_relevel(seed_treatment, "Water", "JA"),
         plot_id = as.factor(plot_id),
         collection_date = as.character(collection_date)) %>% 
  filter(collection_date == "2023-05-02")


# 2. Barplots of the probability of zero counts in the treatments --------------
### By plot treatment
corn_data_clean_for_plotting %>% 
  group_by(plot_treatment) %>% 
  summarise(n = n(),
            n_zero = sum(total_nr_scm_maggots == 0),
            prop_zero = n_zero/n) %>% 
  mutate(plot_treatment_type = str_remove(plot_treatment, "Control_"),
         treatment_control = ifelse(str_detect(plot_treatment, "Control_"), 
                                    "Control", 
                                    "Treatment")) %>% 
  ggplot() + 
  geom_col(aes(x = plot_treatment_type, y = prop_zero, 
               fill = plot_treatment_type, 
               alpha = treatment_control),
           position = position_dodge(0.605),
           color = "black",
           width = 0.6,
           show.legend = F) + 
  scale_x_discrete(labels = c("Fertilizer", "Manure")) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_fill_manual(values = c("#31a354", "#f1a340")) + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  labs(x = "Plot treatment", y = "Proportion of cups without maggots") + 
  annotate(geom = "text", label = "Control", x = c(0.85, 1.85), y = c(0.58, 0.5)) +
  annotate(geom = "text", label = c("+F", "+M"), x = c(1.15, 2.15), y = c(0.62, 0.705)) +
  my_theme + 
  theme(axis.ticks.x = element_blank())

ggsave("./03_Outputs/Figures/Prop_Zero_Plot_Treatment_Barplot.tiff", width = 5, height = 4.5, dpi = 600, device = "tiff")

### By seed treatment
corn_data_clean_for_plotting %>% 
  group_by(seed_treatment) %>% 
  summarise(n = n(),
            n_zero = sum(total_nr_scm_maggots == 0),
            prop_zero = n_zero/n) %>% 
  ggplot() + 
  geom_col(aes(x = seed_treatment, y = prop_zero, fill = seed_treatment),
           color = "black", width = 0.6, show.legend = F) + 
  scale_x_discrete(labels = c("Water", "MeJA")) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_fill_manual(values = c("#2b8cbe", "#e34a33")) + 
  labs(x = "Seed treatment", y = "Proportion of cups without maggots") + 
  my_theme + 
  theme(axis.ticks.x = element_blank())

ggsave("./03_Outputs/Figures/Prop_Zero_Seed_Treatment_Barplot.tiff", width = 3.5, height = 4.5, dpi = 600, device = "tiff")

# 3. Boxplots of the number of maggots in the treatments -----------------------
### By plot treatment
corn_data_clean_for_plotting %>% 
  mutate(plot_treatment_type = str_remove(plot_treatment, "Control_"),
         treatment_control = ifelse(str_detect(plot_treatment, "Control_"), 
                                    "Control", 
                                    "Treatment")) %>%
  ggplot() + 
  geom_boxplot(aes(x = plot_treatment_type, y = total_nr_scm_maggots, 
                   fill = plot_treatment_type, 
                   alpha = treatment_control),
           position = position_dodge(0.7), 
           outlier.alpha = 1,
           width = 0.6,
           show.legend = F) + 
  scale_x_discrete(labels = c("Fertilizer", "Manure")) + 
  scale_y_continuous(limits = c(-2, 50), expand = c(0, 0)) + 
  scale_fill_manual(values = c("#31a354", "#f1a340")) + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  labs(x = "Plot treatment", y = "Number of maggots per cup") + 
  annotate(geom = "text", label = "Control", x = c(0.825, 1.825), y = c(14, 14)) +
  annotate(geom = "text", label = c("+F", "+M"), x = c(1.175, 2.175), y = c(14, 14)) +
  my_theme +
  theme(axis.ticks.x = element_blank())

ggsave("./03_Outputs/Figures/Maggot_Count_Plot_Treatment_Boxplot.tiff", width = 4, height = 4.5, dpi = 600, device = "tiff")

### By seed treatment
corn_data_clean_for_plotting %>% 
  ggplot() + 
  geom_boxplot(aes(x = seed_treatment, y = total_nr_scm_maggots, fill = seed_treatment),
               show.legend = F) + 
  scale_x_discrete(labels = c("Water", "MeJA")) + 
  scale_y_continuous(limits = c(-2, 50), expand = c(0, 0)) + 
  scale_fill_manual(values = c("#2b8cbe", "#e34a33")) + 
  labs(x = "Seed treatment", y = "Number of maggots per cup") + 
  my_theme

ggsave("./03_Outputs/Figures/Maggot_Count_Seed_Treatment_Boxplot.tiff", width = 4, height = 4.5, dpi = 600, device = "tiff")








