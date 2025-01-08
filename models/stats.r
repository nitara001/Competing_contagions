library(igraph)
library(ggplot2)
library(mgcv)
library(marginaleffects)
library(sjPlot)
library(viridis)
library(lme4)
#check for correlations between variables
library(GGally)

relevant_data <- data[, c("network_size", "density", "num_communities", "modularity.x")]

ggpairs(relevant_data)
data<- read.csv("C:\\Users\\s2607536\\Documents\\GitHub\\Competing_contagions\\data.csv")

#mean or median centre variables
data$network_size_centered <- scale(data$network_size, center = TRUE, scale = FALSE)
data$density_centered <- scale(data$density, center = TRUE, scale = FALSE)
data$modularity_centered <- scale(data$modularity.x, center = TRUE, scale = FALSE)
data$num_communities_centered <- data$num_communities - median(data$num_communities, na.rm = TRUE)
data$module_size_variation_centered <- data$module_size_variation - median(data$module_size_variation, na.rm = TRUE)

#compare number of communities---------------------------------------------------------------------------------------------------------
# Linear model
model_linear <- lme4::lmer(
  time ~ type * (modularity_centered + density_centered +
                 network_size_centered + 
                 module_size_variation_centered + num_communities) +
    (1 | species/network), 
  data = data
)
summary(model_linear)
# Ordinal model
data$num_communities <- factor(data$num_communities, ordered = TRUE)
model_ordinal <- lme4::lmer(
  time ~ type * (modularity_centered + density_centered +
                 network_size_centered + 
                 module_size_variation_centered + num_communities) +
    (1 | species/network), 
  data = data
)

# Factor model
data$num_communities <- factor(data$num_communities, ordered = FALSE)
model_factor <- lme4::lmer(
  time ~ type * (modularity_centered + density_centered +
                 network_size_centered + 
                 module_size_variation_centered + num_communities) +
    (1 | species/network), 
  data = data
)

# Compare AICs
aic_results <- data.frame(
  Model = c("Linear", "Ordinal", "Factor"),
  AIC = c(AIC(model_linear), AIC(model_ordinal), AIC(model_factor)),
  Parameters = c(length(fixef(model_linear)), 
                 length(fixef(model_ordinal)), 
                 length(fixef(model_factor)))
)
print(aic_results)
#linear model chosen as although factor has lowest aic, not much lower 
##----------------------------------------------------------------------------------------------------

##try with density as different shapes 
model_quadratic <- lme4::lmer(
  time ~ type * (modularity_centered + network_size_centered +
                 module_size_variation_centered + density_centered + 
                 I(density_centered^2)) +
    (1 | species/network), 
  data = data
)
summary(model_quadratic)
model_cubic <- lme4::lmer(
  time ~ type * (modularity_centered + network_size_centered +
                 module_size_variation_centered + density_centered +
                 I(density_centered^2) + I(density_centered^3)) +
    (1 | species/network), 
  data = data
)
summary(model_cubic)
aic_density_shapes <- data.frame(
  Model = c("Linear", "Quadratic", "Cubic", "Spline"),
  AIC = c(AIC(model_linear), 
          AIC(model_quadratic), 
          AIC(model_cubic), 
          AIC(model_spline)
))
print(aic_density_shapes)

model_spline <- lme4::lmer(
  time ~ type * (modularity_centered + network_size_centered +
                 module_size_variation_centered + bs(density_centered, degree = 3)) +
    (1 | species/network), 
  data = data
)
summary(model_spline)
#=----------------------------------------------------------
#plots
plot_model(
  model_linear,type = "pred", terms = c("modularity_centered", "type"), title = "",
  show.data = TRUE, dot.size = 0.8, dot.alpha = 0.09, line.size = 1.25
) + ylim(0, 800) + xlab("Modularity") + ylab("Time for Spread to Reach 75%") +
  theme_minimal() +
  scale_color_manual(  values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information") ) +
  scale_fill_manual(
    values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information")
  ) +
  theme(
    legend.position = "right",legend.title = element_blank(),
    legend.text = element_text(size = 10),legend.key.size = unit(0.7, "cm"),legend.spacing.x = unit(0.5, "cm")
  )

plot_model(
  model_linear,type = "pred",terms = c("density_centered", "type"),title = "",
  show.data = TRUE,  dot.size = 0.8,
  dot.alpha = 0.09,
  line.size = 1.25
) +
  ylim(0, 800) + xlab("Density") +
  ylab("Time for Spread to Reach 75%") +theme_minimal() +scale_color_manual(
    values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information")
  ) +scale_fill_manual(  values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information")
  ) +
  theme(  legend.position = "right",  legend.title = element_blank(),
    legend.text = element_text(size = 10),legend.key.size = unit(0.7, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  )

plot_model(
  model_linear, type = "pred", terms = c("num_communities", "type"), title = "",show.data = TRUE, dot.size = 0.8, dot.alpha = 0.09,
 line.size = 1.25
) + ylim(0, 800) + xlab("Number of communities") +
  ylab("Time for Spread to Reach 75%") + theme_minimal() +
  scale_color_manual(
    values = c("infection" = "#D2042D", "information" = "#0072B2"),  labels = c("Infection", "Information")) +
  scale_fill_manual(  values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information")
  ) +
  theme(
    legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 10),
 legend.key.size = unit(0.7, "cm"), legend.spacing.x = unit(0.5, "cm") )

plot_model(
  model_linear,
  type = "pred",terms = c("module_size_variation_centered", "type"), title = "", show.data = TRUE, dot.size = 0.8,
  dot.alpha = 0.09, line.size = 1.25) +
  ylim(0, 800) + xlab("Module Size Variation") + ylab("Time for Spread to Reach 75%") +
  theme_minimal() +
  scale_color_manual(
    values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information") ) +
  scale_fill_manual(values = c("infection" = "#D2042D", "information" = "#0072B2"),
    labels = c("Infection", "Information")) +theme(  legend.position = "right", legend.title = element_blank(),legend.text = element_text(size = 10),
    legend.key.size = unit(0.7, "cm"), legend.spacing.x = unit(0.5, "cm"))

plot_model(model_linear, type = "re", title = "Random Effects for Species")
plot_model(
  model_linear,
  type = "re", # Random effects plot
  title = "Random Effects for Species", vline.color = "black", # Set the zero line color
vline.linetype = "dotted", # Set the zero line to dotted colors = c("#FF1493", "#568f58") # Use pink for negative and green for positive effects
) 
##interaction type plots
#  modularity into Low, Medium, High based on quantiles

data$modularity_group <- cut(
  data$modularity_centered,
  breaks = c(-Inf, -0.0586, 0.0383, Inf), # Using 25% and 75% quantiles
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

data$density_group <- cut(
  data$density_centered,
  breaks = c(-Inf, 0.0634, Inf), # Using 50% quantile for Low/High split
  labels = c("Low", "High"),
  include.lowest = TRUE
)

custom_colors <- c("infection" = "#E63B35", "information" = "#3C5488")

ggplot(data, aes(x = modularity_centered, y = time, color = type)) +
  geom_point(alpha = 0.00, size = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, size = 1.5) +  # Increase line thickness
  facet_grid(    density_group ~ modularity_group, 
 scales = "free", 
    labeller = labeller(
      modularity_group = label_value,
      density_group = as_labeller(c("Low" = "Low Density", "High" = "High Density"))
    )
  ) + labs(title = "", x = "Modularity", y = "Time for Spread to Reach 75%", color = "Type"
  ) +
  scale_color_manual(values = custom_colors) + ylim(0, 500) +
  theme_minimal() +
  theme(  strip.text = element_text(size = 12, margin = margin(b = 5)),   legend.position = "bottom", 
    legend.box = "horizontal",   legend.title = element_blank(), 
    legend.text = element_text(size = 10),  legend.key.size = unit(0.6, "cm"),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 10)
  )


#DIFFERNCES IN TIME--------------------------------------------------------------------
mean_times <- data %>%
  group_by(network, type) %>%
  summarise(mean_time = mean(time, na.rm = TRUE)) %>%
  pivot_wider(names_from = type, values_from = mean_time, names_prefix = "time_")
mean_times <- mean_times %>%
  mutate(time_diff = time_infection - time_information)
 <- data %>%
  left_join(mean_times, by = "network")

##time differences by modualritty
proportions_data <- b %>%
  group_by(modularity_group) %>%
  summarize(
    info_faster = mean(time_diff < 0, na.rm = TRUE), # Proportion of networks where information is faster
    infection_faster = mean(time_diff > 0, na.rm = TRUE) # Proportion of networks where infection is faster
  ) %>%
  mutate(
    difference = info_faster - infection_faster # Calculate difference in proportions
  )

ggplot(proportions_data, aes(x = modularity_group, y = difference, fill = "Difference")) +
  geom_bar(stat = "identity", width = 0.6) + # Diverging bar plot
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add a reference line at y = 0
  scale_fill_manual(values = c("purple"), guide = "none") + # Set color for the bars
  labs(
    title = "",
    x = "Modularity Group",
    y = "Difference (Proportion)"
  ) +
  theme_minimal(base_size = 14) + # Minimal theme with base font size
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), # Bold title centered
    axis.title.x = element_text(size = 14, margin = margin(t = 10)), # X-axis label styling
    axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Y-axis label styling
    axis.text = element_text(size = 12) # Axis text styling
  )

##-----------------------------------------------------------------------

custom_colors <- c("infection" = "#E63B35", "information" = "#3C5488")

ggplot(data, aes(x = modularity_centered, y = time, color = type)) +
  geom_point(alpha = 0.00, size = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, size = 1.5) +
  facet_grid(
    density_group ~ modularity_group, 
    scales = "free",  labeller = labeller(   modularity_group = label_value,
      density_group = as_labeller(c("Low" = "Low Density", "High" = "High Density")))
  ) +
  labs( title = "", x = "Modularity", y = "Time for Spread to Reach 75%", color = "Type" ) +
  scale_color_manual(values = custom_colors) + ylim(0, 500) +
  theme_minimal() + theme(
    strip.text = element_text(size = 12, margin = margin(b = 5)), legend.position = "bottom", 
    legend.box = "horizontal",   legend.title = element_blank(), legend.text = element_text(size = 10), 
    legend.key.size = unit(0.6, "cm"),axis.title.x = element_text(size = 14, margin = margin(t = 10)),axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(size = 10), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Keep distinct panel boxes
    strip.text.y = element_text(margin = margin(l = 15, r = 15))  # Extra spacing for density labels
  )


re_plot <- plot_model(
  model_linear,
  type = "re", # Random effects plot
  title = "Random Effects for Species",
  vline.color = "black", 
  vline.linetype = "dotted", 
  colors = c("#FF1493", "#568f58") # Pink for negative, green for positive
)

re_plot[[2]] + theme_minimal()
##---------------------------------------------------------------------------

#num_communities
summary_data <- b %>%
  filter(!is.na(time)) %>% 
  group_by(num_communities, type) %>%
  summarize(
    mean_time = mean(time, na.rm = TRUE), # Mean of raw time values
    sd_time = sd(time, na.rm = TRUE),     # Standard deviation of raw time values
    se_time = sd(time, na.rm = TRUE) / sqrt(n()), # Standard error
    count = n(),
    .groups = "drop"  )

summary_data <- summary_data %>% filter(count > 200)

ggplot(summary_data, aes(x = num_communities, y = mean_time, color = type)) +
  geom_errorbar( aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
    width = 0.2, size = 0.8, position = position_dodge(0.3) ) +
  geom_line(aes(group = type), size = 1, position = position_dodge(0.3)) +
  geom_point(size = 3, position = position_dodge(0.3)) +labs(
    title = "", x = "Number of Communities", y = "Mean Time ± SE",
    color = "Contagion Type" ) +theme_minimal() + theme( legend.position = "right",  plot.title = element_text(hjust = 0.5))

##---------------------------------------------------------------------
summary_data <- b %>%
  group_by(modularity_group) %>%
  summarize(
    mean_diff = mean(time_diff, na.rm = TRUE),        #
    se_diff = sd(time_diff, na.rm = TRUE) / sqrt(n()) # 
  )

ggplot(summary_data, aes(x = modularity_group, y = mean_diff, color = mean_diff > 0)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff), width = 0.2) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_color_manual( values = c("TRUE" = "red", "FALSE" = "blue"),
    labels = c("Information Faster", "Infection Faster")) +
  labs(  title = "", x = "Modularity Group", y = "Mean Time Difference (Infection - Information)",
    color = ""  ) +
  theme_minimal(base_size = 14) + theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12), legend.position = "right")


##---------------------------------------------------------------------------
# Filter data to include only num_communities with more than 200 counts
filtered_summary_data <- b %>%
  group_by(num_communities) %>%
  filter(n() > 300) %>% 
  summarize(
    mean_diff = mean(time_diff, na.rm = TRUE),        # Mean time difference
    se_diff = sd(time_diff, na.rm = TRUE) / sqrt(n()), # Standard error of the mean
    .groups = "drop"
  )

ggplot(filtered_summary_data, aes(x = as.factor(num_communities), y = mean_diff, color = mean_diff > 0)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff), width = 0.2) + # Add error bars
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Reference line at 0
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "blue"),
    labels = c("Information Faster", "Infection Faster")
  ) +
  labs(
    title = "",
    x = "Number of Communities",
    y = "Mean Time Difference (Infection - Information)",
    color = "Outcome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

##-------------------------------------------------------------------------------------------------
##ridge plot
library(ggridges)

filtered_data <- b %>%
  filter(!is.na(num_communities)) # Exclude NAs in num_communities

ggplot(filtered_data, aes(x = time_diff, y = as.factor(num_communities), fill = as.factor(num_communities))) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) + # Ridge plot
  scale_fill_brewer(palette = "Set3", guide = "none") + # Customize colors
  labs(
    title = "",
    x = "Time Difference (Infection - Information)",
    y = "Number of Communities"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12)
  )

##-------------------------------------------------------------------------------
# 3 way interaction
model_linear2 <- lme4::lmer(
  time ~ type * modularity_centered * num_communities +
         density_centered + network_size_centered + 
         module_size_variation_centered +
         (1 | species/network), 
  data = data
)

summary(model_linear2)
