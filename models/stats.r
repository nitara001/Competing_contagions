library(igraph)
library(ggplot2)
library(mgcv)
library(marginaleffects)
#check for correlations between variables

model_1 <- lm(prop_infected ~ modularity + avg_path_length + clustering + density, data = final_resultdf)
model_2 <- lm(conformist.prop_informed ~ modularity + avg_path_length + clustering + density, data = final_resultdf)
coef_model_1 <- coef(model_1)
coef_model_2 <- coef(model_2)

###-----------------------------------------------------

# Plot for modularity
plot_modularity <- ggplot(final_resultdf, aes(x = modularity)) +
  # Smooth line and confidence interval for prop_infected
  geom_smooth(aes(y = prop_infected), method = "loess", color = "#c33a63", se = FALSE, size = 4) +
  # Regression line for prop_infected
  geom_abline(intercept = coef_model_1["(Intercept)"] + coef_model_1["density"] * mean(final_resultdf$density), slope = coef_model_1["modularity"], color = "#c33a63", linetype = "dashed", size = 1.5) +
  # Smooth line and confidence interval for conformist.prop_informed
  geom_smooth(aes(y = conformist.prop_informed), method = "loess", color = "#7895b2", se = FALSE, size = 4) +
  # Regression line for conformist.prop_informed
  geom_abline(intercept = coef_model_2["(Intercept)"] + coef_model_2["density"] * mean(final_resultdf$density), slope = coef_model_2["modularity"], color = "#7895b2", linetype = "dashed", size = 1.5) +
  # Add titles and labels
  labs(title = "",
       y = "Proportion",
       x = "Modularity") +
  # Classic theme with adjustments
  theme_classic() +
  theme(
    axis.text.x = element_text(size = rel(5), color = "black"),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.title.y = element_text(size = rel(3.7), color = "black", margin = margin(r = 25)),
    axis.title.x = element_text(size = rel(3.7), color = "black", margin = margin(t = 25)),
    panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA)  # transparent plot background
  )

ggsave("modularity_plot_with_slopes.svg", plot = plot_modularity, bg = "transparent", width = 10.5, height = 7.5)

# Plot for clustering
plot_clustering <- ggplot(final_resultdf, aes(x = clustering)) +
  # Smooth line and confidence interval for prop_infected
  geom_smooth(aes(y = prop_infected), method = "loess", color = "#c33a63", se = FALSE, size = 4) +
  # Regression line for prop_infected
  geom_abline(intercept = coef_model_1["(Intercept)"] + coef_model_1["density"] * mean(final_resultdf$density), slope = coef_model_1["clustering"], color = "#c33a63", linetype = "dashed", size = 1.5) +
  # Smooth line and confidence interval for conformist.prop_informed
  geom_smooth(aes(y = conformist.prop_informed), method = "loess", color = "#7895b2", se = FALSE, size = 4) +
  # Regression line for conformist.prop_informed
  geom_abline(intercept = coef_model_2["(Intercept)"] + coef_model_2["density"] * mean(final_resultdf$density), slope = coef_model_2["clustering"], color = "#7895b2", linetype = "dashed", size = 1.5) +
  # Add titles and labels
  labs(title = "",
       y = "Proportion",
       x = "Clustering") +
  # Classic theme with adjustments
  theme_classic() +
  theme(
    axis.text.x = element_text(size = rel(5), color = "black"),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.title.y = element_text(size = rel(3.7), color = "black", margin = margin(r = 25)),
    axis.title.x = element_text(size = rel(3.7), color = "black", margin = margin(t = 25)),
    panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA)  # transparent plot background
  )

ggsave("clustering_plot_with_slopes.svg", plot = plot_clustering, bg = "transparent", width = 10.5, height = 7.5)

# Plot for avg_path_length
plot_avg_path_length <- ggplot(final_resultdf, aes(x = avg_path_length)) +
  # Smooth line and confidence interval for prop_infected
  geom_smooth(aes(y = prop_infected), method = "lm", color = "#c33a63", se = FALSE, size = 4) +
  # Regression line for prop_infected
  geom_abline(intercept = coef_model_1["(Intercept)"] + coef_model_1["density"] * mean(final_resultdf$density), slope = coef_model_1["avg_path_length"], color = "#c33a63", linetype = "dashed", size = 1.5) +
  # Smooth line and confidence interval for conformist.prop_informed
  geom_smooth(aes(y = conformist.prop_informed), method = "lm", color = "#7895b2", se = FALSE, size = 4) +
  # Regression line for conformist.prop_informed
  geom_abline(intercept = coef_model_2["(Intercept)"] + coef_model_2["density"] * mean(final_resultdf$density), slope = coef_model_2["avg_path_length"], color = "#7895b2", linetype = "dashed", size = 1.5) +
  # Add titles and labels
  labs(title = "",
       y = "Proportion",
       x = "Average Path Length") +
  # Classic theme with adjustments
  theme_classic() +
  theme(
    axis.text.x = element_text(size = rel(5), color = "black"),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.title.y = element_text(size = rel(3.7), color = "black", margin = margin(r = 25)),
    axis.title.x = element_text(size = rel(3.7), color = "black", margin = margin(t = 25)),
    panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA)  # transparent plot background
  )

ggsave("avg_path_length_plot_with_slopes.svg", plot = plot_avg_path_length, bg = "transparent", width = 10.5, height = 7.5)
