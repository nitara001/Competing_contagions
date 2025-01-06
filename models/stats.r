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
data<- read.csv("C:\\Users\\s2607536\\Downloads\\data.csv")

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


#plots
plot_model(
  model,
  type = "pred",
  terms = c("modularity_centered", "type"),
  title = "modularity",
  show.data = TRUE,
  dot.size = 1.5,
  dot.alpha = 0.1
) + ylim(0,1000)

plot_model(
  model,
  type = "pred",
  terms = c("density_centered", "type"),
  title = "Conditional Effects of Modularity and Type",
  show.data = TRUE,
  dot.size = 1.5,
  dot.alpha = 0.1
)


plot_model(
  model,
  type = "pred",
  terms = c("network_size_centered", "type"),
  title = "network size",
  show.data = TRUE,
  dot.size = 1.5,
  dot.alpha = 0.1
)


plot_model(
  model,
  type = "pred",
  terms = c("module_size_variation_centered", "type"),
  title = "module size variation",
  show.data = TRUE,
  dot.size = 1.5,
  dot.alpha = 0.1
)

plot_model(model, type = "re", title = "Random Effects for Species")


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

library(ggplot2)

ggplot(data, aes(x = modularity_centered, y = time, color = type)) +
  geom_point(alpha = 0.05, size = 1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  facet_grid(density_group ~ modularity_group, scales = "free") +
  labs(
    title = "",
    x = "Modularity (Centered)",
    y = "Time to 75% Infected",
    color = "Type"
  ) +
  ylim(0, 1500) + # Limit the y-axis to 0-1500
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top"
  )

#DIFFERNCES IN TIME
mean_times <- data %>%
  group_by(network, type) %>%
  summarise(mean_time = mean(time, na.rm = TRUE)) %>%
  pivot_wider(names_from = type, values_from = mean_time, names_prefix = "time_")
mean_times <- mean_times %>%
  mutate(time_diff = time_infection - time_information)
