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

#linear model
model<- lme4::lmer(
  time ~ type * (modularity_centered + density_centered 
  + network_size_centered + 
                 module_size_variation_centered + num_communities) + 
    (1 | species/network), 
  data = data
)
summary(model)

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


