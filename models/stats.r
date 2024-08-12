library(igraph)
library(ggplot2)
library(mgcv)
library(marginaleffects)
library(sjPlot)
library(viridis)
#check for correlations between variables

numeric_columns <- sapply(combined_df, is.numeric)
numeric_df <- combined_df[, numeric_columns]
cor(numeric_df, use = "complete.obs")

combined_df$network <- as.factor(combined_df$network)

model_random <- lmer(outbreak_proportion ~ type * (modularity + avg_module_size + avg_path_length + density) + (1| network), data = combined_df)
summary(model_random)
ranef(model_random)
anova(model_random)
plot_model(model_random, type= "re") #how each network varies in their outbreak responses 

mod_plot <- plot_model(model_random, type = "eff",  terms = c("modularity", "type"), ci.lvl = 0.95)
mod_plot<- mod_plot + 
  labs(title= "",
    x = "", 
       y = "", 
       color = "Type") +
  theme_blank()+
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20, margin = margin(t = 23)),
        axis.title.y = element_text(size = 20, margin = margin(r = 23)),
        axis.text = element_text(size = 29,color= "black"))+
  geom_line(size = 2.5) +
  scale_color_manual(values = c("infection" = "#008080","information" = "#d64f80"))+
  scale_fill_manual(values = c("infection" = "#008080", "information" = "#d64f80"))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) 
  
mod_plot

ggsave(filename = "/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/figures/mod_plot.svg", plot = mod_plot, device = "svg", dpi = 300)

module_sizeplot <- plot_model(model_random, type = "eff",  terms = c("avg_module_size", "type"), ci.lvl = 0.95)
module_sizeplot<- module_sizeplot + 
  labs(title= "",
       x = "", 
       y = "", 
       color = "Type",) +
  theme_blank()+
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20, margin = margin(t = 23)),
        axis.title.y = element_text(size = 20, margin = margin(r = 23)),
        axis.text = element_text(size = 29, color= "black"))+
  geom_line(size = 2.5) +
  scale_color_manual(values = c("infection" = "#008080","information" = "#d64f80"))+
  scale_fill_manual(values = c("infection" = "#008080", "information" = "#d64f80"))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) 
module_sizeplot
ggsave(filename = "/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/figures/module_sizeplot.svg", plot = module_sizeplot, device = "svg", dpi = 300)

path_plot <- plot_model(model_random, type = "eff",  terms = c("avg_path_length", "type"), ci.lvl = 0.95)
path_plot<- path_plot + 
  labs(title= "",
       x = "", 
       y = "", 
       color = "Type",) +
  theme_blank()+
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20, margin = margin(t = 23)),
        axis.title.y = element_text(size = 20, margin = margin(r = 23)),
        axis.text = element_text(size = 29, color= "black"))+
  geom_line(size = 2.5) +
  scale_color_manual(values = c("infection" = "#008080","information" = "#d64f80"))+
  scale_fill_manual(values = c("infection" = "#008080", "information" = "#d64f80"))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) 
path_plot
ggsave(filename = "/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/figures/path_plot.svg", plot = path_plot, device = "svg", dpi = 300)

density_plot <- plot_model(model_random, type = "eff",  terms = c("density", "type"), ci.lvl = 0.95)
density_plot<- density_plot + 
  labs(title= "",
       x = "", 
       y = "", 
       color = "Type",) +
  theme_blank()+
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20, margin = margin(t = 23)),
        axis.title.y = element_text(size = 20, margin = margin(r = 23)),
        axis.text = element_text(size = 29, color= "black"))+
  geom_line(size = 2.5) +
  scale_color_manual(values = c("infection" = "#008080","information" = "#d64f80"))+
  scale_fill_manual(values = c("infection" = "#008080", "information" = "#d64f80"))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) 
density_plot
ggsave(filename = "/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/figures/density_plot.svg", plot = path_plot, device = "svg", dpi = 300)

#

#model diagnostics 
#check if residuals are normally distributed
library(performance)
check_model(model_random)

qqnorm(residuals)
qqline(residuals, col = "red")

residuals <- resid(model_random)
fitted_values <- fitted(model_random)

# Plot Residuals vs. Fitted Values
plot(fitted_values, residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red", lty = 2)

# Load necessary library
library(ggplot2)
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#make plot for scaling parameter
gamma <- 1.0
eta <- 0.001
psi <- 10

# Function to calculate the probability of accepting information
probability_accepting_info <- function(proportion, gamma, eta, psi) {
  gamma * ((1 - 2 * eta) / (1 + exp(-psi * (1/gamma) * (proportion - 0.5))) + eta)
}

# Proportion of informed contacts (x-axis)
proportion_informed_contacts <- seq(0, 1, length.out = 500)

# Generate the y values for different gamma values
gamma_values <- c(1.0, 0.7, 0.4)
colors <- c('purple', 'blue', 'yellow')

data <- data.frame()

for (g in gamma_values) {
  probability <- sapply(proportion_informed_contacts, probability_accepting_info, gamma = g, eta = eta, psi = psi)
  df <- data.frame(Proportion = proportion_informed_contacts, Probability = probability, Gamma = as.factor(g))
  data <- rbind(data, df)
}
sig<-ggplot(data, aes(x = Proportion, y = Probability, color = Gamma)) +
  geom_line(size = 1.2) +
  labs(x = "", y = "", title ="") +
  scale_color_manual(values = colors) +
  theme_minimal()  +
  theme(
    legend.text = element_text(size = 18), # Increase legend text size
    legend.title = element_text(size = 18), # Increase legend title size
    axis.text = element_text(size = 18, color= "black"),
    legend.key.size = unit(3, 'lines') # Increase the size of the legend keys
  ) + geom_line(size = 2.5) 

sig
ggsave(filename = "/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/figures/sigmoid.svg", plot = sig, device = "svg", dpi = 300)


## ignore below for report 
model_1 <- lm(prop_infected ~ modularity + avg_path_length + clustering + density, data = final_resultdf)
model_2 <- lm(conformist.prop_informed ~ modularity + avg_path_length + clustering + density, data = final_resultdf)
coef_model_1 <- coef(model_1)
coef_model_2 <- coef(model_2)

###-----------------------------------------------------

# Plot for modularity
plot_modularity <- ggplot(final_resultdf, aes(x = modularity)) +
  geom_smooth(aes(y = prop_infected), method = "loess", color = "#c33a63", se = FALSE, size = 4) +
  geom_abline(intercept = coef_model_1["(Intercept)"] + coef_model_1["density"] * mean(final_resultdf$density), slope = coef_model_1["modularity"], color = "#c33a63", linetype = "dashed", size = 1.5) +
  geom_smooth(aes(y = conformist.prop_informed), method = "loess", color = "#7895b2", se = FALSE, size = 4) +
  geom_abline(intercept = coef_model_2["(Intercept)"] + coef_model_2["density"] * mean(final_resultdf$density), slope = coef_model_2["modularity"], color = "#7895b2", linetype = "dashed", size = 1.5) +
  labs(title = "",
       y = "Proportion",
       x = "Modularity") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = rel(5), color = "black"),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.title.y = element_text(size = rel(3.7), color = "black", margin = margin(r = 25)),
    axis.title.x = element_text(size = rel(3.7), color = "black", margin = margin(t = 25)),
    panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA)  # transparent plot background
  )
plot_modularity
ggsave("modularity_plot_with_slopes.svg", plot = plot_modularity, bg = "transparent", width = 10.5, height = 7.5)

# Plot for clustering
plot_clustering <- ggplot(final_resultdf, aes(x = clustering)) +
  geom_smooth(aes(y = prop_infected), method = "loess", color = "#c33a63", se = FALSE, size = 4) +
  geom_abline(intercept = coef_model_1["(Intercept)"] + coef_model_1["density"] * mean(final_resultdf$density), slope = coef_model_1["clustering"], color = "#c33a63", linetype = "dashed", size = 1.5) +
  geom_smooth(aes(y = conformist.prop_informed), method = "loess", color = "#7895b2", se = FALSE, size = 4) +
  geom_abline(intercept = coef_model_2["(Intercept)"] + coef_model_2["density"] * mean(final_resultdf$density), slope = coef_model_2["clustering"], color = "#7895b2", linetype = "dashed", size = 1.5) +

  labs(title = "",
       y = "Proportion",
       x = "Clustering") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = rel(5), color = "black"),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.title.y = element_text(size = rel(3.7), color = "black", margin = margin(r = 25)),
    axis.title.x = element_text(size = rel(3.7), color = "black", margin = margin(t = 25)),
    panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA)  # transparent plot background
  )

# Plot for avg_path_length
plot_avg_path_length <- ggplot(final_resultdf, aes(x = avg_path_length)) +
  geom_smooth(aes(y = prop_infected), method = "lm", color = "#c33a63", se = FALSE, size = 4) +
  geom_abline(intercept = coef_model_1["(Intercept)"] + coef_model_1["density"] * mean(final_resultdf$density), slope = coef_model_1["avg_path_length"], color = "#c33a63", linetype = "dashed", size = 1.5) +
  geom_smooth(aes(y = conformist.prop_informed), method = "lm", color = "#7895b2", se = FALSE, size = 4) +
  geom_abline(intercept = coef_model_2["(Intercept)"] + coef_model_2["density"] * mean(final_resultdf$density), slope = coef_model_2["avg_path_length"], color = "#7895b2", linetype = "dashed", size = 1.5) +
  labs(title = "",
       y = "Proportion",
       x = "Average Path Length") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = rel(5), color = "black"),
    axis.text.y = element_text(size = rel(5), color = "black"),
    axis.title.y = element_text(size = rel(3.7), color = "black", margin = margin(r = 25)),
    axis.title.x = element_text(size = rel(3.7), color = "black", margin = margin(t = 25)),
    panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
    plot.background = element_rect(fill = "transparent", color = NA)  # transparent plot background
  )

ggplot(final_resultdf, aes(x = avg_module_size)) +
  geom_point(aes(y = prop_infected), color = "#c33a63", size = 2) +
  geom_smooth(aes(y = prop_infected), method = "gam", color = "#c33a63", size = 1.5, fill = "#c33a63", alpha = 0.3, se = TRUE) +
  geom_point(aes(y = conformist.prop_informed), color = "blue", size = 2) +
  geom_smooth(aes(y = conformist.prop_informed), method = "gam", color = "blue", size = 1.5, fill = "blue", alpha = 0.3, se = TRUE) +
  xlim(0,75)

###

summary(mgcv::gam(prop_infected ~ s(modularity) + s(avg_path_length) + s(clustering) + s(avg_module_size)+density, data = final_resultdf_filtered))
summary(mgcv::gam(conformist.prop_informed ~ s(modularity) + s(avg_path_length) + s(clustering) + s(avg_module_size)+density, data = final_resultdf))

model<- lm(response ~ features|network, data= combined_df)
summary(model)
ggplot(combined_df, aes(x= modularity, y= outbreak_proportion, col= type)) +
  geom_point()+
  geom_abline(aes(intercept = 8.520e-02, slope= -1.393e-01, col= 'infection'))+
  geom_abline(aes(intercept = (8.520e-02 +6.472e-02), slope= (-1.393e-01-1.165e-01), col= 'information'))+
  geom_abline(aes(intercept = 8.520e-02, slope= 1.351e-02, col= 'darkgreen'))+
  geom_abline(aes(intercept = (8.520e-02 +6.472e-02), slope= ( 1.351e-02+ 3.017e-03), col= 'darkgreen'))+
  ylim= c(0,0.5)

  plotting random slopes
model_random_slope <- lmer(outbreak_proportion ~ type * (modularity + avg_module_size + avg_path_length + density) +
                             (modularity| network) + (avg_module_size | network) 
                           , data = combined_df)

plot_model(model_random_slope, type= "re")
##
# Extract random effects
random_effects <- ranef(model_random)$network
random_effects_df <- as.data.frame(random_effects)
random_effects_df$network <- rownames(random_effects_df)

# Simplify network names (without adding numbers for duplicates)
random_effects_df$simplified_network <- str_to_title(sub("(_.*)", "", random_effects_df$network))
axis_labels <- random_effects_df$simplified_network
names(axis_labels) <- rownames(random_effects_df)
axis_labels <- make.unique(axis_labels, sep = " ")
plot <- plot_model(model_random, type = "re", show.values = FALSE, value.offset = .3, axis.labels = axis_labels)
plot <-plot+ theme(axis.text = element_text(size = 28),
             axis.text.y = element_text(size = 17),
             axis.title.x = element_text(size = 24),
             axis.title.y = element_text(size = 24) + 
  theme_minimal() 
plot
ggsave(filename = "/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/figures/re_plot.svg", plot = plot, device = "svg", dpi = 300)



gam_model <- gam(response ~ type * modularity + s(avg_module_size) + clustering + avg_path_length + density + s(network, bs = "re"), data = combined_df)
summary(gam_model)
plot.gam(gam_model, select=3)