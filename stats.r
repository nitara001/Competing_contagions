library(igraph)
library(ggplot2)

#infection
summary(lm(Time_to_75_percent_infection ~ Transmission_Probability, data= infection_results_df))
summary(lm(Time_to_75_percent_infection ~ Modularity, data= infection_results_df))
summary(lm(Time_to_75_percent_infection ~ Avg_Module_Size , data =  infection_results_df))
summary(lm(Time_to_75_percent_infection ~ Clustering_coeff , data =  infection_results_df))

summary(lm(Modularity ~ Avg_Module_Size + Clustering_coeff + Network_Size + log(Number_subgroups), data= results_df))

ggplot(results_df, aes(x = Modularity, y = Time_to_75_percent_infection, color = factor(Transmission_Probability))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Transmission_Probability)) +
  labs(title = "Difference Between Modular and Random Networks in Mean Time for Disease/Information Spread",
          x = "Modularity",
       y = "Time to 75% Infection",
       color = "Transmission Probability") +
  theme_minimal()

is.na(results_df)

#information
summary(lm(Time_to_75_percent_informed ~ Modularity, data= results_df))
summary(lm(Time_to_75_percent_informed ~ Avg_Module_Size, data= results_df))
summary(lm(Time_to_75_percent_informed ~ Clustering_coeff, data= results_df))

infection_results_df[sapply(infection_results_df, is.infinite)]<- NA
