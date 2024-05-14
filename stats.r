library(igraph)
library(ggplot2)

#infection

cor.test(infection_results_df$Modularity, infection_results_df$Num_Infected_at_tmax)
cor.test(infection_results_df$Avg_Module_Size, infection_results_df$Num_Infected_at_tmax)
summary(lm(Num_Infected_at_tmax ~ R0, data= infection_results_df))

ggplot(infection_results_df, aes(x = Modularity, y =Num_Infected_at_tmax, color = factor(R0))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = R0)) +
  labs(title = "Difference Between Modular and Random Networks in Mean Time for Disease/Information Spread",
          x = "Modularity",
       y = "Time to 75% Infection",
       color = "Transmission Probability") +
  theme_minimal()

is.na(results_df)

