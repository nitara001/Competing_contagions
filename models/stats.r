library(igraph)
library(ggplot2)
library(mgcv)
library(marginaleffects)
#check for correlations between variables
infection_results_df$Modularity<- as.numeric(infection_results_df$Modularity)
cor.test(infection_results_df$Modularity, infection_results_df$Network_Size)
cor.test(infection_results_df$Avg_Module_Size, infection_results_df$Network_Size)
cor.test(infection_results_df$Avg_Module_Size, infection_results_df$Modularity)
cor.test(infection_results_df$Modularity, infection_results_df$Mean_Degree)
cor.test(infection_results_df$Modularity, infection_results_df$Std_Degree)
#infection
result_df$Modularity<- as.numeric(result_df$Modularity)

model_1= mgcv::gam(Infected_Mean ~ s(Modularity), data = result_df)
plot(model_1, select = 2, shade = TRUE)
abline(h = 0, lty = 'dashed')

###-----------------------------------------------------


p= ggplot(result_df, aes(x = Modularity, y =Prop_infected)) +
  geom_smooth(method = "loess", fill='#982165', color = "#2F4F4F", level=0.75)+
  labs(title = "",
          x = "",
       y = "") +
  theme_classic()+
  theme(axis.text.x = element_text(size = rel(1.8)),
        axis.text.y = element_text(size = rel(1.8)))

p +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )