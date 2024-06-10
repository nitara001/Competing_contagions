library(igraph)
library(ggplot2)
library(mgcv)
library(marginaleffects)
#check for correlations between variables
infection_results_df$Modularity<- as.numeric(infection_results_df$Modularity)
cor.test(infection_results_df$Modularity, infection_results_df$Network_Size)
cor.test(final_resultdf$Avg_Module_Size, final_resultdf$Network_Size)
cor.test(infection_results_df$Avg_Module_Size, infection_results_df$Modularity)
cor.test(infection_results_df$Modularity, infection_results_df$Mean_Degree)
cor.test(infection_results_df$Modularity, infection_results_df$Std_Degree)
#infectione
final_resultdf$modularity<- as.numeric(final_resultdf$modularity)


summary(lm(Prop_infected ~ modularity, data= final_resultdf))
summary(lm(Prop_informed ~ modularity, data= final_resultdf))
model_1= mgcv::gam(Prop_infected ~ s(modularity) + s(avg_module_size), data = final_resultdf)
summary(model_1)
plot(model_1, select = 2, shade = TRUE)
abline(h = 0, lty = 'dashed')

###-----------------------------------------------------


ggplot(final_resultdf, aes(x = modularity)) +
  geom_smooth(aes(y = prop_infected), method = "loess", color = "#2F4F4F", se = FALSE) +  # Disable shading
  geom_smooth(aes(y = proportional.prop_informed), method = "loess", color = "blue", se = FALSE) +  # Disable shading
  labs(title = "", x = "", y = "") +
  theme_classic() +
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

  ggplot(final_resultdf, aes(x = modularity)) +
  geom_smooth(aes(y = difference), method = "loess", color = "#2F4F4F", se = FALSE) +  # Disable shading
  labs(title = "", x = "", y = "") +
  theme_classic() 
