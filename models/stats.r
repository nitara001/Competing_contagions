library(igraph)
library(ggplot2)
library(mgcv)
library(marginaleffects)
#check for correlations between variables

final_resultdf<- read.csv("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\final_results.csv")
#infectione
final_resultdf$modularity<- as.numeric(final_resultdf$modularity)


summary(lm(prop_infected ~ qrel + avg_path_length, data= final_resultdf))
summary(lm(conformist.prop_informed ~ modularity + avg_path_length, data= final_resultdf))
model_1= mgcv::gam(prop_infected ~ s(qrel) + avg_module_size + clustering + avg_path_length, data = final_resultdf)
summary(model_1)
plot(model_1, select = 1, shade = TRUE)

model_2= mgcv::gam(conformist.prop_informed~ s(modularity) + avg_module_size + clustering + avg_path_length, data = final_resultdf)
summary(model_2)
plot(model_1, select = 1, shade = TRUE)


###-----------------------------------------------------


ggplot(final_resultdf, aes(x = modularity)) +
  geom_smooth(aes(y = prop_infected), method = "loess", color = "#2F4F4F", se = FALSE) + 
  geom_smooth(aes(y = proportional.prop_informed), method = "loess", color = "blue", se = FALSE) + 
  geom_smooth(aes(y = conformist.prop_informed), method = "loess", color = "green", se = FALSE) +
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
  geom_smooth(aes(y = difference), method = "loess", color = "#2F4F4F", se = FALSE) +  
  labs(title = "", x = "", y = "") +
  theme_classic() 
