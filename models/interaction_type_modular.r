results<- read.csv("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\results.csv")
networks_with_types<- read.csv("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\sampled_networks.csv")
merged_df <- merge(networks_with_types, results, by = "names")


summary(aov(Modularity ~ interactions, data = merged_df))
summary(aov(Modularity ~ population_type, data = merged_df))


# Or create violin plot
ggplot(merged_df, aes(x = interactions, y = Modularity)) +
  geom_violin() +
  labs(x = "Interaction Type", y = "Modularity") +
  ggtitle("Distribution of Modularity by Interaction Type")
