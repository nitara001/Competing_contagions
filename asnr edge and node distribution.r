##distribution o
asnr<- read.csv("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Data\\asnr-master\\asnr-master\\Python_files\\Network_summary_master_file_10_20_21.csv")

#get ASNR networks
extracted_networks <- list()
for (i in seq_along(all_graphs)) {
  for (n in seq_along(all_graphs[[i]])) {
    graph <- all_graphs[[i]][[n]]
      # Add it to the extracted_networks list
      extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- graph
    }
  }

##edge and node distribution
edge_counts <- list()
node_counts<- list()
edge_dens<- list()

for (network_name in names(extracted_networks)) {
  network <- extracted_networks[[network_name]]
  edge_count <- ecount(network)
  edge_counts[[network_name]] <- edge_count
  node_count <- vcount(network)
  node_counts[[network_name]] <- node_count
  edge_den<- edge_density(network)
  edge_dens[[network_name]]<- edge_den
}

node_counts_vector<- unlist(node_counts)
edge_counts_vector <- unlist(edge_counts)
edge_dens_vector<- unlist(edge_dens)
hist(edge_counts_vector, main = "Distribution of Edge Counts", xlab = "Number of Edges", breaks= 500,  xlim= c(0, 2000))
hist(node_counts_vector, main = "Distribution of Node Counts", xlab = "Number of Nodes", breaks= 50)
hist(edge_dens_vector, main = "Distribution of Edge Densities")
###---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#what networks have attribute data 
# Check if any network has weighted edges
weighted_count <- sum(sapply(extracted_networks, is_weighted))

#what are all the unique attributes
unique_attributes <- unique(unlist(lapply(extracted_networks, list.vertex.attributes)))
print(unique_attributes)

# Check if any network has the "Sex" or "sex" vertex attribute
networks_with_sex <- sum(sapply(extracted_networks, function(network) {
  any(c("Sex", "sex", "SEX") %in% list.vertex.attributes(network))
}))

networks_with_age <- sum(sapply(extracted_networks, function(network) {
  any(c("std_male_age", "male_age", "age(days)", "age.class", "Age (2010)","Age (2009)", "Age","Age class-Gender", "age class", "Age(years)", "Birth year", "age", "age group" ) %in% list.vertex.attributes(network))
}))
print(weighted_count)
print(networks_with_sex)

#which have
networks_with_sex_and_age <- names(Filter(function(network) {
  any(c("Sex", "sex", "SEX") %in% igraph::list.vertex.attributes(network)) &&
    any(c("std_male_age", "male_age", "age(days)", "age.class", "Age (2010)","Age (2009)", "Age","Age class-Gender", "age class", "Age(years)", "Birth year", "age", "age group") %in% igraph::list.vertex.attributes(network))
}, extracted_networks))

print("Networks with both sex and age attributes:")
print(networks_with_sex_and_age)

##--------------------------------------------------------------------------------------------------------------------------------------------------
#captive and wild

table(asnr$population_type)
captive_asnr<- subset(asnr, population_type == "captive")
wild_asnr <- subset(asnr, population_type == "free-ranging")

###--------------------------------------------------------------------------------------------------------------------------------
##interaction type

table(asnr$interaction_type)
