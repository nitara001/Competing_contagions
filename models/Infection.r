library(igraph)
library(foreach)

# Function to calculate transmission probability from R0
calculate_transmission_probability <- function(R0, network) {
  # Calculate transmission probability based on the given R0 value and network structure
  # Average degree of all nodes in the network
  k <- mean(degree(network))
  # Average squared degree of all nodes in the network
  k2 <- mean(degree(network)^2) 
  # Time period over which R0 is calculated- 
  #remember this is now not prop to time of a simulation
  t <- 100
  # Calculate r, the transmission probability for a time period t in that network
  r <- R0 / ((k2 - k) / k)
  # Calculate beta, the per timestep transmission probability per connection
  beta <- 1 - (1 - r)^(1/t)
  return(beta)
}


calculate_modularity_and_avg_module_size <- function(network) {
    community <- cluster_louvain(network, weights= E(network)$weight)
    modularity<- modularity(community)
    avg_module_size <- mean(sizes(community))
    return(list(modularity = modularity, avg_module_size = avg_module_size))
}

calculate_Q_and_avg_module_size <- function(network) {
    community <- cluster_louvain(network, weights = E(network)$weight)
    modularity <- modularity(community)
    avg_module_size <- mean(sizes(community))
    # Calculate L - total edges
    L <- ecount(network)
    # Calculate Lk for each subgroup
    Lk <- numeric(length(community))
    for (i in 1:length(community)) {
        subgraph <- induced_subgraph(network, which(membership(community) == i))
        Lk[i] <- ecount(subgraph)
    }
    
    # Calculate Qmax
    Qmax <- sum(Lk / L * (1 - Lk / L))
    return(list(modularity = modularity, avg_module_size = avg_module_size, Qmax = Qmax))
}

##----------------------------------------------------------------------------------------------------

#get ASNR networks
extracted_networks <- list()
for (i in seq_along(all_graphs)) {
  for (n in seq_along(all_graphs[[i]])) {
    graph <- all_graphs[[i]][[n]]
    if (igraph::vcount(graph) > 10 && igraph::ecount(graph) > 15) {
      # Add it to the extracted_networks list
      extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- graph
    }
  }
}
#str(extracted_networks)
#print(extracted_networks)

############---=------------------------------------------------------------------------------

# Sample networks from larger subset so there is only one of repeated networks
unique_network_names <- unique(gsub("_\\d+$", "", names(extracted_networks)))
sampled_networks <- list()
for (network_name in unique_network_names) {
  matching_networks <- grep(paste0("^", network_name), names(extracted_networks), value = TRUE)
  sampled_network_name <- sample(matching_networks, 1)
  sampled_networks[[sampled_network_name]] <- extracted_networks[[sampled_network_name]]
}


result_df <- data.frame(
  Network = character(0),
  Infected_Mean = numeric(0),
  Prop_Infected= numeric(0),
  Infected_Sd = numeric(0),
  Global_Clustering_Coefficient = numeric(0),
  Network_Size = numeric(0),
  Number_of_Clusters = numeric(0),
  Mean_Degree = numeric(0),
  SD_Degree = numeric(0),
  Average_Path_Length = numeric(0),
  Diameter = numeric(0),
  Modularity = numeric(0),
  Qmax=numeric(0),
  stringsAsFactors = FALSE
)

# Loop over each network
foreach(network_name = names(sampled_networks)[1:min(length(sampled_networks), 57)], .combine = 'c') %do% {
  network <- sampled_networks[[network_name]]

  # Calculate transmission probability
  prob <- calculate_transmission_probability(1.5, network)

  t1s.r.1 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %dopar% {
    res1.r.1 <- do_spr(net = network, 
                        type = "infected", 
                        n_seeds = 1, 
                        loc_seeds = 'R', 
                        s = prob, 
                        tmax = 200, 
                        returnnets = F, 
                        verbose = T, 
                        inform.type = "conformist")
    sum(res1.r.1$infected)
  }

  # Calculate modularity and average module size
  modularity_and_avg_module_size <- calculate_Q_and_avg_module_size(network)
  modularity_score <- modularity_and_avg_module_size$modularity
  avg_module_size <- modularity_and_avg_module_size$avg_module_size
  infected_mean <- mean(t1s.r.1, na.rm = TRUE)
  global_clustering_coefficient <- transitivity(network, type = "global")
  network_size <- vcount(network)
  clusters <- clusters(network) #can choose weak or strong
  num_clusters <- clusters$no #number of connected components of a graph
  degrees <- degree(network)
  mean_degree <- mean(degrees)
  sd_degree <- sd(degrees)
  avg_path_length <- mean_distance(network)
  diameter <- diameter(network)


  # Add the results to the dataframe
  result_df <- rbind(result_df, data.frame(
    Network = network_name,
    Infected_Mean = infected_mean,
    Prop_infected= infected_mean/ vcount(network),
    Infected_Sd= sd(t1s.r.1, na.rm= TRUE),
    Global_Clustering_Coefficient = global_clustering_coefficient,
    Network_Size = network_size,
    Number_of_Clusters = num_clusters,
    Mean_Degree = mean_degree,
    SD_Degree = sd_degree,
    Average_Path_Length = avg_path_length,
    Diameter = diameter,
    Modularity = modularity_score
  ))
}

# Reset row names
rownames(result_df) <- NULL

write.csv(result_df, "C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\results.csv")
