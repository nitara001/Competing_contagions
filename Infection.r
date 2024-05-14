library(igraph)


# Function to calculate transmission probability from R0
calculate_transmission_probability <- function(R0, network) {
  # Calculate transmission probability based on the given R0 value and network structure
  # Average degree of all nodes in the network
  k <- mean(degree(network))
  # Average squared degree of all nodes in the network
  k2 <- mean(degree(network)^2) 
  # Time period over which R0 is calculated- 
  #remember this is now not prop to time of a simulation
  t <- 20
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

calculate_Q_and_avg_module_size(sealions)



#get ASNR networks
extracted_networks <- list()
for (i in seq_along(all_graphs)) {
  for (n in seq_along(all_graphs[[i]])) {
    graph <- all_graphs[[i]][[n]]
    if (igraph::vcount(graph) > 17 && igraph::ecount(graph) > 25) {
      # Add it to the extracted_networks list
      extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- graph
    }
  }
}
str(extracted_networks)
#print(extracted_networks)

############---
# Initialize lists to store results
sim_results <- list()
num_infected_list <- list()
infection_results_list <- list()

# Define the transmission probabilities
R0s <- c(1, 2)

# Iterate over each network
for (network_name in names(extracted_networks)[1:min(length(extracted_networks), 1250)]) { 
  network <- extracted_networks[[network_name]]
  global_clustering_coefficient <- transitivity(network, type = "global")
  network_size <- vcount(network)
  clusters <- clusters(network) #can choose weak or strong
  num_clusters <- clusters$no #number of connected components of a graph
  
  # Calculate transmission probability for each R0
  for (R0 in R0s) {
    prob <- calculate_transmission_probability(R0, network)
    
    # Perform simulation
    simulation_result <- do_spr(net = network, 
                                 type = "infected", 
                                 n_seeds = 1, 
                                 loc_seeds = "R", 
                                 s = prob, 
                                 u = NULL, 
                                 tmax = 100,
                                 thresh = NA, 
                                 infect.limit = NA, 
                                 inform.limit = NA, 
                                 inform.type = "proportional",
                                 min_learn = NULL, 
                                 thr_steep = NULL, 
                                 weighted = FALSE, 
                                 recovery = FALSE, 
                                 recoverytime = NULL, 
                                 recoverprob = NULL, 
                                 death = FALSE, 
                                 deathtime = NULL, 
                                 deathprob = NULL, 
                                 immune = FALSE, 
                                 immunetime = NULL, 
                                 immuneprob = NULL, 
                                 immuneendprob = NULL, 
                                 returnnets = FALSE, 
                                 verbose = TRUE)
    
    # Store simulation result and number of infected at tmax
    sim_results[[paste(network_name, R0, sep = "_")]] <- simulation_result
    num_infected_list[[paste(network_name, R0, sep = "_")]] <- sum(simulation_result$infected)
    
    # Fetch modularity score and average module size for the current network
    modularity_and_avg_module_size <- calculate_modularity_and_avg_module_size(network)
    modularity_score <- modularity_and_avg_module_size$modularity
    avg_module_size <- modularity_and_avg_module_size$avg_module_size
    
    # Append results to the list
    infection_results_list[[length(infection_results_list) + 1]] <- c(network_name, R0, sum(simulation_result$infected), modularity_score, avg_module_size, global_clustering_coefficient, network_size, num_clusters)
  }
}

# Combine all lists into one list
combined_list <- unlist(infection_results_list, recursive = FALSE)
result_matrix <- matrix(combined_list, ncol = length(infection_results_list[[1]]), byrow = TRUE)
infection_results_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
colnames(infection_results_df) <- c("Network", "R0", "Num_Infected_at_tmax", "Modularity", "Avg_Module_Size", "Clustering_coeff", "Network_Size", "Number_subgroups")
head(infection_results_df)


