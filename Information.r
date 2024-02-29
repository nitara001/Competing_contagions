library(igraph)

calculate_modularity_and_avg_module_size <- function(network) {
    community <- fastgreedy.community(network)
    modularity <- modularity(community)
    avg_module_size <- mean(sizes(community))
    return(list(modularity = modularity, avg_module_size = avg_module_size))
}

#get ASNR networks
extracted_networks <- list()

# Iterate over each list element in all_graphs
for (i in seq_along(all_graphs)) {
  # Iterate over each graph within the list element
  for (n in seq_along(all_graphs[[i]])) {
    # Extract the graph and add it to the extracted_networks list
    extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- all_graphs[[i]][[n]]
  }
}

#print(extracted_networks)


# Perform simulations
sim_results <- list()
num_informed_list <- list()

for (i in 1:min(length(extracted_networks), 400)) { 
  network_name <- names(extracted_networks)[i]
  network <- extracted_networks[[i]]

    simulation_result <- do_spr(net = network, 
                                 type = "informed", 
                                 n_seeds = 1, 
                                 loc_seeds = "R", 
                                 s = NULL, 
                                 u = 0.5, 
                                 tmax = 100,
                                 thresh = 0.75, 
                                 infect.limit = NA, 
                                 inform.limit = NA, 
                                 inform.type = "conformist",
                                 min_learn = 0.01, 
                                 thr_steep = 10, 
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
    
    num_informed <- sum(simulation_result$informed)
    sim_results[[paste(network_name, sep = "_")]] <- simulation_result
    num_informed_list[[paste(network_name, R0, sep = "_")]] <- num_informed
  }


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate the time taken for 75% of nodes to become infected


calculate_75_percent_time_informed <- function(simulation_result) {
  target_nodes <- 0.75 * sum(simulation_result$informed) #when 75 percent of total are infected
  cum_infected <- cumsum(simulation_result$informed)
  max(simulation_result$wheninformed[cum_infected >= target_nodes])
}



# Initialize results data frame
results_df <- data.frame(
  Network = character(),
  Time_to_75_percent_informed = numeric(),
  Modularity = numeric(),
  Avg_Module_Size = numeric(),
  Clustering_coeff= numeric(),
  Network_Size = numeric(),
  Number_subgroups= numeric(),
  stringsAsFactors = FALSE
)

# Iterate over each network
for (network_name in names(extracted_networks)[1:min(length(extracted_networks), 400)]) { 
  network <- extracted_networks[[network_name]]
  global_clustering_coefficient <- transitivity(network, type = "global")
  network_size<- vcount(network)
  clusters <- clusters(network); num_clusters <- clusters$no


    prob <- calculate_transmission_probability(R0, network)
    key <- paste(network_name, sep = "_")  # Adjust key to include R0
    simulation_result <- sim_results[[key]]
    num_infected <- num_infected_list[[key]]
    time_to_75_percent <- calculate_75_percent_time_informed(simulation_result)
    
    # Fetch modularity score and average module size for the current network
    modularity_and_avg_module_size <- calculate_modularity_and_avg_module_size(network)
    modularity_score <- modularity_and_avg_module_size$modularity
    avg_module_size <- modularity_and_avg_module_size$avg_module_size
    
    # Add a row to the results data frame
    results_df <- rbind(results_df, data.frame(
      Network = network_name,
      Time_to_75_percent_informed = time_to_75_percent,
      Modularity = modularity_score,
      Avg_Module_Size = avg_module_size,
      Clustering_coeff= global_clustering_coefficient,
      Network_Size= network_size,
      Number_subgroups=num_clusters
    ))
  }

# Remove the first row (initialized empty row)
results_df <- na.omit(results_df)

# Print the results data frame
str(results_df)
