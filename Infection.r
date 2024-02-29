library(igraph)


# Function to calculate transmission probability from R0
calculate_transmission_probability <- function(R0, network) {
  # Calculate transmission probability based on the given R0 value and network structure
  # Average degree of all nodes in the network
  k <- mean(degree(network))
  # Average squared degree of all nodes in the network
  k2 <- mean(degree(network)^2) 
  # Time period over which R0 is calculated
  t <- 100
  # Calculate r, the transmission probability for a time period t in that network
  r <- R0 / ((k2 - k) / k)
  # Calculate beta, the per timestep transmission probability per connection
  beta <- 1 - (1 - r)^(1/t)
  return(beta)
}

calculate_modularity_and_avg_module_size <- function(network) {
    community <- fastgreedy.community(network)
    modularity <- modularity(community)
    avg_module_size <- mean(sizes(community))
    return(list(modularity = modularity, avg_module_size = avg_module_size))
}

#get ASNR networks
extracted_networks <- list()
for (i in seq_along(all_graphs)) {
  for (n in seq_along(all_graphs[[i]])) {
    graph <- all_graphs[[i]][[n]]
    if (igraph::vcount(graph) > 20 && igraph::ecount(graph) > 20) {
      # Add it to the extracted_networks list
      extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- graph
    }
  }
}
str(extracted_networks)
#print(extracted_networks)

############---

# Perform simulations
sim_results <- list()
num_infected_list <- list()
R0s= c(1, 2, 3)

for (i in 1:min(length(extracted_networks), 500)) { 
  network_name <- names(extracted_networks)[i]
  network <- extracted_networks[[i]]
  
  for (R0 in R0s) {

    prob= calculate_transmission_probability(R0, network)

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
    
    num_infected <- sum(simulation_result$infected)
    sim_results[[paste(network_name, R0, sep = "_")]] <- simulation_result
    num_infected_list[[paste(network_name, R0, sep = "_")]] <- num_infected
  }}


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to calculate the time taken for 75% of nodes to become infected

calculate_75_percent_time <- function(simulation_result) {
  target_nodes <- 0.75 * sum(simulation_result$infected) #when 75 percent of total are infected
  cum_infected <- cumsum(simulation_result$infected)
  max(simulation_result$wheninfected[cum_infected >= target_nodes])
}

# store results and network characteristics
infection_results_df <- data.frame(
  Network = character(),
  Transmission_Probability = numeric(),
  Time_to_75_percent_infection = numeric(),
  Modularity = numeric(),
  Avg_Module_Size = numeric(),
  Clustering_coeff= numeric(),
  Network_Size = numeric(),
  Number_subgroups= numeric(),
  stringsAsFactors = FALSE
)

# Iterate over each network
for (network_name in names(extracted_networks)[1:min(length(extracted_networks), 500)]) { 
  network <- extracted_networks[[network_name]]
  global_clustering_coefficient <- transitivity(network, type = "global")
  network_size <- vcount(network)
  clusters <- clusters(network)
  num_clusters <- clusters$no

  # Iterate over each R0
  for (R0 in R0s) {
    prob <- calculate_transmission_probability(R0, network)
    key <- paste(network_name, R0, sep = "_")  # Adjust key to include R0
    simulation_result <- sim_results[[key]]
    num_infected <- num_infected_list[[key]]
    time_to_75_percent <- calculate_75_percent_time(simulation_result)
    
    # Fetch modularity score and average module size for the current network
    modularity_and_avg_module_size <- calculate_modularity_and_avg_module_size(network)
    modularity_score <- modularity_and_avg_module_size$modularity
    avg_module_size <- modularity_and_avg_module_size$avg_module_size
    
    # Create a data frame for the current iteration
    df <- data.frame(
      Network = network_name,
      Transmission_Probability = prob,
      Time_to_75_percent_infection = time_to_75_percent,
      Modularity = modularity_score,
      Avg_Module_Size = avg_module_size,
      Clustering_coeff = global_clustering_coefficient,
      Network_Size = network_size,
      Number_subgroups = num_clusters
    )
    
    # Append the data frame to the results data frame
    infection_results_df <- rbind(infection_results_df, df)
  }
}

# Remove the first row (initialized empty row)
infection_results_df <- infection_results_df[-1, ]

#infection_results_df <- na.omit(infection_results_df)

# Print the results data frame
str(infection_results_df)
####-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# divide transmission probs into 3 
range(infection_results_df$Transmission_Probability)
low_prob_range <- c(0, 0.003)
intermediate_prob_range <- c(0.003, 0.01)
high_prob_range <- c(0.01, 1)


low_prob_df <- infection_results_df[infection_results_df$Transmission_Probability >= low_prob_range[1] & 
                          infection_results_df$Transmission_Probability < low_prob_range[2], ]
intermediate_prob_df <- infection_results_df[infection_results_df$Transmission_Probability >= intermediate_prob_range[1] & 
                                   infection_results_df$Transmission_Probability < intermediate_prob_range[2], ]
high_prob_df <- infection_results_df[infection_results_df$Transmission_Probability >= high_prob_range[1] & 
                            infection_results_df$Transmission_Probability <= high_prob_range[2], ]

nrow(infection_results_df)
nrow(low_prob_df)
nrow(intermediate_prob_df)
nrow(high_prob_df)
