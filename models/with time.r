
#"C:/Users/s2607536/AppData/Local/Programs/R/R-43~1.2/bin/R.exe"
##----------------------------------------------
#set up functions
library(igraph)
library(foreach)
library(dplyr)
library(tidyr)
library(tibble)
library(doParallel)
#lower bhatt coef means more similar distributions
bhatt_coef<-function(vec1,vec2){
  mat1<-as.matrix(vec1[complete.cases(vec1)])
  mat2<-as.matrix(vec2[complete.cases(vec2)])
  mn1<-mean(mat1)
  mn2<-mean(mat2)
  mn_dif<-mn1-mn2
  cov1<-cov(mat1)
  cov2<-cov(mat2)
  p<-(cov1+cov2)/2
  bh<-0.125*t(mn_dif)*p^(-1)*mn_dif+0.5*log(det(p)/sqrt(det(cov1)*det(cov2)))
  return(bh)
}

#function used when optimising u 
get_u_from_s <- function(u_tmp, rnet, infres, inform.type) {
  thresh <- 1.0
  n_seeds <- 1
  loc_seed <- "R"
  min_learn <- 0.001
  thr_steep <- 10
  tmax <- 3000
  netsize <- igraph::vcount(rnet)
  
  # Calculate the threshold point for 75% of the nodes
  thresh_point <- 0.75 * netsize
  
  # Run the simulation 50 times in parallel to find the minimum time when 75% of nodes are informed
  t_tmps <- foreach(i = 1:50, .combine = c, .packages = "igraph") %dopar% {
    res2.r.1 <- do_spr(net = rnet, type = "informed", n_seeds = n_seeds, 
                       loc_seed = loc_seed, inform.type = inform.type, 
                       min_learn = min_learn, thr_steep = thr_steep, 
                       u = u_tmp, tmax = tmax, thresh = thresh, 
                       returnnets = FALSE, verbose = FALSE)
    
    total_informed <- sum(res2.r.1$informed)
  total_nodes <- igraph::vcount(rnet)

  # Check if 75% infection threshold is reached
  if (total_informed >= 0.75 * total_nodes) {
    # Filter times to those greater than zero
    valid_times <- res2.r.1$wheninformed[res2.r.1$wheninformed > 0]

    # Find the minimum time when the cumulative number of infected nodes reaches 75%
    min_time <- sort(valid_times)[cumsum(res2.r.1$informed[res2.r.1$wheninformed > 0]) >= 0.75 * total_nodes][1]

  } else {
    paste("75% informed threshold never reached; only", (total_informed), "informed")
  }}
  
  # Calculate Bhattacharyya coefficient using the new measure
  bh <- bhatt_coef(infres, t_tmps)
  cat(paste(as.vector(bh), "\n"))
  return(as.vector(bh))
}

##---------------------------------------------------------------------------------------------------------------
calculate_modularity_largest_component <- function(network) {
  # Get the largest connected component
  components <- clusters(network, mode = "weak")
  vert_ids <- V(network)[components$membership == which.max(components$csize)]
  largest_component <- induced_subgraph(network, vert_ids)
  
  # Calculate modularity and community sizes
  community <- cluster_louvain(largest_component, weights = E(largest_component)$weight)
  modularity_q <- modularity(community)
  community_sizes <- sizes(community)
  net_size <- vcount(largest_component)
  
  # Calculate ratio of module size to network
  avg_module_size_ratio <- mean(community_sizes) / net_size
  
  return(list(
    modularity = modularity_q, 
    avg_module_size_ratio = avg_module_size_ratio
  ))
}

##---------------------------------------------------------------------------------------------------------------
#extract networks from asnr (need to run read_into_r code)
extracted_networks <- list()
for (i in seq_along(all_graphs)) {
  for (n in seq_along(all_graphs[[i]])) {
    graph <- all_graphs[[i]][[n]]
    if (igraph::vcount(graph) > 10 &&  igraph::is_connected(graph)) {
      # Add it to the extracted_networks list only if the graph is connected
      extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- graph
    }
  }
}

##---------------------------------------------------------------------------------------------------------------

##subset each network so it is just maintains its largest component
largest_components <- list()
for (network_name in names(extracted_networks)) {
  network <- extracted_networks[[network_name]]
  components <- clusters(network, mode = "weak")
  biggest_cluster_id <- which.max(components$csize)
  vert_ids <- V(network)[components$membership == biggest_cluster_id]
  largest_component <- induced_subgraph(network, vert_ids)
  largest_components[[network_name]] <- largest_component
}
##---------------------------------------------------------------------------------------------------------------

# Sample networks from larger subset so there is only one of repeated networks
#unique_network_names <- unique(gsub("_\\d+$", "", names(extracted_networks)))
#sampled_networks <- list()

#for (network_name in unique_network_names) {
  #matching_networks <- grep(paste0("^", network_name), names(extracted_networks), value = TRUE)
  # Exclude networks with 'sexual' in the title
  #matching_networks <- matching_networks[!grepl("sexual", matching_networks, ignore.case = TRUE)]
  #if (length(matching_networks) > 0) {
    #sampled_network_name <- sample(matching_networks, 1)
    #sampled_networks[[sampled_network_name]] <- extracted_networks[[sampled_network_name]]
  #}
#}

##---------------------------------------------------------------------------------------------------------------------------------------------------
# Run simulations

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
clusterExport(cl, c('do_spr', 'get_u_from_s', 'prepseeds', 'calculate_modularity_largest_component'))
t1s.m.1 <- list()
simulation_results<- list()
opt_result <- list()
resultsdfs <- list()
set.seed(123) 

foreach(network_name = names(largest_components)[1:min(length(largest_components), 10)]) %do% {
  network <- largest_components[[network_name]]
  random_network <- erdos.renyi.game(igraph::vcount(network), igraph::ecount(network), type = "gnm")

  t_R0 <- 100
  rs <- 1.5 / ((mean(degree(network)^2) - mean(degree(network))) / mean(degree(network)))
  r_adj <- 1 - (1 - rs)^(1 / t_R0)

  # Run the first simulation for infection
  t1s.r.1 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.1 <- do_spr(net = random_network, type = "infected", n_seeds = 1, loc_seeds = 'R', 
                       s = r_adj, tmax = 3000, returnnets = FALSE, verbose = TRUE, inform.type = "conformist")
    total_infected <- sum(res1.r.1$infected)
    total_nodes <- igraph::vcount(random_network)

    if (total_infected >= 0.75 * total_nodes) {
      valid_times <- res1.r.1$wheninfected[res1.r.1$wheninfected > 0]
      min_time <- sort(valid_times)[cumsum(res1.r.1$infected[res1.r.1$wheninfected > 0]) >= 0.75 * total_nodes][1]

      if (!is.na(min_time) && is.finite(min_time)) {
        min_time
      } else {
        paste("75% infection threshold never reached; only", total_infected, "infected")
      }
    }
  }

  # Optimize to find the best 'u' value
  opt_result <- stats::optimize(get_u_from_s, interval = c(0, 1), rnet = random_network, infres = t1s.r.1, inform.type = 'conformist')
  min_u <- opt_result$minimum

  # Run the second simulation for information spread
  t1s.r.2 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.2 <- do_spr(net = random_network, type = "informed", inform.type = "conformist", min_learn = 0.001, 
                       n_seeds = 1, loc_seeds = 'R', u = min_u, tmax = 3000, returnnets = FALSE, verbose = TRUE)
    total_informed <- sum(res1.r.2$informed)
    total_nodes <- igraph::vcount(random_network)

    if (total_informed >= 0.75 * total_nodes) {
      valid_times <- res1.r.2$wheninformed[res1.r.2$wheninformed > 0]
      min_time <- sort(valid_times)[cumsum(res1.r.2$informed[res1.r.2$wheninformed > 0]) >= 0.75 * total_nodes][1]

      if (!is.na(min_time) && is.finite(min_time)) {
        min_time
      } else {
        paste("75% informed threshold never reached; only", total_informed, "informed")
      }
    }
  }

  # Combined simulations
  t1s.m.1 <- foreach(i = 1:50, .combine = rbind, .packages = "igraph") %do% {
    res1.m.1 <- do_spr(net = network, type = "both", n_seeds = 1, inform.type = "conformist", min_learn = 0.001, 
                       loc_seeds = 'R', s = r_adj, u = min_u, tmax = 3000, returnnets = FALSE, verbose = TRUE)
    total_infected <- sum(res1.m.1$infected)
    total_informed <- sum(res1.m.1$informed)
    total_nodes <- igraph::vcount(network)

    min_time_75_infected <- if (total_infected >= 0.75 * total_nodes) {
      valid_times_infected <- res1.m.1$wheninfected[res1.m.1$wheninfected > 0]
      min_time_75_infected <- sort(valid_times_infected)[cumsum(res1.m.1$infected[res1.m.1$wheninfected > 0]) >= 0.75 * total_nodes][1]
      if (!is.na(min_time_75_infected) && is.finite(min_time_75_infected)) {
        min_time_75_infected
      } else {
        paste("75% infected threshold never reached; only", total_infected, "infected")
      }
    }

    min_time_75_informed <- if (total_informed >= 0.75 * total_nodes) {
      valid_times_informed <- res1.m.1$wheninformed[res1.m.1$wheninformed > 0]
      min_time_75_informed <- sort(valid_times_informed)[cumsum(res1.m.1$informed[res1.m.1$wheninformed > 0]) >= 0.75 * total_nodes][1]
      if (!is.na(min_time_75_informed) && is.finite(min_time_75_informed)) {
        min_time_75_informed
      } else {
        paste("75% informed threshold never reached; only", total_informed, "informed")
      }
    }

    # Collect results
    data.frame(
      network = network_name, 
      min_time_75_infected = min_time_75_infected, 
      min_time_75_informed = min_time_75_informed
    )
  }

  # Calculations for network properties
  modularity_and_avg_module_size <- calculate_modularity_largest_component(network)
  modularity_score <- modularity_and_avg_module_size$modularity
  avg_module_size <- modularity_and_avg_module_size$avg_module_size_ratio
  qrel <- modularity_and_avg_module_size$relative_modularity
  global_clustering_coefficient <- igraph::transitivity(network, type = "global")
  network_size <- igraph::vcount(network)
  density <- igraph::edge_density(network)
  avg_path_length <- igraph::mean_distance(network)

  # Convert variables to strings
  t1s.r.1_str <- toString(t1s.r.1)
  t1s.r.2_str <- toString(t1s.r.2)
  t1s.m.1_str_infected <- toString(t1s.m.1$min_time_75_infected)
  t1s.m.1_str_informed <- toString(t1s.m.1$min_time_75_informed)

  # Prepare the result dataframe
  resultdf <- data.frame(
    t1s.m.1_infected = t1s.m.1_str_infected,
    t1s.m.1_informed = t1s.m.1_str_informed,
    bhs = bhatt_coef(as.numeric(t1s.m.1$min_time_75_infected), as.numeric(t1s.m.1$min_time_75_informed)),
    modularity = modularity_score,
    avg_module_size = avg_module_size,
    clustering = global_clustering_coefficient,
    density = density,
    avg_path_length = avg_path_length,
    network_size = network_size,
    chosenu = min_u,
    t1s_r_1 = t1s.r.1_str,
    t1s_r_2 = t1s.r.2_str,
    r_bhs = bhatt_coef(as.numeric(t1s.r.1), as.numeric(t1s.r.2))
  )
  resultsdfs[[network_name]] <- resultdf
}



new_results <- do.call(rbind, resultsdfs)
new_results$network <- rownames(new_results)
new_results <- new_results[, c("network", setdiff(names(new_results), "network"))]
stopCluster(cl)

infected_df <- new_results %>%
  select(-t1s.m.1_informed) %>%
  mutate(type = "infection") %>%
  separate_rows(t1s.m.1_infected, sep = ", ") %>%
  rename(time = t1s.m.1_infected)

informed_df <- new_results %>%
  select(-t1s.m.1_infected) %>%
  mutate(type = "information") %>%
  separate_rows(t1s.m.1_informed, sep = ", ") %>%
  rename(time = t1s.m.1_informed)

combined_df <- bind_rows(infected_df, informed_df)
View(combined_df)