
#"C:/Users/s2607536/AppData/Local/Programs/R/R-43~1.2/bin/R.exe"
##----------------------------------------------
#set up functions
library(igraph)
library(foreach)
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
  tmax <- 100
  netsize= igraph::V(rnet)
  
  t_tmps <- foreach(i = 1:50, .combine = c, .packages = "igraph") %dopar% {
    
    res2.r.1 <- do_spr(net = rnet, type = "informed", n_seeds = n_seeds, 
                       loc_seed = loc_seed, inform.type = inform.type, 
                       min_learn = min_learn, thr_steep = thr_steep, 
                       u = u_tmp, tmax = tmax, thresh = thresh, 
                       returnnets = FALSE, verbose = FALSE)
                       
      sort(sum(res2.r.1$informed), na.last = TRUE)
    }
  
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
  
  # Calculate ratio of module size to network size
  avg_module_size_ratio <- mean(community_sizes) / net_size
  
  return(list(
    modularity = modularity_q, 
    avg_module_size_ratio = avg_module_size_ratio
  ))
}
  
##---------------------------------------------------------------------------------------------------------------

#extract networks from asnr (need to run read_into_r code first)
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

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c('do_spr', 'get_u_from_s', 'prepseeds', 'calculate_modularity_largest_component'))
t1s.m.1 <- list()
t1s.r.1 <- list()
res1.r.2<- data.frame()
optim_result <- list()
resultsdfs <- list()

# Sample networks from larger subset so there is only one of repeated networks (possibly not doing this anymore)
unique_network_names <- unique(gsub("_\\d+$", "", names(extracted_networks)))
sampled_networks <- list()

for (network_name in unique_network_names) {
  matching_networks <- grep(paste0("^", network_name), names(extracted_networks), value = TRUE)
  # Exclude networks with 'sexual' in the title
  matching_networks <- matching_networks[!grepl("sexual", matching_networks, ignore.case = TRUE)]
  if (length(matching_networks) > 0) {
    sampled_network_name <- sample(matching_networks, 1)
    sampled_networks[[sampled_network_name]] <- extracted_networks[[sampled_network_name]]
  }
}
##---------------------------------------------------------------------------------------------------------------

##subset each network so it is just maintains its largest component(better for sims)
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

##---------------------------------------------------------------------------------------------------------------------------------------------------
# Run simulations
foreach(network_name = names(largest_components)[1:min(length(largest_components),2)]) %do% {
  network <- largest_components[[network_name]]
  random_network <- erdos.renyi.game(igraph::vcount(network), igraph::ecount(network), type = "gnm")
  
  t_R0=100
  rs=1.5/((mean(degree(network)^2)-mean(degree(network)))/mean(degree(network)))
  calculate gamma (r_adj)
  r_adj<-1-(1-rs)^(1/t_R0)
  
  t1s.r.1 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.1 <- do_spr(net = random_network, type = "infected", n_seeds = 1, loc_seeds = 'R', s = r_adj, tmax = 100, returnnets = FALSE, verbose = FALSE, inform.type = "conformist")
    sum(res1.r.1$infected)
  }

opt_result <- stats::optimize(get_u_from_s, interval = c(0, 1), rnet = random_network, infres = t1s.r.1, inform.type = 'conformist')

# Extract the minimum u value corresponding to the minimum bhs value
min_u <- opt_result$minimum

  t1s.r.2 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.2 <- do_spr(net = random_network, type = "informed", inform.type = "conformist", min_learn = 0.001, n_seeds = 1, loc_seeds = 'R', u = min_u, tmax = 100, returnnets = FALSE, verbose = TRUE)
    sum(res1.r.2$informed)
  }
  
  t1s.m.1 <- foreach(i = 1:50, .combine = rbind, .packages = "igraph") %do% {
    res1.m.1 <- do_spr(net = network, type = "both", n_seeds = 1, inform.type = "conformist", min_learn = 0.001, loc_seeds = 'R', s = r_adj, u = min_u, tmax = 100, returnnets = FALSE, verbose = TRUE)
    data.frame(network = network_name, numinfected = sum(res1.m.1$infected), numinformed = sum(res1.m.1$informed))
  }

  modularity_and_avg_module_size <- calculate_modularity_largest_component(network)
  modularity_score <- modularity_and_avg_module_size$modularity
  avg_module_size <- modularity_and_avg_module_size$avg_module_size_ratio
  global_clustering_coefficient <- igraph::transitivity(network, type = "global")
  network_size <- igraph::vcount(network)
  density<- igraph::edge_density(network)
  avg_path_length <- igraph::mean_distance(network)

#collecting data from the simulations
  t1s.r.1_str <- toString(t1s.r.1)
  t1s.r.2_str <- toString(t1s.r.2)
  t1s.m.1_str_infected<- toString(t1s.m.1$numinfected)
  t1s.m.1_str_informed<- toString(t1s.m.1$numinformed)
  resultdf <- data.frame(
    t1s.m.1_infected=  t1s.m.1_str_infected,
    t1s.m.1_informed=  t1s.m.1_str_informed,
    bhs = bhatt_coef(as.vector(t1s.m.1$numinfected), as.vector(t1s.m.1$numinformed)),
    modularity = modularity_score,
    avg_module_size = avg_module_size,
    clustering = global_clustering_coefficient,
    avg_path_length= avg_path_length,
    density= density,
    network_size = network_size,
    chosenu = min_u,
    t1s_r_1 = t1s.r.1_str,
    t1s_r_2= t1s.r.2_str,
    r_bhs = bhatt_coef(as.vector(t1s.r.1), as.vector(t1s.r.2))
  )
  
  resultsdfs[[network_name]] <- resultdf

}

new_results <- do.call(rbind, resultsdfs)
new_results$network <- rownames(new_results)
# Move 'network' column to the first position
new_results <- new_results[, c("network", setdiff(names(new_results), "network"))]
# Clean up parallel backend
stopCluster(cl)

###new dataframe type - so that we arent only using means but collect each result from each sim 
infected_df <- new_results %>%
  select(-t1s.m.1_informed) %>%
  mutate(type = "infection") %>%
  separate_rows(t1s.m.1_infected, sep = ", ") %>%
  rename(response = t1s.m.1_infected)

# Create a data frame for informed type
informed_df <- new_results %>%
  select(-t1s.m.1_infected) %>%
  mutate(type = "information") %>%
  separate_rows(t1s.m.1_informed, sep = ", ") %>%
  rename(response = t1s.m.1_informed)

combined_df <- bind_rows(infected_df, informed_df)
combined_df$network_size<- as.numeric(combined_df$network_size)
combined_df$response<- as.numeric(combined_df$response)
combined_df$outbreak_proportion<- (combined_df$response) / (combined_df$network_size)

