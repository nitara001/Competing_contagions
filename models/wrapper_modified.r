
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


calculate_modularity_largest_component <- function(network) {
  # Get the largest connected component
  components <- clusters(network, mode = "weak")
  biggest_cluster_id <- which.max(components$csize)
  vert_ids <- V(network)[components$membership == biggest_cluster_id]
  largest_component <- induced_subgraph(network, vert_ids)
  # Calculate modularity on the largest component
  community <- cluster_louvain(largest_component, weights = E(largest_component)$weight)
  modularity_q <- modularity(community)
  # Calculate average module size
  avg_module_size <- mean(sizes(community))
  # Get membership 
  mem <- membership(community)
  net_size <- length(mem)
  # Create community membership matrix
  mem_mat <- matrix(0, nr = net_size, nc = net_size)
  for(i in 1:net_size) {
    for(j in 1:net_size) {
      if(mem[i] == mem[j]) {
        mem_mat[i, j] <- 1
      }
    }
  }
  diag(mem_mat) <- 0
  
  # Create graph from the membership matrix
  qmat <- graph.adjacency(mem_mat, mode = "undirected", diag = FALSE)
  # Calculate maximum modularity
  community_qmat <- cluster_louvain(qmat)
  modularity_qmax <- modularity(community_qmat)
  # Calculate relative modularity
  qrel <- modularity_q / modularity_qmax
  return(list(modularity = modularity_q, avg_module_size = avg_module_size, max_modularity = modularity_qmax, relative_modularity = qrel))
}

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

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c('do_spr', 'get_u_from_s', 'prepseeds', 'calculate_modularity_largest_component'))
t1s.m.1 <- list()
t1s.r.1 <- list()
res1.r.2<- data.frame()
optim_result <- list()
resultsdfs <- list()

# Sample networks from larger subset so there is only one of repeated networks
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

##---------------------------------------------------------------------------------------------------------------------------------------------------
# Run simulations
foreach(network_name = names(sampled_networks)[1:min(length(sampled_networks),5)]) %do% {
  network <- sampled_networks[[network_name]]
  random_network <- erdos.renyi.game(igraph::vcount(network), igraph::ecount(network), type = "gnm")


  r_adj= 0.003
 # t_R0=100
 # rs=1.5/((mean(degree(network)^2)-mean(degree(network)))/mean(degree(network)))
			#calculate gamma (r_adj)
			#r_adj<-1-(1-rs)^(1/t_R0)
  
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
  avg_module_size <- modularity_and_avg_module_size$avg_module_size
  qrel<- modularity_and_avg_module_size$relative_modularity
  global_clustering_coefficient <- igraph::transitivity(network, type = "global")
  network_size <- igraph::vcount(network)
  degrees <- igraph::degree(network)
  mean_degree <- mean(degrees)
  sd_degree <- sd(degrees)
  avg_path_length <- igraph::mean_distance(network)

#collecting data from the simulations
  t1s.r.1_str <- toString(t1s.r.1)
  t1s.r.2_str <- toString(t1s.r.2)
  t1s.m.1_str_infected<- toString(t1s.m.1$numinfected)
  t1s.m.1_str_informed<- toString(t1s.m.1$numinformed)
  resultdf <- data.frame(
    infected.mean = mean(t1s.m.1$numinfected, na.rm = TRUE),
    infected.sd = sd(t1s.m.1$numinfected, na.rm = TRUE),
    prop_infected=  mean(t1s.m.1$numinfected, na.rm = TRUE) / vcount(network),
    conformist.informed.mean = mean(t1s.m.1$numinformed, na.rm = TRUE),
    conformist.informed.sd = sd(t1s.m.1$numinformed, na.rm = TRUE),
    conformist.prop_informed=  mean(t1s.m.1$numinformed, na.rm = TRUE) / vcount(network),
    t1s.m.1_infected=  t1s.m.1_str_infected,
    t1s.m.1_informed=  t1s.m.1_str_informed,
    bhs = bhatt_coef(as.vector(t1s.m.1$numinfected), as.vector(t1s.m.1$numinformed)),
    modularity = modularity_score,
    qrel= qrel,
    avg_module_size = avg_module_size,
    clustering = global_clustering_coefficient,
    mean_degree = mean_degree,
    sd_degree = sd_degree,
    avg_path_length= avg_path_length,
    network_size = network_size,
    chosenu = min_u,
    r_infected_mean = mean(t1s.r.1, na.rm = TRUE),
    r_infected_sd = sd(t1s.r.1, na.rm = TRUE),
    r_informed_mean = mean(t1s.r.2, na.rm = TRUE),
    r_informed_sd = sd(t1s.r.2, na.rm = TRUE),
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

