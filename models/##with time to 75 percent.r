
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c('do_spr', 'get_u_from_s', 'prepseeds', 'calculate_modularity_largest_component'))
t1s.m.1 <- list()
t1s.r.1 <- list()
res1.r.2<- data.frame()
optim_result <- list()
resultsdfs <- list()


##---------------------------------------------------------------------------------------------------------------------------------------------------
# Run simulations
foreach(network_name = names(sampled_networks)[1:min(length(sampled_networks))]) %do% {
  network <- sampled_networks[[network_name]]
  random_network <- erdos.renyi.game(igraph::vcount(network), igraph::ecount(network), type = "gnm")


  #r_adj= 0.003
  t_R0=100
  rs=1.5/((mean(degree(network)^2)-mean(degree(network)))/mean(degree(network)))
			#calculate gamma (r_adj)
			r_adj<-1-(1-rs)^(1/t_R0)
  
  t1s.r.1 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.1 <- do_spr(net = random_network, type = "infected", n_seeds = 1, loc_seeds = 'R', s = r_adj, tmax = 200, returnnets = FALSE, verbose = FALSE, inform.type = "conformist")
    sum(res1.r.1$wheninfected)}

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

