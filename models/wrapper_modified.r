source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\Original_BEAS\\spreadfunctions_2.R")
source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\Original_BEAS\\generate_networks.R")

#"C:/Users/s2607536/AppData/Local/Programs/R/R-43~1.2/bin/R.exe"
##----------------------------------------------

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
get_u_from_s <- function(u_tmp, rnet, infres) {
  thresh <- 1.0
  n_seeds <- 1
  loc_seed <- "R"
  inform.type <-"conformist"
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

calculate_modularity_and_avg_module_size <- function(network) {
    community <- cluster_louvain(network, weights= E(network)$weight)
    modularity<- modularity(community)
    avg_module_size <- mean(sizes(community))
    return(list(modularity = modularity, avg_module_size = avg_module_size))
}

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


# Sample networks from larger subset so there is only one of repeated networks
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c('do_spr', 'get_u_from_s', 'prepseeds', 'calculate_modularity_and_avg_module_size'))
t1s.m.1 <- list()
t1s.r.1 <- list()
optim_result <- list()
resultsdfs <- list()

# Sample networks from larger subset so there is only one of repeated networks
unique_network_names <- unique(gsub("_\\d+$", "", names(extracted_networks)))
sampled_networks <- list()
for (network_name in unique_network_names) {
  matching_networks <- grep(paste0("^", network_name), names(extracted_networks), value = TRUE)
  sampled_network_name <- sample(matching_networks, 1)
  sampled_networks[[sampled_network_name]] <- extracted_networks[[sampled_network_name]]
}

# Run simulations
foreach(network_name = names(sampled_networks)[1:min(length(sampled_networks), 57)]) %do% {
  network <- sampled_networks[[network_name]]
  random_network <- erdos.renyi.game(igraph::vcount(network), igraph::ecount(network), type = "gnm")
  
r_adj= 0.005
  
  t1s.r.1 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.1 <- do_spr(net = random_network, type = "infected", n_seeds = 1, loc_seeds = 'R', s = r_adj, tmax = 100, returnnets = FALSE, verbose = FALSE, inform.type = "conformist")
    sum(res1.r.1$infected)
  }

opt_result <- stats::optimize(get_u_from_s, interval = c(0, 1), rnet = random_network, infres = t1s.r.1)

# Extract the minimum u value corresponding to the minimum bhs value
min_u <- opt_result$minimum

  t1s.r.2 <- foreach(i = 1:50, .combine = c, .packages = "igraph") %do% {
    res1.r.2 <- do_spr(net = random_network, type = "informed", inform.type = "conformist", min_learn = 0.0001, n_seeds = 1, loc_seeds = 'R', u = min_u, tmax = 100, returnnets = FALSE, verbose = TRUE)
    sum(res1.r.2$informed)
  }
  
  t1s.m.1 <- foreach(i = 1:50, .combine = rbind, .packages = "igraph") %dopar% {
    res1.m.1 <- do_spr(net = network, type = "both", n_seeds = 1, inform.type = "conformist", min_learn = 0.0001, loc_seeds = 'R', s = r_adj, u = min_u, tmax = 100, returnnets = FALSE, verbose = TRUE)
    data.frame(network = network_name, numinfected = sum(res1.m.1$infected), numinformed = sum(res1.m.1$informed))
  }

  modularity_and_avg_module_size <- calculate_modularity_and_avg_module_size(network)
  modularity_score <- modularity_and_avg_module_size$modularity
  avg_module_size <- modularity_and_avg_module_size$avg_module_size
  
  resultdf <- data.frame(
    infected.mean = mean(t1s.m.1$numinfected, na.rm = TRUE),
    infected.sd = sd(t1s.m.1$numinfected, na.rm = TRUE),
    Prop_infected=  mean(t1s.m.1$numinfected, na.rm = TRUE) / vcount(network),
    informed.mean = mean(t1s.m.1$numinformed, na.rm = TRUE),
    informed.sd = sd(t1s.m.1$numinformed, na.rm = TRUE),
    Prop_informed=  mean(t1s.m.1$numinformed, na.rm = TRUE) / vcount(network),
    bhs = bhatt_coef(as.vector(t1s.m.1$numinfected), as.vector(t1s.m.1$numinformed)),
    modularity = modularity_score,
    avg_module_size = avg_module_size,
    chosenu = min_u,
    r_infected_mean = mean(t1s.r.1, na.rm = TRUE),
    r_infected_sd = sd(t1s.r.1, na.rm = TRUE),
    r_informed_mean = mean(t1s.r.2, na.rm = TRUE),
    r_informed_sd = sd(t1s.r.2, na.rm = TRUE),
    r_bhs = bhatt_coef(as.vector(t1s.r.1), as.vector(t1s.r.2))
  )
  
  resultsdfs[[network_name]] <- resultdf
}

final_resultdf <- do.call(rbind, resultsdfs)

# Clean up parallel backend
stopCluster(cl)

# Save results
write.csv(final_resultdf, "C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\final_results.csv")