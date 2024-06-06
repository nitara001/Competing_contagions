source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\spreadfunctions_2.R")
source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\generate_networks.R")

#"C:/Users/s2607536/AppData/Local/Programs/R/R-43~1.2/bin/R.exe"
##----------------------------------------------
tmax=200
nseeds=1
locseeds="R"
nreps=20

## apply functions
library(igraph)
library("nloptr")
library(foreach)
library(doParaellel)
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
  tmax <- 500
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
unique_network_names <- unique(gsub("_\\d+$", "", names(extracted_networks)))
sampled_networks <- list()
for (network_name in unique_network_names) {
  matching_networks <- grep(paste0("^", network_name), names(extracted_networks), value = TRUE)
  sampled_network_name <- sample(matching_networks, 1)
  sampled_networks[[sampled_network_name]] <- extracted_networks[[sampled_network_name]]
}

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set up parallel processing
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c('do_spr', 'get_u_from_s', 'prepseeds', 'calculate_modularity_and_avg_module_size'))
t1s.m.1 <- list()
t1s.r.1 <- list()
opt_results <- list()
resultsdfs <- list()
curr_u_list <- list() 
# Loop over extracted networks and corresponding t1s.r.1 values
foreach(network_name = names(sampled_networks)[1:min(length(sampled_networks), 3)], .combine = 'c') %do% {
  network <- sampled_networks[[network_name]]

  random_network<- erdos.renyi.game(igraph::vcount(network), igraph::ecount(network), type= "gnm") #make random network for each asnr network
  t_R0 <- 100
  # Calculate r
  rs <- 1 / ((mean(degree(random_network)^2) - mean(degree(random_network))) / mean(degree(random_network)))
  # Calculate (r_adj)
  r_adj <- 1 - (1 - rs)^(1 / t_R0)
  # Generate distribution of infection times for this network, at this r0
  t1s.r.1<- foreach(i = 1:50, .combine = c, .packages = "igraph") %dopar% {

    res1.r.1 <- do_spr(net = random_network, type = "infected", n_seeds = 1, loc_seeds = 'R', s = r_adj, tmax = 500, returnnets = F, verbose = T, inform.type = "conformist")
    sort(sum(res1.r.1$infected), na.last = TRUE)
  }
  
  # Perform optimization for each network and corresponding t1s.r.1 values
opt_results=stats::optimize(get_u_from_s, interval = c(0, 1), rnet = network, infres = t1s.r.1)
  
  all_opt_results[[network_name]] <- opt_results
    
    # get curr_u for the current network and store it
    curr_u <- opt_results$minimum
    curr_u_list[[network_name]] <- curr_u

#do information spread at that curr_u
			t1s.r.2=foreach(i=1:50,.combine=c,.packages="igraph")%dopar%{
				res1.r.2<-do_spr(net=random_network,type="informed",inform.type="conformist",min_learn=0.001 ,n_seeds=1,loc_seeds='R',u=curr_u,thresh=1.0 ,tmax=500,returnnets=F,verbose=T)
  				sort(sum(res1.r.2$informed,na.last=TRUE))
			}
      
    #simulate spread of disease and information 50 times in the real netwrok
				t1s.m.1=foreach(i=1:50,.combine=rbind,.packages="igraph")%dopar%{
					res1.m.1=do_spr(net=network,type="both",n_seeds=1,inform.type="conformist",min_learn=0.001,loc_seeds='R',s=r_adj,u=curr_u,tmax=500,returnnets=F,verbose=T)

						data.frame(network=network_name,numinfected=sort(sum(res1.m.1$infected,na.last=TRUE)),
						sort(sum(res1.m.1$informed,na.last=TRUE)))
					
   # Get some network characteristics
  modularity_and_avg_module_size <- calculate_modularity_and_avg_module_size(network)
  modularity_score <- modularity_and_avg_module_size$modularity
  avg_module_size <- modularity_and_avg_module_size$avg_module_size
  

  resultdf <- data.frame(
    infected.mean = mean(sum(t1s.m.1$infected, na.rm = TRUE)),
    infected.sd = sd(sum(t1s.m.1$infected, na.rm = TRUE)),
    informed.mean = mean(sum(t1s.m.1$informed, na.rm = TRUE)),
    informed.sd = sd(sum(t1s.m.1$informed, na.rm = TRUE)),
    bhs = as.vector(bhatt_coef(sum(t1s.m.1$infected), sum(t1s.m.1$informed))),
    modularity = modularity_score,
    avg_module_size = avg_module_size,
    umean = mean(opt_results$minimum),
    usd = sd(opt_results$minimum),
    chosenu = curr_u,
    r_infected_mean = mean(t1s.r.1, na.rm = TRUE),
    r_infected_sd = sd(t1s.r.1, na.rm = TRUE),
    r_informed_mean = mean(t1s.r.2, na.rm = TRUE),
    r_informed_sd = sd(t1s.r.2, na.rm = TRUE),
    r_bhs = as.vector(bhatt_coef(t1s.r.1, t1s.r.2))
  )

  resultsdfs[[network_name]] <- resultdf}}


# Combine resultdfs into a single data frame
final_resultdf <- do.call(rbind, resultsdfs)
stopcluster(cl)