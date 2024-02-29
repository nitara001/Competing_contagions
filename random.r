library(igraph)
# Generate the random network
random_sealions <- erdos.renyi.game(1007, 131351, type = "gnm")

# Initialize lists to store results
num_infected_list <- vector("list", length = 50)
time_to_75_percent_list <- vector("numeric", length = 50)

# Number of simulations
num_simulations <- 50

for (i in 1:num_simulations) {
  # Perform simulation with a different starting individual each time
  simulation_result <- do_spr(
    net = random_sealions,
    type = "infected",
    n_seeds = 1,
    loc_seeds = "R", # Select a random starting individual
    s = 0.0005,
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
    returnnets = FALSE)

  num_infected <- sum(simulation_result$infected)
  
  # Calculate time to 75% infection
  time_to_75_percent <- calculate_75_percent_time(simulation_result, vcount(random_sealions))
  
  # Store time to 75% infection in the list
  time_to_75_percent_list[i] <- time_to_75_percent

}
  
calculate_75_percent_time <- function(simulation_result, total_nodes) {
  infected_nodes <- simulation_result$infected
  cum_infected <- cumsum(infected_nodes)

  target_nodes <- 0.75 * total_nodes
  time_to_75_percent <- which.max(cum_infected >= target_nodes)
  
  return(time_to_75_percent)
}


###----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Step 2: Simulate Information Spread

gamma_values <- c(0.01, 0.05, 0.1)  

# Initialize list to store results
time_to_75_percent_informed_list <- list()

for (gamma in gamma_values) {
  results_for_gamma <- numeric(num_simulations)  # Initialize results for this gamma
  
  for (i in 1:num_simulations) {
    info_simulation_result <- do_spr(
      net = random_sealions,
      type = "informed",
      n_seeds = 1,
      loc_seeds = "R",  # Select a random starting individual
      s = NULL,  # Use null for information spread
      u = gamma,  # Use different values of gamma
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
      returnnets = FALSE
    )
  
    time_to_75_percent_informed <- calculate_75_percent_time_informed(info_simulation_result, vcount(random_sealions))
    results_for_gamma[i] <- time_to_75_percent_informed
  }
  
  # Store results for this gamma
  time_to_75_percent_informed_list[[as.character(gamma)]] <- results_for_gamma / gamma
}

print(time_to_75_percent_informed_list)

  # Record time to inform 75% of nodes for each value of gamma

calculate_75_percent_time_informed <- function(simulation_result, total_nodes) {
  informed_nodes <- simulation_result$informed
  cum_informed <- cumsum(informed_nodes)

  target_nodes <- 0.75 * total_nodes
  time_to_75_percent <- which.max(cum_informed >= target_nodes)
  
  return(time_to_75_percent)
}


########################--------------------------------------------------------------------------------------------------------------------------------------------------
source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\generate_networks.R")
source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\spreadfunctions_2.R")
library(igraph)
library(foreach)

random_sealions <- erdos.renyi.game(1007, 131351, type = "gnm")

#from https://stats.stackexchange.com/questions/78849/measure-for-separability
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
get_u_from_s=function(u_tmp,rnet,infres,thresh=0.75,n_seeds=1,loc_seed="R",inform.type=c("proportional","conformist"),min_learn=0.01,thr_steep=10,tmax=3500){
	thresh_point=thresh*length(V(rnet))
	#cat(paste(u_tmp,"\n"))
	t_tmps=foreach(i=1:50,.combine=c,.packages="igraph")%dopar%{

		source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\spreadfunctions_2.r")
		source("C:\\Users\\s2607536\\OneDrive - University of Edinburgh\\Code\\BEAS Code implementation\\wrapper.r")
    netsize= igraph::V(rnet)
		
    res2.r.1<-do_spr(net=rnet,type="informed",n_seeds=n_seeds,loc_seed=loc_seed,inform.type=inform.type,min_learn=min_learn,thr_steep=thr_steep,u=u_tmp,tmax=tmax,thresh=thresh,returnnets=F,verbose=F)
  			if(is.na(sort(res2.r.1$wheninformed,na.last=TRUE)[thresh_point])==FALSE){
   				sort(res2.r.1$wheninformed,na.last=TRUE)[thresh_point]
  			}
	}
	bh<-bhatt_coef(infres,t_tmps)
	cat(paste(as.vector(bh),"\n"))
	return(as.vector(bh))
}

tmax=3500
nseeds=1
locseeds="R"
nreps=20

thresh=0.75
thresh_points=c(0.5,0.75,0.9)
thresh_point=thresh*V(random_sealions)

#conformist learning steepness
minlearn=0.001;thrsteep=10
R0s=-c(2,3)
r_adj_values <- numeric(length(R0s))
starttime = Sys.time()

# Calculate r_adj for diseases
# Set timesteps 
t_R0 = 100

# Initialize vectors to store results
t1s.r.1 <- numeric(length(R0s))
r_adj_values <- numeric(length(R0s))

## this needs to be run 50x for a distribution at some point 
foreach(rep=c(1:nreps))%do%{
# Loop over R0 values
for (R0 in seq_along(R0s)) {
    # Calculate r_adj for diseases
    # Set timesteps 
    t_R0 = 100

    # Calculate r
    rs = R0 / ((mean(degree(random_sealions)^2) - mean(degree(random_sealions))) / mean(degree(random_sealions)))

    r_adj = 1 - (1 - rs)^(1/t_R0)
    r_adj_values[i] <- r_adj 


    res1.r.1 <- do_spr(
        net = random_sealions,
        type = "infected",
        n_seeds = 1,
        loc_seeds = "R", 
        s = r_adj_values, 
        u = NULL, 
        tmax = 3500,
        thresh =thresh,
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
        returnnets = FALSE
    )
    
    # Store the result for this R0 value
    tls.r.1= sort(res1.r.1$wheninfected,na.last=TRUE)[thresh_point]
}}

#tls.r.1= na.omit( tls.r.1)

##allus
library(foreach)

allus= foreach(i=1:5)%do%{optimise(get_u_from_s,c(0,1),rnet=random_sealions,infres=t1s.r.1,thresh=thresh,n_seeds=nseeds,loc_seed=locseeds,inform.type="conformist",min_learn=minlearn,thr_steep=thrsteep)}
