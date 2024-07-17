#Function to simulate transmission of information and disease in networks
#By Evans and Silk (2021)
#Parameters:
#net: network in which to simulate spread. Can be a network previously output by this function to continue spread.
#type: type of spread to simulate
#n_seeds: Number of starting seeds
#loc_seeds: Start seeds at random (R) or in a random cluster (CL) if n_seeds>1. Can also choose individuals of high betweenness (B) or using a vector of vertex number
#s: the per timestep disease transmission probability per connection
#u: scaling parameter for information transfer
#tmax: stop after a certain number of timesteps
#thresh: stop after spreading process reaches this number of nodes
#infect.limit: allow only this number of nodes to be infected in a timestep
#inform.limit: allow only this number of nodes to be informed in a timestep.
#inform.type: type of learning rule to use
#min_learn: parameter controlling baseline of conformist learning rule
#thr_steep: parameter controlling steepness of conformist learning rule
#weighted: take association strengths into account
#recovery: nodes can recover from infection
#recoverytime: how long before nodes can start to recover from infection
#recoverprob: probability of a node recovering from infection after a certain amount of time
#death: infection can result in death
#deathtime: how long does a node have to be infected before they might die
#deathprob: probability of death after being infected for a certain amount of time
#immune: nodes can become immune to infection
#immunetime: minimum time immunity lasts
#immuneprob: probability of becoming immune after recovering
#immuneendprob: probability of losing immunity after a certain amount of time
#returnnets: return a list of network objects if true, otherwise return a dataframe of when individuals became infected/informed
#verbose: show progress bar, give messages.
library(igraph)
do_spr<-function(net,type=c("both","infected","informed"),n_seeds=1,loc_seeds=c("R","CL"),s,u,tmax,thresh=NA,infect.limit=NA,inform.limit=NA,inform.type=c("proportional","conformist"),min_learn=0.001,thr_steep=10,weighted=F,recovery=F,recoverytime=1,recoverprob=1,death=F,deathtime=1,deathprob=0.5,immune=F,immunetime=1,immuneprob=0.5,immuneendprob=0.5,returnnets=T,verbose=T){
  N=length(V(net))
  
  if(!is.na(thresh)){
    thresh_point=thresh*N
  }else{
    thresh_point=N
  }
  
  
  
  #Set up seeds
  
  if(type%in%c("informed","both")&is.null(V(net)$informed)){
    net=prepseeds(net,"informed",n_seeds,loc_seeds)
  }
  if(type%in%c("infected","both")&is.null(V(net)$infected)){
    net=prepseeds(net,"infected",n_seeds,loc_seeds)
  }
  
  if(is.null(V(net)$alive)){
    V(net)$alive=T
  }
  if(is.null(V(net)$whendead)){	
    V(net)$whendead=NA
  }
  if(is.null(V(net)$immune)){
    V(net)$immune=F
  }
  
  if(recovery|immune&(is.null(V(net)$recovered)|!is.null(V(net)$whenrecovered))){
    V(net)$recovered=F
    V(net)$whenrecovered=NA	
  }	
  
  #store when an individual became infected
  
  if(is.null(V(net)$timestep)){
    V(net)$timestep=0
  }
  
  if(returnnets){
    
    nets=list(net)
  }
  
  
  if(verbose){
    pb <- txtProgressBar((V(net)$timestep[1]),(V(net)$timestep[1]+tmax),style=3)
  }
  for(t in (V(net)$timestep[1])+1:(V(net)$timestep[1]+tmax)){
    if(verbose){
      setTxtProgressBar(pb, t)
    }
    if(!is.null(V(net)$infected)){
      
      if(weighted){
        
        tmat=as_adjacency_matrix(net,attr="weight",sparse=FALSE)*as.numeric(V(net)$infected&V(net)$alive)
        
      }else{
        
        tmat=as_adjacency_matrix(net,sparse=FALSE)*as.numeric(V(net)$infected&V(net)$alive)
        
      }
      
      trisk<-colSums(tmat)
      
      tvec<-rbinom(n=N,size=1,prob=1-(1-s)^trisk)
      if(!is.na(infect.limit)){
        limitinfect=sample(which(tvec==1&!V(net)$infected&V(net)$alive&!V(net)$immune),min(c(infect.limit,sum(tvec==1&!V(net)$infected&V(net)$alive&!V(net)$immune))),replace=F)
        tvec=rep(0,length(tvec))
        tvec[limitinfect]=1
      }
      V(net)$wheninfected[tvec==1&!V(net)$infected&V(net)$alive&!V(net)$immune]=t
      V(net)$infected[tvec==1&!V(net)$infected&V(net)$alive&!V(net)$immune]=T
      
      if(recovery){
        #can recover after a certain time
        canrecover=t-V(net)$wheninfected>recoverytime&V(net)$infected&V(net)$alive
        recovers=rbinom(sum(canrecover),1,recoverprob)
        V(net)$infected[canrecover]=!recovers>0
        V(net)$recovered[canrecover]=recovers>0
        V(net)$whenrecovered[canrecover&V(net)$recovered]=t
        if(immune){
          #can become immune
          canbeimmune=V(net)$whenrecovered==t&!is.na(V(net)$whenrecovered)
          beimmune=rbinom(sum(canbeimmune),1,immuneprob)
          V(net)$immune[canbeimmune]=beimmune>0
          
          #can be immune for a certain time then lose it at a certain probability
          canendimmune=t-V(net)$whenrecovered>immunetime&V(net)$recovered&V(net)$alive&V(net)$immune
          endimmune=rbinom(sum(canendimmune),1,immuneprob)
          
          V(net)$immune[canendimmune]=!endimmune>0
        }
      }
      
      if(death){
        #can recover after a certain time
        candie=t-V(net)$wheninfected>deathtime&V(net)$infected&V(net)$alive
        deaths=rbinom(sum(candie),1,recoverprob)
        V(net)$alive[candie]=!deaths>0
        V(net)$whendead[candie&!V(net)$alive]=t
      }
    }

    if(!is.null(V(net)$informed)&inform.type=="proportional"){
      if(weighted){
        tmat=as_adjacency_matrix(net,attr="weight",sparse=FALSE)*as.numeric(V(net)$informed&V(net)$alive)
      }else{
        tmat=as_adjacency_matrix(net,sparse=FALSE)*as.numeric(V(net)$informed&V(net)$alive)
      }
	  trisk<-colSums(tmat)/(colSums(as_adjacency_matrix(net,sparse=FALSE)))
      tvec<-rbinom(n=N,size=1,prob=u*trisk)
      V(net)$wheninformed[tvec==1&!V(net)$informed&V(net)$alive]=t
      V(net)$informed[tvec==1&!V(net)$informed&V(net)$alive]=T
	V(net)$trisk=u*trisk
	  
    }else if(!is.null(V(net)$informed)&inform.type=="conformist"){
      if(weighted){
        tmat2=as_adjacency_matrix(net,attr="weight",sparse=FALSE)
	  tmat=tmat2*as.numeric(V(net)$informed&V(net)$alive)
      }else{
        tmat2=as_adjacency_matrix(net,sparse=FALSE)
	  tmat=tmat2*as.numeric(V(net)$informed&V(net)$alive)
      }
	deg=(colSums(tmat2))
	deg[deg<0]=0
    t_x<-colSums(tmat)
	t_x[deg>0]=t_x[deg>0]/deg[deg>0]
      trisk<-u*(1-2*min_learn)/(1+exp(-thr_steep*(1/u)*(t_x-0.5)))+min_learn
	trisk=trisk*(t_x>0)
      tvec<-rbinom(n=N,size=1,prob=trisk)
      V(net)$wheninformed[tvec==1&!V(net)$informed&V(net)$alive]=t
      V(net)$informed[tvec==1&!V(net)$informed&V(net)$alive]=T
	V(net)$t_x=t_x
	V(net)$trisk=trisk
    }
    if(returnnets){
      V(net)$timestep=t
      nets=c(nets,list(net))
    }
    if(!is.null(V(net)$informed)&!is.null(V(net)$infected)){
      if(sum(V(net)$informed)>=thresh_point&sum(V(net)$infected)>=thresh_point){
        break()
      }
    }else if(!is.null(V(net)$informed)){
      if(sum(V(net)$informed)>=thresh_point){
        break()
      }
    }else if(!is.null(V(net)$infected)){
      if(sum(V(net)$infected)>=thresh_point){
        break()
      }
    }
  }
  if(verbose){
    cat("\n")
  }
  if(returnnets){
    return(nets)
  }else{
    df=as_data_frame(net,"vertices")
    return(df[,colnames(df)!="timestep"])
  }
}

prepseeds=function(net,type,n_seeds=1,loc_seeds){
  N=length(V(net))
  
  
  #generic function for setting up seeds
  vertex_attr(net,type)=F
  if (is.numeric(loc_seeds)){
    vertex_attr(net,type,loc_seeds)=T
  }else if(loc_seeds=="R"){
    vertex_attr(net,type,sample(1:N,n_seeds,replace=F))=T
  }else if (loc_seeds=="B"){
    vertex_attr(net,type,order(betweenness(net),decreasing=T)[1:n_seeds])=T
  }else if (loc_seeds=="CL"){
    vertex_attr(net,type,sample(1:N,1,replace=F))=T
    if(n_seeds>1){
      connected=neighbors(net,vertex_attr(net,type))
      vertex_attr(net,type,sample(which(V(net)%in%connected),min(length(connected),n_seeds-1),replace=F))=T
    }
  }
  
  vertex_attr(net,paste0("when",type))=NA
  vertex_attr(net,paste0("when",type),which(vertex_attr(net,type)))=0
  return(net)
}
