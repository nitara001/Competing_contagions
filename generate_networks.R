generate_random=function(N=40,e_prob=0.05,directed=FALSE,loops=FALSE,plot=F){
  library(igraph)
  A<-igraph::erdos.renyi.game(n=N,p.or.m=e_prob,type="gnp",directed=directed,loops=loops)
  if(plot==TRUE){
    plot(A,vertex.label=NA,vertex.size=8)
  }
  return(A)
}

#Creates modular networks with the same number of edges an existing networks
#There are two existing implementations, one which calculates intermediate modularities based on a block model design (style=bm), and the other that calculates intermediate modularities based on current community structure (style=cg). The former is currently better supported

generate_modular2<-function(in_net,N,n_groups=10,max_dens,Qrel=0.5,style=c("bm","cg"),cutpoint=2000,updates=FALSE,plot=FALSE,returnqrel=F){
  
  net_size=N
  
  pop_dat<-rep(1:n_groups,each=(net_size/n_groups))
  
  fqmat<-matrix(0,nr=net_size,nc=net_size)
  for(i in 1:net_size){
    for(j in 1:net_size){
      if(pop_dat[i]==pop_dat[j]){fqmat[i,j]<-1}
    }
  }
  diag(fqmat)<-0
  
  #This value must be close to or greater than the maximum graph density
  t_ed<-edge_density(graph.adjacency(fqmat,mode="undirected"))
  
  inet<-in_net
  net<-as_adjacency_matrix(inet,sparse=FALSE)
  fqnet<-graph.adjacency(fqmat,mode="undirected")
  
  if(style=="bm"){
    q<-modularity(inet,membership=pop_dat)
    qmax<-modularity(fqnet,membership=pop_dat)
  }
  
  if(style=="cg"){
    q<-modularity(inet,membership=cluster_louvain(inet)$membership)
    mem<-cluster_louvain(inet)$membership
    mem_mat<-matrix(0,nr=net_size,nc=net_size)
    for(i in 1:net_size){
      for(j in 1:net_size){
        if(mem[i]==mem[j]){mem_mat[i,j]<-1}
      }
    }
    diag(mem_mat)<-0
    qmat<-graph.adjacency(mem_mat,mode="undirected")
    qmax<-modularity(qmat,membership=cluster_louvain(qmat)$membership)
  }
  
  qrel<-q/qmax
  
  c<-1
  
  while(qrel<Qrel&c<=cutpoint){
    tmp_rows<-which(net==1,arr.ind=T)
    tmp_samp<-sample(1:nrow(tmp_rows),1)
    tmp1<-tmp_rows[tmp_samp,1]
    tmp2<-tmp_rows[tmp_samp,2]
    if(tmp1!=tmp2){
    	tval<-net[tmp1,tmp2]
      net[tmp1,tmp2]<-fqmat[tmp1,tmp2]
      net[tmp2,tmp1]<-fqmat[tmp2,tmp1]
      tval2<-net[tmp1,tmp2]
      tdiff<-tval-tval2
      if(tdiff<0){
        tmp3<-which(net==1&fqmat==0,arr.ind=TRUE)
	  if(nrow(tmp3)>0){
          tmp4<-tmp3[sample(1:nrow(tmp3),1),]
          net[tmp4[1],tmp4[2]]<-fqmat[tmp4[1],tmp4[2]]
          net[tmp4[2],tmp4[1]]<-fqmat[tmp4[2],tmp4[1]]
	  }
      }
      if(tdiff>0){
        tmp3<-which(net==0&fqmat==1,arr.ind=TRUE)
	  if(nrow(tmp3)>0){
        	tmp4<-tmp3[sample(1:nrow(tmp3),1),]
        	net[tmp4[1],tmp4[2]]<-fqmat[tmp4[1],tmp4[2]]
        	net[tmp4[2],tmp4[1]]<-fqmat[tmp4[2],tmp4[1]]
	  }
      }
    }
    
    tnet<-graph.adjacency(net,mode="undirected")
	V(tnet)$pop_dat=pop_dat
    
    #print(gsize(tnet))
    
    if(style=="bm"){
      q<-modularity(tnet,membership=pop_dat)
      qmax<-modularity(fqnet,membership=pop_dat)
    }
    
    if(style=="cg"){
      q<-modularity(tnet,membership=cluster_louvain(tnet)$membership)
      mem<-cluster_louvain(tnet)$membership
      mem_mat<-matrix(0,nr=net_size,nc=net_size)
      for(i in 1:net_size){
        for(j in 1:net_size){
          if(mem[i]==mem[j]){mem_mat[i,j]<-1}
        }
      }
      diag(mem_mat)<-0
      qmat<-graph.adjacency(mem_mat,mode="undirected")
      qmax<-modularity(qmat,membership=cluster_louvain(qmat)$membership)
    }
    
    qrel<-q/qmax
    if(updates==TRUE){if(c%%10==1){print(qrel)}}
    
    c<-c+1
    
  }
  
  if(plot==TRUE){plot(tnet,vertex.label=NA,vertex.size=6)}
  if(returnqrel){
	return(list(net=tnet,qrel=qrel))
	}else{
  		A<-tnet
  		return(A)
	}
  
} #end function