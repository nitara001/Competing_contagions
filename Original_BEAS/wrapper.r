source("/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/Original_BEAS/generate_networks.R")
source("/Users/nitarawijayatilake/Documents/GitHub/Competing_contagions/Original_BEAS/spreadfunctions_2.R")

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

		source("spreadfunctions_2.R")
		source("wrapper.R")
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
netsize= 200
thresh=0.75
thresh_points=c(0.5,0.75,0.9)

#define parameters
densities=seq(0.02,0.1,0.01) #density
R0s<-c(1.1,1.25,1.5,2,3) # R0s
allngroups=c(5,8,10,20,25,40) # n groups
allqrel=c(0.4,0.6,0.8) # modularities.

#Set measuring point for similarity in epidemics
thresh_point=thresh*netsize

#conformist learning steepness
minlearn=0.001;thrsteep=10

#for each density
foreach(density=densities)%do%{
	#20 reps
	foreach(rep=c(1:nreps))%do%{
		#generate a random network of certain density
		start_net=generate_random(N=netsize,e_prob=density)	
		
		#for each R0	
		foreach(R0=R0s)%do%{
			
			#if output file already exists we can skip this			if(file.exists(file.path(outputdir,paste(density,rep,R0,".csv",sep="_"),fsep="\\"))){
				#cat(paste("skip",density,rep,R0,"\n"))

				return()
			}
			starttime=Sys.time()
			#calculate r_adj for diseases
			#set timsteps 
			t_R0=100
			
			
			#Caclculate r
			rs=1.5/((mean(degree(start_net)^2)-mean(degree(start_net)))/mean(degree(start_net)))
			#calculate gamma (r_adj)
			r_adj<-1-(1-rs)^(1/t_R0)
			
			#generate distribution of infection times for this network, at this r0
			t1s.r.1=foreach(i=1:50,.combine=c,.packages="igraph")%do%{
				res1.r.1<-do_spr(net=start_net,type="infected",n_seeds=nseeds,loc_seeds=locseeds,s=r_adj,thresh=0.75,tmax=tmax,returnnets=F,verbose=F)
  				sort(res1.r.1$wheninfected,na.last=TRUE)[thresh_point]
			}
			#generate potential u
			allus=foreach(i=1:5)%do%{optimise(get_u_from_s,c(0,1),rnet=start_net,infres=t1s.r.1,thresh=thresh,n_seeds=nseeds,loc_seed=locseeds,inform.type="conformist",min_learn=minlearn,thr_steep=thrsteep)}
			
			#choose best u		
			allus=do.call(rbind,lapply(allus,function(x){data.frame(u=x$minimum,bhs=x$objective)}))
			curr_u=allus$u[allus$bhs==min(allus$bhs)][1]
			
			#do information spread through random network using this u
			t1s.r.2=foreach(i=1:50,.combine=c,.packages="igraph")%do%{
				res1.r.2<-do_spr(net=start_net,type="informed",inform.type="conformist",min_learn=minlearn,thr_steep=thrsteep,n_seeds=nseeds,loc_seeds=locseeds,u=curr_u,thresh=0.75,tmax=tmax,returnnets=F,verbose=F)
  				sort(res1.r.2$wheninformed,na.last=TRUE)[thresh_point]
			}
			
			#all combos of group n and modularity
			allpars=expand.grid(ngroups=allngroups,Qrel=allqrel)
			
			#foreach combo of groupsize
			allresultdf=foreach(j = 1:nrow(allpars),.combine=rbind)%do%{
				cat(paste(density,rep,R0,allpars$ngroups[j],allpars$Qrel[j],"\n"))
				
				#Rewire network
				m1=generate_modular2(in_net=start_net,N=netsize,n_groups=allpars$ngroups[j],max_dens=0.1,Qrel=allpars$Qrel[j],style="bm",cutpoint=20000,updates=F,plot=FALSE,returnqrel=T)
				m=m1$net
				returnedqrel=m1$qrel#store target modularity and actual modularity 
				#give vertices names
				V(m)$names=1:netsize
				
				#simulate spread of disease and information 50 times
				t1s.m.1=foreach(i=1:50,.combine=rbind,.packages="igraph")%do%{
					res1.m.1=do_spr(net=m,type="both",n_seeds=nseeds,inform.type="conformist",min_learn=minlearn,thr_steep=thrsteep,loc_seeds=locseeds,s=r_adj,u=curr_u,tmax=tmax,returnnets=F,verbose=F)
					#get when spreads reached threshold points
					do.call(rbind,lapply(thresh_points,function(x){
						data.frame(thresh_point=x,wheninfected=sort(res1.m.1$wheninfected,na.last=TRUE)[x*nrow(res1.m.1)],
						wheninformed=sort(res1.m.1$wheninformed,na.last=TRUE)[x*nrow(res1.m.1)])
					}))
				}
				#store network dfs

				resultdf=do.call(rbind,lapply(thresh_points,function(x){
					data.frame(thresh_point=x,
					infected.mean=mean(t1s.m.1$wheninfected[t1s.m.1$thresh_point==x],na.rm=T),
					infected.sd=sd(t1s.m.1$wheninfected[t1s.m.1$thresh_point==x],na.rm=T),
					informed.mean=mean(t1s.m.1$wheninformed[t1s.m.1$thresh_point==x],na.rm=T),
					informed.sd=sd(t1s.m.1$wheninformed[t1s.m.1$thresh_point==x],na.rm=T),
					bhs=as.vector(bhatt_coef(t1s.m.1$wheninfected[t1s.m.1$thresh_point==x],t1s.m.1$wheninformed[t1s.m.1$thresh_point==x]))
				)}))
			
				data.frame(density=density,rep=rep,ngroups=allpars$ngroups[j],
					TargetQrel=allpars$Qrel[j],ReturnQrel=returnedqrel,
					R0=R0,s=r_adj,umean=mean(allus$u),usd=sd(allus$u),chosenu=curr_u,
					r_infected_mean=mean(t1s.r.1,na.rm=T),r_infected_sd=sd(t1s.r.1,na.rm=T),
					r_informed_mean=mean(t1s.r.2,na.rm=T),r_informed_sd=sd(t1s.r.2,na.rm=T),
					r_bhs=as.vector(bhatt_coef(t1s.r.1,t1s.r.2)),
					resultdf,starttime=starttime,endtime=Sys.time())
			
			
			}
			cat(paste("export",density,rep,R0,"\n"))

			write.csv(allresultdf,file.path(outputdir,paste(density,rep,R0,".csv",sep="_"),fsep="\\"),row.names=F)

			
		}	
	}			
