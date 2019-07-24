########################################################################################
# close_loop_sims.R
#
# Closed-loop simulation of alternative harvest policies for Kuskokwim Chinook
# July 5, 2019
# Author: B. Connors (DFO)
#        
########################################################################################

# --- Set common conditions for simulations  --------------------------------------------
num.sims = 500 # number of Monte Carlo trials
ny = 50 # number of years in forward simulation
pm.yr <- ny-20
for.error <- 0.27 
OU <- 0.1 

# --- Create array to store outcomes ----------------------------------------------------
harvest_goal <- seq(1000,100000,length.out=40)
egfloor <- seq(1,100000,length.out=40)
sim.outcomes <- array(NA,dim=c(length(egfloor),9, length(harvest_goal),num.sims))
sim.outcomes.spw.time <- array(NA,dim=c(ny,13,length(egfloor),length(harvest_goal),num.sims))

# --- Stationary Ricker SR dynamics ----------------------------------------------------

	# set structural form of SR relationship
	SR_rel <-  "Ricker" 
	dir.SR <- "F"
	SR_devs <- array(1,dim=c(ny,2,13))
	
	# run simulations
	ptm <- proc.time()
	for (w in 1:length(harvest_goal)){
	  for (k in 1:length(egfloor)){
	    for (l in 1: num.sims){
  			draw <- sample(10000,1)
  			alpha <- process.iteration(samps[draw,])$alpha
  			if(SR_rel == "Beverton-Holt"){alpha <- alpha* 0.7424242}
  			beta <- process.iteration(samps[draw,])$beta
  			Ro <- log(alpha)/beta
  			vcov.matrix <- process.iteration(samps[draw,])$Sigma_R
  			mat <- process.iteration(samps[draw,])$pis
  			Rec <- process.iteration(samps[draw,])$R
  			Spw <- process.iteration(samps[draw,])$S
  			lst.resid <- process.iteration(samps[draw,])$last_resid
  			phi <- process.iteration(samps[draw,])$phi
  			sub <- ifelse(harvest_goal[w]<45000,harvest_goal[w],45000)
  			com <- ifelse(harvest_goal[w]<45000,0,harvest_goal[w]-45000)
  			expan <- 1/(rnorm(1,0.56,0.05))
	
				out <- process(ny,Ro,vcov.matrix,phi,mat,alpha,beta,sub,com,egfloor[k],pm.yr,
							   for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV,period,dir.SR,SR_devs,expan)
				sim.outcomes[k,,w,l] <- out$PMs
				sim.outcomes.spw.time[,,k,w,l] <- out$S
			}
		}	
	}	
(proc.time() - ptm)/60

saveRDS(sim.outcomes,"outputs/base_sims.ricker")  
saveRDS(sim.outcomes.spw.time,"outputs/base_sims_projections.ricker")    

# --- Time-varying Ricker SR dynamics ----------------------------------------------------

	# set structural form of SR relationship
	SR_rel <-  "Ricker" 
	dir.SR <- "T"
	SR_devs  <- SR_para_devs(apply(samps[,1:13],c(2),quantile,probs=c(0.5),na.rm=T),
							 apply(samps[,14:26],c(2),quantile,probs=c(0.5),na.rm=T),
							 5,35,ny)
	
	# run simulations
	ptm <- proc.time()
	for (w in 1:length(harvest_goal)){
	  for (k in 1:length(egfloor)){
	    for (l in 1: num.sims){
  			draw <- sample(10000,1)
  			alpha <- process.iteration(samps[draw,])$alpha
  			if(SR_rel == "Beverton-Holt"){alpha <- alpha* 0.7424242}
  			Ro <- log(alpha)/beta
  			beta <- process.iteration(samps[draw,])$beta
  			vcov.matrix <- process.iteration(samps[draw,])$Sigma_R
  			mat <- process.iteration(samps[draw,])$pis
  			Rec <- process.iteration(samps[draw,])$R
  			Spw <- process.iteration(samps[draw,])$S
  			lst.resid <- process.iteration(samps[draw,])$last_resid
  			phi <- process.iteration(samps[draw,])$phi
  			sub <- ifelse(harvest_goal[w]<45000,harvest_goal[w],45000)
  			com <- ifelse(harvest_goal[w]<45000,0,harvest_goal[w]-45000)
  			expan <- 1/(rnorm(1,0.56,0.05))
	
				out <- process(ny,Ro,vcov.matrix,phi,mat,alpha,beta,sub,com,egfloor[k],pm.yr,
							   for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV,period,dir.SR,SR_devs,expan)
				sim.outcomes[k,,w,l] <- out$PMs
				sim.outcomes.spw.time[,,k,w,l] <- out$S
			}
		}	
	}	
(proc.time() - ptm)/60

saveRDS(sim.outcomes,"outputs/base_sims.rickerTV")  
saveRDS(sim.outcomes.spw.time,"outputs/base_sims_projections.rickerTV")  

# --- Time-varying Beverton-Holt SR dynamics ----------------------------------------------------

	# set structural form of SR relationship
	SR_rel <-  "Beverton-Holt" 
	BH.alpha.CV <- 0.4
	period <- 12
	dir.SR <- "F"
	SR_devs <- array(1,dim=c(ny,2,13))
	
	# run simulations!
	ptm <- proc.time()
	for (w in 1:length(harvest_goal)){
	  for (k in 1:length(egfloor)){
	    for (l in 1: num.sims){
  			draw <- sample(10000,1)
  			alpha <- process.iteration(samps[draw,])$alpha
  			if(SR_rel == "Beverton-Holt"){alpha <- alpha*0.7424242}
  			Ro <- log(alpha)/beta
  			beta <- process.iteration(samps[draw,])$beta
  			vcov.matrix <- process.iteration(samps[draw,])$Sigma_R
  			mat <- process.iteration(samps[draw,])$pis
  			Rec <- process.iteration(samps[draw,])$R
  			Spw <- process.iteration(samps[draw,])$S
  			lst.resid <- process.iteration(samps[draw,])$last_resid
  			phi <- process.iteration(samps[draw,])$phi
  			sub <- ifelse(harvest_goal[w]<45000,harvest_goal[w],45000)
  			com <- ifelse(harvest_goal[w]<45000,0,harvest_goal[w]-45000)
  			expan <- 1/(rnorm(1,0.56,0.05))
			
				out <- process(ny,Ro,vcov.matrix,phi,mat,alpha,beta,sub,com,egfloor[k],pm.yr,
							   for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV,period,dir.SR,SR_devs,expan)
				sim.outcomes[k,,w,l] <- out$PMs
				sim.outcomes.spw.time[,,k,w,l] <- out$S
			}
		}	
	}	
(proc.time() - ptm)/60

saveRDS(sim.outcomes,"outputs/base_sims.bhTV")  
saveRDS(sim.outcomes.spw.time,"outputs/base_sims_projections.bhTV")  