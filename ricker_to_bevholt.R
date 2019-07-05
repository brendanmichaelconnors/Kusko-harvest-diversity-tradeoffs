#------------------------------------------------------------------------------#
# This code seeks to determine the magnitude of increase in alpha in the BH SR relationship required to match equilibrium outcomes in the Kusko that are equivalent to a Ricker SR formulation 
#------------------------------------------------------------------------------#
#### Base Ricker
alpha <- colMedians(samps)[1:13]
beta <- colMedians(samps)[14:26]
Ro <- log(alpha)/beta
rick_s_eq <- sum(Ro)

com <- 0
for.error <- 0.27 # from ben (0.27)
OU <- 0.1 #(0.2 collie and peterman high end YK chum)

# set conditions for simulation
ny = 50
pm.yr <- 30

SR_rel <-  "Beverton-Holt" #"Beverton-Holt"
BH.alpha.CV <- 0.4
period <- 12
dir.SR <- "F"
num.sims<-100
alpha.multiplier <- seq(0.5,1.5,length.out=100)
outcomes.BH <- array(NA,dim=c(length(alpha.multiplier),9, num.sims))
harvest_goal = 0
egfloor =0
# run simulations!
ptm <- proc.time()
for (w in 1:length(alpha.multiplier)){
		
	for (l in 1: num.sims){
		draw <- sample(10000,1)
		alpha <- process.iteration(samps[draw,])$alpha * alpha.multiplier[w]
		beta <- process.iteration(samps[draw,])$beta
		vcov.matrix <- process.iteration(samps[draw,])$Sigma_R
		mat <- process.iteration(samps[draw,])$pi
		Rec <- process.iteration(samps[draw,])$R
		Spw <- process.iteration(samps[draw,])$S
		lst.resid <- process.iteration(samps[draw,])$last_resid
		phi <- process.iteration(samps[draw,])$phi
		
		sub <- ifelse(harvest_goal[1]<45000,harvest_goal[1],45000)
		com <- ifelse(harvest_goal[1]<45000,0,harvest_goal[1]-45000)
		
		out <- process(ny,Ro,vcov.matrix,phi,mat,alpha,beta,sub,com,egfloor[1],pm.yr,for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV,period,dir.SR,SR_devs)
		outcomes.BH[w,,l] <- out$PMs
		}
	}	
	
bh_s_seq <- apply(outcomes.BH,c(1,2),quantile,probs=c(0.5),na.rm=T)[,1]

alpha.multiplier[which(abs(bh_s_seq-rick_s_eq)==min(abs(bh_s_seq-rick_s_eq)))]

#### at end of day you need to multiply alpha by 0.742 to generate similar equilibrium outcomes between BH and Ricker formulations