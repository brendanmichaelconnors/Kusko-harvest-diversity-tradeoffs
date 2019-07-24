########################################################################################
# functions.R
#
# Functions for closed loop simulations and assocaited analyses 
# July 5, 2019
# Author: B. Connors (DFO) and B. Staton (CRITFC)
#        
########################################################################################


	
#------------------------------------------------------------------------------#
# Single-stock simulation function with alternative structural forms
#------------------------------------------------------------------------------#
# ny <- the number of years
# phi <- the expected correlation through time
# rec_sd <- process error in SD units
# mat <- stock-specific maturation schedules
# alpha <- sub-stock productivity (not in log space)
# beta <- sub-stock density depedence 
# U <- exploitation rate
# OU <- outcome uncertainty (CV)
# Spw <- spawner abundance in year 1:4
# Rec <- recruitment from years 1:4
# SR_rel <- structural form of the SR relationship ("Ricker" or "Beverton-Holt")
# BH.alpha.CV <- magnitude (amplitude) of environmental forcing on alpha if SR_rel = "Beverton-Holt"
# period <- period of enviro forcing cycle if SR_rel = "Beverton-Holt"
# alpha.scalar <- scalar for BH alpha

sr_sim = function(ny,phi,rec_sd, mat,alpha,beta,U,OU,Spw,Rec,SR_rel,BH.alpha.CV,period,alpha.scalar){
  ns <- 1
  if(SR_rel == "Beverton-Holt"){alpha <- alpha* alpha.scalar}
  
	# create vectors of time varying alpha
	if (SR_rel == "Beverton-Holt"){ 
		beta.tim <- (alpha/beta)*exp(-1)
		alpha.time <- matrix(NA,ny,length(alpha))
		for (t in 1:ny){
			alpha.time[t,] <- sin(2*pi*(t/period))*((alpha + (alpha * BH.alpha.CV)) - alpha) + alpha
		}
	}
	
	#Create recruitment deviations  
	epi = rnorm(ny, 0, rec_sd)

	#Build time series of Spawners (S), abundance of returning spawners pre-harvest
	# (N), and the component of the residual that is correlated throught time (v)
	R = t(matrix(0,ns,ny));  R[1:8,] <- Rec
	S = R ; S[8,] <- Spw
	v = R; v[,]=0
	N = array(0,dim=c(ny,4,ns))
	Ntot = R; Ntot[,]=0
	H = Ntot; S = Ntot
	predR = Ntot
	
	v[7,] <- Rec[1]
	
	N[4:7,1,]=R[4:7-(3),] * mat[1]
	N[5:7,2,]=R[5:7-(4),] * mat[2]
	N[6:7,3,]=R[6:7-(5),] * mat[3]
	N[7,4,]=R[7-(6),] * mat[4]
		
	# Loop through years of simulation	
	for(i in (7+1):ny){
		N[i,1,] = R[i-(4),] * mat[1]
		N[i,2,] = R[i-(5),] * mat[2]
		N[i,3,] = R[i-(6),] * mat[3]
		N[i,4,] = R[i-(7),] * mat[4]
		Ntot[i,] = sum(N[i,,])

		# apply harvest
		outcome_error <- (1+rnorm(1,0,OU))
		H[i,] =  U*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
		S_exp = Ntot[i,]-H[i,]
		S_exp[S_exp<0] = 0
		S_exp[S_exp<50] = 0
		S[i,] = S_exp

		# predict recruitment
		if (SR_rel == "Ricker"){
			R[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i])
			predR[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,])
			v[i,] = log(R[i,])-log(predR[i,])
			v[v[,]=='NaN'] <- 0
		      }
	
		if (SR_rel == "Beverton-Holt"){
			R[i,] = alpha.time[i,]*S[i,]/(1+(alpha.time[i,]/beta.tim[])*S[i,])*exp(phi*v[i-1,]+epi[i])
			predR[i,] = alpha.time[i,]*S[i,]/(1+(alpha.time[i,]/beta.tim[])*S[i,])
			v[i,] = log(R[i,])-log(predR[i,])
			v[v[,]=='NaN'] <- 0
		     }
	  }
 
	S[S[,]=='NaN'] <- 0
	Ntot[Ntot[,]=='NaN'] <- 0

	list(S=S[,],R=R[,], N=Ntot[,],H=H[,])
	
}

