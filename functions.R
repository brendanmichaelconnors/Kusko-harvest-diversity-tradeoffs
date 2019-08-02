########################################################################################
# functions.R
#
# Functions for closed loop simulations and associated analyses 
# July 5, 2019
# Author: B. Connors (DFO) and B. Staton (CRITFC)
#        
########################################################################################

#------------------------------------------------------------------------------#
# Subsistence harvest control rule function
#------------------------------------------------------------------------------#
# sub < - subsistence requirement
# egfloor <- escapement goal
# run <- run size (i.e., pre-harvest recruitment) 
# com <- maximum commercial harvest
# for.error <- forecast error 
# OU <- outcome uncertainty
sub_hcr = function(sub, com, egfloor, run,for.error){
	run.est <- run * rlnorm(1, 0, for.error); if(is.na(run.est)==TRUE){run.est <- run};if(is.na(run)==TRUE){run <- 0}
	if(run.est - egfloor <= 0){hr = 0}
	if(run.est > egfloor){
		if((run.est - egfloor) <= (sub)){hr = (run - egfloor)/run}
		if((run.est - egfloor) > (sub)){
			if((run.est - egfloor) > (sub + com)){hr = (sub + com)/run}
			if((run.est - egfloor) <= (sub + com)){hr = (sub + (run - egfloor - sub))/run}
			}
		}
	if(hr < 0){hr=0}
	if(hr >1 ){hr=1}
	return(hr)
	}

#------------------------------------------------------------------------------#
# Status function (to estimate whether stock is overfished or predicted to go 
#  extinct at a given harvest rate, over the long-term)
#------------------------------------------------------------------------------#
# U <- harvest rate
# a <- productivity (Ricker a parameter)
# b <- density dependence (Ricker beta parameter)
SC.eq <- function(U,a,b){
	a <- log(a)
	S.eq <- max(0,(a-(-log(1-U)))/b)
	C.eq <- max(0,((a-(-log(1-U)))/b)*exp(a-b*((a-(-log(1-U)))/b))-((a-(-log(1-U)))/b))
	OF <- ifelse(U>0.5*a-0.07*a^2,1,0)
	EX <- ifelse(S.eq==0,1,0)
	return(c(S.eq,C.eq,OF,EX))
	}
	
#------------------------------------------------------------------------------#
# Multi-stock simulation function with alternative structural forms
#------------------------------------------------------------------------------#
# ny <- the number of years
# vcov.matrix <- process error variance-covarance matrix
# phi <- the expected correlation through time
# mat <- maturation schedule
# alpha <- sub-stock productivity (NOT in log space)
# beta <- sub-stock density depedence 
# sub <- subsistence requirement
# com <- maximum commercial harvest
# egfloor <- escapement goal
# pm.yr <- year of simulation that pms start to be calculated over
# for.error <- forcast error (CV)
# OU <- outcome uncertainty (CV)
# Rec <- estimated recruitments from last years of empirical data 
# Spw <- estimated spawers from last years of empirical data
# lst.resid <- estimated recruitment deviation from last year of empirical data
# SR_rel <- structural form of the SR relationship ("Ricker" or "Beverton-Holt")
# BH.alpha.CV <- magnitude (amplitude) of environmental forcing on alpha if SR_rel = "Beverton-Holt"
# period <- period of enviro forcing cycle if SR_rel = "Beverton-Holt"
# dir.SR <- flag for directional change in SR parameters ("Y" or "N")
# SR.devs <- deviations in SR parameters by time step if dir.SR == "Y"
# expan <- expantion of system to account for unmodelled spawning populations

process = function(ny,vcov.matrix,phi,mat,alpha,beta,sub,com,egfloor,pm.yr,for.error,OU,Rec,Spw,lst.resid,SR_rel,BH.alpha.CV,period,dir.SR,SR_devs,expan){
	ns = length(alpha) #number of sub-stocks
	for.error = for.error
	OU = OU
	m.alpha <- alpha
	m.beta <- beta

	# create vectors of time varying alpha
	if (SR_rel == "Beverton-Holt"){ 
		beta.tim <- (alpha/beta)*exp(-1)
		alpha.time <- matrix(NA,ny,length(alpha))
		for (t in 1:ny){
			alpha.time[t,] <- sin(2*pi*(t/period))*((alpha + (alpha * BH.alpha.CV)) - alpha) + alpha
		}
	}
	
	#Create recruitment deviations that are correlated among stocks 
	epi = rmvnorm(ny, sigma= vcov.matrix)

	#Build time series of Spawners (S), abundance of returning spawners pre-harvest
	# (N), and the component of the residual that is correlated throught time (v)
	R = t(matrix(0,ns,ny))
	S = R * (1-0)
	v = R; v[,]=0
	R[1:3,]=Rec
	N = array(0,dim=c(ny,4,ns))
	Ntot = R; Ntot[,]=0
	H = Ntot; S = Ntot
	S[4:7,] = Spw
	predR = Ntot
	
	# populate first few years with realized states
	R[4,] = alpha[]*S[4,]*exp(-beta[]*S[4,]+(phi*lst.resid)+epi[4,])
	predR[4,] = alpha[]*S[4,]*exp(-beta[]*S[4,])
	v[4,] = log(R[4,])-log(predR[4,])
	v[v[,]=='NaN'] <- 0
	
	for(i in 5:7){
	  R[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i,])
	  predR[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,])
	  v[i,] = log(R[i,])-log(predR[i,])
	  v[v[,]=='NaN'] <- 0	
	}
	
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
		Ntot[i,] = colSums(N[i,,])

		# apply harvest control rule
		run.size <- sum(Ntot[i,])
		if(is.na(run.size)==TRUE){run.size <- 0}
		if(run.size > 999000) {run.size <- 1000000}
		HR.all = sub_hcr(sub,com,egfloor, run.size,for.error)
		HR_adj = 1
		realized.HR <- (HR.all*HR_adj); realized.HR[realized.HR < 0] <- 0; realized.HR[realized.HR > 1] <-1
		outcome_error <- (1+rnorm(1,0,OU))
		H[i,] =  realized.HR*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
		S_exp = Ntot[i,]-H[i,]
		S_exp[S_exp<0] = 0
		S_exp[S_exp<50] = 0
		S[i,] = S_exp
		if (dir.SR == "T") {
			alpha <- m.alpha* SR_devs[i,1,]
			beta <- m.beta* SR_devs[i,2,]
				}

		# predict recruitment
		if (SR_rel == "Ricker"){
			R[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,]+phi*v[i-1,]+epi[i,])
			predR[i,] = alpha[]*S[i,]*exp(-beta[]*S[i,])
			v[i,] = log(R[i,])-log(predR[i,])
			v[v[,]=='NaN'] <- 0
		      }
	
		if (SR_rel == "Beverton-Holt"){
			R[i,] = alpha.time[i,]*S[i,]/(1+(alpha.time[i,]/beta.tim[])*S[i,])*exp(phi*v[i-1,]+epi[i,])
			predR[i,] = alpha.time[i,]*S[i,]/(1+(alpha.time[i,]/beta.tim[])*S[i,])
			v[i,] = log(R[i,])-log(predR[i,])
			v[v[,]=='NaN'] <- 0
		     }
	  }
 
	# Output
	# Performance measures:
	#	1: escapement
	#	2: harvest
	#	3: harvest rate
	#	4: predicted overfished
	#	5: predicted trending towards extinction
	#	6: empirical extinction
	#	7: proportion of years failed to meet subsistence goal
	#	8: CV in harvest
	#	9: proportion of tributary goals met
	
	pms <- matrix(NA,1,9) 

	S[S[,]=='NaN'] <- 0
	Ntot[Ntot[,]=='NaN'] <- 0
	over <- matrix(NA,length(alpha))
	ext <- matrix(NA,length(alpha))
	ext.emp <-ext
	trib.gl <-ext
	harvest_rate <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])[,1]
	harvest_rate[harvest_rate>1] <- 1
	harvest_rates <- (H[pm.yr:ny,]/Ntot[pm.yr:ny,])
	harvest_rates[harvest_rates>1] <- 1
	harvest_rate[harvest_rate=='NaN']<-1
	harvest_rates[harvest_rates=='NaN']<-1
	Smax <- round((m.alpha/m.beta)/m.alpha,digits=0)  
	ln.alpha <- log(m.alpha)
	Smsy <- round((ln.alpha*(0.5-0.07* ln.alpha))/m.beta)
	for(j in 1:length(alpha)){
		over[j] <- SC.eq(median(harvest_rates[,j]),alpha[j],beta[j])[3]
		ext[j] <- SC.eq(median(harvest_rates[,j]),alpha[j],beta[j])[4]
		ext.emp[j] <- ifelse(median(S[(ny-pm.yr):ny,j]) < ((log(alpha)/beta)*0.05)[j],1,0) # less than 5% of unfished biomass/abundance
		trib.gl[j] <- ifelse(median(S[(ny-pm.yr):ny,j]) >= (Smsy[j]),1,0)
		}

	pms[,1] <- (sum(S[pm.yr:ny,])/(ny - pm.yr +1)) * expan
	pms[,2] <- (sum(H[pm.yr:ny,])/(ny - pm.yr +1)) * expan
	pms[,3] <- median(harvest_rate)
	pms[,4] <- sum(over)/length(alpha)
	pms[,5] <- sum(ext)/length(alpha)
	pms[,6] <- sum(ext.emp)/length(alpha)
	pms[,7] <- sum(rowSums(H[pm.yr:ny,]) < (sub*0.90))/(ny - pm.yr +1)
	pms[,8] <- sd(H[pm.yr:ny,])/mean(H[pm.yr:ny,]) 
	pms[,9] <- sum(trib.gl)/length(alpha) 
	
	if (SR_rel == "Beverton-Holt"){
	  list(S=S[,],R=R[,], N=Ntot[,],H=H[,],BH_alpha = alpha.time, BH_beta = beta.tim, PMs=pms)}
	else{
	  list(S=S[,],R=R[,], N=Ntot[,],H=H[,],PMs=pms)}
	
}

#------------------------------------------------------------------------------#
# Functions to process posterior samples
#------------------------------------------------------------------------------#
# Kusko stocks are in this order (# in brackets is ID on map in ms)
 # 1.) Aniak (4)
 # 2.) George (9)
 # 3.) Holitna (7)
 # 4.) Holokuk (5)
 # 5.) Kisaralik (2)
 # 6.) Kogrukluk (8)
 # 7.) Kwethluk (1)
 # 8.) Oskawalik (6)
 # 9.) Pitka (13)
 # 10.) Swift (11)
 # 11.) Takotna (12)
 # 12.) Tatlawiksuk (10)
 # 13.) Tuluksak (3)

process.iteration = function(samp) {
  # 1.) extract names
  nms = names(samp)
  A = 4
  ns = 13
  # 2.) extract elements according to the names and put them into the appropriate data structure
  
  # parameters
  alpha = unname(samp[substr(nms, 1, 5) == "alpha"])
  beta = unname(samp[substr(nms, 1, 4) == "beta"])
  last_resid = unname(samp[substr(nms, 1, 10) == "last_resid"])
  phi = unname(samp["phi"])
  Sigma_R = matrix(samp[substr(nms, 1, 7) == "Sigma_R"], ns, ns)
  pis = c(as.numeric(samp["pi_1"]), as.numeric(samp["pi_2"]), as.numeric(samp["pi_3"]), as.numeric(samp["pi_4"]))
  
  # states
  S = matrix(samp[substr(nms, 1, 2) == "S_"], A, ns)
  R = matrix(samp[substr(nms, 1, 2) == "R_"], A - 1, ns)
  
  # 3.) create output list
  output = list(
    alpha = alpha,
    beta = beta,
    phi = phi,
    last_resid = last_resid,
    Sigma_R = Sigma_R,
    S = S,
    R = R,
    pis = pis
  )

  # 4.) return output
  return(output)

}

process.iteration.rs = function(samp) {
  # 1.) extract names
  nms = names(samp)
  s_years = 40
  r_years = 43
  ns = 13
  # 2.) extract elements according to the names and put them into the appropriate data structure
  
  # states
  S = matrix(samp[substr(nms, 1, 2) == "S_"], s_years, ns)
  R = matrix(samp[substr(nms, 1, 2) == "R_"], r_years, ns)
  
  # 2.) create output list
  output = list(
    S = S,
    R = R
  )

  # 3.) return output
  return(output)

}

#------------------------------------------------------------------------------#
# Generate deviations in alpha and beta for directional change simulations
#------------------------------------------------------------------------------#
# alpha <- median posterior etimate of alpha
# beta <- median posterior etimate of beta
# step.num <- number of time steps transition occurs over
# start.step <- time step transition begins
# ny <- number of years simulations runs over
SR_para_devs = function(alpha, beta, step.num, start.step,ny){
	id <- seq(1:13)

	end.states.order <- c( # who is switching with who? 
	4,	#1
	7,	#2
	13,	#3
	1,	#4
	5,	#5
	10,	#6
	2,	#7
	9,	#8
	8,	#9
	6,	#10
	12,	#11
	11,	#12
	3	#13
		)

	alpha.end <- alpha[end.states.order]
	beta.end <- beta[end.states.order]
	
	df <- data.frame(
	  alpha=vector(mode="numeric", length=2),
	  beta=vector(mode="numeric", length=2),
	  steps=vector(mode="numeric", length=2))
	
	input <- list(df,df,df,df,df,df,df,df,df,df,df,df,df)
	
	for(j in 1:13) {                        
	  input[[j]]$steps<- c(1,step.num)
	  input[[j]]$alpha<- c(alpha[j]/alpha[j], alpha.end[j]/alpha[j])
	  input[[j]]$beta<- c(beta[j]/beta[j], beta.end[j]/beta[j])		
	}
	  output <- array(NA,dim=c(step.num,2,13))
	  for (i in 1:length(input)) {
	    # in future; process 'n1' 2-at-a-time (1,2; 2,3; 3,4)
	    # not just 2-1,calling approx() each time and removing duplicates
	    n1 <- input[[i]]$steps[2]-input[[i]]$steps[1]+1
		if(i!= step.num){if(n1>1) {
			      ans <- approx(
			        x=input[[i]]$alpha,
			        y=input[[i]]$beta,
			        method="linear",
			        n=n1
			      )
			      if(ans$x[1] < 1){ans$x <- rev(ans$x)
			      	ans$y <- rev(ans$y)}
			      output[,,i] <- cbind(ans$x,ans$y)
			    }
		   	}
		if(i== step.num) {output[,,i] <- matrix(1, step.num,2)}   	
	  }
	
	SR_devs <- array(1,dim=c(ny,2,13))
	SR_devs[start.step:(start.step+step.num-1),,] <- output
	
	for(w in (start.step+step.num):ny){
		SR_devs[w,,] <- output[step.num,,]
	}

	return(SR_devs)	
}

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
    beta.tim <- (alpha/beta)*exp(-1)
    alpha.time <- matrix(NA,ny,length(alpha))
    for (t in 1:ny){
      alpha.time[t,] <- sin(2*pi*(t/period))*((alpha + (alpha * BH.alpha.CV)) - alpha) + alpha
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
    har_rate <- sub_hcr(sub, 0, egfloor, Ntot[i,],0)
    #har_rate = U
    outcome_error <- (1+rnorm(1,0,OU))
    H[i,] =  har_rate*Ntot[i,]*ifelse(outcome_error<0, 0, outcome_error) 
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
  
  list(S=S[,],R=R[,], N=Ntot[,],H=H[,],BH_prod = alpha.time)
  
}

