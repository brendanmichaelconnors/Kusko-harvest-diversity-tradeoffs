########################################################################################
# figures.R
#
# Code to generate equilibrium trade-off, time varying populaiton diversity and 
#	harvest control rule figures 
# July 5, 2019
# Author: B. Connors (DFO)
#        
########################################################################################

# --- Colors and stock IDs  -----------------------------------------------------------
MyColour <- viridis(13)
MyTextColour <- c("white","white","white","white","white","white","white","white","white","black","black","black","black")
stock_id <- c(4,9,7,5,2,8,1,6,13,11,12,10,3)

# --- Calculate outcomes across a range of mixed-stock harvest rates  ------------------
wide.posteriors <- array(NA,dim=c(10000,13,2),dimnames=list(NULL,seq(1,13),c("alphas","betas")))
wide.posteriors[,,1]<-samps[1:10000,1:13]
wide.posteriors[,,2]<-samps[1:10000,14:26]

alphas <-13
betas <- 13

# estimate equilibrium yield profile and corresponding risk
t3 <- array(NA,dim=c(length(U),4,10000))

U <- seq(0,1,0.01)

for(w in 1:10000){
	draw <- sample(10000,1)
	aa<-seq(1, alphas,1); bb <- seq(1, alphas,1); for (k in 1: alphas){aa[k] <- log(wide.posteriors[draw,k,1]); bb[k] <- wide.posteriors[draw,k,2] }
	for (i in 1:length(U)) {
  		t1 <- matrix(nrow=length(aa),ncol=4)
  		for (j in 1:length(aa)){
   			 t1[j,] <- SC.eq(U[i],aa[j],bb[j])
   		 }
  		t3[i,,w] <- apply(t1,2,sum); t3[i,3:4,w] <- t3[i,3:4,w]/length(aa)
  	}
}

t3[,1:2,] <- t3[,1:2,]*1.75
t3.median <- apply(t3,1:2,quantile,probs=c(0.5),na.rm=T)
t3.upper <-apply(t3,1:2,quantile,probs=c(0.9),na.rm=T)
t3.lower <- apply(t3,1:2,quantile,probs=c(0.1),na.rm=T)
t3.median[101,2]<-0
t3.upper[101,2]<-0
t3.lower[101,2]<-0

yyy <- as.numeric(t3.median[,3])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
over.med <- predict(xx)
over.med[over.med <0] =0 
over.med[over.med >1] =1 
over.med[1:8]=0

yyy <- as.numeric(t3.upper[,3])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
over.up <- predict(xx)
over.up[over.up <0] =0 
over.up[over.up >1] =1 
over.up[97:101] =1 
over.up[1:8]=0

yyy <- as.numeric(t3.lower[,3])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
over.low <- predict(xx)
over.low[over.low <0] =0 
over.low[over.low >1] =1 
over.low[97:101] =1 
over.low[1:8]=0

yyy <- as.numeric(t3.median[,4])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
ext.med <- predict(xx)
ext.med[ext.med <0] =0 
ext.med[ext.med >1] =1 
#ext.med[1:8]=0

yyy <- as.numeric(t3.upper[,4])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
ext.up <- predict(xx)
ext.up[ext.up <0] =0 
ext.up[ext.up >1] =1 
#ext.up[97:101] =1 
#ext.up[1:8]=0

yyy <- as.numeric(t3.lower[,4])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
ext.low <- predict(xx)
ext.low[ext.low <0] =0 
ext.low[ext.low >1] =1 
ext.low[1:8]=0

# --- Figures  -----------------------------------------------------------------

	#------------------------------------------------------------------------------#
	# (a) productivity vs size plot and (b) equilibrium tradeoffs
	#------------------------------------------------------------------------------#
	
	m.alpha <- apply(wide.posteriors[,,1],c(2),quantile,probs=c(0.5),na.rm=T)
	alpha.upper <- apply(wide.posteriors[,,1],c(2),quantile,probs=c(0.95),na.rm=T)
	alpha.lower<- apply(wide.posteriors[,,1],c(2),quantile,probs=c(0.05),na.rm=T)
	a <-length(m.alpha); b <- a
	size <- (log(wide.posteriors[,,1])/wide.posteriors[,,2])/1000
	
	m.size <- apply(size,c(2),quantile,probs=c(0.5),na.rm=T)
	size.lower <- apply(size,c(2),quantile,probs=c(0.05),na.rm=T)
	size.upper <- apply(size,c(2),quantile,probs=c(0.95),na.rm=T)
	
	jpeg("figures/fig_3.July2019.jpeg",width=7.75, height=3.5, units="in",res=800)
	
		#dev.new(width=7.75, height=3.5,new=FALSE)
		par(mfrow=c(1,2),bty="o", mar=c(3,3,1,3),oma=c(1,1,1,1))#set dimensions to plots
		
		# --- panel a --------------------------------------------------------------
		plot(m.size,m.alpha,yaxt="n",col="white",xlim=c(0,25),ylim=c(1,8))
		axis(2,las=2)
		
		for(i in 1:13){
			xlcl = size.lower[i]; xucl= size.upper[i]; y= m.alpha[i]; ylcl= alpha.lower[i]; yucl =alpha.upper[i];x= m.size[i]
			lines(x=c(xlcl,xucl),y=c(y,y),col=MyColour[i]);lines(x=c(x,x),y=c(ylcl,yucl),col= MyColour[i])
		}
		
		points(m.size,m.alpha, cex=1.5, pch=16,col=MyColour)
		text(m.size,m.alpha, stock_id,cex=0.5,col= MyTextColour)
		axis(2,labels=F)
		
		text(0.5,7.8,"(a)")
		mtext("Productivity ",2,line=2.75,cex=1)
		mtext("(recruits/spawner)",2,line=1.75,cex=1)
		mtext("Equilibrium size (000s)",1,line=2,cex=1)
		box(col="grey")
		# panel b
		plot(U*100, t3.upper[,2]/1000,type="l",yaxt="n",ylab="",xlab="",lwd=1,lty=2,ylim=c(0,110))
		axis(2,las=2,col="blue",col.axis="blue")
		
		polygon(c(U*100,rev(U*100)),c((t3.upper[,2]/1000),rev(t3.lower[,2]/1000)),col="#0000FF25",border=NA)
		points(U*100, (t3.median[,2]/1000),type="l",col="blue",lwd=3)
		points(U*100, (t3.upper[,2]/1000),type="l",col="blue",lwd=1,lty=2)
		points(U*100, (t3.lower[,2]/1000),type="l",col="blue",lwd=1,lty=2)
		text(3,105,"(b)")
		
		#panel b
		par(new=TRUE)
		plot(U*100, ext.med,type="l",yaxt="n",xaxt="n",ylab="",xlab="",lwd=3,col="red")
		polygon(c(U*100,rev(U*100)),c((ext.up),rev(ext.low)),col="#FF000025",border=NA)
		points(U*100, ext.med,type="l",col="red",lwd=3,lty=1)
		points(U*100, ext.up,type="l",col="red",lwd=1,lty=2)
		points(U*100, ext.low,type="l",col="red",lwd=1,lty=2)
		axis(4,las=2,labels=c(0,20,40,60,80,100),at=c(0,0.2,0.4,0.6,0.8,1),col.ticks="red",col.axis="red")
		box(col="grey")
		
		mtext("Harvest (000s)",2,line=2.25,col="blue")
		mtext("Harvest rate (%)",1,line=2.25)
		mtext("Populations extirpated (%) ",4,line=2.5,col="red")
	
	dev.off()
	
	#------------------------------------------------------------------------------#
	# Time varying population diversity 
	#------------------------------------------------------------------------------#
	
	alpha <- apply(samps[,1:13],c(2),quantile,probs=c(0.5),na.rm=T)
	beta <- apply(samps[,14:26],c(2),quantile,probs=c(0.5),na.rm=T)
	size <- (log(alpha)/beta)/1000
	id <- seq(1:13)
	
	end.states.order <- c(
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
	size.end <- (log(alpha.end)/beta.end)/1000
	
	jpeg("figures/TV_alpha_beta.Mar42019.jpeg",width=3.5, height=3.25, units="in",res=400)
		#dev.new(width=3.5, height=3.25,new=FALSE)
		par(mfrow=c(1,1),bty="o", mar=c(1,2,2,1),oma=c(2,2.5,0,0))#set dimensions to plots
		
		plot(size,alpha,pch=16,cex=2,col= MyColour,xlim=c(0,20),ylim=c(1,6.5),yaxt="n")
		axis(2,at=c(1,2,3,4,5,6,7),las=2)
		Arrows(size,alpha, size.end,alpha.end,lwd=1, arr.type="triangle",col="light grey",arr.adj=1,arr.tpe="triangle")
		points(size,alpha,pch=16,cex=2,col= MyColour)
		mtext("Productivity ",2,line=3.5,cex=1)
		mtext("(recruits/spawner)",2,line=2.5,cex=1)
		mtext("Equilibrium size (000s)",1,line=2,cex=1)
		box(col="grey")
	dev.off()