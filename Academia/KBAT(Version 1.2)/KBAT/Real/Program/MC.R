
p.tot1 <- as.matrix(p.tot1)
MAP <- as.matrix(MAP)

ACF <- function(x,lags)
{
 k <- acf(x[,],lags,"correlation",F)
 miss <- sum(k$acf,na.rm=F)
 L <- J <- 1:(lags+1)
 k1 <- abs(outer(L,J,'-'))
 k2 <- matrix(k$acf[,,1][k1+1],(lags+1))
 if(is.na(miss))
 {
	k2 <- matrix(0.9,ncol=(lags+1),nrow=(lags+1))
	diag(k2) <- 1
 }

 return(k2)
}

genMC <- function(win.p,acf.m,B)
{
 size.L <- length(win.p)
# U<-matrix(0,B,size.L)
# for (j in 1:B)
#	U[j,]<-runif(size.L,0,1)

 U <- matrix(runif(B*size.L),B,size.L)

 window.TP <- rbind(win.p,U)
 c.window <- ncol(window.TP)
 G <- t(chol(acf.m[1:c.window,1:c.window]))
 window.TP[2:(B+1),] <- t(1-matrix(pnorm(G%*%t(matrix(qnorm(1-window.TP[2:(B+1),]),B,))),,B))

 return(window.TP)
}

trun.pv <- function(win.ep,leave.i,tuo)
{
 trunc.i <- I(win.ep<=tuo)*1
 trunc.i[,leave.i] <- 1
 trunc.TP <- win.ep*trunc.i
 trunc.TP[trunc.TP==0] <- 1

 return(trunc <- list(trunc.i=trunc.i,trunc.TP=trunc.TP))
}

weight.fun <- function(trunc.i,w.map)
{
 w <- apply(trunc.i,1,function(x){x*w.map})
 w <- apply(w,2,function(x){x/sum(x)})
 w[is.na(w)] <- 0

 return(w)
}

est.fun <- function(trunc.TP,W,sel,B)
{
 MPM.rank <- PPM.rank <- WPDPPM.rank <- WLDPPM.rank <- WPDLDPPM.rank <- WPD1PPM.rank <- WPDLD1PPM.rank <- 0
 c.window <- ncol(trunc.TP)

 time <- ifelse(B>1000,ceiling(B/1000),1)
 a <- 1
 repeat
 {
	bn <- ifelse(a==time,B%%1000,1000)
	st <- 2+1000*(a-1)
	end <- ifelse((a==time & bn!=0),(st-1)+bn,(st-1)+1000)
	trunc.TP1 <- rbind(trunc.TP[1,],as.matrix(trunc.TP[st:end,]))

	if(any(sel=="PPM"))
	{
		I1 <- matrix(1,c.window,1)
		PPM <- -2*log(trunc.TP1) %*% I1
		PPM.rank <- PPM.rank+length(which(PPM[1]<=PPM[2:(end-st+2)]))
	}
	if(any(sel=="WPPM-PD"))
	{
		ifelse(c.window==1,WPDPPM <- -2*log(trunc.TP1) %*% matrix(1,1,1),
					 WPDPPM <- apply((-2*log(trunc.TP1) * W$W.PD[c(1,st:end),]),1,sum))
		WPDPPM.rank <- WPDPPM.rank+length(which(WPDPPM[1]<=WPDPPM[2:(end-st+2)]))
	}
	if(any(sel=="WPPM-LD"))
	{
		ifelse(c.window==1,WLDPPM <- -2*log(trunc.TP1) %*% matrix(1,1,1),
					 WLDPPM <- apply((-2*log(trunc.TP1) * W$W.LD[c(1,st:end),]),1,sum))
		WLDPPM.rank <- WLDPPM.rank+length(which(WLDPPM[1]<=WLDPPM[2:(end-st+2)]))
	}
	if(any(sel=="WPPM-PDLD"))
	{
		ifelse(c.window==1,WPDLDPPM <- -2*log(trunc.TP1) %*% matrix(1,1,1),
					 WPDLDPPM <- apply((-2*log(trunc.TP1) * W$W.PDLD[c(1,st:end),]),1,sum))
		WPDLDPPM.rank <- WPDLDPPM.rank+length(which(WPDLDPPM[1]<=WPDLDPPM[2:(end-st+2)]))
	}
	if(any(sel=="KBAT-PD"))
	{
		ifelse(c.window==1,WPD1PPM <- -2*log(trunc.TP1) %*% matrix(1,1,1),
					 WPD1PPM <- apply((-2*log(trunc.TP1) * W$K.PD[c(1,st:end),]),1,sum))
		WPD1PPM.rank <- WPD1PPM.rank+length(which(WPD1PPM[1]<=WPD1PPM[2:(end-st+2)]))
	}
	if(any(sel=="KBAT-PDLD"))
	{
		ifelse(c.window==1,WPDLD1PPM <- -2*log(trunc.TP1) %*% matrix(1,1,1),
					 WPDLD1PPM <- apply((-2*log(trunc.TP1) * W$K.PDLD[c(1,st:end),]),1,sum))
		WPDLD1PPM.rank <- WPDLD1PPM.rank+length(which(WPDLD1PPM[1]<=WPDLD1PPM[2:(end-st+2)]))
	}

	a <- a+1
	if(a>time)
		break
 }

 if(any(sel=="MPM"))
 {
	MPM <- apply(trunc.TP,1,min)
	MPM.rank <- length(which(MPM[1]>=MPM[2:(B+1)]))
 }

 stat.pv <- c(MPM.rank,PPM.rank,WPDPPM.rank,WLDPPM.rank,WPDLDPPM.rank,WPD1PPM.rank,WPDLD1PPM.rank)/B
 names(stat.pv) <- c("MPM","PPM","WPPM-PD","WPPM-LD","WPPM-PDLD","KBAT-PD","KBAT-PDLD") 

 stat.pv <- stat.pv[names(stat.pv)%in%sel]

 return(stat.pv)
}


all <- NULL
for(tuo in Tuo.size)
{
	for(l in op)
	{
		for(i in St:End)
		{
			if(p.tot1[i]<bonf)
			{
				all <- rbind(all,rep(p.tot1[i],length(Stat)))
				next
			}

			if(borw2==1)
			{
				m <- l
				window.H <- ifelse(i-m<St,St,i-m)
				window.T <- ifelse(i+m>End,End,i+m)
				win.m <- MAP[window.H:window.T]
				win.p <- p.tot1[window.H:window.T]
				size.k <- min((End-St),(2*m))
				leave.i <- ifelse(i-m<St,i-St+1,m+1)
			}
			else
			{
				h <- l
				win.s <- which(abs(MAP[i]-MAP)<=h & MAP>=MAP[St] & MAP<=MAP[End])
				win.m <- MAP[win.s]
				win.p <- p.tot1[win.s]
				size.k <- min((End-St),length(win.m))
				leave.i <- which(win.m==MAP[i])
			}
#			set.seed(i)
			k2 <- ACF(matrix(p.tot1[St:End]),size.k)
			window.TP <- genMC(win.p,k2,B)
			trunc <- trun.pv(window.TP,leave.i,tuo)
			trunc.i <- trunc$trunc.i
			trunc.TP <- trunc$trunc.TP

			if(sum(trunc.i[1,])==1)
			{
				all <- rbind(all,rep(p.tot1[i],length(Stat)))
				next
			}

			if(borw=="Window")
			{
				h1 <- max(abs(MAP[i]-MAP[c(window.H,window.T)]))
				h2 <- ifelse(abs(MAP[window.H]-MAP[i])>=abs(MAP[window.T]-MAP[i]),window.H,window.T)
				h3 <- ifelse(h2>i,ifelse(h2==End,0.5,(MAP[h2+1]-MAP[h2])/2),
					ifelse(h2==St,0.5,(MAP[h2]-MAP[h2-1])/2))
				h4 <- h3+h1
			}

			window.map <- abs(MAP[i]-win.m)
			window.map1 <- 1/(1+window.map)
			window.map2 <- window.map/ifelse(borw=="Window",h4,h)
			window.map3 <- (3/4*(1-(window.map2)^2))/sum(3/4*(1-(window.map2)^2))
			tmp <- ifelse(LDdata!=0,ifelse(borw=="Window",window.ld<-LD[i,window.H:window.T],window.ld<-LD[i,which(abs(MAP[i]-MAP)<=h)]),window.ld<-0)

			tmp <- ifelse(any(Stat=="WPPM-PD"),W.PD<-t(weight.fun(trunc.i,window.map1)),W.PD<-0)
			tmp <- ifelse(any(Stat=="WPPM-LD"),W.LD<-t(weight.fun(trunc.i,window.ld)),W.LD<-0)
			tmp <- ifelse(any(Stat=="WPPM-PDLD"),W.PDLD<-t(weight.fun(trunc.i,window.map1*window.ld)),W.PDLD<-0)
			tmp <- ifelse(any(Stat=="KBAT-PD"),K.PD<-t(weight.fun(trunc.i,window.map3)),K.PD<-0)
			tmp <- ifelse(any(Stat=="KBAT-PDLD"),K.PDLD<-t(weight.fun(trunc.i,window.map3*window.ld)),K.PDLD<-0)
			window.w <- list(W.PD=W.PD,W.LD=W.LD,W.PDLD=W.PDLD,K.PD=K.PD,K.PDLD=K.PDLD)
			pv <- est.fun(trunc.TP,window.w,Stat,B)

			check <- which(pv==0)
			bonf.num <- ceiling(1/bonf)
			if(length(check)!=0)
			{
				if(bonf.num>B)
				{
					a <- ifelse(bonf.num-B>10000,ceiling((bonf.num-B)/10000),1)
					check.pv <- rep(0,length(check))

					repeat
					{
						b <- ifelse(a==1,(bonf.num-B)%%10000,10000)

						window.TP <- genMC(win.p,k2,b)
						trunc <- trun.pv(window.TP,leave.i,tuo)
						trunc.i <- trunc$trunc.i
						trunc.TP <- trunc$trunc.TP

						tmp <- ifelse(any(Stat=="WPPM-PD"),W.PD<-t(weight.fun(trunc.i,window.map1)),W.PD<-0)
						tmp <- ifelse(any(Stat=="WPPM-LD"),W.LD<-t(weight.fun(trunc.i,window.ld)),W.LD<-0)
						tmp <- ifelse(any(Stat=="WPPM-PDLD"),W.PDLD<-t(weight.fun(trunc.i,window.map1*window.ld)),W.PDLD<-0)
						tmp <- ifelse(any(Stat=="KBAT-PD"),K.PD<-t(weight.fun(trunc.i,window.map3)),K.PD<-0)
						tmp <- ifelse(any(Stat=="KBAT-PDLD"),K.PDLD<-t(weight.fun(trunc.i,window.map3*window.ld)),K.PDLD<-0)
						window.w <- list(W.PD=W.PD,W.LD=W.LD,W.PDLD=W.PDLD,K.PD=K.PD,K.PDLD=K.PDLD)

						tmp <- ifelse(sum(is.element(Stat,"SLM"))==1,stat<-Stat[check+1],stat<-Stat[check])
						check.pv <- check.pv+b*est.fun(trunc.TP,window.w,stat,b)

						a <- a-1
						if(a==0)
							break
					}
					pv[check] <- (pv[check]*B+check.pv)/bonf.num
				}
				pv[pv==0] <- bonf
			}

			if(any(Stat=="SLM"))
				all <- rbind(all,c(p.tot1[i],pv))
			else
				all <- rbind(all,pv)
		}
	}
}

plottable <- paste(KBAT.path,"\\Program\\Plot.R",sep="")
source(plottable)




