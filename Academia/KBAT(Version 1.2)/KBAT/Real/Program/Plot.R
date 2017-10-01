

p2<-all
all<--log10(all)



MAP<-MAP[St:End,]

titlborw<-ifelse(borw=="Window","m",borw)
h<-length(op)
k<-length(Stat)
n<-length(Tuo.size)

par(mfrow=c(4,2))
for(i in 1:h)
{
	for(m in 1:n)
	{
		for(j in 1:k)
		{
			titlborw<-ifelse(borw=="Window","m",borw)
			plot(MAP,all[((End-St+1)*((i+(h*(m-1)))-1)+1):((End-St+1)*(i+(h*(m-1)))),j],ylim=c(0,max(all)),xlim=c(min(MAP),max(MAP)),xlab=Xlabe, 
			     ylab = "-log(p value)",type = "n", axes =T)
			title(main=paste("Theta=",Tuo.size[m],", ",Stat[j],", ",titlborw,"=",op[i],sep=""))
			lines(MAP,all[((End-St+1)*((i+(h*(m-1)))-1)+1):((End-St+1)*(i+(h*(m-1)))),j],col=4,lty=1)
		}

		if(k<8)
		{
			for(y in 1:abs(8-k))
				plot(c(1,2),c(1,2),type="n",axes=F,xlab="",ylab="")
		}
		savePlot(filename=paste(PLOT.out,"_",titlborw,"_",op[i],"_Theta_",Tuo.size[m],".jpeg",sep=""),type="jpg")
	}
}



PATH.power<-paste(PATHout,"\\Output.txt",sep="")

times2<-0
for(m in 1:n)
{
	for(i in 1:h)
	{
		b1<-matrix("-",1,73)
		if(times2==0)
			write.table(b1,PATH.power,quote=F,row.names=F,col.names=F,append=F)
		else
			write.table(b1,PATH.power,quote=F,row.names=F,col.names=F,append=T)

		titlborw<-ifelse(borw=="Window","m",borw)
		c1<-paste("Table: P-values of all statistics for each marker while the" ,titlborw, "is",op[i],"and truncation threshold is",Tuo.size[m],sep=" ")
		write.table(c1,PATH.power,quote=F,row.names=F,col.names=F,append=T)

		b<-matrix("-",1,73)
		write.table(b,PATH.power,quote=F,row.names=F,col.names=F,append=T)
		a<-matrix(0,1,(length(Stat)+1))
		a[1,]<-c("Marker",Stat)
		a[1,]<-as.matrix(format(a[1,],width=15,justify="right"))
		write.table(a,PATH.power,quote=F,row.names=F,col.names=F,na = "NA",append=T)
		b<-matrix("-",1,73)
		write.table(b,PATH.power,quote=F,row.names=F,col.names=F,append=T)
		c<-matrix(0,Num,(length(Stat)+1))
		c[1:(End-St+1),1]<-c(1:(End-St+1))
		c[1:(End-St+1),1]<-as.matrix(format(c[1:(End-St+1),1],width=15,justify="right"))
		data<-NULL
		data<-as.matrix(p2[((End-St+1)*((i+(h*(m-1)))-1)+1):((End-St+1)*(i+(h*(m-1)))),])
		c[1:(End-St+1),2:(length(Stat)+1)]<-format(data,nsmall=4,justify="right",width=15)
		write.table(c,PATH.power,quote=F,row.names=F,col.names=F,append=T)

		b<-matrix("-",1,73)
		write.table(b,PATH.power,quote=F,row.names=F,col.names=F,append=T)
		b<-matrix(" ",1,73)
		write.table(b,PATH.power,quote=F,row.names=F,col.names=F,append=T)
		times2<-times2+1
	}
}

msgover <- paste("The computation of KBAT is finished.",sep="")
tkmessageBox(title="KBAT message",message=msgover,icon="info",type="ok")


