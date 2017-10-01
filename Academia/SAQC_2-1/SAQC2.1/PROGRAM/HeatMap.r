##---------------------- Interface parameter -----------------------------------
popu_file=tclvalue(QIfigure_inputpath)
temp_output_QIfigure=tclvalue(QIfigure_outputpath)
if(popu_file!="Example")
{
	if(file.exists(popu_file)==F){stop("SAQC error message: The pathway of input is incorrect!!")}
	if(file.exists(popu_file)==T & file.exists(temp_output_QIfigure)==F){stop("SAQC error message: The pathway of output is incorrect!!")}  ## YTLin
	if(file.exists(temp_output_QIfigure)==F){stop("SAQC error message: The pathway of output is incorrect!!")}  ## YTLin
	output_QIfigure = temp_output_QIfigure  ## YTLin
}
if(popu_file=="Example")
{
	if(file.exists(temp_output_QIfigure)==T){output_QIfigure=paste(temp_output_QIfigure,"/Test_QIfigure",sep="")}
	if(file.exists(temp_output_QIfigure)==F){output_QIfigure=paste(pathway,"OUTPUT/Test_QIfigure",sep="")}
	dir.create(output_QIfigure,showWarnings=F)
	popu_file=paste(pathway,"EXAMPLE/Test_QIfigure",sep="")
}
output_popu=paste(output_QIfigure,"/HeatMap and Polygon figure",sep="")
dir.create(output_popu,showWarnings=F)

## Plot HeatMap, Polygon or both
QIfigure=c(tclvalue(QIfigure.Value1),tclvalue(QIfigure.Value2),tclvalue(QIfigure.Value3))
QIfigure=as.numeric(QIfigure)
if(sum(QIfigure)==0)
{stop("Please make sure that the option of plotting is correct!!")}

QIsource=c(tclvalue(QI_resouce.val))  ## 1: SAQC database, 2: User-provided database
QIsource=as.numeric(QIsource)
ipopu=as.numeric(tclvalue(QIfigure_groupValue))
ichip=as.numeric(tclvalue(QIfigure_chipValue))
Quantile<-c(as.numeric(tclvalue(QIfigure_quantile.Value1)),as.numeric(tclvalue(QIfigure_quantile.Value2)),as.numeric(tclvalue(QIfigure_quantile.Value3)),as.numeric(tclvalue(QIfigure_quantile.Value4)))

if(QIfigure[2] == 1){  ## Plot Polygon: yes
	if(QIsource==1) ## 1: SAQC database   		# Upperbound of QI
	{
	  #popu=Asia,Africa,Europe,Taiwan,Combine
	  pathway_DB=paste(pathway,"DATABASE",sep="")
	  chipname=c("Affy 100K","Affy 500K")
	  pathway_chip=paste(pathway_DB,"/",chipname[ichip],sep="")
	  popuname=c("Asian","YRI","CEU","Taiwan","Combined")
	  pathway_popu=paste(pathway_chip,"/",popuname[ipopu],sep="")
	  QIchip1_file=paste(pathway_popu,"/Upper_bound/Chip1.txt",sep="")
	  QIchip2_file=paste(pathway_popu,"/Upper_bound/Chip2.txt",sep="")
	  QImerge_file=paste(pathway_popu,"/Upper_bound/Merge.txt",sep="")
	  QIchip1=read.table(QIchip1_file,head=T,sep="\t")[-c(3,6),]
	  QIchip2=read.table(QIchip2_file,head=T,sep="\t")[-c(3,6),]
	  QImerge=read.table(QImerge_file,head=T,sep="\t")[-c(3,6),]

	  wholeUB=NULL;
	  #95%  [[1]][[1]][[2]]
	  #975% [[1]][[2]]
	  #99%  [[2]]
	  for(k in 1:3) #for 95%, 97.5%, 99% respectively
	  {
		UB=c(QIchip1[,k+1],QIchip2[,k+1],QImerge[,k+1])
		UB=UB[c(3,7,11,4,8,12,1,5,9,2,6,10)] #Median QI1, Median QI2, WinsQI1, WinsQI2
		UB=matrix(UB,4,3, byrow=T) #column: Chip1, Chip2, Merge; row: Median QI1, Median QI2, WinsQI1, WinsQI2
		wholeUB=list(wholeUB,UB)
	  }
	}

	if(QIsource==2)  ## 2: User-provided database  ## user defined
	{

	  UBQI95=c(tclvalue(MQI1_chip1.q95.val),tclvalue(MQI1_chip2.q95.val),tclvalue(MQI1_merge.q95.val),
	  tclvalue(MQI2_chip1.q95.val),tclvalue(MQI2_chip2.q95.val),tclvalue(MQI2_merge.q95.val),
	  tclvalue(WQI1_chip1.q95.val),tclvalue(WQI1_chip2.q95.val),tclvalue(WQI1_merge.q95.val),
	  tclvalue(WQI2_chip1.q95.val),tclvalue(WQI2_chip2.q95.val),tclvalue(WQI2_merge.q95.val))
	  UBQI975=c(tclvalue(MQI1_chip1.q975.val),tclvalue(MQI1_chip2.q975.val),tclvalue(MQI1_merge.q975.val),
	  tclvalue(MQI2_chip1.q975.val),tclvalue(MQI2_chip2.q975.val),tclvalue(MQI2_merge.q975.val),
	  tclvalue(WQI1_chip1.q975.val),tclvalue(WQI1_chip2.q975.val),tclvalue(WQI1_merge.q975.val),
	  tclvalue(WQI2_chip1.q975.val),tclvalue(WQI2_chip2.q975.val),tclvalue(WQI2_merge.q975.val))
	  UBQI99=c(tclvalue(MQI1_chip1.q99.val),tclvalue(MQI1_chip2.q99.val),tclvalue(MQI1_merge.q99.val),
	  tclvalue(MQI2_chip1.q99.val),tclvalue(MQI2_chip2.q99.val),tclvalue(MQI2_merge.q99.val),
	  tclvalue(WQI1_chip1.q99.val),tclvalue(WQI1_chip2.q99.val),tclvalue(WQI1_merge.q99.val),
	  tclvalue(WQI2_chip1.q99.val),tclvalue(WQI2_chip2.q99.val),tclvalue(WQI2_merge.q99.val))

	  wholeUB=NULL;
	  for(k in 1:3)
	  {
		if(k==1){UB=UBQI95}
		if(k==2){UB=UBQI975}
		if(k==3){UB=UBQI99}
		UB=as.numeric(UB)
		UB=matrix(UB,4,3, byrow=T)
		wholeUB=list(wholeUB,UB)
	  }
	}
	Q95=wholeUB[[1]][[1]][[2]] #row: Median QI1, Median QI2, WinsQI1, WinsQI2; column: Chip1, Chip2, Merge
	Q975=wholeUB[[1]][[2]]
	Q99=wholeUB[[2]]
	UB_MQI1=cbind(Q95[1,],Q975[1,],Q99[1,]) #row: Chip1, Chip2, Merge; column: 95%, 97.5%, 99%
	UB_MQI2=cbind(Q95[2,],Q975[2,],Q99[2,])
	UB_WQI1=cbind(Q95[3,],Q975[3,],Q99[3,])
	UB_WQI2=cbind(Q95[4,],Q975[4,],Q99[4,])

	if(QIsource==1) ## 1: SAQC database  ## if the quantile option is empty.
	{
	  iQ=c(tclvalue(QIfigure_quantile.Value1),tclvalue(QIfigure_quantile.Value2),tclvalue(QIfigure_quantile.Value3))
	  iQ=as.numeric(iQ)
	  for(j in 1:3)
	  {
		if(iQ[j]==0)
		{
		  UB_MQI1[,j]=NA;UB_MQI2[,j]=NA
		  UB_WQI1[,j]=NA;UB_WQI2[,j]=NA
		}
	  }
	}
} #end of QIfigure[2] == 1
##------------------------------------------------------------------------------
libname=.packages(all.available = TRUE)
if(length(nchar(libname[libname=="TeachingDemos"]))==0)
{install.packages("TeachingDemos", repos="http://cran.csie.ntu.edu.tw")}
if(length(nchar(libname[libname=="tkrplot"]))==0)
{install.packages("tkrplot", repos="http://cran.csie.ntu.edu.tw")}
library(TeachingDemos)
##-------------------------- call subfunction ----------------------------------
pathway_program=paste(pathway,"/PROGRAM",sep="")
HeatMap=paste(pathway_program,"/HeatMap_subfunction.r",sep="")
source(HeatMap)
##==============================================================================

QI_file=list.files(popu_file,full.names=T)
QI_name=list.files(popu_file,full.names=F)
QI_name=strsplit(QI_name,".txt");QI_name=unlist(QI_name)
QI_name[1]=paste(QI_name[1],"     ",sep="")
QI_name[2]=paste(QI_name[2],"     ",sep="")

#chip_name=c("Merge","Chip2","Chip1")
#chip_name=paste(chip_name,"               ",sep="")
chip_name = c(format("Merge", width = 20), format("Chip2", width = 22), format("Chip1", width = 22))  
MQI1=read.table(QI_file[1],head=T,sep="\t")
MQI2=read.table(QI_file[2],head=T,sep="\t")
WQI1=read.table(QI_file[3],head=T,sep="\t")
WQI2=read.table(QI_file[4],head=T,sep="\t")
index.2chip = "Merge" %in% colnames(MQI1)

if(index.2chip){
	index.column = c(1, 3, 4, 2)
} else {
	index.column = c(1, 2, 2, 2)  ## Fill the other two columns with same data to avoid the only-one-chip problem
	if(colnames(MQI1)[2] %in% c("Hind", "Nsp")){
		chip_name = rep(format("Chip1", width = 22), 3)
		index.chip = 1
	} else {
		chip_name = rep(format("Chip2", width = 22), 3)
		index.chip = 2
	}
	if(QIfigure[2]){
		UB_MQI1=matrix(rep(UB_MQI1[index.chip,], 3), ncol = 3, byrow = T) #row: Chip1, Chip2, Merge; column: 95%, 97.5%, 99%
		UB_MQI2=matrix(rep(UB_MQI2[index.chip,], 3), ncol = 3, byrow = T)
		UB_WQI1=matrix(rep(UB_WQI1[index.chip,], 3), ncol = 3, byrow = T)
		UB_WQI2=matrix(rep(UB_WQI2[index.chip,], 3), ncol = 3, byrow = T)
	}
}
MQI1 = MQI1[, index.column]
MQI2 = MQI2[, index.column]
WQI1 = WQI1[, index.column]
WQI2 = WQI2[, index.column]

#chip_name=format(colnames(MQI1)[4:2])
#chip_name=c(format(colnames(MQI1)[4], width = 20), format(colnames(MQI1)[3], width = 23), format(colnames(MQI1)[2], width = 23))
#chip_name=c(format(colnames(MQI1)[4], width = 20), format("Sty", width = 23), format("Nsp", width = 23))

RangeMQI1=range(as.numeric(as.matrix(MQI1[,-1])));
RangeMQI2=range(as.numeric(as.matrix(MQI2[,-1])));
RangeWQI1=range(as.numeric(as.matrix(WQI1[,-1])));
RangeWQI2=range(as.numeric(as.matrix(WQI2[,-1])));
RangeQI=rbind(RangeMQI1,RangeMQI2,RangeWQI1,RangeWQI2)
MeanQI=apply(RangeQI,1,mean)
MaxQI=ceiling(RangeQI[,2])
QI_range=cbind(0,RangeQI,MeanQI,MaxQI)
colnames(QI_range)[1]="inf_min" 
colnames(QI_range)[2]="Min" 
colnames(QI_range)[3]="Max" 
colnames(QI_range)[4]="Median" 
colnames(QI_range)[5]="sup_max" 
QI_range=QI_range[,c(1,2,4,3,5)] 
#z=z[order(z[,4],decreasing=T),]

##------------------------------------------------------------------------------
#for polygon plot
Mintable=NULL;Maxtable=NULL;
for(k in 1:4)
{
  #k=1
  if(k==1){QItable=MQI1}
  if(k==2){QItable=MQI2}
  if(k==3){QItable=WQI1}
  if(k==4){QItable=WQI2}
  QItable=QItable[,-1]
  wholeRangeQI=NULL;
  for(j in 1:3)    #Chip1, Chip2, Merge
  {
    RangeQI=range(QItable[,j])
    wholeRangeQI=rbind(wholeRangeQI,RangeQI)
  }
  Mintable=rbind(Mintable,t(wholeRangeQI[,1])) #row: MedianQI1, MedianQI2, WinsQI1, WinsQI2; column: Chip1, Chip2, Merge
  Maxtable=rbind(Maxtable,t(wholeRangeQI[,2]))
}
##       chip1 chip2 merge
##  MQI1
##  MQI2
##  WQI1
##  WQI2

##------------------------------------------------------------------------------
#for HeatMap
z=if(index.2chip){
	c(as.numeric(as.matrix(MQI1[,-1])),as.numeric(as.matrix(MQI2[,-1])),as.numeric(as.matrix(WQI1[,-1])),as.numeric(as.matrix(WQI2[,-1])))
} else c(MQI1[,2], MQI2[,2], WQI1[,2], WQI2[,2])
MeanZ=mean(z);StdZ=sd(z);
QImin=min(z);
QImax=MeanZ+6*StdZ
RangeZ=c(QImin,QImax)
d=diff(RangeZ)/2
r=d/20
inf_Zmin=0
sup_Zmin=QImin+d
Zmin=MeanZ
#Zmin=QImin
#Zmin=0

inf_Zmax=QImax-d
sup_Zmax=QImax+d
Zmax=MeanZ+6*StdZ #QImax
nWidth=14

if(QIfigure[1]==1){ #HeatMap
	HeatMap <- function(QI_Selection=QI_name[4],Min_QI=QImin,Max_QI=QImax, N_color=10,Sort_By_Array=chip_name[1]){ #,Print_option=Print)#,Reset_option=Reset
		#QIstat=QI_name[4];MinQI=QImin;MaxQI=QImax;Ncolor=10;ChipSelection="Merge";;Print="Yes"
		QI_selection=QI_Selection
		cols <-terrain.colors(N_color)
		#if(Reset_option==T){MinQI=Zmin;MaxQI=Zmax}
		if(QI_selection==QI_name[1]){QI=MQI1}
		if(QI_selection==QI_name[2]){QI=MQI2}
		if(QI_selection==QI_name[3]){QI=WQI1}
		if(QI_selection==QI_name[4]){QI=WQI2}
		# iRank<- switch(ChipSelect,2="Merge",3="Chip2",4="Chip1")
		if(Sort_By_Array==chip_name[1]){iRank=4}
		if(Sort_By_Array==chip_name[2]){iRank=3}
		if(Sort_By_Array==chip_name[3]){iRank=2}
		z=QI[order(QI[,iRank],decreasing=T),]
		indname=z[,1]
		z=z[,-1];z=as.matrix(z);z=matrix(as.numeric(z),nrow(z),ncol(z))
		nind=nrow(z)
		## For convience, Hsin-Chi copied one-chip data to three ata column instead of changing programs
		## and then just set different ylim to avoid showing the duplicate copied one-chip result
		y.lim = if(index.2chip){
			3
		} else 1
		y=c(1:3)
		x=c(1:nind)
		par(mai=c(0.75,0.35,0.25,0.35))
		##==========================  HeatMap figure ============================
		##==========================  HeatMap ===================================
		##  replace the extremely large data points with (mean+6*sd) 
		##	replace the extremely small data points with (min)
		print(Min_QI);print(Max_QI)
		tempz=z
		print(range(tempz))
		tempz[tempz<Min_QI]=Min_QI;  ##  QImax=MeanZ+6*StdZ
		tempz[tempz>Max_QI]=Max_QI   ##	QImin=min(z);
		print(Min_QI)
		print(range(tempz))
		m1 <<- Min_QI
		masaki <<- tempz
		#plot(1,1,ylim=c(0.5,3.5),xlim=c(0.5,(nind+0.5)),xlab="",ylab="",axes=F)
		dspace=0.5
		#image(x, y, as.matrix(tempz[,1]), col = cols, 
		image(x, y, tempz, col = cols, 
			#ylim=c(0.5,3.5+dspace),
			ylim=c(0.5, y.lim + 0.5 + dspace),
			#ylim=c(0.5, 1.5 + dspace),
			xlim=c(0.5,(nind+0.5)),zlim=c(Min_QI,Max_QI),axes=F,ylab="",xlab="")
		##-------------------------------------------
		if(nind<61)
		{
			Xlist=(rep(x,3));Ylist=c(rep(1,nind),rep(2,nind),rep(3,nind))
			a=as.numeric(z);a=round(a*100)/100
			text(Xlist,Ylist,labels=a,srt = 60,col="blue",cex=0.75)
		}
		##-------------------------------------------
		block=10;n=ceiling(nind/block);Xlab=c(1:n)*block
		Xlab[length(Xlab)]=nind;
		x=c(0,x);y=c(0,y)
		abline(v=c(x+0.5),lty=2,col="skyblue");abline(h=c(y+0.5),lty=1,col="blue");abline(v=c(Xlab+0.5),lty=1,col="blue");
		#axis(2,at=c(1,2,3),labels=c("Chip1","Chip2","Merge"),tick=F,font=4,line=-0.5,col.axis="blue")
		axis(2,at=1:y.lim,labels=if(index.2chip){ c("Chip1","Chip2","Merge") } else unlist(strsplit(Sort_By_Array, " "))[1],
			tick=F,font=4,line=-0.5,col.axis="blue")
		#axis(2,at=c(1,2,3),labels=colnames(QI)[-1],tick=F,font=4,line=-0.5,col.axis="blue")
		x=x[-1]
		axis(1,at=x,labels=indname,col.axis="blue",tick=F,las=2,line=-2.25,cex.axis=0.75,font=4)
		##------------------------------------------------------------------
		#rect(0.5,3.5,nind+0.5,3.5+dspace,col="aliceblue")
		rect(0.5, y.lim + 0.5,nind+0.5,y.lim + 0.5 +dspace,col="aliceblue")
		xleft=0.5+nind*0.25;xright=0.5+nind*0.75
		#ybottom=3.5+dspace*0.375
		ybottom = y.lim + 0.5 + dspace*0.375
		#ytop=3.5+dspace*0.75
		ytop = y.lim + 0.5 + dspace*0.75
		xspace=(xright-xleft)/N_color
		xblock=xleft+(0:N_color)*xspace
		for(i in 1:N_color){rect(xblock[i],ybottom,xblock[i+1],ytop,col=cols[N_color+1-i],border=F)}
		Z=range(tempz);Zset=Z[1]+c(0:6)*diff(Z)/6;Zset=round(Zset*1000)/1000;Zset=Zset[7:1]
		print(Z);print(Zset)
		U=c(xleft,xright);Uset=U[1]+c(0:6)*diff(U)/6;Uset=round(Uset*1000)/1000;
		#text(Uset,rep(3.5+dspace*0.3,length(Uset)),labels=Zset,col="blue",font=4)
		text(Uset,rep(y.lim + 0.5 + dspace*0.3,length(Uset)),labels=Zset,col="blue",font=4)

		box(col="blue")
	}

	##-----------------------figure output----------------------------------
	##=======================================================================
	#temp_QI_range=QI_range[4,];inf_Zmin=temp_QI_range[1];sup_Zmin=temp_QI_range[3];Zmin=temp_QI_range[2] 
	#inf_Zmax=temp_QI_range[3];sup_Zmax=temp_QI_range[5];Zmax=temp_QI_range[4]         
	HeatMaplist <- list(
		list(QI_Selection=list('radiobuttons', values=c(QI_name),init=QI_name[4]),
		#Sort_By_Array=list('radiobuttons', values=chip_name,init=chip_name[1]),
		Sort_By_Array=list('radiobuttons', values=unique(chip_name),init=chip_name[1]),
		Color_Scale=list(
		Min_QI=list('slider',from=inf_Zmin, to=sup_Zmin, init=Zmin, resolution=r,width=nWidth),
		Max_QI=list('slider',from=inf_Zmax, to=sup_Zmax, init=Zmax, resolution=r,width=nWidth),
		N_color=list('slider',from=2, to=40,init=18,resolution=1,width=nWidth)))
	)
	heatmapinfo=Mytkexamp_heatmap(HeatMap,HeatMaplist, vscale=1.5, hscale=2.75,plotloc='top',wait=F,Tkexample_title="QI HeatMap plot")        ##<------let wait option to be True or False
}
##==============================================================================
##==============================================================================

##----------------------------------------------------------------------------------------------------------------
if(QIfigure[2]==1) #Polygon
{
	Polygon <- function(QI_Selection=QI_name[4],Sort_By_Array=chip_name[1],Chip1_QI=1,Chip2_QI=1,Merge_QI=1){ #,Print_option=Print)
		#QIstat=QI_name[2];zmin=QImin;zmax=QImax;ncol=10;ChipSelection="Merge";PlotingSelection="Polygon";Chip1QIline=1;Chip2QIline=1;MergeQIline=1
		QI_selection=QI_Selection
		if(QI_selection==QI_name[1]){QI=MQI1;UBQI=UB_MQI1}
		if(QI_selection==QI_name[2]){QI=MQI2;UBQI=UB_MQI2}
		if(QI_selection==QI_name[3]){QI=WQI1;UBQI=UB_WQI1}
		if(QI_selection==QI_name[4]){QI=WQI2;UBQI=UB_WQI2}
		# iRank<- switch(ChipSelect,2="Merge",3="Chip2",4="Chip1")
		if(Sort_By_Array==chip_name[1]){iRank=4}
		if(Sort_By_Array==chip_name[2]){iRank=3}
		if(Sort_By_Array==chip_name[3]){iRank=2}
		z=QI[order(QI[,iRank],decreasing=T),]
		indname=z[,1]
		z=z[,-1];z=as.matrix(z);z=matrix(as.numeric(z),nrow(z),ncol(z))
		nind=nrow(z)
		y.lim = if(index.2chip){
			3
		} else 1
		y=c(1:3)
		x=c(1:nind)
		par(mai=c(0.75,0.35,0.25,0.35))
		##==========================  Polygoen figure ============================
		if(index.2chip){
			QIline=c(Chip1_QI,Chip2_QI,Merge_QI)
		} else if(index.chip == 1){
			QIline=rep(Chip1_QI, 3)
		} else QIline=rep(Chip2_QI, 3)
		PolyGonPlot_interactive(z,UBQI,indname,QIline)
		#print(QIline)
		##------------------------------------------------------------------------
	}
    # UBinit_Merge=UB_WQI2[1,1]
    # UBinit_Chip2=UB_WQI2[2,1]
    # UBinit_Chip1=UB_WQI2[3,1]
	
	UBinit_Merge=UB_WQI2[3,1] #row: Chip1, Chip2, Merge; column: 95%, 97.5%, 99%
    UBinit_Chip2=UB_WQI2[2,1]
    UBinit_Chip1=UB_WQI2[1,1]
	
	unit = 0.05
	resolution.unit = 0.05
	Polygonlist <- list(
		list(QI_Selection=list('radiobuttons', values=c(QI_name),init=QI_name[4]),
			#Sort_By_Array=list('radiobuttons', values=chip_name,init=chip_name[3]),
			#Sort_By_Array=list('radiobuttons', values=unique(chip_name),init=chip_name[3]),
			Sort_By_Array=list('radiobuttons', values=unique(chip_name),init=chip_name[1]),
			#Print=list('checkbox',stat=print_stat),
			Ref._Line=if(index.2chip){
				# list(Merge_QI=list('slider',from=min(Mintable[,3]), to=max(Maxtable[,3]),init=UBinit_Merge,  resolution=0.05,width=nWidth),
				# Chip2_QI=list('slider',from=min(Mintable[,2]), to=max(Maxtable[,2]), init=UBinit_Chip2, resolution=0.05,width=nWidth),
				# Chip1_QI=list('slider',from=min(Mintable[,1]), to=max(Maxtable[,1]), init=UBinit_Chip1, resolution=0.05,width=nWidth))
				list(Merge_QI=list('slider',from=(min(Mintable[,3]) %/% unit)*unit, to=max(Maxtable[,3], UB_WQI2[3, 3]),init=UBinit_Merge,  resolution=resolution.unit,width=nWidth),
				Chip2_QI=list('slider',from=(min(Mintable[,2]) %/% unit)*unit, to=max(Maxtable[,2], UB_WQI2[2, 3]), init=UBinit_Chip2, resolution=resolution.unit,width=nWidth),
				Chip1_QI=list('slider',from=(min(Mintable[,1]) %/% unit)*unit, to=max(Maxtable[,1], UB_WQI2[1, 3]), init=UBinit_Chip1, resolution=resolution.unit,width=nWidth))
				#Chip1_QI=list('slider',from=0.376, to=max(Maxtable[,1], UB_WQI2[1, 3]), init=UBinit_Chip1, resolution=0.05,width=nWidth))
				#Maxtable
				##       chip1 chip2 merge
				##  MQI1
				##  MQI2
				##  WQI1
				##  WQI2
			} else if(index.chip == 1){
				#list(Chip1_QI=list('slider',from=min(Mintable[,1]), to=max(Maxtable[,1], UB_WQI2[1, 3]), init=UBinit_Chip1, resolution=0.05,width=nWidth))
				list(Chip1_QI=list('slider',from=(min(Mintable[,1]) %/% unit)*unit, to=max(Maxtable[,1], UB_WQI2[1, 3]), init=UBinit_Chip1, resolution=resolution.unit,width=nWidth))
			} else {
				#list(Chip2_QI=list('slider',from=min(Mintable[,2]), to=max(Maxtable[,2], UB_WQI2[2, 3]), init=UBinit_Chip2, resolution=0.05,width=nWidth))
				list(Chip2_QI=list('slider',from=(min(Mintable[,2]) %/% unit)*unit, to=max(Maxtable[,2], UB_WQI2[2, 3]), init=UBinit_Chip2, resolution=resolution.unit,width=nWidth))
			}
		)
	)

  #polygoninfo=tkexamp(Polygon,Polygonlist, vscale=1.5, hscale=2.75,plotloc='top',wait=F)
  polygoninfo=Mytkexamp_polygon(Polygon,Polygonlist, vscale=1.5, hscale=2.75,plotloc='top',wait=F,Tkexample_title="QI Polygon plot")
}

