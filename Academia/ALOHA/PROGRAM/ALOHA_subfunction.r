library(gtools)
sort.append<-function(x,base)
{
	library(gtools)
	sort.base<-mixedsort(paste(x[,base],1:length(x[,base])))
	tmp<-strsplit(sort.base," ")
	tmp1<-matrix(as.integer(unlist(tmp)),2,)
	y<-x[tmp1[2,],]
	return(y)
}
AI = function(AF, AFref, alpha, n)
{
	alphaset = 1 - c(alpha/(3*n), alpha/(6*n), alpha/(3*n)) 
	Z = qnorm(alphaset)
	sdset = AFref[, c(5, 8, 11)]  
	sdset = sapply(1:3, function(i) Z[i]*sdset[, i])
	meanset = AFref[, c(4, 7, 10)] 
	is.na(meanset[is.na(sdset)]) = T
	mean.diff = abs(meanset - AF[, 4]) 
	diff_mean_sd = mean.diff - sdset
	index.CI = diff_mean_sd > 0 
	diff_mean_sd[index.CI] = 1 
	diff_mean_sd[!index.CI] = 0
	index.AI = as.integer(apply(diff_mean_sd, 1, function(x) all(x == 1, na.rm = T)))
	homo_index = meanset[, c(1, 3)] - AF[, 4] 
	homo_index[[2]] = -homo_index[[2]]
	index.AI.2 = homo_index > 0
	index.AI.2[is.na(index.AI.2)] = 1
	AF.new = cbind(AF, AIindex = index.AI*index.AI.2[,1]*index.AI.2[,2])
	return(AF.new)
}
LOH = function(AF, AFref, alpha, n)
{
	alphaset = 1 - c(alpha/(2*n))
	Z = qnorm(alphaset)
	sdset = AFref[, 8]
	sdset = Z*sdset
	meanset = AFref[, 7]
	is.na(meanset[is.na(sdset)]) = T
	mean.diff = abs(meanset - AF[, 4])
	diff_mean_sd = mean.diff - sdset
	index.CI = diff_mean_sd > 0
	diff_mean_sd[index.CI] = 1 
	diff_mean_sd[!index.CI] = 0
	diff_mean_sd[is.na(index.CI)] = 0
	AF = cbind(AF, LOHindex = diff_mean_sd)
	return(AF)
}
WindowSlidingProp = function(a, w)
{
	nsnp = length(a)
	d = w %/% 2
	A = sapply(1:(2*d + 1), function(i){
		index1 = d-i+1 
		index2 = i-d-1 
		if(index1 < 0) index1 = 0
		if(index2 < 0) index2 = 0
		c(rep(NA, index1), a[(1+index2):(nsnp-index1)], rep(NA, index2))
	})
	ave_index = apply(A, 1, function(x) mean(x, na.rm = T))
	return(ave_index)
}
Double_pp_checking = function(Dpp){
	pp=as.numeric(Dpp[,2]);Dpp[,2]=pp
	temp_pp=table(Dpp[,2]);temp_pp=temp_pp[temp_pp>1];temp_pp=names(temp_pp)
	npp=length(temp_pp)
	if(npp==0){whole_double_pp=Dpp}
	if(npp>0){
		single_SNP=Dpp[Dpp[,2]!=temp_pp,]
		whole_double_pp=single_SNP;
		for(k in 1:npp){
			double_pp=temp_pp[k]
			temp_chr=Dpp[Dpp[,2]==double_pp,]
			num_double_pp=nrow(temp_chr)
			add_position=(c(1:num_double_pp)-1)*1
			pp=as.numeric(temp_chr[,2])
			pp=pp+add_position
			temp_chr[,2]=pp
			whole_double_pp=rbind(whole_double_pp,temp_chr)
		}
	}
	whole_double_pp=sort.append(whole_double_pp,2)
	return(whole_double_pp)
}
AFref_checking=function(AFpath,NR_AFref){
	cat("Allele frequecy checking process is preparing.\n")
	chr_AFref_file=list.files(AFpath,full.names=T)
	chr_AFref_name=list.files(AFpath,full.names=F)
	nchr=length(chr_AFref_file)
	if(nchr!=24){stop("ALOHA error message: The number of AF reference is incorrect, it should be 24 (22+X+Y) chromosomes.")}
	a=c(1:9);a=paste(0,a,sep="");b=c(a,10:24);b=paste("Chr_",b,".txt",sep="")
	if(sum(chr_AFref_name==b)!=24)
	{stop("ALOHA error message: The file name of AF reference is incorrect.Please check the file name.")}
	whole_empty_chr=NULL;
	for(j in 1:nchr){
		chr=read.table(chr_AFref_file[j],head=T,sep="\t")
		nCol=ncol(chr)
		if(nCol!=12)
		{stop("ALOHA error message: The number of column in ",chr_AFref_name[j]," is incorrect.Please check the number of column.")}
		Coltag=c("Probe_set","Phy_position","Chiptype","num_AA","mean_AA","std_AA","num_AB","mean_AB","std_AB","num_BB","mean_BB","std_BB")
		Colchr=colnames(chr)
		if(sum(Colchr==Coltag)!=12)
		{stop("ALOHA error message: The column name of ",chr_AFref_name[j]," is incorrect.Please check the colname name.")}
		empty_chr=chr[,c(1,4,7,10)]
		empty_chr=cbind(empty_chr,apply(empty_chr[,-1],1,sum))
		empty_chr_probe_set=empty_chr[empty_chr[,5]==0,1]
		if(length(empty_chr_probe_set)>0){
			empty_chr=cbind(j,chr[chr[,1]%in% empty_chr_probe_set,])
			whole_empty_chr=rbind(whole_empty_chr,empty_chr)
			chr=chr[!chr[,1]%in% empty_chr_probe_set,]
		}
		NR_AFref_chr=paste(NR_AFref,"/",chr_AFref_name[j],sep="")
		write.table(chr,NR_AFref_chr,col.names=T,row.names=F,quote=F,sep="\t")
	}
	colnames(whole_empty_chr)[1]="Chr"
	whole_empty_chr=whole_empty_chr[,c(1:4)]
	cat("Allele frequecy checking process is finished.\n")  
	return(whole_empty_chr)
}
PPcum = function(A)
{
	whole_d = whole_PP = NULL; d = 0
	for(ichr in 1:23)
	{
		pp = A[A[,1] == ichr, 2]
		cumPP = pp + d
		whole_PP = c(whole_PP, cumPP)
		d = max(whole_PP)
		whole_d = c(whole_d, d)
	}
	A[, 2] = whole_PP
	A = list(A, whole_d)
	return(A)
}
AF_figure = function(AF_table)
{
	pp = AF_table[,1]; geno = AF_table[,2]; AF = AF_table[,3];
	geno = as.character(geno)
	geno[geno == "AA"] = 3; geno[geno == "AB"] = 2
	geno[geno == "BB"] = 1; geno[geno == "NoCall"] = 0
	geno = as.numeric(geno)
	AF_table[, 2] = geno
	Yrange = c(0, 1)
	Xrange = range(pp)
	plot(1, 1, type = "n", ylab = "", xlab = "", xlim = Xrange, ylim = Yrange, axes = F)
	points(pp, AF, pch = ".", col = "blue")
	AI = AF_table[AF_table[, 4] == 1,]; pp = AI[, 1]; AF = AI[, 3]
	points(pp, AF, pch = ".", col = "red", cex = 1.75)
	Ylabel = c(Yrange[1], mean(Yrange), Yrange[2])
	Ylabel = round(Ylabel, 2)
	return(Ylabel)
}
AF_hist_clustering = function(B,rspan)
{
    af = B[,2]
    nsnp = length(af)
    nbreak = round(nsnp)*0.25
    dfset = 6
    a = hist(af,breaks = nbreak,plot = F)
    a = sm.spline(a$mids, a$counts, df = dfset)
    a_table = cbind(a[2]$x,a[3]$ysmth)
    nbreak = nrow(a_table)
    nspan = round(nbreak*rspan)    
    dummy_tail = rep(c(1,0),nspan)
    dummy_tail = matrix(dummy_tail,nspan,2,byrow = T)
    a_table = rbind(a_table,dummy_tail)
    maxmin <- msExtrema(a_table,span = nspan)
    local_max = maxmin$index.max*1
    a_table = a_table*local_max[,2]
    a_table = a_table[1:nbreak,]
    a_critical = a_table[a_table[,1]!=0,1]
    ncritical = length(a_critical)
    a_critical = matrix(rep(a_critical,nsnp),ncritical,nsnp)
    a_critical = t(a_critical)
    b = cbind(B,a_critical)
    pp = b[,1];b = b[,-1]
    for(k in 1:ncritical){b[,k+1]=b[,1]-b[,k+1]}
    u = b[,-1]
    u = abs(u)
    minU = cbind(u,apply(u,1,min))
    for(k in 1:ncritical){minU[,k]=minU[,ncritical+1]-minU[,k]}
    u = minU[,-(ncritical+1)]
    u[u==0]=1;u[u!=1]=0
    u = as.matrix(u)
    u = matrix(u,nsnp,ncritical)
    for(k in 1:ncritical){u[,k]=u[,k]*k}
    u = apply(u,1,sum)
    B = cbind(B,u)
    return(B)
}
figure_public_setting = function(X, Y, label_name)
{
    cexsize = 1.5
    labelsize = 1.75
    titlesize = 2
    box() 
    X = c(0, X)
    abline(v = X, col = "black", lty = 2)
    Xchr = X + c(diff(X)/2, 0) 
    chrset = seq(1, 23, 2)
    Xchr = cbind(chrset, Xchr[chrset])
    axis(1, at = Xchr[, 2], labels = Xchr[, 1], line = -0.5, cex.axis = cexsize, tick = F)
    axis(1, at = max(X), labels = "Chromosome", line = -0.5, cex.axis = labelsize, tick = F, col.axis = "blue", hadj = 0.2)
    axis(3, at = max(X)/2, labels = label_name[1], line = -0.5, cex.axis = titlesize, tick = F, col.axis = "darkgreen")
    Ylabel = round(Y, 2)
    axis(2, at = Y, labels = Ylabel, line = -4, cex.axis = cexsize, tick = F, las = 2)
    axis(2, at = Y[2], labels = label_name[2], line = 0.25, cex.axis = labelsize, tick = F)
}
AILOH_figure = function(A, def_set){
	pp=A[,1]
	case_index=A[,2]
	control_index=A[,4]
	S=c(A[,2],A[,4])
	yregion=range(S,na.rm=T)
	rb=0.025
	wb=rb*diff(yregion)
	xregion=range(pp,na.rm=T)
	plot(1,1,type="n",xlim=xregion,ylim=yregion,axes=F,xlab="",ylab="")
	points(pp,case_index,type="l",col="pink")
	points(pp,control_index,type="l",col="skyblue")
	for(ichr in 1:23){
		chr=A[A[,5]==ichr,];chr_def=def_set[ichr]
		temp_pp=chr[,1]
		temp_case_index=chr[,2];temp_control_index=chr[,4];
		case_smth=NULL;
		control_smth=NULL;
		if(!is.na(sum(temp_case_index))){
			a<-sm.spline(temp_pp, temp_case_index, df=chr_def) 
			lines(a, lty=1, col = "red",lwd=1)
			case_smth=unlist(a[3]); 
		}
		if(!is.na(sum(temp_control_index))){
			b<-sm.spline(temp_pp, temp_control_index, df=chr_def)
			lines(b, lty=1, col = "blue",lwd=1)
			control_smth=unlist(b[3]);
		}
		diff_smth=control_smth-case_smth
		if(length(diff_smth)==0){
			smoth_pp=temp_pp
			smth=cbind(smoth_pp,diff_smth,max(yregion))
			Xright=smth[,1];Xleft=smth[,1]
			Yupper=smth[,2];Ylower=smth[,2]-wb
			segments(Xright,Ylower,Xleft,Yupper,col="green")
		}
		if(length(diff_smth)>0){
			smoth_pp=unlist(a[2])
			smth=cbind(smoth_pp,diff_smth,max(yregion))
			positive_smth=smth[smth[,2]>0,]      
			negative_smth=smth[smth[,2]<0,]      
			positive_smth=as.matrix(positive_smth)
			negative_smth=as.matrix(negative_smth)
			if(ncol(positive_smth)==3){
				Xright=positive_smth[,1];Xleft=positive_smth[,1]
				Yupper=positive_smth[,3];Ylower=positive_smth[,3]-wb
				segments(Xright,Ylower,Xleft,Yupper,col="blue")
			}
			if(ncol(negative_smth)==3){
				Xright=negative_smth[,1];Xleft=negative_smth[,1]
				Yupper=negative_smth[,3];Ylower=negative_smth[,3]-wb
				segments(Xright,Ylower,Xleft,Yupper,col="red")
			}
		}
	}
	yregion=c(yregion[1],mean(yregion),yregion[2])
	return(yregion)
}
BPMG.chr <- function(Index_list, N_list, IndivName_list, Col_list, chr, len_list, Outpath){
	chrname <- ifelse(chr<10, paste("0", chr, sep=""), chr)
	C <- read.table(paste(Outpath, "\\temp\\CforBP_chr", chrname,".txt",sep=""), header = F, sep="\t",
				 colClasses = rep("numeric", 2))
	G <- read.table(paste(Outpath, "\\temp\\GforBP_chr", chrname,".txt",sep=""), header = F, sep="\t",
				 colClasses = rep("numeric", 2))
	Gxlim <- c(floor(quantile(G[,1], probs = 0.01, na.rm = T)*100)/100, ceiling(quantile(G[,1], probs = 0.99, na.rm = T)*100)/100)
	Gylim<- c(floor(quantile(G[,2], probs = 0.01, na.rm = T)*100)/100, ceiling(quantile(G[,2], probs = 0.99, na.rm = T)*100)/100)
	Cxlim<- c(floor(min(C[,1], na.rm = T)*1.001*100)/100, ceiling(max(C[,1], na.rm = T)*1.001*100)/100)
	Cylim<- c(floor(min(C[,2], na.rm = T)*1.001*100)/100, ceiling(max(C[,2], na.rm = T)*1.001*100)/100)
	Gxtick <- unique(round(seq(from = Gxlim[1], to = Gxlim[2], len = 5), 2))
	Gytick <- unique(round(seq(from = Gylim[1], to = Gylim[2], len = 5), 2))
	Cxtick <- unique(round(seq(from = Cxlim[1], to = Cxlim[2], len = 5), 2))
	Cytick <- unique(round(seq(from = Cylim[1], to = Cylim[2], len = 5), 2))
	Cx.mid <- mean(Cxlim);				Cy.mid <- mean(Cylim);
	Cx.Q1 <- (Cxlim[1]*3+Cx.mid)/4;		Cy.Q1 <- (Cylim[1]*3+Cy.mid)/4;
	Cx.Q3 <- (Cxlim[2]*3+Cx.mid)/4;		Cy.Q3 <- (Cylim[2]*3+Cy.mid)/4;
	tiff(filename  =  paste(Outpath, "\\AFBP_chr", chrname,".tiff",sep=""), width  =  1024, height  =  768)
	par(mar = c(5-1, 4+8, 4+0.5, 2+3))
	plot(G[,1],G[,2],xlim = Gxlim,ylim = Gylim,pch="~",col="green",axes = F,xlab="",ylab="",main="",xaxt="n",yaxt="n")
	axis(side = 1, cex.axis =  1, at = Gxtick, labels = Gxtick, col.axis="gray60", col.tick="gray60")
	axis(side = 2, cex.axis =  1, at = Gytick, labels = Gytick, col.axis="gray60", col.tick="gray60")
	mtext("Marker (component 1)", line = 1, side = 1, col="black", padj = 1, cex = 1.5)
	mtext("Marker (component 2)", line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
	par(xpd = NA)
	text(Gxlim[1]*1.4, Gylim[2]*1.25, paste("Chr ", chrname, sep=""), font = 2, col="blue", cex = 3)
	legend(x="bottomleft", legend = len_list, col = Col_list, pch = 16, bty="o", box.col="gray",
		cex = 1.5, inset  =  c(-0.2, 0), title="AF biplot")
	par(new = T)
	plot(C[,1],C[,2],type="n",xlim = Cxlim,ylim = Cylim,axes = T,xlab="",ylab="",main="",xaxt="n",yaxt="n")
	axis(side = 3, cex.axis =  1, at = Cxtick, labels = Cxtick, col.axis="gray60", col.tick="gray60")
	axis(side = 4, cex.axis =  1, at = Cytick, labels = Cytick, col.axis="gray60", col.tick="gray60")
	mtext("Individual (component 1)", line = 3, side = 3, col="black", padj = 1, cex = 1.5)
	mtext("Individual (component 2)", line = 1.5, side = 4, col="black", padj = 1, cex = 1.5)
	cumN_list <- c(1, cumsum(N_list))
	for (i in 1:length(N_list)){
		{if(i==1)	posi <- cumN_list[i]: cumN_list[i+1]
		else{		posi <- (cumN_list[i]+1): cumN_list[i+1]}
		}
		xx <- C[posi,1]
		yy <- C[posi,2]
		nowname <- IndivName_list[posi]
		arrows(Cx.mid, Cy.mid, C[posi,1], C[posi,2], length = 0, col = Col_list[i])
		cri <- xx >= Cx.Q1 & xx <= Cx.Q3
		if(length(which(cri)==T)!=0)	text(xx[cri], yy[cri], nowname[cri], col = Col_list[i], cex = 1);
		cri <- xx > Cx.Q3
		if(length(which(cri)==T)!=0)	text(xx[cri], yy[cri], nowname[cri], col = Col_list[i], cex = 1, adj = 1);
		cri <- xx < Cx.Q1
		if(length(which(cri)==T)!=0)	text(xx[cri], yy[cri], nowname[cri], col = Col_list[i], cex = 1, adj = 0);
	}
	dev.off()
} 
BPMG.WG <- function(Index_list, N_list, IndivName_list, Col_list, len_list, Outpath){
	tiff(paste(Outpath, "\\GAFBP.tiff", sep=""), width = 1024, height = 768)
	nf <- layout(matrix(1:25,5, 5, byrow = F), TRUE)
	for(chr in 1:23){
		chrname <- ifelse(chr<10, paste("0", chr, sep=""), chr)
		C <- read.table(paste(Outpath, "\\temp\\CforBP_chr", chrname,".txt",sep=""), header = F, sep="\t",
					 colClasses = rep("numeric", 2))
		G <- read.table(paste(Outpath, "\\temp\\GforBP_chr", chrname,".txt",sep=""), header = F, sep="\t",
					 colClasses = rep("numeric", 2))
		Gxlim <- c(floor(quantile(G[,1], probs = 0.01, na.rm = T)*100)/100, ceiling(quantile(G[,1], probs = 0.99, na.rm = T)*100)/100)
		Gylim<- c(floor(quantile(G[,2], probs = 0.01, na.rm = T)*100)/100, ceiling(quantile(G[,2], probs = 0.99, na.rm = T)*100)/100)
		Cxlim<- c(floor(min(C[,1], na.rm = T)*1.001*100)/100, ceiling(max(C[,1], na.rm = T)*1.001*100)/100)
		Cylim<- c(floor(min(C[,2], na.rm = T)*1.001*100)/100, ceiling(max(C[,2], na.rm = T)*1.001*100)/100)
		Gxtick <- unique(round(seq(from = Gxlim[1], to = Gxlim[2], len = 5), 2))
		Gytick <- unique(round(seq(from = Gylim[1], to = Gylim[2], len = 5), 2))
		Cxtick <- unique(round(seq(from = Cxlim[1], to = Cxlim[2], len = 5), 2))
		Cytick <- unique(round(seq(from = Cylim[1], to = Cylim[2], len = 5), 2))
		Cx.mid <- mean(Cxlim);				Cy.mid <- mean(Cylim);
		Cx.Q1 <- (Cxlim[1]*3+Cx.mid)/4;		Cy.Q1 <- (Cylim[1]*3+Cy.mid)/4;
		Cx.Q3 <- (Cxlim[2]*3+Cx.mid)/4;		Cy.Q3 <- (Cylim[2]*3+Cy.mid)/4;
		par(mar = c(5, 4+0.5, 2+2, 2+1))
		plot(G[,1],G[,2],xlim = Gxlim,ylim = Gylim,pch=".",col="green",axes = F,xlab="",ylab="",main="",xaxt="n",yaxt="n")
		axis(side = 1, cex.axis =  1, at = Gxtick, labels = Gxtick, col.axis="gray60", col.tick="gray60")
		axis(side = 2, cex.axis =  1, at = Gytick, labels = Gytick, col.axis="gray60", col.tick="gray60")
		mtext("Marker_C1", line = 1, side = 1, col="black", padj = 1, cex = 0.9)
		mtext("Marker_C2", line = 3, side = 2, col="black", padj = 1, cex = 0.9)
		par(xpd = NA)
		text(Gxlim[1]*1.3, Gylim[2]*2.1, paste("Chr ", chrname, sep=""), font = 2, col="blue", cex = 2.3)
		par(new = T)
		plot(C[,1],C[,2],type="n",xlim = Cxlim,ylim = Cylim,axes = T,xlab="",ylab="",main="",xaxt="n",yaxt="n")
		axis(side = 3, cex.axis =  1, at = Cxtick, labels = Cxtick, col.axis="gray60", col.tick="gray60")
		axis(side = 4, cex.axis =  1, at = Cytick, labels = Cytick, col.axis="gray60", col.tick="gray60")
		mtext("Ind_C1", line = 3, side = 3, col="black", padj = 1, cex = 0.9)
		mtext("Ind_C2", line = 1, side = 4, col="black", padj = 1, cex = 0.9)
		cumN_list <- c(1, cumsum(N_list))
		for (i in 1:length(N_list)){
			{if(i==1)	posi <- cumN_list[i]: cumN_list[i+1]
			else{		posi <- (cumN_list[i]+1): cumN_list[i+1]}
			}
			xx <- C[posi,1]
			yy <- C[posi,2]
			arrows(Cx.mid, Cy.mid, C[posi,1], C[posi,2], length = 0, col = Col_list[i])
			cri <- xx >= Cx.Q1 & xx <= Cx.Q3
			if(length(which(cri)==T)!=0)	text(xx[cri], yy[cri], i, col = Col_list[i], cex = 0.7);
			cri <- xx > Cx.Q3
			if(length(which(cri)==T)!=0)	text(xx[cri], yy[cri], i, col = Col_list[i], cex = 0.7, adj = 1);
			cri <- xx < Cx.Q1
			if(length(which(cri)==T)!=0)	text(xx[cri], yy[cri], i, col = Col_list[i], cex = 0.7, adj = 0);
		}
	} 
	plot(1:2,type="n",xlab="",ylab="",main="",xaxt="n",yaxt="n", bty="n")
	mtext("AF biplot", side = 1, padj = 0, adj = 0, line = 1, font = 2, cex = 2)
	plot(1:length(AllIndivNum)*1.5,type="n",xlab="",ylab="",main="",xaxt="n",yaxt="n", bty="n")
	legend(x = 1, y = length(AllIndivNum)*1.5, legend = len_list, col = Col_list, pch = 16, bty="n", cex = 2)
	dev.off()
} 
Biplot_NR_manage = function(path, A, ind_name, chip_info)
{
	marker = A[[1]]
	ind = A[[2]]
	singular_value = A[[3]]
	DBline = c("===========================================================")
	SGline = c("-----------------------------------------------------------")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = F, sep = "\t")
	write.table(DBline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("===  1. Proportion of explained variation (PEV)  ===", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(DBline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(paste("Variable", "Singular_value", "PEV_(%)", "Cumulative_PEV_(%)", sep = "\t"), file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(SGline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	per_sv = singular_value/sum(singular_value)
	n = length(per_sv)
	sv_list = c(1:n)
	whole_sv = cbind(paste("SV", sv_list, sep = ""), "  ", 
		format(sprintf("% .4f", singular_value), justify = "right"),
		format(sprintf("% .4f", per_sv), justify = "right"),
		format(sprintf("% .4f", cumsum(per_sv)), justify = "right")
	)
	write.table(whole_sv, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	DBline = c("=================================")
	SGline = c("------------------------------------------")
	write.table(DBline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("===  2. Individual component  ===", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(DBline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(paste("Name", "Ind_C1", "Ind_C2", sep = "\t"), file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(SGline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	whole_IndComp = cbind(ind_name, 
							format(sprintf("% 2.4f", ind[, 1]), justify = "right"), 
							format(sprintf("% 2.4f", ind[, 2]), justify = "right"))
	write.table(whole_IndComp, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	DBline = c("=================================")
	SGline = c("-------------------------------------------------------------------")
	write.table(DBline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("===  3. Marker component  ===", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(DBline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	chip_info[, 2] = format(sprintf("% 3.6f", chip_info[, 2]/1000000), justify = "right")
	whole_marker = cbind(chip_info[, 1:2], format(sprintf("% 1.4f", marker[, 1]), justify = "right"),
							format(sprintf("% 1.4f", marker[, 2]), justify = "right"))
	write.table(paste("Probe_set", "Phy_posi_(Mb)", "Mar_C1", "Mar_C2", sep = "\t"), file = path, 
					col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(SGline, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table(whole_marker, file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
	write.table("", file = path, col.names = F, row.names = F, quote = F, append = T, sep = "\t")
}
data_format_checking = function(inputfile, temp_group_choice, Logfile, SampleList)
{
	checking_message = NULL
	popuname = list.files(inputfile)
	error_message = NULL
	if(temp_group_choice == 1){
		if(length(popuname) < 1) stop("There are not enough population for calculating")
		if(any(c("Case", "Control") %in% popuname)) stop("The input path have 'Case' and 'Control' folder", sep = "")
	}
	if(temp_group_choice == 2){
		if(!all(c("Case", "Control") %in% popuname)) stop("ALOHA error message: Directory Case and/or Control is not found.")
	}
	popufile = list.files(inputfile, full.names = T) 
	popuname = list.files(inputfile, full.names = F)
	for(i in 1:length(popufile)){ 
		checking_message = paste("   \n" ,sep=""); cat(checking_message, file = Logfile, append = T, sep = "")
		checking_message = paste("   ~The data format checking of ", popuname[i], " is preparing.\n", sep = "")
		cat(checking_message, file = Logfile, append = T, sep = "")
		indfile = list.files(popufile[i], full.names = T) 
		indname = list.files(popufile[i], full.names = F) 
		nind = length(indfile)
		if(nind == 0)
		{
			error_message = paste("There are no individual files within directory ", popuname[i], sep = "")
			stop(error_message)
		}
		SLtable_file = paste(SampleList, "/", popuname[i], ".txt", sep = "")
		digital_Indname = max(nchar(indname))
		Colnames = c("Obs", "Ind_ID   ", "Gender", "CR")
		Colnames = matrix(Colnames, 1, 4)
		write.table(Colnames, SLtable_file, row.names = F, quote = F, col.names = F, append = F, sep = "\t")
		Sample_list_table = NULL
		Sample_list_table = cbind(c(1:nind), indname)
		whole_gender = NULL
		whole_CR = NULL
		for(j in 1:nind){ 
			checking_message = paste("         ~The data format checking of ", indname[j], " is preparing.\n" ,sep="")
			cat(checking_message, file = Logfile, append = T, sep = "")
			chrfile = list.files(indfile[j], full.names = T)
			chrname = list.files(indfile[i], full.names = F)
			nchr = length(chrfile)
			if(nchr != 23){
				error_message = paste("ALOHA error message: The individual's chromosomal information is uncompleted", sep = "")
				stop(error_message)
			} else {
				whole_geno = NULL
				for(ichr in 1:nchr){ 
					if(file.exists(chrfile[ichr]) == F){
						error_message = paste("ALOHA error message: The file of ", chrname[ichr], " of ", indname[i], "is not exist.", sep = "")
						stop(error_message)
					} else { 
						chr = read.table(chrfile[ichr], head = T, sep = "\t")
						ColName = colnames(chr)
						chr = as.matrix(chr)
						if(dim(chr)[1] == 0){
							error_message = paste("ALOHA error message: No information in chromosome ", ichr, ". ", sep = "")
							stop(error_message)
						}
						if(dim(chr)[2] != 6){
							error_message = paste("ALOHA error message: The column's infomation are incorrect. ", sep = "")
							stop(error_message)
						}
						if(ColName[1] != "Probe_set"){
							error_message = paste("ALOHA error message: The column of 'Probe_set' are incorrect. ", sep = "")
							stop(error_message)
						}
						if(ColName[2] != "Chr"){
							error_message = paste("ALOHA error message: The column of 'Chr' are incorrect. ", sep = "")
							stop(error_message)
						}
						if(ColName[3] != "Phy_position"){
							error_message = paste("ALOHA error message: The column of 'Phy_position' are incorrect. ", sep = "")
							stop(error_message)
						}
						if(ColName[4] != "Genotype"){
							error_message = paste("ALOHA error message: The column of 'Genotype' are incorrect. ", sep = "")
							stop(error_message)
						}
						if(ColName[5] != "Chiptype"){
							error_message = paste("ALOHA error message: The column of 'Chiptype' are incorrect. ", sep = "")
							stop(error_message)
						}
						if(ColName[6] != "AF"){
							error_message = paste("ALOHA error message: The column of 'AF' are incorrect. ", sep = "")
							stop(error_message)
						}
					}
					whole_geno = c(whole_geno, chr[, 4])
					if(ichr == nchr){sex_geno = chr[, 4]}
				} 
				genotable = table(whole_geno)
				genotable = data.frame(t(as.matrix(genotable)))
				CR = sum(genotable$AA + genotable$AB + genotable$BB)/sum(genotable)
				sexgenotable = table(sex_geno)
				sexgenotable = data.frame(t(as.matrix(sexgenotable)))
				heterRatio = sum(sexgenotable$AB)/sum(sexgenotable)
				if(heterRatio < 0.1){
					gender = "M"
				} else gender = "F"
				whole_CR = c(whole_CR, CR)
				whole_gender = c(whole_gender, gender)
			} 
			checking_message = paste("         ~The data format checking of ", indname[j], " is finished.\n" , sep = "")
			cat(checking_message, file = Logfile, append = T, sep = "")
		} 
		Sample_list_table = cbind(Sample_list_table, whole_gender, whole_CR)
		write.table(Sample_list_table, SLtable_file, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
		checking_message = paste("   ~The data format checking of ", popuname[i], " is finished.\n" , sep = "")
		cat(checking_message, file = Logfile, append = T, sep = "")
	} 
	return_message = paste("The input data is suitable for ALOHA calculating.\n", sep = "")
	cat(return_message)
}
All_sample_AILOH_figure=function(A,PP_bound,title_name)
{
  library(pspline) 
   gender=A[1,];gender=as.character(gender[-c(1,2)])
   A=A[-1,]      
   def.par <- par(no.readonly = TRUE) 
   nchr=length(PP_bound)
   pp_info=A[,c(1,2)]
   temp_pp=as.numeric(pp_info[,2])/1000000
   PP_bound = PP_bound/1000000
   ind=A[,-c(1,2)]
   ind=as.matrix(ind);
   Xwidth=2 
   nind=ncol(ind)
   Xlim=c(0,max(temp_pp,na.rm=T))
   Ylim=c(0,nind)+0.5   
   kk=10 
   nf <- layout(matrix(c(1,3,2,4),2,2,byrow=TRUE),c(kk,1),c(kk,1))
   layout.show(nf)
   par(mar=c(1,3,4,1))
   plot(1,1,type="n",xlim=Xlim,ylim=Ylim,axes=F,ylab="",xlab="")
   for(i in 1:nind)   
   {
      index=as.numeric(ind[,i])
      yupper=i+0.5;ylower=i-0.5; 
      index=cbind(temp_pp,index,(nind+2)-yupper,nind-ylower)         
      NA_index=index[index[,2]==2,]
      positive_index=index[index[,2]==1,]
      negative_index=index[index[,2]==0,]
      NA_index=as.matrix(NA_index);
      if(dim(NA_index)[2]!=4){NA_index=t(NA_index)}
      positive_index=as.matrix(positive_index);
      if(dim(positive_index)[2]!=4){positive_index=t(positive_index)}
      negative_index=as.matrix(negative_index);
      if(dim(negative_index)[2]!=4){negative_index=t(negative_index)}            
      if(nrow(NA_index)>0)
      {
        xleft=NA_index[,1];xright=NA_index[,1]   
        yupper=NA_index[,3];ylower=NA_index[,4]             
        segments(xleft,ylower,xright,yupper,col="green")
      }      
      if(nrow(negative_index)>0)
      {
        xleft=negative_index[,1];xright=negative_index[,1]   
        yupper=negative_index[,3];ylower=negative_index[,4]             
        segments(xleft,ylower,xright,yupper,col="blue")
      }
      if(nrow(positive_index)>0)
      {
        xleft=positive_index[,1];xright=positive_index[,1]   
        yupper=positive_index[,3];ylower=positive_index[,4]             
        segments(xleft,ylower,xright,yupper,col="red")
      }
   }
   Yind=c(1:nind)+0.5;Ylab=c(1:nind)
   Yind=(nind+2)-Yind 
   colset="skyblue"
   if(length(Yind)>30)
   {
      nblock=10 
      n=length(Yind)
      Dind=floor(n/10)
      Indlist=c(1:nblock)*Dind
      Yind=Yind[Indlist]
      Ylab=Ylab[Indlist]
   } 
   abline(h=Yind-1,col=colset,lty=2,lwd=1.25)
   abline(v=c(0,PP_bound),col=colset,lty=2,lwd=1.25**Xwidth)
   cexsize=2
   labelsize=2.5
   titlesize=3
   xloc=c(diff(PP_bound),0)/2
   xloc=xloc+PP_bound
   xloc=c(PP_bound[1]/2,xloc[-nchr])
   chr_list=cbind(c(1:nchr),xloc)   
   chr_list=chr_list[c(1:ceiling(nchr/2))*2-1,]
   YlabcexLine=-5
   axis(3,at=chr_list[,2],labels=chr_list[,1],line=YlabcexLine+1,cex.axis=cexsize,tick=F,font.axis=4)
   axis(3,at=0,labels="Chromosome",line=YlabcexLine+2.5,font.axis=2,cex.axis=labelsize*0.75,tick=F,col.axis="blue",hadj=0.2)     
   axis(2,at=Yind-0.5,labels=Ylab,line=YlabcexLine,cex.axis=cexsize,tick=F,las=2,font.axis=4)
   axis(2,at=mean(Ylim),labels="Obs",line=-0.5,cex.axis=labelsize,tick=F,col.axis="blue",font.axis=2)
   axis(3,at=max(temp_pp)/2,labels=title_name,line=0.25,cex.axis=titlesize*1.2,tick=F,col.axis="darkgreen",font.axis=4)   
   box(col="black") 
   par(mar=c(1,3,0,1))
   ind_SNP=ind
   ind_SNP=matrix(as.numeric(ind_SNP),nrow(ind_SNP),ncol(ind_SNP))
   ind_SNP[ind_SNP==2]=NA
   SNP_nind=apply(ind_SNP,1,function(x){sum(x,na.rm=T)})
   ind_SNP[is.na(ind_SNP)==F]=1
   ind_SNP[is.na(ind_SNP)==T]=0 
   SNP_all=apply(ind_SNP,1,function(x){sum(x,na.rm=T)})        
   prop_snp=SNP_nind/SNP_all
   prop_snp[is.na(prop_snp)==T]=NA   
   prop_snp=as.matrix(prop_snp);prop_snp=as.numeric(prop_snp) 
   temp_Ylim=c(0,ceiling(max(prop_snp,na.rm=T)*100)/100)
   plot(1,1,type="n",xlim=Xlim,ylim=temp_Ylim,ylab="",xlab="",axes=F)
   points(temp_pp,prop_snp,type="l",col="purple",lwd=2)
   abline(v=c(0,PP_bound),col=colset,lty=2,lwd=1.25*Xwidth)
   n=4;d=temp_Ylim[2]/n 
   Ylab=c((0:4)*d)
   abline(h=Ylab,col=colset,lty=2,lwd=1.25)
   axis(2,at=Ylab,labels=Ylab*100,line=YlabcexLine,cex.axis=cexsize,tick=F,las=2,font.axis=4)    
   axis(2,at=median(Ylab),labels="Aberration (%)",line=YlabcexLine+5,col.axis="blue",cex.axis=labelsize*0.65,tick=F,font.axis=2)
   box(col="black") 
   par(mar=c(1,0,4,1))
   ind_SNP=ind
   ind_SNP=matrix(as.numeric(ind_SNP),nrow(ind_SNP),ncol(ind_SNP))
   ind_SNP[ind_SNP==2]=NA
   Ind_nsnp=apply(ind_SNP,2,function(x){sum(x,na.rm=T)})
   ind_SNP[is.na(ind_SNP)==F]=1
   ind_SNP[is.na(ind_SNP)==T]=0 
   Ind_all=apply(ind_SNP,2,function(x){sum(x,na.rm=T)})        
   prop_ind=Ind_nsnp/Ind_all
   prop_ind[is.na(prop_ind)==T]=NA   
   X=cbind(prop_ind,prop_ind)
   X=as.numeric(t(X)) 
   Yupper=nind-c(1:nind)+1+0.5
   Ylower=nind-c(1:nind)+1-0.5
   Y=cbind(Yupper,Ylower)
   genderY=cbind(gender,Y)
   Y=as.numeric(t(Y)) 
   temp_Xlim=c(0,ceiling(max(prop_ind,na.rm=T)*100)/100)
   genderY=cbind(genderY,temp_Xlim[1],temp_Xlim[2])
   plot(1,1,type="n",xlim=temp_Xlim,ylim=Ylim,ylab="",xlab="",main="",axes=F)
   nind=nrow(genderY)
   for(i in 1:nind)
   {
      temp_ind=genderY[i,]
      gender=temp_ind[1]
      loc=as.numeric(temp_ind[-1])
      if(gender=="M"){gendercol="skyblue"}
      if(gender=="F"){gendercol="pink"}      
      rect(loc[3],loc[2],loc[4],loc[1],col=gendercol,border=NA)
   }
   X=c(0,0,X,X[length(X)]) 
   Y=c(0.5,Y[1],Y,0.5) 
   polygon(X,Y,col="purple",density=45)
   abline(h=Yind-1,col=colset,lty=2,lwd=1.25)
   n=2;d=max(X)/n;Xlab=c(0:2)*d
   Xlab=round(Xlab*100)/100 
   axis(1,at=Xlab,labels=Xlab*100,line=YlabcexLine+2.5,cex.axis=cexsize,tick=F,font.axis=4)
   axis(1,at=median(Xlab),labels="Aberration (%)",line=YlabcexLine+5,col.axis="blue",cex.axis=labelsize*0.65,tick=F,font.axis=2)   
   abline(v=Xlab,col=colset,lty=2,lwd=1.25) 
   box(col="black") 
   par(mar=c(1,0,3,1))
   plot(1,1,type="n",xlim=c(0,10),ylim=c(0,10),axes=F,ylab="",xlab="")
   r=1 
   rect(0,0,10*r,5,col="skyblue",border="skyblue")
   rect(0,5,10*r,10,col="pink",border="pink")
   Xtext=5
   text(Xtext,2.5,labels="Male",cex=2,font=4)   
   text(Xtext,7.5,labels="Female",cex=2,font=4)      
   par(def.par)
   return()
}
simu_AILOH=function(nSNP,nind,p,r)
{
  set.seed(200)
  a=runif(nSNP*nind)
  a[a>p]=1;a[!a>p]=0
  a=matrix(a,nSNP,nind)
  if(r<1)
  {
    NAsnp=nSNP*r;NAsnp=c(1:(nSNP-NAsnp))+NAsnp
    a[NAsnp,]=2
  }
  a=cbind((1:nSNP),(1:nSNP)*1000000,a)
  return(a)
}
