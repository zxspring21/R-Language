
cat("Please wait a while, KBAT is running...\n")

pkg.http <- "http://cran.opensourceresources.org/"
pkgs <- .packages(all.available=TRUE)
if(all(pkgs!="genetics")) install.packages("genetics",repos=pkg.http)

require(genetics)

file <- dir(PATH1)

PVALUE.path <- paste(PATH1,"\\pv.txt",sep="")
MAP.path <- paste(PATH1,"\\map.txt",sep="")
LD.path <- paste(PATH1,"\\ld.txt",sep="")
Geno.path <- paste(PATH1,"\\geno.txt",sep="\\")
STAT.out <- paste(PATHout,"\\Output.txt",sep="")
PLOT.out <- paste(PATHout,"\\Plot",sep="")

tmp <- ifelse(length(grep("pv",file))!=0,p.tot1<-read.table(file=PVALUE.path),p.tot1<-0)
tmp <- ifelse(length(grep("map",file))!=0,MAP<-read.table(file=MAP.path),MAP<-0)
tmp <- ifelse(LDdata==1,LDa<-read.table(file=LD.path,sep=","),
	 ifelse(LDdata==2,LDa<-read.table(file=Geno.path),LDa<-0))

Genotype <- function(x,Num)
{
 if(ncol(x)>Num)
 {
	for(i in 1:(ncol(x)/2))
	{
		name <- (paste("g",i,sep=""))
		assign(name,genotype(x[,i*2-1],x[,i*2]))
		if(i==1)
			gt <- g1
		else
			gt <- data.frame(gt,get(name))
	}
	g <- makeGenotypes(gt)
 }
 else
	g <- makeGenotypes(x)

 return(g)
}
tmp <- ifelse(LDdata==2,LDa <- Genotype(LDa,Num),LDa)

ld <- function(data,LDdata)
{
 if(LDdata==1)
 {
	data1 <- as.matrix(data)
	maxL <- max(data1[,1],data1[,2])
	AA <- matrix(0,maxL,maxL)
	AA[(data1[,1]-1)*maxL+data1[,2]] <- data1[,3]
	LD <- abs(AA)
 }
 else
 {
	ld <- LD(data)

	attributes(ld)
	newld <- as.matrix(ld$"D'")
	newt <- as.matrix(newld)
	newt[lower.tri(newt)] <- 0
	diag(newt) <- 1
	new2 <- t(newt)
	diag(new2) <- 0
	LD <- abs(as.matrix(new2+newt))
 }
 
 return(LD)
}

tmp <- ifelse(LDdata!=0,LD <- ld(LDa,LDdata),LDa)

prog<-paste(KBAT.path,"\\Program\\MC.R",sep="") 
source(prog)
