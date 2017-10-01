"%w/o%" <- function(x,y) x[!x %in% y]

#################
### function  ###
#################
libname = .packages(all.available = TRUE)
if(!"locfit" %in% libname) install.packages("locfit", repos="http://cran.csie.ntu.edu.tw")
library(locfit)
Size2Prop <- function(SNPselectype, propsize, LengthChr){
  prop <- propsize
  if(SNPselectype == 1){
    if(propsize > LengthChr) propsize <- LengthChr
    prop <- propsize/LengthChr
  }
  prop
}

Read.User <- function(i, j, chromosome){
### 2010-12-27
### 20130117
### Data format can be fixed (Probe_set, Chromosome, Physical_position, Genotype)
### or can be provided by user
  Data <- read.table(paste(population[i,], Path[[i]][j], sep=""), sep = "\t", as.is = TRUE, stringsAsFactors=FALSE, header = T)
  #Data <- data.frame(Data[, 1], NA, Data[, 2:4])
  #Data[, 1] <- as.character(Data[, 1])
  if (!is.na(RsCol.column[i])){
    Data <- Data[, c(SnpCol.column[i], RsCol.column[i], ChrCol.column[i], PosiCol.column[i], CallCol.column[i])]
  }else{
    Data <- cbind(Data[, SnpCol.column[i]], rep("", nrow(Data)), Data[, c(ChrCol.column[i], PosiCol.column[i], CallCol.column[i])])
  }
  colnames(Data) <- c("Probe", "RS_number", "chr", "pos", "call")
  Data <- Data[which(Data[, 3] %in% chromosome), ]
  Data
}

Read.Affy <- function(i, j, Chromosome){
### 2010-12-25
### need path , AFFY6.SNP
### i    number- 1:number.popu
### j    number- 1:IDnumber
### Chr  vector- 1:23
### cancel the first two blank rows (skip = 2 in read.table) in the array data 20111226
### 20130109
### change into an individual a file instead of one chromosome of each individual a file
### 20130117
### Annotation file is not provided
  Data <- NULL
  Chr <- Chromosome
#  if(any(Chromosome %in% 23)) Chr <- c(Chr, 24)
  #for(chr in Chr){
  #  Data <- rbind(Data,
  #    cbind(read.table(paste(population[i,], path[[i]][[j]][chr], sep=""), header = T, sep="\t", stringsAsFactors=FALSE)[c(SnpCol.column[i], CallCol.column[i])], ifelse(chr==24, 23, chr)))
  #}
  Data <- read.table(paste(population[i,], Path[[i]][j], sep=""), sep = "\t", as.is = TRUE, stringsAsFactors = FALSE, header = T)
  Data <- Data[(Data[, ChrCol.column[i]] %in% Chr), ]
  #InAnnoSNP <- intersect(Data[, 1], AFFY6.SNP[, 1])
  #Data <- data.frame(Data[match(InAnnoSNP, Data[, 1]), ],
  #  AFFY6.SNP[match(InAnnoSNP, AFFY6.SNP[, 1]), 2:3])
  if (!is.na(RsCol.column[i])){
    Data <- Data[, c(SnpCol.column[i], PosiCol.column[i], ChrCol.column[i], CallCol.column[i], RsCol.column[i])]
  }else{
    Data <- cbind(Data[, c(SnpCol.column[i], PosiCol.column[i], ChrCol.column[i], CallCol.column[i])], rep("", nrow(Data)))
  }
  colnames(Data) <- c("Probe", "pos", "chr", "call", "RS_number")
  Data[, "pos"] <- as.numeric(Data[, "pos"])
  #Data[, c(1, 4, 3, 5, 2)]
  Data
}

Read.500K <- function(i, j, Chromosome, all = FALSE, chip = Chip){
### 2010-12-29
### need path
### i    number- 1:number.popu
### j    number- 1:IDnumber
### Chr  vector- 1:23
### cancel the first two blank rows (skip = 2 in read.table) in the array data 20111226
  Data <- NULL
#  if(any(Chromosome %in% 23)) Chr <- c(Chr, 24)
  for(a in chip){
    if (!is.na(RsCol.column[i])){
      Data <- rbind(Data,
        read.table(paste(population[i,], path[[i]][[j]][a], sep=""), header = T, sep="\t", stringsAsFactors = FALSE)[c(SnpCol.column[i], RsCol.column[i], ChrCol.column[i], PosiCol.column[i], CallCol.column[i])])
    }else{
      Data_temp <- read.table(paste(population[i,], path[[i]][[j]][a], sep=""), header = T, sep="\t", stringsAsFactors = FALSE)[c(SnpCol.column[i], ChrCol.column[i], PosiCol.column[i], CallCol.column[i])]
      Data_temp <- cbind(Data_temp[, 1], rep("", nrow(Data_temp)), Data_temp[, 2:4])
      Data <- rbind(Data, Data_temp)
    }
  }
  if (length(which(Data[, 3] == "X")) > 0) Data[which(Data[, 3] == "X"), 3] <- 23
  if (length(which(Data[, 3] == "")) > 0) Data[which(Data[, 3] == ""), 3] <- NA
  Data[, 3] <- as.numeric(Data[, 3])
  if(!all) Data <- Data[which(Data[, 3] %in% Chromosome), ]
  colnames(Data) <- c("Probe", "RS_number", "chr", "pos", "call")
  Data
}


ReadData <- function(type, i, j, chromosome, ...){
  switch(type,
        "1" = Read.500K(i, j, chromosome, ...), 
        "2" = Read.500K(i, j, chromosome, ...),
        "3" = Read.Affy(i, j, chromosome),
        "4" = Read.User(i, j, chromosome)
  )
}



OrgRawData <- function(Data, chr.num){
### 2010-09-16 for user-provided
### 2010-12-25 use name
  Data <- Data[which(Data[, "Probe"] %in% interSNP), ]
  Data <- Data[Data[, "chr"] == chr.num, ]
  Data <- Data[order(Data[, "pos"]), ]
  Data
}

if(FALSE){
CombRawData <- function(i, chr.num, type = paste(chiptype), ...){
### 2010-12-27 add chiptype
  chr <- sprintf("%02.f", chr.num)    
  DATA <- ReadData(type, i, 1, chr.num, ...)
  Data <- OrgRawData(DATA, chr.num)[, c("Probe", "RS_number", "pos", "call")]
  #dim(Data) 
  if(IDnumber[i] > 1){
    for(j in 2:IDnumber[i]){      
      DATA <- ReadData(type, i, j, chr.num, ...)
      Data <- cbind(Data, OrgRawData(DATA, chr.num)[, "call"])
      rm(DATA); gc()
    }
  }
  colnames(Data) <- c("Probe", "RS_number", "Pos", paste(i, Label[[i]], sep="_"))
  Data
}
}

CombRawData <- function(i, chr.num, type = paste(chiptype), ...){
### 2010-12-27 add chiptype
### 2011-01-07 change cbind to  assign  
### 2011-05-15 
  chr <- sprintf("%02.f", chr.num)    
  DATA <- ReadData(type, i, 1, chr.num, ...)
  DATA <- OrgRawData(DATA, chr.num)
  Data <- as.data.frame(matrix(NA, nrow=nrow(DATA), ncol=3+length(IDnumber[i])))
  Data[,1:4] <- DATA[, c("Probe", "RS_number", "pos", "call")]  
  if(IDnumber[i] > 1){
    for(j in 2:IDnumber[i]){
      DATA <- ReadData(type, i, j, chr.num, ...)
      Data[, 3+j] <- OrgRawData(DATA, chr.num)[, "call"]
      rm(DATA); gc()
    }
  }
  colnames(Data) <- c("Probe", "RS_number", "Pos", paste(i, Label[[i]], sep="_"))
  Data
}


Coding <- function(SNPdata, Type=c("LOH", "RV"), Sequencing = FALSE, Rare.threshold = NULL, MAF.weight = NULL){
### 2010-09-24
### 2012-12-10
### 2012-12-25 Change "BA" call into "AB" call
### 2013-06-11 return coding and MAF
### 2014-06-20 save rare variants major homozygousity as 1 and add MAF/Rare.threshold as weight for minor
  if(MAF.weight == "Y") snp.weight <- rep(1, length(SNPdata))
  SNPdata[which(SNPdata == "BA")] <- "AB"
  if(Type == "LOH"){
    
    TABLE <- table(factor(SNPdata, c("AA", "AB", "BB")))
    number.A <- TABLE[1]*2+ TABLE[2]
    number.B <- TABLE[3]*2+ TABLE[2]
    minor.prop <- min(number.A/(number.A+number.B), number.B/(number.A+number.B))
    if (all(TABLE == 0)){
        Data <- rep(NA, length(SNPdata))
        minor.prop <- 0
    }else{}

    if (!Sequencing){
      Data <- as.numeric((SNPdata=="AA")|(SNPdata=="BB"))
    } else {
      Data <- rep(0, length(SNPdata))
      Data[is.na(SNPdata)] <- NA

      if (minor.prop >= 0 & minor.prop <= Rare.threshold){
        # For rare SNP in sequencing data, we are interesting in minor homozygous and heterzygous call
        Data[which(SNPdata == "AA" | SNPdata == "BB")] <- 1
        # For rare SNP in sequencing data, major homozygous call is non-informative
        if(MAF.weight == "Y") snp.weight[which(SNPdata == ifelse(number.A > number.B, "AA", "BB"))] <- minor.prop/Rare.threshold
      } else {
        # For common variants in sequencing data, we are interesting in homozygous call
        Data[which(SNPdata == "AA" | SNPdata == "BB")] <- 1
      }
    }
  }
  if(Type == "RV"){
    TABLE <- table(factor(SNPdata, c("AA", "AB", "BB")))
    number.A <- TABLE[1]*2+ TABLE[2]
    number.B <- TABLE[3]*2+ TABLE[2]
    minor.prop <- min(number.A/(number.A+number.B), number.B/(number.A+number.B))
    Data <- rep(0, each=length(SNPdata))
    if(minor.prop > 0 & minor.prop <= 0.01){
      if(number.A > number.B){
        Data[which(SNPdata=="BB"|SNPdata=="AB")] <- 1
        Data[which(SNPdata=="AA")] <- 0
		  }
      else{
        Data[which(SNPdata=="AA"|SNPdata=="AB")] <- 1
        Data[which(SNPdata=="BB")] <- 0
      }
    }
  } 
  #Data
  if(MAF.weight == "Y"){
	list(coding = Data, maf = minor.prop, weight = snp.weight)
  }else{
    list(coding = Data, maf = minor.prop)
  }
}

Localfit <- function(SNPname, RSnumber, Coading, Position, alpha, deg=Deg, weights = rep(1, length(Coading))){
  ## 2010-09-16 based on LOHfit
  ## 2010-12-27
  Coading <- as.numeric(Coading) ###
  Position <- as.numeric(Position)/10^6
  lf <- locfit(Coading ~ Position, alpha = alpha, deg=deg, weights = weights);
  LOH.est <- predict(lf, Position)
  LOH.est[LOH.est > 1] <- 1
  LOH.est[LOH.est < 0] <- 0
  #data.frame(SNPname, as.character(round(Position, 6)), as.character(round(LOH.est, digit=4)))
  #data.frame(SNPname, round(Position, 6), round(LOH.est, digit=4))
  cbind(SNPname, RSnumber, Position, round(LOH.est, digit=4))
}

CombData <- function(i, chr.num, chip12=NULL){
### 2010-09-20
### 2010-12-28
### 2010-12-30 for 500K, affy 6.0,
### 2010-01-03  for 500K merge two or only 1
### 2010-01-07 change cbind to assign
### 2010-01-1 change for User
### 2012-12-24 save .RData to speed up
### 2013-08-15 change for draw biplot only
  chr <- sprintf("%02.f", chr.num)
  #if(chiptype %in% 3:4)
  #  DATA <- read.csv(paste(pathnumLOH, PopGroup, i, "/", Label[[i]][1], "/Num_Popu",i,"_", Label[[i]][1],"_Chr", chr,".csv", sep=""))
  #if(chiptype %in% 1:2)
  #  DATA <- read.csv(paste(pathnumLOH, PopGroup, i, "/", Label[[i]][1], "/Num_Popu",i,"_", Label[[i]][1], "_Chip", paste(chip12, collapse=""), "_Chr", chr,".csv", sep=""))
  if (!plot.only){
    load(paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][1], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(chip12, collapse = ""), sep = "")), "_Chr", chr, ".RData", sep = ""))
  }else{
    load(paste(inpath, "/Num_Popu", i, "_", Label[[i]][1], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(chip12, collapse = ""), sep = "")), "_Chr", chr, ".RData", sep = ""))
  }
  DATA <- partdata
  Data <- as.data.frame(matrix(NA, nrow=nrow(DATA), ncol= 3 + length(IDnumber[i])))
  #partdata <- DATA
  #if (!file.exists(paste(pathnumLOH, "temp", sep = ""))) dir.create(paste(pathnumLOH, "temp", sep = ""), showWarnings = F)
  #if (!file.exists(paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][1], ifelse(chiptype %in% 3:4, "", paste("_Chip", chip12, sep = "")), "_Chr", chr, ".RData", sep = ""))) save(partdata, file = paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][1], ifelse(chiptype %in% 3:4, "", paste("_Chip", chip12, sep = "")), "_Chr", chr, ".RData", sep = ""))
  if (is.na(RsCol.column[i])){
  #if(chiptype %in% 4){
    Data[, 1] <- DATA[, 2]; Data[, 3:4] <- DATA[, 3:4]
  }else{
    Data[, 1:4] <- DATA
  }
  if(IDnumber[i] > 1){
    for(j in 2:IDnumber[i]){
      #if(chiptype %in% 3:4)
      #  save.dir <- paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], "/Num_Popu",i,"_", Label[[i]][j],"_Chr", chr,".csv",sep="");
      #if(chiptype %in% 1:2)
      #  save.dir <- paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], "/Num_Popu",i,"_", Label[[i]][j], "_Chip", paste(chip12, collapse=""), "_Chr", chr,".csv", sep="");
      ##Data[, 3+j] <- read.csv(save.dir)[, 4]
      #partdata <- read.csv(save.dir)
      if (!plot.only){
        load(paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(chip12, collapse = ""), sep = "")), "_Chr", chr, ".RData", sep = ""))
      }else{
        load(paste(inpath, "/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(chip12, collapse = ""), sep = "")), "_Chr", chr, ".RData", sep = ""))
      }
      Data[, 3+j] <- partdata[, 4]
      #if (!file.exists(paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", chip12, sep = "")), "_Chr", chr, ".RData", sep = ""))) save(partdata, file = paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", chip12, sep = "")), "_Chr", chr, ".RData", sep = ""))
    }
  }
  colnames(Data) <- c("Probe", "RS_number", "Pos", paste(i, Label[[i]], sep="_"))
  Data
}
# CombData(1, 1)


nonzero <- function(pos, threshold){
  ### pick regions with more than 5 snps nearby.
  ### chr.data <- comb.data[chromo==names.chromo[k],];
  #s <- seq(min(pos),max(pos),0.5);
  pos <- as.numeric(pos)
  nz <- sapply(1:length(pos), function(k) sum(abs(pos-pos[k]) < 0.5))
  nz > threshold
}

KL <- function(lambda0, lambda1, r, type=c("max", "sum")){
###uses sum or max of KL distances 
  lambda0 <- apply(as.matrix(lambda0), 1, function(x){min(x, 0.99)});
  lambda0 <- apply(as.matrix(lambda0), 1, function(x){max(x, 0.01)});
  lambda1 <- apply(as.matrix(lambda1), 1, function(x){min(x, 0.99)});
  lambda1 <- apply(as.matrix(lambda1), 1, function(x){max(x, 0.01)});
  a1 <- (lambda0/lambda1);
  a2 <- (1-lambda0)/(1-lambda1);
  a3 <- log2(a1)*lambda0+log2(a2)*(1-lambda0);
  a4 <- a1>1;
  #print(a3*r*a4)
  switch(type,
    max = max(a3*r*a4),
    sum = sum(a3*r*a4))
}

rank.test <- function(fit.0, fit.1, alpha0, snp, deg=Deg){
### 20100916, 20100924
  mu.norm <- apply(cbind(fit.0, fit.1), 1, function(k) mean(as.numeric(k)));
  m11 <- locfit(mu.norm ~ snp, alpha = alpha0, deg=deg);
  m1 <- predict(m11, snp);
  if(MemoryEnough){
    yy <- sapply(1:ncol(fit.1), function(j) KL(as.numeric(fit.1[,j]), m1, nonzero(snp, 5), "max"))
    xx <- sapply(1:ncol(fit.0), function(j) KL(as.numeric(fit.0[,j]), m1, nonzero(snp, 5), "max"))
  }
  if(!MemoryEnough){
    xx <- rep(0,ncol(fit.0));
    yy <- rep(0,ncol(fit.1));
    for(j in 1:ncol(fit.1)){
      yy[j] <- KL(as.numeric(fit.1[, j]), m1, nonzero(snp, 5), "max"); #can change function
    }
    for(j in 1:ncol(fit.0)){
      xx[j] <- KL(as.numeric(fit.1[, j]), m1, nonzero(snp, 5), "max") #can change function
    }
  }
  wilcox.test(xx, yy, alternative="less")
}
# rank.test(Data[[1]][, -(1:2)], Data[[2]][, -(1:2)], Size2Prop(length(Data[[1]][, 2]), propsize), Data[[1]][, 2])


  

RankTest <- function(chr.num, interSNP, ...){
### 2010-09-01
### 2010-12-28
### 2011-05-15 :  inter_SNP
  Data <- list(NULL)
  for(i in 1:number.popu) {Data[[i]] <- CombData(i, chr.num, ...)}
  inter_SNP <- as.character(Data[[1]][which(Data[[1]][, 1] %in% interSNP), 1])
  for(i in 1:number.popu){
    Data[[i]] <- Data[[i]][match(inter_SNP, Data[[i]][, 1]), ]
  }
  Position <- Data[[1]][, 3]
  rank.test(Data[[1]][, -(1:3)], Data[[2]][, -(1:3)], Size2Prop(SNPselectype, propsize, length(Position)), Data[[1]][, 3])$p.value;
}
#RankTest(1, interSNP)

####################################
#####    Quantitle function    #####
####################################
QLOHcontrol <- function(num.group, chr.num, popu, Data, QLOHselect=QLOHselect){
  ## only use control data
  ## 2011-01-03 for 500K
  ## 2011-01-04 for Chr23
  Quantile <- function(x){quantile(as.numeric(x), QLOHselect)}
  if(num.group==1){
    QLOH <- matrix(QLOHselect, nrow=nrow(Data), ncol=ifelse(chr.num == 23, 2, 1))  
    if(chr.num == 23)
      colnames(QLOH) <- c(paste(QLOHselect*100, "%_QI_F", sep=""),paste(QLOHselect*100, "%_QI_M", sep=""))
    else
      colnames(QLOH) <- paste(QLOHselect*100, "%_QI", sep="")
  }
  if(num.group == 2){
    if(chr.num == 23){
      female <- colnames(Data) %in% paste(popu, femaleID[[popu]], sep="_")
      male <- !female; male[1:3] <- FALSE
      QLOH <- cbind(apply(as.matrix(Data[, female], byrow=TRUE, ncol=sum(female)), 1, Quantile),
        apply(as.matrix(Data[, male], byrow=TRUE, ncol=sum(male)), 1, Quantile))
      colnames(QLOH) <- c(paste(QLOHselect*100, "%_QI_F", sep=""),paste(QLOHselect*100, "%_QI_M", sep=""))
    }
    if(chr.num != 23){            
      QLOH <- as.matrix(apply(as.matrix(Data[, -(1:3)], byrow=TRUE, ncol=IDnumber[[popu]]), 1, Quantile), ncol=1)
      colnames(QLOH) <- paste(QLOHselect*100, "%_QI", sep="")
    }    
  }
  QLOH  
}
############
svd.GB <- function(A){
### from LC
## set "Inf" and "-Inf" to 0 in inverse_diagEV 20111226
  evA=t(A)%*%A
  ##-----------------------------------------
  evA_value=eigen(evA)$values;evA_value=round(evA_value*1000000)/1000000
  V=eigen(evA)$vectors;#V=round(V*1000)/1000
  ##-----------------------------------------
  transposeV <- t(V)
  evA_value_sqrt <- sqrt(evA_value);
  diagEV=diag(evA_value_sqrt);
  inverse_diagEV <- diag(1/evA_value_sqrt);#diagEV=round(diagEV*1000)/1000
  inverse_diagEV[is.infinite(inverse_diagEV)]=0
  U12=A%*%V[,1:2]
  U12[,1]=U12[,1]*inverse_diagEV[1,1]
  U12[,2]=U12[,2]*inverse_diagEV[2,2]
  list(d=sqrt(evA_value),u=U12,v=V)
}

svd.CG <- function(svd.GB, alpha){
## for svd.GB from LC 20100902
    n <- length(svd.GB$d)
    lambda <- diag(svd.GB$d)
    if(alpha==1){
      C <- svd.GB$u%*%lambda
 	    G <- svd.GB$v
  	}
 	  else if(alpha==0){
	      C <- svd.GB$u*sqrt(n)
        G <- t(lambda%*%t(svd.GB$v))/sqrt(n)
 	  }    
    list(C = C, G = G)
}
# svd.CG(as.data.frame(MM)[1:5,], 0)

marker.color<-function(Length){
  num <- round(Length/15)
  rep(terrain.colors(16)[15:1], c(rep(num, 14), Length - num*14))
}

demo.color<-function(){
  ch.col <- terrain.colors(16)[15:1]
  plot(1:30,type="n",axes=F)
  for(i in 1:15){
    rect(8,8+(i-1),14,9+(i-1),col=ch.col[i])
  }
  text(rep(14,3),c(8,15.5,23), c(1, round(length(CG$C[, 1])/2), length(CG$C[, 1])),pos=4)
}

#########
demo.newcolor.pq <- function(P, Q, position){
## from LC
  ch.col <- terrain.colors(16)[15:1]
  plot(1:30,type="n",axes=F)
  #############################################
  start.text <- (floor( min(position)*100 ))/100
  end.text<- (floor( max(position)*100 ))/100
  number.position <- length(position)
  p.text<- (floor( P*100 ))/100
  q.text<- (floor( Q*100 ))/100
  #############################################
  CUT.point<-c()
  length.center<- 0
  CUTpoint.length <- (max(position)-min(position)-(Q-P))/15
  #############################################
  for(i in 1:15){
    CUT.point[i]<-min(position)+ ( (max(position)-min(position)-(Q-P) )/15)*i + length.center
    if(CUT.point[i]>P & CUT.point[i]<Q){
      length.center<- Q-P
      CUT.point[i]<- CUT.point[i]+ length.center
    }
  }
  #############################################
  b <- 0
  k <- 0
  for(i in 1:15){
    if(P<=CUT.point[i] & k==0){
	  	if(i==1){
        length.temp <- CUT.point[i]-min(position)
        length.up <- CUT.point[i]-P
        length.down <- P- min(position)
        prop.up <- (CUT.point[i]-P)/length.temp
        prop.down <- (P- min(position))/length.temp
        rect(8,6+(i-1)+b,14,6+ prop.down+ (i-1)+ b,col=ch.col[i])
        coord.y.pend <- 6+ prop.down+ (i-1)+ b
        prop.center <- length.center/ CUTpoint.length
        if(prop.center<0.3){
          b<- prop.center
          rect(8,6+ prop.down+(i-1)+b,14,6+ prop.down+ prop.up+ (i-1)+ b,col=ch.col[i])
          points(11,6+ prop.down+ (i-1)+ prop.center/2,pch=1,cex=2, col="red")
        }
        else{
          points(11,6+ prop.down+ (i-1)+ prop.center/2,pch=1,cex=2, col="red")
          # line from y: 8+ prop.down+ (i-1)+ b  to 8+ prop.down+ (i-1)+ b+ prop.center
          b<- prop.center
          lines(c(11,11),c(6+ prop.down+ (i-1),6+ prop.down+(i-1)+b), lty=2, col="red")
          rect(8,6+ prop.down+(i-1)+b,14,6+ prop.down+ prop.up+ (i-1)+ b,col=ch.col[i])
        }
        coord.y.qstart<- 6+ prop.down+(i-1)+b
         coord.y.end<- 6+ prop.down+ prop.up+ (i-1)+ b
      }
      else{
        length.temp <- CUT.point[i]- CUT.point[i-1]
					 length.up <- CUT.point[i]-P
					 length.down<- P- CUT.point[i-1]
					 prop.up<- (CUT.point[i]-P)/length.temp
					 prop.down<- (P- CUT.point[i-1])/length.temp
					 rect(8,6+(i-1)+b,14,6+ prop.down+ (i-1)+ b,col=ch.col[i])
					 coord.y.pend<- 6+ prop.down+ (i-1)+ b
					 prop.center<- length.center/ CUTpoint.length
					 if(prop.center<0.3){
					  b<- prop.center
					  rect(8,6+ prop.down+(i-1)+b,14,6+ prop.down+ prop.up+ (i-1)+ b,col=ch.col[i])
					  points(11,6+ prop.down+ (i-1)+ prop.center/2,pch=1,cex=2, col="red")
     					 }
					 else
					 {
					  points(11,6+ prop.down+ (i-1)+ prop.center/2,pch=1,cex=2, col="red")
					  # line from y: 8+ prop.down+ (i-1)+ b  to 8+ prop.down+ (i-1)+ b+ prop.center
					  b<- prop.center
					  lines(c(11,11),c(6+ prop.down+ (i-1),6+ prop.down+(i-1)+b), lty=2, col="red")
					  rect(8,6+ prop.down+(i-1)+b,14,6+ prop.down+ prop.up+ (i-1)+ b,col=ch.col[i])
					 }
					 
					 coord.y.qstart<- 6+ prop.down+(i-1)+b
					 coord.y.end<- 6+ prop.down+ prop.up+ (i-1)+ b
					}
				 	k<-k+1
    				 }
				 else
				 {
				  length.temp<- CUT.point[i]- CUT.point[i-1]
				  prop.up<- (CUT.point[i]-P)/length.temp
				  prop.down<- (P- CUT.point[i-1])/length.temp
				  rect(8,6+(i-1)+b,14,7+(i-1)+b,,col=ch.col[i])
				  coord.y.end<- 6+ prop.down+ prop.up+ (i-1)+ b
				 }			

				}
			 ##### Text #####
			 if(chr.num==21)
			 {
			  text(14, 4.7, sprintf("%.2f",start.text), pos=4)
			  text(14, 4.2, "(1)", pos=4, cex=0.8)
			  text(14, coord.y.end+0.8, sprintf("%.2f",end.text), pos=4)
			  text(14, coord.y.end+0.3, paste("(", length(position), ")", sep=""), pos=4, cex=0.8)
			  if(sum(position<=P)==1)
			  {}
			  else
			  {
			   text(14, coord.y.pend-0.3, sprintf("%.2f",p.text), pos=4)
			   text(14, coord.y.pend-0.8, paste("(", sum(position<=P), ")", sep=""), pos=4, cex=0.8)
			  }
			  text(14, coord.y.qstart+0.8, sprintf("%.2f",q.text), pos=4)
			  text(14, coord.y.qstart+0.3, paste("(", sum(position<=P)+1, ")", sep=""), pos=4, cex=0.8)

			 }
      else{
          text(14, 5.7, sprintf("%.2f",start.text), pos=4)
          text(14, 5.3, "(1)", pos=4, cex=0.8)
          text(14, coord.y.end+0.8, sprintf("%.2f",end.text), pos=4)
          text(14, coord.y.end+0.3, paste("(", length(position), ")", sep=""), pos=4, cex=0.8)
          text(14, coord.y.pend-0.3, sprintf("%.2f",p.text), pos=4)
          text(14, coord.y.pend-0.8, paste("(", sum(position<=P), ")", sep=""), pos=4, cex=0.8)
          text(14, coord.y.qstart+0.8, sprintf("%.2f",q.text), pos=4)
          text(14, coord.y.qstart+0.3, paste("(", sum(position<=P)+1, ")", sep=""), pos=4, cex=0.8)
      }
}

demo.newcolor.q <- function(Q, position){
### from LC
  ch.col<-terrain.colors(16)[15:1]
  plot(1:30,type="n",axes=F)
  #############################################
  end.text <- (floor( max(position)*100 ))/100
  number.position<- length(position)
  q.text <- (floor( Q*100 ))/100
  #############################################
  CUT.point <- c()
  CUTpoint.length <- ( max(position)-min(position) )/15
	#############################################
  CUT.point[1]<- CUT.point[1]+ min(position)
  for(i in 2:15){
			  CUT.point[i] <- CUTpoint.length*i+ min(position)
  }
  #############################################
  b <- 0
  for(i in 1:15){
    if(i==1){
      points(11,6+ (i-1)+ b,pch=1,cex=2, col="red")
      prop.center<- Q/CUTpoint.length
      b<- prop.center
      lines(c(11,11),c(6+ (i-1),6+ (i-1)+b), lty=2, col="red")
      rect(8,6+(i-1)+b,14,7+ (i-1)+ b,col=ch.col[i])
      coord.y.qstart<- 6+(i-1)+b
    }
    else{
      rect(8,6+(i-1)+b,14,7+(i-1)+b,,col=ch.col[i])
      coord.y.end<- 7+ (i-1)+ b
    }
  }
  ##### Text #####
  text(14, coord.y.qstart-0.3, sprintf("%.2f",q.text), pos=4)
  text(14, coord.y.qstart-0.8, "(1)", pos=4, cex=0.8)
  text(14, coord.y.end+0.8, sprintf("%.2f",end.text), pos=4)
  text(14, coord.y.end+0.3, paste("(", length(position), ")", sep=""), pos=4, cex=0.8)
}


BIPLOT <- function(C, G, row.col=1, col.col=4, ind_name_list=NULL, row_name_list=NULL, title=NULL, All=FALSE,...){
  range_C1 <- c(min(C[,1]), mean(C[,1],trim=0.05), max(C[,1]) );
  Crangex <- max(range_C1[2]-range_C1[1], range_C1[3]-range_C1[2]);
  new.Cxlim <- c(range_C1[2]-Crangex, range_C1[2]+Crangex);
  range_C2 <- c(min(C[,2]), mean(C[,2],trim=0.05), max(C[,2]) );
  Crangey <- max(range_C2[2]-range_C2[1], range_C2[3]-range_C2[2]);
  new.Cylim <-c(range_C2[2]-Crangey, range_C2[2]+Crangey);
  ###  
  plot(C[,1], C[,2], type="n", main="", xlim=new.Cxlim,
    ylim=new.Cylim, xlab="", ylab="", cex.axis=ifelse(All, 0.9, 0.8), ...);
  mtext(ifelse(All, "Marker_C1", "Marker (Component 1)"),side=1,line=2, cex=ifelse(All, 0.8, 1));
  mtext(ifelse(All, "Marker_C2", "Marker (Component 2)"),side=2,line=2, cex=ifelse(All, 0.8, 1));
  title(title, line=ifelse(All, 3, 4))
#  text(C[,1], C[,2], "~", cex=0.7, col=row.col);
  text(C[,1], C[,2], row_name_list, col=row.col, cex=0.7);
  text.color <- c(1, round(length(C[,1])/2),length(C[,1]))
  par(new=T);
  # G[,1]=G[,1]-mean(G[,1]); G[,2]=G[,2]-mean(G[,2])
  range_1<- c(min(G[,1]), mean(G[,1],trim=0.05), max(G[,1]));
  Grangex<- max(range_1[2]-range_1[1], range_1[3]-range_1[2]);
  new.Gxlim<-c(range_1[2]-Grangex, range_1[2]+Grangex);
  
  range_2<- c(min(G[,2]), mean(G[,2],trim=0.05), max(G[,2]) );
  Grangey <- max(range_2[2]-range_2[1], range_2[3]-range_2[2]);
  new.Gylim<-c(range_2[2]-Grangey, range_2[2]+Grangey);
  
  plot(G[,1],G[,2],type="n",xaxt="n",yaxt="n",xlab="",ylab="", xlim=new.Gxlim, ylim=new.Gylim);
  text(G[,1]*1.02, G[,2]*1.02, ind_name_list, col=col.col);
#  arrows(mean(new.Gxlim), mean(new.Gylim), G[,1], G[,2], length=ifelse(All, .05, .1), col=col.col); before 2010-09-23
  arrows(mean(new.Gxlim),mean(new.Gylim), G[,1], G[,2], length=ifelse(All, .05, .1), col=col.col);  
  axis(3, cex.axis=ifelse(All, 0.9, 0.8));
  axis(4, cex.axis=ifelse(All, 0.9, 0.8));
  mtext(ifelse(All, "Individual_C1", "Individual (Component 1)"), side=3, line=2, cex=ifelse(All, 0.8, 1))
  mtext(ifelse(All, "Individual_C2", "Individual (Component 2)"), side=4, line=ifelse(All, 2, 3), cex=ifelse(All, 0.8, 1))
}
