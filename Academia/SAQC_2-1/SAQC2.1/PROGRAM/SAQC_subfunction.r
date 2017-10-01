##==============================================================================
##------------------------  Quality Index of Chip ------------------------------
##==============================================================================
QI_chip=function(QI)
{
  #QI=Chip1
  QI[,3]=QI[,3]-QI[,2]
  for(k in 4:6){QI[,k]=(QI[,3])*QI[,k]}
  ind_QI=apply(QI,2,sum)
  ind_QI[4]=ind_QI[4]/ind_QI[3];
  ind_QI[5]=ind_QI[5]/ind_QI[3];
  ind_QI[6]=ind_QI[6]/ind_QI[3];
  ind_QI[2]=sprintf("%.f",ind_QI[2]);
  ind_QI[3]=sprintf("%.f",ind_QI[3]);
  ind_QI[4]=sprintf("%.4f",ind_QI[4]);
  ind_QI[5]=sprintf("%.4f",ind_QI[5]);
  ind_QI[6]=sprintf("%.4f",ind_QI[6]);
  ind_QI=matrix(ind_QI,1,6)
  return(ind_QI)
}
################################################################################
##                        input file checking                                 ##
################################################################################
infile_checking=function(infile_path, inputformat, Ind_Array)
{
	#infile_path=pathway_input;inputformat=1
	#infile_path=paste("F:/working_area/SAQC/EXAMPLE/Test_Affy100K",sep="");inputformat=2 
	cat("   Input files checking process.\n") 
	if(inputformat==1)## Genotype/Intensity
	{
		cat("Input data is the Genotype/Intensity based.\n") 
		GIfile=list.files(infile_path,full.names=T)
		GIname=list.files(infile_path,full.names=F)
		GIname=mixedsort(GIname)
		if(!("IndGeno" %in% GIname)){stop("SAQC error message: IndGeno folder cannot be found.\n")}
		if(Ind_Array == 1 | Ind_Array == 2){
			if(!("IndPI" %in% GIname)){stop("SAQC error message: IndPI folder cannot be found.\n")}
			GT_name=list.files(GIfile[2],full.names=F);nGT=length(GT_name)
			PI_name=list.files(GIfile[3],full.names=F);nPI=length(PI_name)
			if(nGT!=nPI){stop("SAQC error message: The number of files are different between IndGeno and IndPI.\n")}
		}
	}
  if(inputformat==2)## Allele frequency
  {
    cat("Input data is the adjust allele frequency based.\n") 
    indfile=list.files(infile_path,full.names=T)
    indname=list.files(infile_path,full.names=F)
    nind=length(indname)
    #whole_gender=NULL;whole_CR=NULL;
    
    cat("IndID","\t","gender","\t","CallRate"," \n") 
    for(k in 1:nind)
    {
      #k=1 
      chrfile=list.files(indfile[k],full.names=T);
      chrname=list.files(indfile[k],full.names=F);
      nchr=length(chrname)
      chrname=strsplit(chrname,".txt");chrname=unlist(chrname)
      chrname=unlist(strsplit(chrname,"Chr_"))[(1:nchr)*2]
      chrname=as.numeric(chrname)
      #cat(indname[k],": ",chrname,"\n")
      if(nchr!=23)
      {
         error_message=paste("SAQC error message: The individual's chromosomal information is uncompleted.",sep="")
         stop(error_message)
      }
      
        if(nchr==23)
        {
          whole_geno=NULL;
          for(ichr in 1:nchr)
          {
            #ichr=1
            if(file.exists(chrfile[ichr])==F)
            {
                error_message=paste("SAQC error message: The file of ",chrname[ichr]," of ",indname[k],"does not exist.",sep="")
                stop(error_message)
            }
            if(file.exists(chrfile[ichr])==T)
            {
              chr=read.table(chrfile[ichr],head=T,sep="\t")
              ColName=colnames(chr)
              chr=as.matrix(chr)
              if(dim(chr)[1]==0){error_message=paste("SAQC error message: No information in chromosome ",ichr,". ",sep="");stop(error_message)}
              # if(dim(chr)[2]!=5){error_message=paste("SAQC error message: The column's infomation are incorrect. ",sep="");stop(error_message)}
              # if(ColName[1]!="Probe_set"){error_message=paste("SAQC error message: The column of 'Probe_set' are incorrect. ",sep="");stop(error_message)}
              # if(ColName[2]!="Phy_position"){error_message=paste("SAQC error message: The column of 'Phy_position' are incorrect. ",sep="");stop(error_message)}
              # if(ColName[3]!="Chiptype"){error_message=paste("SAQC error message: The column of 'Chiptype' are incorrect. ",sep="");stop(error_message)}
              # if(ColName[4]!="Genotype"){error_message=paste("SAQC error message: The column of 'Genotype' are incorrect. ",sep="");stop(error_message)}
              # if(ColName[5]!="adjustAF"){error_message=paste("SAQC error message: The column of 'adjustAF' are incorrect. ",sep="");stop(error_message)}
              # Probeset=chr[,1];PhyPos=chr[,3]; Geno=chr[,4]; AF=chr[,6]
                ##----------------------------------------------------------------------------
                #   Probe_set     Phy_position  Chiptype  Genotype     adjustAF
                #  SNP_A-1677174   836727         2         AA           0.916
                #  SNP_A-1718890   2224111        1         BB           0.048
                #  SNP_A-1678466   2915399        2         BB           0.215
               #if(nrow(chr)==0){stop(error_message)}
            }
            whole_geno=c(whole_geno,chr[,4])
            if(ichr==nchr){sex_geno=chr[,4]}
          }
          genotable=table(whole_geno); #print(names(genotable)); 
		  names(genotable) = gsub(" ", "", names(genotable)); #print(names(genotable))
		  genotable=data.frame(t(as.matrix(genotable))); #print(names(genotable));
		   
		  # print(genotable$AA);print(genotable$AB);print(genotable$BB);print(sum(genotable))
          CR=sum(genotable$AA+genotable$AB+genotable$BB)/sum(genotable)
          sexgenotable=table(sex_geno);sexgenotable=data.frame(t(as.matrix(sexgenotable)))
          heterRatio=sum(sexgenotable$AB)/sum(sexgenotable)
          if(heterRatio<0.1){gender="M"}
          if(heterRatio>=0.1){gender="F"}
          #whole_CR=c(whole_CR,CR);whole_gender=c(whole_gender,gender)
          cat(indname[k],"\t", gender,"\t",CR," \n") 
        }    
    }
  }
  return("The input files are suitable for SAQC process")
}
##----------------- User provided the AFref checking function ------------------
AFref_checking=function(AFpath)
{
  #AFpath=pathway_AF;
  #AFpath=paste("F:/working_area/AFref_database/Asia_temp/SAQC-AF_ref",sep="") 
  chr_AFref_file=list.files(AFpath,full.names=T)
  chr_AFref_name=list.files(AFpath,full.names=F)
  nchr=length(chr_AFref_file)
  if(nchr!=24){stop("SAQC error message: The number of AF reference is incorrect, it should be 24 (22+X+Y) chromosomes.")}
  a=c(1:9);a=paste(0,a,sep="");b=c(a,10:24);b=paste("Chr_",b,".txt",sep="")
  if(sum(chr_AFref_name==b)!=24)
  {stop("SAQC error message: The file name of AF reference is incorrect.Please check the file name.")}
  whole_null_geno=NULL; 
  for(j in 1:nchr)
  {
    #j=1 
    chr=read.delim(chr_AFref_file[j],head=T,sep="\t", stringsAsFactors = F)
    nCol=ncol(chr)
    if(nCol!=12)
    {stop("SAQC error message: The number of column in ",chr_AFref_name[j]," is incorrect.Please check the number of column.")}
    ##--------------------------------------------------------------------------
    Coltag=c("Probe_set","Phy_position","Chiptype","num_AA","mean_AA","std_AA","num_AB","mean_AB","std_AB","num_BB","mean_BB","std_BB")
    Colchr=colnames(chr)
    if(sum(Colchr==Coltag)!=12)
    {stop("SAQC error message: The column name of ",chr_AFref_name[j]," is incorrect.Please check the column name.")}
    ##--------------------------------------------------------------------------
    ##  empty SNP checking
    num_Geno=chr[,c(1,3,4,7,10)]
    total_geno=apply(num_Geno[,-c(1,2)],1,sum)
    num_Geno=cbind(num_Geno,total_geno)
    null_geno=num_Geno[num_Geno[,6]==0,]     
    if(nrow(null_geno)>0){null_geno=cbind(j,null_geno[,c(1,2)])}    
    colnames(null_geno)[1]="Chr" 
    probeset=null_geno[,2] 
    chr=chr[!chr[,1] %in% probeset,]
    write.table(chr,chr_AFref_file[j],col.names=T,row.names=F,quote=F,sep="\t") 
    whole_null_geno=rbind(whole_null_geno,null_geno)
  }
  
  return(whole_null_geno)
}
##==============================================================================



################################################################################
##                        public sub-function area                            ##
################################################################################
##--------------------------public tool-----------------------------------------
sort.append<-function(x,base)
{
	library(gtools)
	sort.base<-mixedsort(paste(x[,base],1:length(x[,base])))
	tmp<-strsplit(sort.base," ")
	tmp1<-matrix(as.integer(unlist(tmp)),2,)
	y<-x[tmp1[2,],]
	return(y)
}
#-------------------------------------------------------------------------------
##==============================================================================
################################################################################
##                              Zone 1                                        ##
##                calculate the AF from intensity data                        ##
################################################################################
##----------------------- sub funtion-------------------------------------------
##該程式主要是從14個quartets上面求得intensity的值
##return的部份是IRAS

MPAMindiv = function(AFdata){

#==========================#
#===  feature extration ===#   
#==========================#    
# PAs: perfect match probes for the sense strand of allele A  #    # PAas: perfect match probes for the antisense strand of allele A  #
# PBs: perfect match probes for the sense strand of allele B  #    # PBas: perfect match probes for the antisense strand of allele B  #
# MAs:     mis-match probes for the sense strand of allele A  #    # MAas:     mis-match probes for the antisense strand of allele A  #
# MBs:     mis-match probes for the sense strand of allele B  #    # MBas:     mis-match probes for the antisense strand of allele B  #
	pp = seq(1, 56, by=8) + 1 # The first column: SNP ID
		PAs=AFdata[, pp];        	PBs=AFdata[, pp + 1];	MAs=AFdata[, pp + 2];	MBs=AFdata[, pp + 3];        
		PAas=AFdata[, pp + 4];	PBas=AFdata[, pp + 5];	MAas=AFdata[, pp + 6];	MBas=AFdata[, pp + 7];        
		
	##	Sense strand
	Ns = matrix(0, nrow=nrow(MAs), ncol=ncol(MAs));
		Ns[!is.na(MAs)] = Ns[!is.na(MAs)] + 1;		MAs[is.na(MAs)] = 0;
		Ns[!is.na(MBs)] = Ns[!is.na(MBs)] + 1;		MBs[is.na(MBs)] = 0;
	mPAs = PAs;	mPAs[is.na(mPAs)] = 0;
	mPBs = PBs;	mPBs[is.na(mPBs)] = 0;	
	mPAs = (mPAs-(MAs+MBs)/Ns);						mPAs[mPAs<0] = 0;
	mPBs = (mPBs-(MAs+MBs)/Ns);						mPBs[mPBs<0] = 0;

	posi = ((mPAs==0) & (mPBs==0) );      mPAs[posi] = NA;          mPBs[posi] = NA;       
	RAS_s = mPAs/(mPAs+mPBs);
			
	##	Anti-sense strand
	Nas = matrix(0, nrow=nrow(MAas), ncol=ncol(MAas));
		Nas[!is.na(MAas)] = Nas[!is.na(MAas)] + 1;	MAas[is.na(MAas)] = 0;
		Nas[!is.na(MBas)] = Nas[!is.na(MBas)] + 1;	MBas[is.na(MBas)] = 0;
	mPAas = PAas;	mPAas[is.na(mPAas)] = 0;
	mPBas = PBas;	mPBas[is.na(mPBas)] = 0;	
	mPAas = (PAas-(MAas+MBas)/Nas);                   mPAas[mPAas<0] = 0;
	mPBas = (PBas-(MAas+MBas)/Nas);                   mPBas[mPBas<0] = 0;                      
	
	posi = ((mPAas==0) & (mPBas==0) );    mPAas[posi] = NA;         mPBas[posi] = NA;
	RAS_as = mPAas/(mPAas+mPBas);      
	
	##	RASav
	RAS_s = apply(RAS_s, 1, median, na.rm=T);		Ns = (!is.na(RAS_s))*1;		RAS_s[is.na(RAS_s)] = 0
	RAS_as = apply(RAS_as, 1, median, na.rm=T);	Nas = (!is.na(RAS_as))*1;	RAS_as[is.na(RAS_as)] = 0
	RASav = (RAS_s + RAS_as)/(Ns + Nas)

	##
	##	HI
	##
	NA_s = (!is.na(PAs))*1;		NA_as = (!is.na(PAas))*1
	NB_s = (!is.na(PBs))*1;		NB_as = (!is.na(PBas))*1
	PAs[is.na(PAs)] = 0;			PAas[is.na(PAas)] = 0;	
	PBs[is.na(PBs)] = 0;			PBas[is.na(PBas)] = 0;
	HI = PAs + PAas + PBs + PBas;		HI[is.na(HI)] = 0
		HI = HI[, 1] + HI[, 2] + HI[, 3] + HI[, 4] + HI[, 5] + HI[, 6] + HI[, 7]
	N = NA_s + NA_as + NB_s + NB_as
		N = N[, 1] + N[, 2] + N[, 3] + N[, 4] + N[, 5] + N[, 6] + N[, 7]
	HI = log2(HI/N)
	
	return(list(RASav=RASav, HI=HI))
} # MPAMindiv

##############################################################################################
#######                  Calculate AF of Genotype/intensity based data
##############################################################################################

Feature_RData_ft = function(inpath, source_name, Ind_Array, NA_str, ParaSet, X_label = 23, chrlist, Call_level=c("AA", "AB", "BB", "NoCall"), outpath, chiptype_choice, onechip_choice, logfile){
#	Ind_Array = 1: AFFY_100K, 2: AFFY_500K, 3: AFFY_Array6, 4: AFFY_Axiom, 5: Illumina
#	Ind_Datatype: 1. Genotype/Intensity-based, 2. AF-based, 3. RData-based, 4. CEL-based
	dir.create(paste(outpath, "TEMP/", source_name, sep = ""), showWarnings = F)
	
	outpath_indiv = paste(outpath, "/TEMP/", source_name,"/TEMP_Indiv", sep="")
	dir.create(outpath_indiv, showWarnings = F)
	
	source_name_G = paste(outpath, "TEMP/", source_name, "/", source_name, "_GAF", sep = "")
	dir.create(source_name_G, showWarnings = F)
	
	source_name_I = paste(outpath, "TEMP/", source_name, "/", source_name, "_ind", sep = "")
	dir.create(source_name_I, showWarnings = F)

# source_name = nowGname
# NA_str = NA_SNP
# ParaSet = ParaSet_SNP

	# dir.create(paste(outpath, "/TEMP_Indiv/", source_name, sep=""), showWarnings=F, recursive=T)

	##################################################################################################
	##
	##			Genotype 
	##
	##################################################################################################
	Gpath = paste(inpath, "IndGeno", sep="/");
	# cat(Gpath, sep = "\n")
	filename_g = mixedsort(dir(Gpath))
	i_len = ifelse(Ind_Array>=3, length(filename_g), length(filename_g)/2)
	if(i_len==0)		print("*** Error message: There is no files in \"IndGeno\" directory.")	##	IndGeno must have files
	if(i_len%%1!=0)	print("*** Error message: The number of genotype files are not complete for each individual.")	
	samplelist = data.frame(Index=1:i_len, Sample_ID=NA, Gender=NA, stringsAsFactors=F)

	##
	##	Parameter setting for APT (apt-cel-extract) output format	
	##
		##	ParaSet_Data = c(SkipRow_Data, SNPcol_Data, Chrcol_Data, Posicol_Data, Callcol_Data, PIAcol_Data, PIBcol_Data);
	SkipRow = ParaSet[1];		SNPcol = ParaSet[2];		Chrcol = ParaSet[3];			Posicol = ParaSet[4];		
	CallAcol = ParaSet[5];		CallBcol = ParaSet[6];		PIAcol = ParaSet[7];			PIBcol = ParaSet[8];			BAFcol = ParaSet[9]

	if(Ind_Array==1 | Ind_Array==2)	SkipRow = ifelse(is.na(SkipRow), 2, SkipRow)
	
	##
	##	colClasses and col.names setting
	##				
	Gdata1 = read.table(paste(Gpath, filename_g[1], sep="/"), skip=SkipRow, header=F, sep="\t", nrow=1, na.strings=NA_str, colClasses="character")
	g_ColClasses = g_Colnames = rep("NULL", ncol(Gdata1))
	g_ColClasses[SNPcol] = "character";		g_Colnames[SNPcol] = "SNP_ID"
	if(Ind_Array!=4)		g_ColClasses[Chrcol] = "character";		g_Colnames[Chrcol] = "Chr"
	if(Ind_Array!=4)		g_ColClasses[Posicol] = "integer";		g_Colnames[Posicol] = "PhyPosi"
	g_ColClasses[CallAcol] = "character";		g_Colnames[CallAcol] = "Genotype"

	# if(is.na(Chrcol) & is.na(Posicol)){
		# stopifnot(file.exists(Annopath))
		Anno = read.table(dir(paste(inpath, "Annotation", sep = "/"), full.names = T), header=F, skip=1, col.names=c("SNP_ID", "Chr", "PhyPosi"), colClasses=c("character", "character", "integer"), sep="\t", as.is=T)

		##	Chanage chr label 
		if("X" %in% unique(Anno$Chr)){		Anno$Chr[Anno$Chr=="X"] = X_label} # X_label = 23
		if("XY" %in% unique(Anno$Chr)){		Anno$Chr[Anno$Chr=="XY"] = X_label}
		# if("Y" %in% unique(Anno$Chr)){		Anno$Chr[Anno$Chr=="Y"] = 24}

		##	Remove SNPs on MT chr or unknown chr
		Anno = Anno[(Anno$Chr!="MT") & (Anno$Chr!="M") & (!is.na(Anno$Chr)) & (Anno$Chr!="") & (Anno$Chr!="Y") & (Anno$Chr!=24), ]
		Anno$Chr = as.integer(Anno$Chr)
		
		##	Remove control SNPs
		Anno = Anno[substring(Anno$SNP_ID, 1, 4)!="AFFX", ]
	# }
	
	cat("			==========                Raw AF calculation                ==========\n")
	cat("   ~ The calculation of unadjust allele frequency is preparing.\n", file=logfile, append=T)
	##################################################################################################
	##	 	1: AFFY_100K, 2: AFFY_500K
	##
	if(Ind_Array==1 | Ind_Array==2){		
		rm("CallBcol", "PIAcol", "PIBcol", "BAFcol")
		p_ColClasses = c("NULL", "character", rep("numeric", 56))
		p_Colnames = c("NULL", "SNP_ID", "Q1_PA_s", "Q1_PB_s", "Q1_MA_s", "Q1_MB_s", "Q1_PA_as", "Q1_PB_as", "Q1_MA_as", "Q1_MB_as", "Q2_PA_s", "Q2_PB_s", "Q2_MA_s", "Q2_MB_s", "Q2_PA_as", "Q2_PB_as", "Q2_MA_as", "Q2_MB_as", "Q3_PA_s", "Q3_PB_s", "Q3_MA_s", "Q3_MB_s", "Q3_PA_as", "Q3_PB_as", "Q3_MA_as", "Q3_MB_as", "Q4_PA_s", "Q4_PB_s", "Q4_MA_s", "Q4_MB_s", "Q4_PA_as", "Q4_PB_as", "Q4_MA_as", "Q4_MB_as", "Q5_PA_s", "Q5_PB_s", "Q5_MA_s", "Q5_MB_s", "Q5_PA_as", "Q5_PB_as", "Q5_MA_as", "Q5_MB_as", "Q6_PA_s", "Q6_PB_s", "Q6_MA_s", "Q6_MB_s", "Q6_PA_as", "Q6_PB_as", "Q6_MA_as", "Q6_MB_as", "Q7_PA_s", "Q7_PB_s", "Q7_MA_s", "Q7_MB_s", "Q7_PA_as", "Q7_PB_as", "Q7_MA_as", "Q7_MB_as")	
			
			Ppath = paste(inpath, "IndPI", sep="/");		filename_p = dir(Ppath)
			i_len_p = ifelse(Ind_Array>=3, length(filename_p), length(filename_p)/2)
		
			
			g_ColClasses[c(Chrcol, Posicol)] = c("character", "integer")
			if(chiptype_choice == "2chip"){
				for(i in 1:i_len){
					cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is prepareing.\n",file=logfile,append=T)
					if(i==1){
						cat("\n*** The RData transformation of ", nrow(samplelist), " samples of ", source_name, " data is starting.\n", sep="")
						cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
					}
					if(i%%10==0)	cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
					
					samplelist[i, "Sample_ID"] = strsplit(filename_g[2*i-1], "_|.txt|.CEL")[[1]][1]				
					
					##
					##	Import data
					##	
					##	Genotype data
					Gdata1 = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=c("", NA_str), fill = T)	
					Gdata1$Chiptype = "chip1"
					##
					##	Check column content
					##
					cri1 = !any(grep("SNP_A-", Gdata1$SNP_ID))
					cri2 = !all(c(1:22, "X") %in% Gdata1$Chr)
					cri3 = !is.integer(Gdata1$PhyPosi)
					cri4 = !all(c("AA", "AB", "BB") %in% Gdata1$Genotype)
					if(any(c(cri1, cri2, cri3, cri4)==TRUE)){
						Gdata1 = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow-1, header=T, sep="\t", as.is=T, na.strings=NA_str, nrow=1000, fill = T)	

						g_ColClasses_mdy = g_Colnames_mdy = rep("NULL", ncol(Gdata1))
						for(jj in 1:ncol(Gdata1)){
							if(any(grep("SNP_A-", Gdata1[, jj]))){													g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "SNP_ID"
							}else if(length(grep("CHR", toupper(colnames(Gdata1)[jj])))!=0){								g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Chr"
							}else if(is.integer(Gdata1[, jj]) & (length(grep("POSI", toupper(colnames(Gdata1)[jj])))!=0)){	g_ColClasses_mdy[jj] = "integer"
							g_Colnames_mdy[jj] = "PhyPosi"
							}else if(("AA" %in% Gdata1[, jj]) | ("AB" %in% Gdata1[, jj]) | ("BB" %in% Gdata1[, jj])){g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Genotype"}
						}
						
						Gdata1 = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)
						Gdata1$Chiptype = "chip1"
						Gdata2 = read.table(paste(Gpath, filename_g[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)
						Gdata2$Chiptype = "chip2"
					}else{
						Gdata2 = read.table(paste(Gpath, filename_g[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=NA_str, fill = T)	
						Gdata2$Chiptype = "chip2"
					}
					Gdata2$Chiptype = "chip2"
					Gdata = rbind(Gdata1, Gdata2);		rm("Gdata1", "Gdata2")

					##	Intensity data
					Pdata1 = read.table(paste(Ppath, filename_p[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=p_ColClasses, col.names=p_Colnames, as.is=T, na.strings="null", fill = T)	
					Pdata2 = read.table(paste(Ppath, filename_p[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=p_ColClasses, col.names=p_Colnames, as.is=T, na.strings="null", fill = T)	
					Pdata = rbind(Pdata1, Pdata2);		rm("Pdata1", "Pdata2");		gc()
					
					ft_out = MPAMindiv(Pdata)
					Pdata = data.frame(SNP_ID=Pdata[, 1], AF=ft_out$RASav, stringsAsFactors=F);		rm("ft_out");		gc()

			

					##
					##	Individual data
					##
					pp = match(c("Chr", "PhyPosi"), colnames(Gdata))
					pp = pp[!is.na(pp)]
					if(length(pp) > 0){
						Out = merge(Anno, Gdata[, -pp], by="SNP_ID", all.x=T, sort=F)
					}else{
						Out = merge(Anno, Gdata, by="SNP_ID", all.x=T, sort=F)
					}
					Out = merge(Out, Pdata, by="SNP_ID", all.x=T, sort=F)
					
					Out$QC_rm_index = 0					
					posi = !is.element(Out$SNP_ID, Gdata$SNP_ID)
					Out$QC_rm_index[posi] = 1		
					Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "AF", "Genotype", "Chiptype", "QC_rm_index")]
	
					save(Out, file=paste(outpath_indiv, "/", samplelist[i, "Sample_ID"], ".RData", sep=""));
					cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is finished.\n",file=logfile,append=T)
					cat("\n",file=logfile,append=T)
				} # for(i in 1:i_len){
			}else if(chiptype_choice=="1chip"){
				if(onechip_choice =="chip1"){
					for(i in 1:i_len){
						cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is prepareing.\n",file=logfile,append=T)
						if(i==1){
							cat("\n*** The RData transformation of ", nrow(samplelist), " samples of ", source_name, " data is starting.\n", sep="")
							cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
						}
						if(i%%10==0)	cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
						
						samplelist[i, "Sample_ID"] = strsplit(filename_g[2*i-1], "_|.txt|.CEL")[[1]][1]				
						
						##
						##	Import data
						##	
						##	Genotype data
						Gdata = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=c("", NA_str), fill = T)	
						Gdata$Chiptype = "chip1"
						##
						##	Check column content
						##
						cri1 = !any(grep("SNP_A-", Gdata$SNP_ID))
						cri2 = !all(c(1:22, "X") %in% Gdata$Chr)
						cri3 = !is.integer(Gdata$PhyPosi)
						cri4 = !all(c("AA", "AB", "BB") %in% Gdata$Genotype)
						if(any(c(cri1, cri2, cri3, cri4)==TRUE)){
							Gdata = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow-1, header=T, sep="\t", as.is=T, na.strings=NA_str, nrow=1000, fill = T)	

							g_ColClasses_mdy = g_Colnames_mdy = rep("NULL", ncol(Gdata))
							for(jj in 1:ncol(Gdata)){
								if(any(grep("SNP_A-", Gdata[, jj]))){													g_ColClasses_mdy[jj] = "character"
								g_Colnames_mdy[jj] = "SNP_ID"
								}else if(length(grep("CHR", toupper(colnames(Gdata)[jj])))!=0){								g_ColClasses_mdy[jj] = "character"
								g_Colnames_mdy[jj] = "Chr"
								}else if(is.integer(Gdata[, jj]) & (length(grep("POSI", toupper(colnames(Gdata)[jj])))!=0)){	g_ColClasses_mdy[jj] = "integer"
								g_Colnames_mdy[jj] = "PhyPosi"
								}else if(("AA" %in% Gdata[, jj]) | ("AB" %in% Gdata[, jj]) | ("BB" %in% Gdata[, jj])){g_ColClasses_mdy[jj] = "character"
								g_Colnames_mdy[jj] = "Genotype"}
							}
							
							Gdata = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)
							Gdata$Chiptype = "chip1"
						}
						# Gdata2 = read.table(paste(Gpath, filename_g[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=NA_str)	
						# Gdata2$Chiptype = "chip2"
						# Gdata = Gdata1;		rm("Gdata1")

						##	Intensity data
						Pdata = read.table(paste(Ppath, filename_p[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=p_ColClasses, col.names=p_Colnames, as.is=T, na.strings="null", fill = T)	
						# Pdata2 = read.table(paste(Ppath, filename_p[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=p_ColClasses, col.names=p_Colnames, as.is=T, na.strings="null")	
						# Pdata = rbind(Pdata1, Pdata2);		rm("Pdata1", "Pdata2");		gc()
						
						ft_out = MPAMindiv(Pdata)
						Pdata = data.frame(SNP_ID=Pdata[, 1], AF=ft_out$RASav, stringsAsFactors=F);		rm("ft_out");		gc()
						

						##
						##	Individual data
						##
						pp = match(c("Chr", "PhyPosi"), colnames(Gdata))
						pp = pp[!is.na(pp)]
						if(length(pp) > 0){
							Out = merge(Anno, Gdata[, -pp], by="SNP_ID", all.x=T, sort=F)
						}else{
							Out = merge(Anno, Gdata, by="SNP_ID", all.x=T, sort=F)
						}
						Out = merge(Out, Pdata, by="SNP_ID", all.x=T, sort=F)
						Out$QC_rm_index = 0					
						posi = !is.element(Out$SNP_ID, Gdata$SNP_ID)
						Out$QC_rm_index[posi] = 1		
						Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "AF", "Genotype", "Chiptype", "QC_rm_index")]

						save(Out, file=paste(outpath_indiv, "/", samplelist[i, "Sample_ID"], ".RData", sep=""));  		
						cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is finished.\n",file=logfile,append=T)
						cat("\n",file=logfile,append=T)						
					} # for(i in 1:i_len){
				}else if(onechip_choice == "chip2"){
					for(i in 1:i_len){
						cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is prepareing.\n",file=logfile,append=T)
						if(i==1){
							cat("\n*** The RData transformation of ", nrow(samplelist), " samples of ", source_name, " data is starting.\n", sep="")
							cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
						}
						if(i%%10==0)	cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
						
						samplelist[i, "Sample_ID"] = strsplit(filename_g[2*i-1], "_|.txt|.CEL")[[1]][1]				
						
						##
						##	Import data
						##	
						##	Genotype data
						# Gdata1 = read.table(paste(Gpath, filename_g[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=c("", NA_str))	
							# Gdata1$Chiptype = "chip1"
						Gdata = read.table(paste(Gpath, filename_g[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=NA_str, fill = T)	
							Gdata$Chiptype = "chip2"
						##
						##	Check column content
						##
						cri1 = !any(grep("SNP_A-", Gdata$SNP_ID))
						cri2 = !all(c(1:22, "X") %in% Gdata$Chr)
						cri3 = !is.integer(Gdata$PhyPosi)
						cri4 = !all(c("AA", "AB", "BB") %in% Gdata$Genotype)
						if(any(c(cri1, cri2, cri3, cri4)==TRUE)){
							Gdata = read.table(paste(Gpath, filename_g[2*i - 0], sep="/"), skip=SkipRow-1, header=T, sep="\t", as.is=T, na.strings=NA_str, nrow=, fill = T)	

							g_ColClasses_mdy = g_Colnames_mdy = rep("NULL", ncol(Gdata))
							for(jj in 1:ncol(Gdata)){
								if(any(grep("SNP_A-", Gdata[, jj]))){													g_ColClasses_mdy[jj] = "character"
								g_Colnames_mdy[jj] = "SNP_ID"
								}else if(length(grep("CHR", toupper(colnames(Gdata)[jj])))!=0){								g_ColClasses_mdy[jj] = "character"
								g_Colnames_mdy[jj] = "Chr"
								}else if(is.integer(Gdata[, jj]) & (length(grep("POSI", toupper(colnames(Gdata)[jj])))!=0)){	g_ColClasses_mdy[jj] = "integer"
								g_Colnames_mdy[jj] = "PhyPosi"
								}else if(("AA" %in% Gdata[, jj]) | ("AB" %in% Gdata[, jj]) | ("BB" %in% Gdata[, jj])){g_ColClasses_mdy[jj] = "character"
								g_Colnames_mdy[jj] = "Genotype"}
							}
							
							Gdata = read.table(paste(Gpath, filename_g[2*i - 0], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)
							Gdata$Chiptype = "chip1"
						}
						# Gdata = rbind(Gdata1, Gdata2);		rm("Gdata1", "Gdata2")

						##	Intensity data
						# Pdata1 = read.table(paste(Ppath, filename_p[2*i - 1], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=p_ColClasses, col.names=p_Colnames, as.is=T, na.strings="null")	
						Pdata = data.frame(SNP_ID=Pdata[, 1], AF=ft_out$RASav, stringsAsFactors=F);		rm("ft_out");		gc()
						

						##
						##	Individual data
						##
						pp = match(c("Chr", "PhyPosi"), colnames(Gdata))
						pp = pp[!is.na(pp)]
						if(length(pp) > 0){
							Out = merge(Anno, Gdata[, -pp], by="SNP_ID", all.x=T, sort=F)
						}else{
							Out = merge(Anno, Gdata, by="SNP_ID", all.x=T, sort=F)
						}
						Out = merge(Out, Pdata, by="SNP_ID", all.x=T, sort=F)
						Out$QC_rm_index = 0					
						posi = !is.element(Out$SNP_ID, Gdata$SNP_ID)
						Out$QC_rm_index[posi] = 1		
						Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "AF", "Genotype", "Chiptype", "QC_rm_index")]

						save(Out, file=paste(outpath_indiv, "/", samplelist[i, "Sample_ID"], ".RData", sep=""));  	
						cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is finished.\n",file=logfile,append=T)
						cat("\n",file=logfile,append=T)								
					} # for(i in 1:i_len){
				}
			}
			
			
		##################################################################################################
		## 		3: AFFY_Array6 
		##
		}else if(Ind_Array==3){				
			
			g_ColClasses[PIAcol] = "numeric";		g_Colnames[PIAcol] = "PI_A"
			g_ColClasses[PIBcol] = "numeric";		g_Colnames[PIBcol] = "PI_B"			
			for(i in 1:i_len){
				cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is prepareing.\n",file=logfile,append=T)
				if(i==1){
					cat("\n*** The RData transformation of ", nrow(samplelist), " samples of ", source_name, " data is starting.\n", sep="")
					cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
				}
				if(i%%10==0)	cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
	
				samplelist[i, "Sample_ID"] = strsplit(filename_g[i], "_|.txt|.CEL")[[1]][1]				
	
				##
				##	Import data
				##		
				data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=NA_str, fill = T)										
				##
				##	Check column content
				##
				cri1 = !any(grep("SNP_A-", data$SNP_ID))
				cri2 = !all(c(1:22, "X") %in% data$Chr)
				cri3 = !is.integer(data$PhyPosi)
				cri4 = !all(c("AA", "AB", "BB") %in% data$Genotype)
				cri5 = !is.numeric(data$PI_A)
				cri6 = !is.numeric(data$PI_B);
				if(any(c(cri1, cri2, cri3, cri4, cri5, cri6)==TRUE)){
					data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow-1, header=T, sep="\t", as.is=T, na.strings=NA_str, nrow=1000, fill = T)	

					g_ColClasses_mdy = g_Colnames_mdy = rep("NULL", ncol(data))
					for(jj in 1:ncol(data)){
						if(any(grep("SNP_A-", data[, jj]))){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "SNP_ID"
						}else if(length(grep("CHR", toupper(colnames(data)[jj])))!=0){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Chr"
						}else if(is.integer(data[, jj]) & (length(grep("POSI", toupper(colnames(data)[jj])))!=0)){
							g_ColClasses_mdy[jj] = "integer"
							g_Colnames_mdy[jj] = "PhyPosi"
						}else if(("AA" %in% data[, jj]) | ("AB" %in% data[, jj]) | ("BB" %in% data[, jj])){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Genotype"
						}else if(is.numeric(data[, jj]) & (length(grep("A", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "numeric"
							g_Colnames_mdy[jj] = "PI_A"
						}else if(is.numeric(data[, jj]) & (length(grep("B", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "numeric"
							g_Colnames_mdy[jj] = "PI_B"
						}
					}
					data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)
				}
				
				##
				##	Match onto Annotation
				##
				pp = match(c("Chr", "PhyPosi"), colnames(data))
				pp = pp[!is.na(pp)]
				if(length(pp) > 0){
					Out = merge(Anno, data[, -pp], by="SNP_ID", all.x=T, sort=F)
				}else{
					Out = merge(Anno, data, by="SNP_ID", all.x=T, sort=F)
				}

				##
				##	Individual data
				##
				Out$Chiptype = "chip1"
				Out$AF = Out$PI_A/(Out$PI_A + Out$PI_B)
				Out$QC_rm_index = 0					
				posi = !is.element(Out$SNP_ID, data$SNP_ID)
				Out$QC_rm_index[posi] = 1

				
				Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "AF", "Genotype", "Chiptype", "QC_rm_index")]
				# ##	Chr label
				# if("X" %in% unique(Out$Chr)){		Out$Chr[Out$Chr=="X"] = X_label} # X_label = 23
				# if("XY" %in% unique(Out$Chr)){		Out$Chr[Out$Chr=="XY"] = X_label}
				# # if("Y" %in% unique(Out$Chr)){		Out$Chr[Out$Chr=="Y"] = Y_label} # Y_label = 25					
				# # if(24 %in% unique(Out$Chr)){		Out$Chr[Out$Chr==24] = Y_label} # Y_label = 25					
				save(Out, file=paste(outpath_indiv, "/", samplelist[i, "Sample_ID"], ".RData", sep=""));
				cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is finished.\n",file=logfile,append=T)
				cat("\n",file=logfile,append=T)
			} # for(i in 1:i_len){
			
		##################################################################################################
		##		4: AFFY_Axiom (from GTC)
		##
		}else if(Ind_Array==4){				
			
			g_ColClasses[PIAcol] = "numeric";		g_Colnames[PIAcol] = "Log2Ratio"
			g_ColClasses[PIBcol] = "numeric";		g_Colnames[PIBcol] = "Strength"	

			
			for(i in 1:i_len){
				cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is prepareing.\n",file=logfile,append=T)
				if(i==1){				
					cat("\n*** The RData transformation of ", nrow(samplelist), " samples of ", source_name, " data is starting.\n", sep="")
					cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
				}
				if(i%%10==0)	cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")

				samplelist[i, "Sample_ID"] = strsplit(filename_g[i], "_|.txt|.CEL")[[1]][1]				

				##
				##	Import data
				##		
				data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=NA_str, fill = T)
				
				##
				##	Check column content
				##
				cri1 = !((any(grep("AX-", data$SNP_ID)) | any(grep("AFFX-", data$SNP_ID)) | any(grep("rs", data$SNP_ID)) | any(grep("b36", data$SNP_ID))))
				cri4 = !all(c("AA", "AB", "BB") %in% data$Genotype)
				cri5 = !is.numeric(data$Log2Ratio)
				cri6 = !is.numeric(data$Strength);
				if(any(c(cri1, cri4, cri5, cri6)==TRUE)){
					data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow-1, header=T, sep="\t", as.is=T, na.strings=NA_str, nrow=1000, fill = T)	

					g_ColClasses_mdy = g_Colnames_mdy = rep("NULL", ncol(data))
					for(jj in 1:ncol(data)){
						if(any(grep("AX-", data[, jj])) | any(grep("AFFX-", data[, jj])) | any(grep("rs", data[, jj])) | any(grep("b36", data[, jj]))){						
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "SNP_ID"
						}else if(length(grep("CHR", toupper(colnames(data)[jj])))!=0){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Chr"
						}else if(is.integer(data[, jj]) & (length(grep("POSI", toupper(colnames(data)[jj])))!=0)){
							g_ColClasses_mdy[jj] = "integer"
							g_Colnames_mdy[jj] = "PhyPosi"
						}else if(("AA" %in% data[, jj]) | ("AB" %in% data[, jj]) | ("BB" %in% data[, jj])){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Genotype"
						}else if(is.numeric(data[, jj]) & (length(grep("Ratio", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "numeric"
							g_Colnames_mdy[jj] = "Log2Ratio"
						}else if(is.numeric(data[, jj]) & (length(grep("Strength", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "numeric"
							g_Colnames_mdy[jj] = "Strength"	
						}
					}
					data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)									
				}
				

				##
				##	Intensity calculation when data is from GTC
				##
				##			A=2^((T+0.5*R)) and B=2^((T-0.5*R)), where Strength (T) and Log Ratio (R)
				LR = data[, "Log2Ratio"]
				ST = data[, "Strength"]	
				
				data$PI_A = 2^((ST + 0.5*LR))
				data$PI_B = 2^((ST - 0.5*LR))

				##
				##	Match onto Annotation
				##
				pp = match(c("Chr", "PhyPosi"), colnames(data))
				pp = pp[!is.na(pp)]
				if(length(pp) > 0){
					Out = merge(Anno, data[, -pp], by="SNP_ID", all.x=T, sort=F)
				}else{
					Out = merge(Anno, data, by="SNP_ID", all.x=T, sort=F)
				}
				Out$AF = Out$PI_A/(Out$PI_A + Out$PI_B)
				Out$Chiptype = "chip1"
				Out$QC_rm_index = 0					
				posi = !is.element(Out$SNP_ID, data$SNP_ID)
				Out$QC_rm_index[posi] = 1

				
				Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "AF", "Genotype", "Chiptype", "QC_rm_index")]

				save(Out, file=paste(outpath_indiv, "/", samplelist[i, "Sample_ID"], ".RData", sep=""));
				cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is finished.\n",file=logfile,append=T)
				cat("\n",file=logfile,append=T)
			} # for(i in 1:i_len){
			

		##################################################################################################
		## 		5: Illumina
		##
		}else if(Ind_Array==5){				

			if(!is.na(CallBcol)){
															g_Colnames[CallAcol] = "Genotype1"
				g_ColClasses[CallBcol] = "character";		g_Colnames[CallBcol] = "Genotype2"
			}			
			g_ColClasses[PIAcol] = "numeric";		g_Colnames[PIAcol] = "PI_A"
			g_ColClasses[PIBcol] = "numeric";		g_Colnames[PIBcol] = "PI_B"
			g_ColClasses[BAFcol] = "numeric";		g_Colnames[BAFcol] = "BAF"
			
			for(i in 1:i_len){
				cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is prepareing.\n",file=logfile,append=T)
				if(i==1){
					cat("\n*** The RData transformation of ", nrow(samplelist), " samples of ", source_name, " data is starting.\n", sep="")
					cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")
				}
				if(i%%10==0)	cat("*** The RData transformation of ", i, "th samples of ", source_name, " data is starting.\n", sep="")

				samplelist[i, "Sample_ID"] = strsplit(filename_g[i], "_|.txt|.CEL")[[1]][1]				

				##
				##	Import data
				##		
				if(i==1)	g_ColClasses[c(Chrcol, Posicol)] = c("character", "integer")
				data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses, col.names=g_Colnames, as.is=T, na.strings=NA_str, fill = T)
				
				##
				##	Check column content
				##
				cri1 = !any(grep("rs", data$SNP_ID))
				cri2 = !all(c(1:22, "X") %in% data$Chr)
				cri3 = !is.integer(data$PhyPosi)
				cri4 = !all(c("AA", "AB", "BB") %in% data$Genotype)
				cri5 = !is.numeric(data$PI_A)
				cri6 = !is.numeric(data$PI_B);
				
				if(!is.na(CallBcol)){
					cri7 = !all(c("A", "B") %in% data$Genotype1)
					cri8 = !all(c("A", "B") %in% data$Genotype2)
				}else{
					cri7 = FALSE
					cri8 = FALSE
				}			
				if(!is.na(BAFcol)){
					cri9 = !is.numeric(data$BAF)
				}else{
					cri9 = FALSE
				}
				
				if(any(c(cri1, cri2, cri3, cri4, cri5, cri6, cri7, cri8, cri9)==TRUE)){
					data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow-1, header=T, sep="\t", as.is=T, na.strings=NA_str, nrow=1000, fill = T)	

					g_ColClasses_mdy = g_Colnames_mdy = rep("NULL", ncol(data))
					for(jj in 1:ncol(data)){
						if(any(grep("rs", data[, jj]))){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "SNP_ID"
						}else if(length(grep("CHR", toupper(colnames(data)[jj])))!=0){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Chr"
						}else if(is.integer(data[, jj]) & (length(grep("POSI", toupper(colnames(data)[jj])))!=0)){
							g_ColClasses_mdy[jj] = "integer"
							g_Colnames_mdy[jj] = "PhyPosi"
						}else if(("AA" %in% data[, jj]) | ("AB" %in% data[, jj]) | ("BB" %in% data[, jj])){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Genotype"
						}else if((("A" %in% data[, jj]) | ("B" %in% data[, jj])) & (length(grep("1...AB", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Genotype1"
						}else if((("A" %in% data[, jj]) | ("B" %in% data[, jj])) & (length(grep("2...AB", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "character"
							g_Colnames_mdy[jj] = "Genotype2"
						}else if(is.numeric(data[, jj]) & (length(grep("Freq", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "numeric"
							g_Colnames_mdy[jj] = "BAF"
						}else if(is.integer(data[, jj]) & (length(grep("X", colnames(data)[jj]))!=0) & (length(grep("Raw", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "integer"
							g_Colnames_mdy[jj] = "PI_A"	
						}else if(is.integer(data[, jj]) & (length(grep("Y", colnames(data)[jj]))!=0) & (length(grep("Raw", colnames(data)[jj]))!=0)){
							g_ColClasses_mdy[jj] = "integer"
							g_Colnames_mdy[jj] = "PI_B"
						}
					}
					
					data = read.table(paste(Gpath, filename_g[i], sep="/"), skip=SkipRow, header=F, sep="\t", colClasses=g_ColClasses_mdy, col.names=g_Colnames_mdy, as.is=T, na.strings=NA_str, fill = T)										
				}
				
				
				##	Combine genotype
				if(!is.na(CallBcol)){
															g_Colnames[CallAcol] = "Genotype1"
					g_ColClasses[CallBcol] = "character";		g_Colnames[CallBcol] = "Genotype2"					
					data[, "Genotype"] = paste(data[, "Genotype1"], data[, "Genotype2"], sep="")
					##	Replace the genotype with NA in any allele
						posi = is.na(data[, "Genotype1"]) | is.na(data[, "Genotype2"]);			data[posi, "Genotype"] = "NoCall"
					data = data[, -match(c("Genotype1", "Genotype2"), colnames(data))]
				}

				##
				##	Individual data
				##		
				data$Chiptype = "chip1"
				if(!is.na(BAFcol)){
					data$AF = 1 - data[, "BAF"]
				}else{
					data$AF = data$PI_A/(data$PI_A + data$PI_B)
				}
				data$Call = data[, "Genotype"]
				
				pp = match(c("Chr", "PhyPosi"), colnames(data))
				pp = pp[!is.na(pp)]
				if(length(pp) > 0){
					Out = merge(Anno, data[, -pp], by="SNP_ID", all.x=T, sort=F)
				}else{
					Out = merge(Anno, data, by="SNP_ID", all.x=T, sort=F)
				}
				
				Out$QC_rm_index = 0					
				posi = !is.element(Out$SNP_ID, data$SNP_ID)
				Out$QC_rm_index[posi] = 1

				
				Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "AF", "Genotype", "Chiptype", "QC_rm_index")]

				save(Out, file=paste(outpath_indiv, "/", samplelist[i, "Sample_ID"], ".RData", sep=""));
				cat(paste("     ~  The unadjust allele frequency of the individual ",i,sep="")," is finished.\n",file=logfile,append=T)
				cat("\n",file=logfile,append=T)
				
			} # for(i in 1:i_len){
			
		} # end for "if(Ind_Array==1 | Ind_Array==2){"		
		
	cat("*** The RData transformation of ", i_len, " samples of ", source_name, " data is finished.\n", sep="")
		
	return(invisible(NULL))
	
} # Feature_RData_ft



Feature_Stat_Fast_ft = function(source_name, Ind_Array, X_label = 23, chrlist, Call_level=c("AA", "AB", "BB", "NoCall"), inpath, outpath){
	#	Ind_Array = 1: AFFY_100K, 2: AFFY_500K, 3: AFFY_Array6, 4: AFFY_Axiom, 5: Illumina
	#	Ind_Datatype = 1. Genotype/Intensity-based, 2. AF-based, 3. RData-based, 4. CEL-based
	# source_name = nowGname
	dir.create(paste(outpath, "TEMP/", source_name, sep = ""), showWarnings = F)
	source_name_G = paste(outpath, "TEMP/", source_name, "/", source_name, "_GAF", sep = "")
	source_name_I = paste(outpath, "TEMP/", source_name, "/", source_name, "_ind", sep = "")
	dir.create(source_name_G, showWarnings = F)
	dir.create(source_name_I, showWarnings = F)
	######c############################################################################################
	##
	##			Data pre-processing 
	##
	##################################################################################################	
	RData_inpath = paste(outpath, "TEMP/", source_name, "/TEMP_Indiv", sep="")
	# dir.create(paste(outpath, "/TEMP/", source_name, sep=""), showWarnings=F, recursive=T)

	filename_g = mixedsort(dir(RData_inpath));			
	i_len = length(filename_g)
	samplelist = data.frame(Index=1:i_len, Sample_ID=NA, Gender=NA, stringsAsFactors=F)		
	
	Anno = read.table(dir(paste(inpath, "Annotation", sep = "/"), full.names = T), header=F, skip=1, col.names=c("SNP_ID", "Chr", "PhyPosi"), colClasses=c("character", "character", "integer"), sep="\t", as.is=T)

	##	Chanage chr label 
	if("X" %in% unique(Anno$Chr)){		Anno$Chr[Anno$Chr=="X"] = X_label} # X_label = 23
	if("XY" %in% unique(Anno$Chr)){		Anno$Chr[Anno$Chr=="XY"] = X_label}
	# if("Y" %in% unique(Anno$Chr)){		Anno$Chr[Anno$Chr=="Y"] = 24}

	##	Remove SNPs on MT chr or unknown chr
	Anno = Anno[(Anno$Chr!="MT") & (Anno$Chr!="M") & (!is.na(Anno$Chr)) & (Anno$Chr!="") & (Anno$Chr!="Y") & (Anno$Chr!=24), ]
	Anno$Chr = as.integer(Anno$Chr)
	
	##	Remove control SNPs
	Anno = Anno[substring(Anno$SNP_ID, 1, 4)!="AFFX", ]
	##
	##	Chromosomal 
	##
	for(chr in chrlist){
		chrname = ifelse(chr<10, paste("0", chr, sep=""), chr)

		##
		##	Individual 
		##
		for(i in 1:i_len){

			if(i==1){
				cat("\n*** Chr ", chrname, ": Chrosomal data combination of ", nrow(samplelist), " samples of ", source_name, " is starting.\n", sep="")
				cat("*** Chr ", chrname, ": Chrosomal data combination of ", i, "th samples of ", source_name, " is starting.\n", sep="")
				
				##	Pre-allocate
				Anno_chr = Anno[(Anno$Chr==chr) & (!is.na(Anno$Chr)), ]		
				Anno_chr = Anno_chr[order(Anno_chr$PhyPosi), ]
			}
			if(i%%10==0)	cat("*** Chr ", chrname, ": Chrosomal data combination of ", i, "th samples of ", source_name, " is starting.\n", sep="")
			
			##	Load individual RData: Out
			load(file=paste(RData_inpath, "/", filename_g[i], sep=""))
			if(chr==1){
					##	Chr label
					if("X" %in% unique(Out$Chr)){		Out$Chr[Out$Chr=="X"] = X_label} # X_label = 23
					if("XY" %in% unique(Out$Chr)){		Out$Chr[Out$Chr=="XY"] = X_label}
					# if("Y" %in% unique(Out$Chr)){		Out$Chr[Out$Chr=="Y"] = Y_label} # Y_label = 25					
					# if(24 %in% unique(Out$Chr)){		Out$Chr[Out$Chr==24] = Y_label} # Y_label = 25					
					
					# if(!(("QC_rm_index" %in% colnames(Out))))	Out$QC_rm_index = 0
					
					##
					##	Save RData
					##
					# if(Ind_Array<=2){		Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "UAF", "PI_A", "PI_B", "HI_raw", "Genotype", "ChipIndex", "QC_rm_index", "HI_QN")]
					# }else{					Out = Out[, c("SNP_ID", "Chr", "PhyPosi", "UAF", "PI_A", "PI_B", "Genotype", "ChipIndex", "QC_rm_index", "HI_QN")]	}
					# save(Out, file=paste(RData_outpath, "/", filename_g[i], sep=""));  			
			}
			sample_ID = strsplit(filename_g[i], ".RData")[[1]][1]
			##
			##	Remove SNPs on MT chr or unknown chr
			##
			chrposi = as.numeric(na.omit(match(Anno_chr$SNP_ID, Out$SNP_ID[Out$QC_rm_index==0])))
			chip = Out[chrposi, c("SNP_ID", "PhyPosi", "Chiptype")]
			data_chr = Out[chrposi, ]
			
			##	Pre-allocate
			if(i==1)		AF = geno = matrix(NA, nrow=nrow(chip), ncol=nrow(samplelist))		

			##
			##	Chr data
			##
			AF[, i] = data_chr[, "AF"]
			geno[, i] = data_chr[, "Genotype"]

			##
			##	Gender detection
			##			
			if(chr==X_label){
				samplelist[i, "Sample_ID"] = sample_ID

				gg = table(factor(geno[, i], levels=Call_level))
				cc = gg[2]/sum(gg[1:3])
				samplelist[i, "Gender"] = ifelse(cc <= 0.1, "M", "F")
			}
			
			ind_chr = data.frame(chip[,c("SNP_ID", "PhyPosi", "Chiptype")], Genotype=geno[, i], AF = AF[, i])
			colnames(ind_chr)[1:2] = c("Probe_set", "Phy_position")
			dir.create(paste(source_name_I, "/", sample_ID, sep = ""), showWarnings = F)
			write.table(ind_chr, file = paste(source_name_I, "/", sample_ID, "/Chr_", chrname, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
			
		} # for(i in 1:i_len){

		##
		##	Export file
		##
		dir.create(paste(source_name_G, "/Chr_", chrname, sep = ""), showWarnings = F)
		save(AF, file=paste(source_name_G, "/Chr_", chrname, "/AF.RData", sep=""));  
		save(geno, file=paste(source_name_G, "/Chr_", chrname, "/geno.RData", sep=""));		
		colnames(chip) = c("Probe_set", "Phy_position", "Chiptype");
		save(chip, file=paste(source_name_G, "/Chr_", chrname, "/chip.RData", sep=""));   					
		cat("*** Chr ", chrname, ": Chrosomal data combination of ", nrow(samplelist), " samples of ", source_name, " is finished.\n", sep="")
		gc()
	} # for(chr in chrlist){
	dir.create(paste(source_name_G, "/Chr_24", sep = ""), showWarnings = F)
	male_index = (samplelist[,"Gender"]=="M")
	load(file=paste(source_name_G, "/Chr_23/AF.RData", sep=""))
	load(file=paste(source_name_G, "/Chr_23/geno.RData", sep=""))
	AF = as.matrix(AF[,which(male_index)])
	geno = as.matrix(geno[,which(male_index)])
	save(AF, file=paste(source_name_G, "/Chr_24/AF.RData", sep=""))
	save(geno, file=paste(source_name_G, "/Chr_24/geno.RData", sep=""))
	load(file=paste(source_name_G, "/Chr_23/AF.RData", sep=""))
	load(file=paste(source_name_G, "/Chr_23/geno.RData", sep=""))
	AF = as.matrix(AF[,which(!male_index)])
	geno = as.matrix(geno[,which(!male_index)])
	save(AF, file=paste(source_name_G, "/Chr_23/AF.RData", sep=""))
	save(geno, file=paste(source_name_G, "/Chr_23/geno.RData", sep=""))
	file.copy(from = paste(source_name_G, "/Chr_23/chip.RData", sep=""), to = paste(source_name_G, "/Chr_24/chip.RData", sep=""), overwrite = T)
	rm("AF", "chip", "geno", "ind_chr");		gc()

	
	##
	##	Export samplelist
	##	
	# samplelist = samplelist[-3]
	# colnames(samplelist) = c("Obs", "Ind_ID")
	# save(samplelist, file=paste(outpath, "/TEMP/", source_name, "/samplelist.RData", sep="")); 
	write.table(samplelist, paste(outpath, "/TEMP/", source_name, "/samplelist.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

	return(invisible(NULL))

} #	Feature_Stat_Fast_ft



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
##                        Data description                                    ##
################################################################################


YesNo_option=function(choice)
{
    # if(choice==3)
      # return("NA")
    # if(choice==1)
      # return("No")
    # if(choice==2)
      # return("Yes")
	if(choice == 1){
		return("No")
	} else if(choice == 2){
		return("Yes")
	} else return("NA")
	
}

based_option=function(IG_based)
{
    #IG_based=c(1,1)
     IG_based=sum(IG_based%*%c(1,2))
    if(IG_based==0)
      return(" ")
    if(IG_based==1)
      return("Intensity-based")
    if(IG_based==2)
      return("Genotype-based")
    if(IG_based==3)
      return("Intensity-based and Genotype-based")
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

################################################################################
##                        subfunction for CPA                                 ##
################################################################################

##-----------------------sub function - kh--------------------------------------
kh_ku=function(AF,geno)
{
	geno=as.matrix(geno)
	geno[geno!="AB"]=0;  geno[geno=="AB"]=1;
	geno=matrix(as.numeric(geno),dim(geno)[1],dim(geno)[2])
	heterN=apply(geno,1,sum);

	IRAS1=AF;IRAS2=1-AF;
	IRAS=IRAS1/IRAS2
	IRAS=geno*IRAS;IRAS[IRAS=="NaN"]=0
	IRAS=apply(IRAS,1,sum)
	kh=IRAS/heterN

	IRAS1=geno*IRAS1;  IRAS1=apply(IRAS1,1,sum);IRAS1=IRAS1/heterN
	IRAS2=geno*IRAS2;  IRAS2=apply(IRAS2,1,sum);IRAS2=IRAS2/heterN
	IRAS=IRAS1/IRAS2
	ku=kh+(heterN/(heterN-1))*(IRAS-kh)
	CPA=cbind(kh, ku, heterN)
	return(CPA)
}

##------------------------------------------------------------------------------
################################################################################
##                    NoCall Ratio and figure                                 ##
################################################################################
Call_Ratio=function(AF_file)
{
  #AF_file=pathway_adjAF
  ind_file=list.files(AF_file,full.names=T)
  ind_name=list.files(AF_file,full.names=F)
  nind=length(ind_name)
  whole_ratio=NULL;
  for(k in 1:nind)
  {
	#k=1
	chr_file=list.files(ind_file[k],full.names=T)
	nchr=length(chr_file)
	whole_geno=NULL;
	for(j in 1:nchr)
	{
	  #j=1
	  chr=read.table(chr_file[j],head=T,sep="\t")
	  geno=as.character(as.matrix(chr[,4]))
	  GenoName=names(table(geno)) #why use table instead of directly use == "NoCall"?
	  #because HC re-format the column, so it will be "      NoCall", which contains some spaces.
	  #use the name of table seem not bad
	  if(length(GenoName)==3){GenoName=c(GenoName,"NoCall")}
	  geno[geno!=GenoName[4]]=0;
	  geno[geno==GenoName[4]]=1
	  geno=as.numeric(geno)
	  chip=as.character(as.matrix(chr[,3]))
	  chipName=names(table(chip))
	  chip[chip==chipName[1]]=1;
	  chip[chip==chipName[2]]=2;
	  geno=as.numeric(geno)
	  chip=as.numeric(chip)
	  geno=cbind(geno,chip)
	  whole_geno=rbind(whole_geno,geno)
	}
	merge_ncr=sum(whole_geno[,1])/length(whole_geno[,1])
	chip1=whole_geno[whole_geno[,2]==1,];chip1_ncr=sum(chip1[,1])/length(chip1[,2])
	chip2=whole_geno[whole_geno[,2]==2,];chip2_ncr=sum(chip2[,1])/length(chip2[,2])
	merge_ncr=1-merge_ncr;chip1_ncr=1-chip1_ncr;chip2_ncr=1-chip2_ncr;##<------Call Rate
	merge_ncr=sprintf("%.4f",merge_ncr)
	chip1_ncr=sprintf("%.4f",chip1_ncr)
	chip2_ncr=sprintf("%.4f",chip2_ncr)
	nocall_ratio=c(ind_name[k],chip1_ncr,chip2_ncr,merge_ncr) #this is call rate
	whole_ratio=rbind(whole_ratio,nocall_ratio)
  }
  whole_ratio=as.matrix(whole_ratio)
  whole_ratio=matrix(whole_ratio,nind,4)
  whole_ratio=data.frame(whole_ratio)

  colnames(whole_ratio)[1]="Ind_ID"
  colnames(whole_ratio)[2]="CR_Chip1"
  colnames(whole_ratio)[3]="CR_Chip2"
  colnames(whole_ratio)[4]="CR_Merge"
  return(whole_ratio)
}
##----------------------------
##
##  Thinks about that as the NoCall ratio is greater than 0.1
##  How could you display it??
##
NoCall_plot=function(QI,option,Ylab)
{
   #QI=NCr;option="2chip";Ylab=c("Hind","Xba","Merge")
   maintitle="NoCall Rate"
   QI=QI[,c(3:5)]
   chrname=c(1:(ncol(QI)))
	z=QI;
	nind=nrow(QI);
	nchr=ncol(QI)/3
	ynchr=nchr
	chrname=chrname[1:nchr]
	z=as.numeric(as.matrix(z))
	MAXz=ceiling(max(z)*10)/10
	R=1/MAXz
	x=c(1:nind);y=c(1:(3*nchr))
	normalz=as.numeric(as.matrix(QI[c(1:nind),]))
	##--------------------------------------------------------
	z=matrix(z,nind,3*nchr)
	n=18  #n=16
	def.par <- par(no.readonly = TRUE)
	nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(7,1))
	layout.show(nf)
	par(mar=c(5,5,5,0.25))
	if(option=="1chip"){
		Yn=1
		p.choice = 1
	} else {
		Yn=3
		p.choice = 1:3
	}
	image(x,y,z,zlim=range(z),col=terrain.colors(n, alpha = 1),xlim=c(0.5,(nind+0.5)),ylim=c(0.5,Yn*ynchr+0.5),axes=F,cex.lab=1.4,oldstyle =T,ylab=" ",xlab=" ",main=maintitle,cex.main=3,font.main=4)
	for(k in p.choice)
	{
	  #k=1
	  temp_QI=QI[,k];#temp_QI[temp_QI>=0.1]=1;temp_QI[temp_QI<0.1]=0
	  temp_QI=temp_QI*50
	  ##--------------------------
	  #Y=z[,k]*10+0.5+(k-1)
	  Y=z[,k]*R+0.5+(k-1)
	  X=c(x,x[length(x)]+1)-0.5;
	  Y=c(Y,Y[length(Y)])
	  #for(j in 1:nind){rect(X[j],k-0.5,X[j+1],k+0.5,density=temp_QI[j],col="gray")}
	  ##--------------------------
	  polyX=matrix(0,1,length(X)*2);
	  polyX[c(1:length(X))*2]=X
	  polyX[c(1:length(X))*2-1]=X
	  polyY=matrix(0,1,length(Y)*2);
	  polyY[c(1:length(Y))*2]=Y
	  polyY[c(1:length(Y))*2-1]=Y
	  polyY=c(k+0.5,polyY[1:(length(Y)*2-2)],k+0.5)
	  polygon(polyX,polyY,col="white",border="black",lwd=2)
	  #points(X,Y,type="s",lwd=2)
	}
	for(yy in 0:(nchr*3)){points(c(0.5,nind+0.5),c(yy-0.5,yy-0.5),type="l",col="skyblue",lty=2)}
	for(k in 1:length(ind_name)){abline(v=k-0.5,col="skyblue",lty=3)}
	  points(c(0.5,nind+0.5),c(nchr+0.5,nchr+0.5),type="l")
	  points(c(0.5,nind+0.5),c(nchr*2+0.5,nchr*2+0.5),type="l")
	  chrMARK=c(1:nchr)
	  cexsize=1.75
	  ymove=0.5
	  axis(2,at=c(nchr*1,nchr*2,nchr*3),labels=Ylab,font.axis=4,tick=F,line=-0.25+ymove,col.axis="blue",cex.axis=cexsize)
	  axis(1,at=nind/2+0.5,labels=c("Sample"),tick=F,line=1.5,col.axis="darkgreen",cex.axis=cexsize+0.75,font.axis=4)
	  labsize=1.75
	  nbroke=10;
	  ntemp=floor(nind/nbroke);
	if(nind>=30)
	{
	  xind=c(1:nbroke)*ntemp
	  if((nind-max(xind))>=ntemp/2){xind=c(xind,nind)}
	  if((nind-max(xind))<ntemp/2){xind=c(xind[-length(xind)],nind)}
	}
	if(nind<30)
	{
	  xind=c(1:nind)
	}
	 axis(1,at=xind,labels=xind,tick=F,line=-0.25,cex.axis=labsize)
	 abline(v=xind+0.5,col="blue",lwd=1,lty=2)
	 abline(v=c(0,nind)+0.5,col="black",lwd=2)
	 abline(h=c(0,1,2,3)+0.5,col="black",lwd=2)

	 for(j in 1:10)
	 {
		d=0.1		   
	   if(option=="1chip"){
			Ytag=(1 - 0.5 + (j-1)*d)
			len.y = 1
		} else {
			Ytag=(c(1,2,3)-0.5+(j-1)*d)
			len.y = 3
		}
	   abline(h=Ytag,col="blue",lwd=1.5,lty=2)
	   if(j%%2==1){axis(2,at=Ytag,labels=rep(j-1,len.y)*d*MAXz,col="blue",tick=F,font.axis=4,cex.axis=cexsize-0.5,line=-3,las=2)}
	   if(j%%2==0){axis(4,at=Ytag,labels=rep(j-1,len.y)*d*MAXz,col="blue",tick=F,font.axis=4,cex.axis=cexsize-0.5,line=-3,las=2)}
	 }
	axis(4,at=c(1,2,3)+0.5,labels=rep(sprintf("%.2f",MAXz),3),col.axis="Red",tick=F,font.axis=4,cex.axis=cexsize-0.5,line=-3,las=2)
	##----------------------------------------------------------------------------
	par(mar=c(5,0.25,5,1))
	plot(1,1,xlim=c(0,3),ylim=c(0,n+2),type="n",axes=F,ylab="",xlab="")
	leftside=0.25;
	rightside=leftside+1
	d=c(1:n)
	rect(leftside,d,rightside,d+1,col = terrain.colors(n),border = NA)
	rect(leftside,d,rightside,d+1,border ="black")
	xtext=rightside+0.1
	block=4;
	r=diff(range(z))/block#r=(0.1)/block;
	cpoint=min(z)+((1:(block+1))-1)*r
	cpoint=round(cpoint*1000)/1000
	cpoint=sprintf("%1.3f",cpoint)
	text(xtext,1+(0),labels=cpoint[1],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(4.5),labels=cpoint[2],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(9),labels=cpoint[3],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(13.5),labels=cpoint[4],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(18),labels=cpoint[5],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	par(def.par)#- reset to default
}
##---------------------Call Rate figure-------------------------------------
Call_plot=function(QI,option,Ylab)
{
	#QI=CR;option="2chip";Ylab=c("Hind","Xba","Merge")
	maintitle="Call Rate"
	QI=QI[,c(3:5)] #only remain call rate
	#QI=runif(150,0,0.5);QI=matrix(QI,50,3)
	# chrname=c(1:(ncol(QI)))  ## YTLin
	z=QI; #call rate
	nind=nrow(QI);
	nchr=ncol(QI)/3
	ynchr=nchr
	# chrname=chrname[1:nchr]  ## YTLin
	z=as.numeric(as.matrix(z))
	MAXz=ceiling(max(z)*10)/10
	R=1/MAXz
	x=c(1:nind);y=c(1:(3*nchr))
	normalz=as.numeric(as.matrix(QI))
	##--------------------------------------------------------
	z=matrix(z,nind,3*nchr)
	n=18  #n=16
	def.par <- par(no.readonly = TRUE)
	nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(7,1))
	layout.show(nf)
	par(mar=c(5,5,5,0.25))
	# if(option=="1chip"){Yn=1};if(option=="2chip"){Yn=3}
	if(option=="1chip"){
		Yn=1
		p.choice = 1
	} else {
		Yn=3
		p.choice = 1:3
	}
	image(x,y,z,zlim=range(z),col=terrain.colors(n, alpha = 1)[n:1],xlim=c(0.5,(nind+0.5)),ylim=c(0.5,Yn*ynchr+0.5),axes=F,cex.lab=1.4,oldstyle =T,ylab=" ",xlab=" ",main=maintitle,cex.main=3,font.main=4)
	for(k in p.choice) #use polygon to represent call rate
	{
		#k=1
		# temp_QI=QI[,k];#temp_QI[temp_QI>=0.1]=1;temp_QI[temp_QI<0.1]=0  ## YTLin
		# temp_QI=temp_QI*50  ## YTLin
		##--------------------------
		#Y=z[,k]*10+0.5+(k-1)
		Y=z[,k]*R+0.5+(k-1)
		X=c(x,x[length(x)]+1)-0.5;
		Y=c(Y,Y[length(Y)])
		#for(j in 1:nind){rect(X[j],k-0.5,X[j+1],k+0.5,density=temp_QI[j],col="gray")}
		##--------------------------
		polyX=matrix(0,1,length(X)*2);
		polyX[c(1:length(X))*2]=X
		polyX[c(1:length(X))*2-1]=X
		polyY=matrix(0,1,length(Y)*2);
		polyY[c(1:length(Y))*2]=Y
		polyY[c(1:length(Y))*2-1]=Y
		polyY=c(k+0.5,polyY[1:(length(Y)*2-2)],k+0.5)
		polygon(polyX,polyY,col="white",border="black",lwd=2)
		#points(X,Y,type="s",lwd=2)
	}
	for(yy in 0:(nchr*3)){points(c(0.5,nind+0.5),c(yy-0.5,yy-0.5),type="l",col="skyblue",lty=2)}
	for(k in 1:length(ind_name)){abline(v=k-0.5,col="skyblue",lty=3)}
	points(c(0.5,nind+0.5),c(nchr+0.5,nchr+0.5),type="l")
	points(c(0.5,nind+0.5),c(nchr*2+0.5,nchr*2+0.5),type="l")
	chrMARK=c(1:nchr)
	cexsize=1.75
	ymove=0.5
	axis(2,at=c(nchr*1,nchr*2,nchr*3),labels=Ylab,font.axis=4,tick=F,line=-0.25+ymove,col.axis="blue",cex.axis=cexsize)
	axis(1,at=nind/2+0.5,labels=c("Sample"),tick=F,line=1.5,col.axis="darkgreen",cex.axis=cexsize+0.75,font.axis=4)
	labsize=1.75
	nbroke=10;
	ntemp=floor(nind/nbroke);
	if(nind>=30) #if there are too many samples
	{
		xind=c(1:nbroke)*ntemp
		if((nind-max(xind))>=ntemp/2){xind=c(xind,nind)}
		if((nind-max(xind))<ntemp/2){xind=c(xind[-length(xind)],nind)}
	}
	if(nind<30)
	{
		xind=c(1:nind)
	}
	axis(1,at=xind,labels=xind,tick=F,line=-0.25,cex.axis=labsize)
	abline(v=xind+0.5,col="blue",lwd=1,lty=2)
	abline(v=c(0,nind)+0.5,col="black",lwd=2)
	abline(h=c(0,1,2,3)+0.5,col="black",lwd=2)

	for(j in 1:10)
	{
		d=0.1
		if(option=="1chip"){
			Ytag=(1 - 0.5 + (j-1)*d)
			len.y = 1
		} else {
			Ytag=(c(1,2,3)-0.5+(j-1)*d)
			len.y = 3
		}
		abline(h=Ytag,col="blue",lwd=1.5,lty=2)
		if(j%%2==1){axis(2,at=Ytag,labels=rep(j-1,len.y)*d*MAXz,col="blue",tick=F,font.axis=4,cex.axis=cexsize-0.5,line=-3,las=2)}
		if(j%%2==0){axis(4,at=Ytag,labels=rep(j-1,len.y)*d*MAXz,col="blue",tick=F,font.axis=4,cex.axis=cexsize-0.5,line=-3,las=2)}
	}
	##----------------------------------------------------------------------------
	par(mar=c(5,0.25,5,1)) #right palette
	plot(1,1,xlim=c(0,3),ylim=c(0,n+2),type="n",axes=F,ylab="",xlab="")
	leftside=0.25;
	rightside=leftside+1
	d=c(1:n)
	rect(leftside,d,rightside,d+1,col = terrain.colors(n)[n:1],border = NA)
	rect(leftside,d,rightside,d+1,border ="black")
	xtext=rightside+0.1
	block=4;
	r=diff(range(z))/block#r=(0.1)/block;
	cpoint=min(z)+((1:(block+1))-1)*r
	cpoint=round(cpoint*1000)/1000
	cpoint=sprintf("%1.3f",cpoint)
	text(xtext,1+(0),labels=cpoint[1],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(4.5),labels=cpoint[2],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(9),labels=cpoint[3],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(13.5),labels=cpoint[4],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	text(xtext,1+(18),labels=cpoint[5],adj=c(0,0.5),font=4,cex=cexsize-0.5)
	par(def.par)#- reset to default
}





################################################################################
##                          AF - figure                                       ##
################################################################################


AFpp_figure=function(AF,pp,k,i)
{
	#k=ichr;i=1
	font_axis=2;cex_axis=1.5
	n=round(max(pp)/100)
	xtick=c((0:n)*100,max(pp))
	xtick=round(xtick*10)/10
	if(i==1) #intensity-based
	{
		denAF=density(AF);
		denAF.x=denAF$x;
		denAF.y=denAF$y;
		denAF.y=denAF.y-denAF.y[length(denAF.y)]
		plot(pp,AF,pch=1,cex=1,xlim=c(0,max(pp)),ylim=c(0,1),col="blue",axes=F,ylab="",xlab="",main="")
		axis(1,at=xtick,tick = F,cex.axis=cex_axis-0.2,line=-1)
		par(new=T)
		plot(1:10,type="n",ylim=c(0,1),xlim=c(0,3),ylab="",xlab="",axes=F);
		axis(3,at=c(0,3),labels=c(0,1),cex=0.45,line=-0.75,tick=F)
		polygon(denAF.y,denAF.x,col="pink",density=40,angle=90)
		lines(denAF.y,denAF.x,col="purple",cex=1)
		X=c(0,3)
	}
	if(i==2) #genotype-based
	{
		plot(max(pp),max(AF),xlim=c(0,max(pp)),ylim=c(0,1),axes=F,ylab="",xlab="",main="")
		axis(1,at=xtick,tick = F,cex.axis=cex_axis-0.2,line=-1)
		genoAF=table(AF);
		rBB=genoAF[1]/sum(genoAF[1:3]);rAB=genoAF[2]/sum(genoAF[1:3]);rAA=genoAF[3]/sum(genoAF[1:3])
		genoR=cbind(round(rAA*1000)/1000,round(rAB*1000)/1000,round(rBB*1000)/1000);
		r=genoR*max(pp);yline=c(1,0.5,0);d=0.035*4

		rect(0,yline[1]-d,r[1],yline[1]+d,col="pink")
		rect(0,yline[2]-d,r[2],yline[2]+d,col="pink")
		rect(0,yline[3]-d,r[3],yline[3]+d,col="pink")

		text(r[1],yline[1]-1.75*d,labels=sprintf("%.3f", genoR[1]),cex=1.2)
		text(r[2],yline[2]-1.75*d,labels=sprintf("%.3f", genoR[2]),cex=1.2)
		text(r[3],yline[3]+1.75*d,labels=sprintf("%.3f", genoR[3]),cex=1.2)
		points(pp,AF,pch="|",cex=0.75,col="blue")
		X=c(0,range(pp))
	}

	box(lty ="solid")
	axis(3,at=max(X)/2,labels=c(paste("Chromosome=",k,sep="")),tick = F,cex.axis=1.5,font.axis=2)
	axis(1,at=max(X)/2,labels=c("Physical position (Mb)"),tick = F,cex.axis=cex_axis,font.axis=font_axis,col.axis="darkgreen",line=0.5)
	axis(2,at=c(0,0.5,1),tick = F,las=2,cex.axis=cex_axis-0.2,line=-0.85)
	axis(2,at=0.5,labels=c("Allele frequency"),tick = F,cex.axis=cex_axis,font.axis=font_axis,col.axis="darkgreen",line=1.5)
	return(0)
}

AF_plotting=function(pathway_output,chr_file,chromosome_choice,temp_ind_name,p,GenoAF)
{
	#pathway_output=pathway_Intensity_AFplot;p=1;GenoAF=2
	nchr=length(chr_file)
	if(p==1)## whole chip
	{
		dir.create(paste(pathway_output,"/Merge",sep=""),showWarnings=F);
		tiff(filename=paste(pathway_output,"/Merge/",temp_ind_name,".tiff",sep=""),width = 1200*3, height = 800*3, pointsize = 14, res = 216, compression = "lzw")
		chipname="Merge";
	}
	if(p==2)## Nsp chip
	{
		dir.create(paste(pathway_output,"/Chip 1",sep=""),showWarnings=F);
		tiff(filename=paste(pathway_output,"/Chip 1/",temp_ind_name,".tiff",sep=""),width = 1200*3, height = 800*3, pointsize = 14, res = 216, compression = "lzw")
		chipname="Chip 1"; #Nsp
	}
	if(p==3)## Sty chip
	{
		dir.create(paste(pathway_output,"/Chip 2",sep=""),showWarnings=F);
		tiff(filename=paste(pathway_output,"/Chip 2/",temp_ind_name,".tiff",sep=""),width = 1200*3, height = 800*3, pointsize = 14, res = 216, compression = "lzw")
		chipname="Chip 2"; #Sty
	}
	if(nchr>1)
	{
		if(nchr>1 & nchr<5){a=2}
		if(nchr>4 & nchr<10){a=3}
		if(nchr>9 & nchr<17){a=4}
		if(nchr>16){a=5}
		par(mfcol=c(a,a),mai=c(0.5,0.5,0.5,0.2),oma=c(0.1,1,0.1,0.1))
	}
	for(i in 1:nchr)
	{
		#i=11
		ichr=chromosome_choice[i]
		chr=read.table(chr_file[i],fill=T,header=T,sep="")
		chiptype_name=table(chr[,3])
		if(p==2){chr=chr[chr$Chiptype=="chip1",]}
		if(p==3){chr=chr[chr$Chiptype=="chip2",]}

		if(nrow(chr)!=0)
		{
			chr=chr[is.na(chr[,5])==F,]
			ppaf=cbind(chr[,2],chr[,5])
			ppaf=ppaf[is.na(ppaf[,2])==F,]
			pp=ppaf[,1]/1000000;

			if(GenoAF==1) ##AF-physical position
			{
				AF=ppaf[,2]; AF=AF[is.na(AF)==F]
				AFpp_figure(AF,pp,ichr,1)
			}
			if(GenoAF==2) ##Geno-physical position
			{
				Geno=chr[,4];  ##Genotypr
				Geno=as.character(Geno)
				Geno[Geno=="AA"]=1;Geno[Geno=="AB"]=0.5
				Geno[Geno=="BB"]=0;Geno[Geno=="NoCall"]=2
				AF=as.numeric(Geno)
				AFpp_figure(AF,pp,ichr,2)
			}
		}
	}
	dev.off()
	return(0)
}
  
  
  
################################################################################
##                          AFref calculation                                 ##
################################################################################
  AF_ref=function(AF,GE)
  {
    #AF=AF_temp;GE=GE_AA
    GE[GE==0]=NA
    a=AF*GE;
    num_GE=apply(GE,1,function(x){sum(x,na.rm=T)});
    #naAF=(is.na(AF)*1-1)*(-1)*GE;nonNAAF=apply(naAF,1,sum);num_GE=nonNAAF##<<---special critior for have genotype without allele frequency
    sum_AF=apply(a,1,function(x){sum(x,na.rm=T)});mean_AF=sum_AF/num_GE
    std_AF=apply(a*a,1,function(x){sum(x,na.rm=T)});std_AF=std_AF-num_GE*(mean_AF)^2;std_AF=std_AF/(num_GE-1)
    std_AF[std_AF*100000000<1]=0  ##<--------------noise
    whole_display=cbind(num_GE,mean_AF,sqrt(std_AF))
	whole_display[is.nan(whole_display)] = NA
    return(whole_display)
  }
  

################################################################################
##                          AFref calculation                                 ##
################################################################################

##111111111111111111111111111111111111111111111111111111111111111111111111111111
##------------------------------  Quality Index 1  -----------------------------
##==============================================================================
Quality_Index_1=function(a)
{
	#a=chr
	#Probe_set Genotype Chiptype Phy_position Genotype.1 adjustAF   mean_AA     std_AA   mean_AB     std_AB    mean_BB     std_BB
	AF=as.numeric(a[,6])
	## standardize => Mahalanobis distance
	Q1AA=((AF-as.numeric(a[,7]))/as.numeric(a[,8]))^2;
	Q1AB=((AF-as.numeric(a[,9]))/as.numeric(a[,10]))^2;
	Q1BB=((AF-as.numeric(a[,11]))/as.numeric(a[,12]))^2;
	QI=cbind(Q1AA,Q1AB,Q1BB);

	#QI[is.na(QI)==T]=0
	a=cbind(a[,1:6],Q1AA,Q1AB,Q1BB)
	GT=as.character(a[,2])
	#GT=as.integer(factor(GT, levels = c(3, 2, 1), labels = c(1:3)))
	# GTAA=as.numeric(GT==3)  ##AA=3
	# GTAB=as.numeric(GT==2)  ##AB=2
	# GTBB=as.numeric(GT==1)  ##BB=1  ##NoCall=0
	GTAA = GTAB = GTBB = NA
	## GTAA=1: genotype=AA, o.w. GTAA=NA, etc
	GTAA[GT == 3] = GTAB[GT == 2] = GTBB[GT == 1] = 1
	GT=cbind(GTAA,GTAB,GTBB)

	QI=QI*GT
	#QI[QI=="NaN"]=0;QI[QI=="Inf"]=0;
	# QI=apply(QI,1,sum) 	
	num.infinite = rowSums(is.finite(QI))	# possible value = 0, 1 
	index.infinite = which(num.infinite == 0)
	QI = rowSums(QI, na.rm = T)
	QI[index.infinite] = NA		## The non-inofmative SNPs (no-available-QI-value in all genotpyes or NoCall) are given NA
	#stopifnot(all(num.infinite != 0))
	a=cbind(a[,1:6],QI) # 1. Probe_set, 2. Genotype, 3. Chiptype, 4. Phy_position, 5. Genotype.1, 6. adjustAF,
	#a[a[,2]==0,7]="NA"
	colnames(a)[1]="Probe_set";
	colnames(a)[2]="Genotype";
	colnames(a)[3]="Chiptype"
	colnames(a)[4]="Phy_position";
	colnames(a)[5]="Genotype";
	colnames(a)[6]="AF";
	colnames(a)[7]="QI";
	a=data.frame(a, stringsAsFactors = F)
	# a[,1]=as.numeric(a[,1])
	a[,2]=as.character(a[,2])
	a[,3]=as.character(a[,3])
	a[,4]=as.numeric(a[,4])
	# a=data.frame(a, stringsAsFactors = F)
	return(a)
}
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##222222222222222222222222222222222222222222222222222222222222222222222222222222
##------------------------------  Quality Index 2  -----------------------------
##==============================================================================

Quality_Index_2=function(a)
{
	#a=chr
    #1. Probe_set, 2. Genotype, 3. Chiptype, 4. Phy_position, 5. Genotype.1, 6. adjustAF, 7. mean_AA, 8. std_AA, 9. mean_AB, 10. std_AB, 11. mean_BB, 12.std_BB
	# z=qnorm(0.95, mean = 0, sd = 1)
	##===================================================================
	##		Difference from QI1: without considering genotype
	##===================================================================
	##	Judge which genotype that SNP closest to 
	distance=cbind(as.numeric(a[,6])-as.numeric(a[,7]),as.numeric(a[,6])-as.numeric(a[,9]),as.numeric(a[,6])-as.numeric(a[,11]))
	distance=abs(distance);
	distance[is.na(distance)==T]=5
	RangeLevel=apply(distance,1,function(x){x==min(x)})
	RangeLevel= t(RangeLevel)*1
	RangeLevel[RangeLevel == 0] = NA
	
	level=a[,c(6:12)]	#adjustAF   mean_AA     std_AA   mean_AB     std_AB    mean_BB     std_BB
	level=as.matrix(level)
	level=matrix(as.numeric(level),nrow(level),ncol(level))
	## standardize => Mahalanobis distance
	AA=((level[,1]-level[,2])/level[,3])^2
	AB=((level[,1]-level[,4])/level[,5])^2
	BB=((level[,1]-level[,6])/level[,7])^2
	QI=cbind(AA,AB,BB);
	
	# QI[is.na(QI)==T]=0
	QI=QI*RangeLevel
	# QI[QI=="NaN"]=0;
	# QI[QI=="Inf"]=0;
	
	#QI[QI=="NaN"]=0;QI[QI=="Inf"]=0;
	# QI=apply(QI,1,sum) 
	num.infinite = rowSums(is.finite(QI))
	index.infinite = which(num.infinite == 0)
	QI = rowSums(QI, na.rm = T)
	QI[index.infinite] = NA
	
	a=cbind(a[,1:6],QI)
	colnames(a)[1]="Probe_set";
	colnames(a)[2]="Genotype";
	colnames(a)[3]="Chiptype"
	colnames(a)[4]="Phy_position";
	colnames(a)[5]="Genotype";
	colnames(a)[6]="AF";
	colnames(a)[7]="QI";
	a=data.frame(a, stringsAsFactors = F)
	# a[,1]=as.numeric(a[,1])
	a[,2]=as.character(a[,2])
	a[,3]=as.character(a[,3])
	a[,4]=as.numeric(a[,4])
	# a=data.frame(a)
	return(a)
}
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputmanage=function(result)
{
  result[,3]=signif(as.numeric(result[,3]),4)
  result[,6]=signif(as.numeric(result[,6]),4)
  result[,9]=signif(as.numeric(result[,9]),4)
  result[,3]=sprintf("%.4f",result[,3])
  result[,6]=sprintf("%.4f",result[,6])
  result[,9]=sprintf("%.4f",result[,9])
  result[,2]=sprintf("%5.f",result[,2])
  result[,5]=sprintf("%5.f",result[,5])
  result[,8]=sprintf("%5.f",result[,8])
  return(result)
}


QIstat=function(QI,qtrimmean,qmedian)
{
  #QI=QImerge;qtrimmean=0.95
  # QI[QI=="NA"]=NA
  QI=as.matrix(QI)#[,c(6,7,8)])
  QI=matrix(as.numeric(QI),nrow(QI),ncol(QI))
  medianQI=apply(QI,2,function(x){quantile(x,qmedian,na.rm=T)})
  ##	One-sided winsorized mean 
  quantile_mean=apply(QI,2,function(x){quantile(x,qtrimmean,na.rm =T)})
  QI[QI[,1]>quantile_mean[1],1]=quantile_mean[1]
  QI[QI[,2]>quantile_mean[2],2]=quantile_mean[2]
  QI[QI[,3]>quantile_mean[3],3]=quantile_mean[3]
  meanQI=apply(QI,2,function(x){mean(x,na.rm=T)})
  QI=c(meanQI,medianQI)
  QI=matrix(QI,1,6)
  return(QI)
}


################################################################################
###                  Graphical Results area                                  ###
################################################################################

lnorm_data_select=function(chr_QI,a)
{
  #chr_QI=z;a=0.95
  temp_chr_QI<-(chr_QI-mean(chr_QI))/sd(chr_QI)
  d_temp_chr_QI<-density(temp_chr_QI)

  threshold<-min(d_temp_chr_QI$x);temp_chr_QI_unlnorm=temp_chr_QI-threshold;#<--------warnings

  temp_chr_QI_lognormal<-log(temp_chr_QI_unlnorm,base=exp(1))
  m_log_QI<-mean(temp_chr_QI_lognormal)
  sd_log_QI<-sd(temp_chr_QI_lognormal)

  zMAX=qlnorm(a,m_log_QI,sd_log_QI)
  zMAX=(exp(zMAX)+threshold)*sd(chr_QI)+mean(chr_QI);
 # temp_chr_QI_lognormal[temp_chr_QI_lognormal>zMAX]=zMAX
 # temp_chr_QI=(exp(temp_chr_QI_lognormal)+threshold)*sd(chr_QI)+mean(chr_QI);
  return(zMAX)
}

##==============================================================================
##=================== Image plotting for chip base =============================
##==============================================================================

HeatMap_Plot=function(QI,maintitle,option,Ylab,QIstat)
{
	## QIstat: QI from database
	#QI=A;
	#maintitle=QI_name[k];option=chiptype_choice;Ylab=chipname_list;QIstat=DB_QI
	chrname=c(1:(ncol(QI))) ##	y-axis: Merge, Chip2, Chip1
	z=QI;
	nind=nrow(QI)
	nchr=ncol(QI)/3 ## Current version doesn't support chromosomal plot, so nchr = 1
	ynchr=nchr
	chrname=chrname[1:nchr]
	z=as.numeric(as.matrix(z))
	x=c(1:nind);y=c(1:(3*nchr))
	normalz=as.numeric(as.matrix(QI[c(1:nind),]))
	ave_z=QIstat[1];sd_z=QIstat[2]
	#zMAX=lnorm_data_select(normalz,0.9)
	sigma=ave_z+c(0:6)*sd_z	#	For the palette, will text 6 labels on it
	##--------avoid outline value-----------------------------
	z[z<sigma[1]]=sigma[1];		#	Replace the extremely small value (lower than mean) with mean+1*sd
	z[z>sigma[7]]=sigma[7];		#	Replace the extremely large value (larger than mean+6*sd) with mean+6*sd
	##--------------------------------------------------------
	#	For the palette, will text 6 labels on it
	block=6;
	r=(max(z)-min(z))/block;
	cpoint=((1:(block+1))-1)*r+min(z)
	cpoint=sigma
	cpoint=round(cpoint*100)/100
	##--------------------------------------------------------
	z=matrix(z,nind,3*nchr)
	n=18  #n=16
	def.par <- par(no.readonly = TRUE)
	nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(7,1))
	layout.show(nf)
	par(mar=c(5,5,5,0.25))
	if(option=="2chip"){nY=3}	# Merge, Chip2, Chip1
	if(option=="1chip"){nY=1}	# Only draw one-chip result (Still draw three figures but )
	# zlim=c(sigma[1],sigma[7]): mean+1*sd, mean+6*sd
	# The way of "image" command is drawing from 0.5, 1.5, ...
	#		ylim=c(0.5,nY*ynchr+0.5): 0.5-1.5: the first chip, 1.5-2.5: the second chip, 2.5-3.5: Merge
	image(x,y,z,zlim=c(sigma[1],sigma[7]),col=terrain.colors(n, alpha = 1),xlim=c(0.5,(nind+0.5)),ylim=c(0.5,nY*ynchr+0.5),axes=F,cex.lab=1.4,oldstyle =T,ylab=" ",xlab=" ",main=maintitle,cex.main=3,font.main=4)

	#for(yy in 0:(nchr*3)){points(c(0.5,nind+0.5),c(yy-0.5,yy-0.5),type="l",col="skyblue",lty=2)}
	#for(k in 1:length(ind_name)){abline(v=k-0.5,col="skyblue",lty=3)}
	points(c(0.5,nind+0.5),c(nchr+0.5,nchr+0.5),type="l")	 	#	i.e. abline(h=nchr+0.5)
	points(c(0.5,nind+0.5),c(nchr*2+0.5,nchr*2+0.5),type="l")	#	i.e. abline(h=nchr*2+0.5)
	cexsize=1.75;ymove=0.5
	axis(2,at=c(nchr*1,nchr*2,nchr*3),labels=Ylab,font.axis=4,tick=F,line=-0.25+ymove,col.axis="blue",cex.axis=cexsize)
	axis(1,at=nind/2+0.5,labels=c("Sample"),tick=F,line=1.5,col.axis="darkgreen",cex.axis=cexsize+0.75,font.axis=4)
	labsize=1.75
	nbroke=10

	##	If nind > 30, then the labels of samples will be texted every 10 samples
	ntemp=floor(nind/nbroke)
	if(nind>=30)
	{
		xind=c(1:nbroke)*ntemp
		if((nind-max(xind))>=ntemp/2){xind=c(xind,nind)}
		if((nind-max(xind))<ntemp/2){xind=c(xind[-length(xind)],nind)}
	}
	if(nind<30){xind=c(1:nind)}
	abline(v=xind+0.5,col="blue",lwd=1,lty=2)
	axis(1,at=xind,labels=xind,tick=F,line=-0.25,cex.axis=labsize)
	box(col="blue",lwd=2,lty=1)
	##----------------------------------------------------------------------------
	##
	##	Palette
	##
	par(mar=c(5,0.25,5,1))
	plot(1,1,xlim=c(0,3),ylim=c(0,n+2),type="n",axes=F,ylab="",xlab="")
	leftside=0.25;
	rightside=leftside+1
	d=c(1:n)
	rect(leftside,d,rightside,d+1,col = terrain.colors(n),border = NA)
	rect(leftside,d,rightside,d+1,border ="black")
	cpoint=round(sigma*100)/100
	xtext=rightside+0.1
	tagsize=0.25
	cpoint=sprintf("%1.2f",cpoint)
	##	The six text labels on the right side of palette
	text(xtext,1+(0+1)/2,labels=cpoint[1],adj=c(0,0.5),cex=0.85+tagsize)
	text(xtext,1+(3),labels=cpoint[2],adj=c(0,0.5),cex=0.85+tagsize)
	text(xtext,1+(6),labels=cpoint[3],adj=c(0,0.5),cex=0.85+tagsize)
	text(xtext,1+(9),labels=cpoint[4],adj=c(0,0.5),cex=0.85+tagsize)
	text(xtext,1+(12),labels=cpoint[5],adj=c(0,0.5),cex=0.85+tagsize)
	text(xtext,1+(15),labels=cpoint[6],adj=c(0,0.5),cex=0.85+tagsize)
	text(xtext,1+(17+18)/2,labels=cpoint[7],adj=c(0,0.5),cex=0.85+tagsize)
	par(def.par)#- reset to default
	##----------------------------------------------------------------------
	return(0)
}


Polygon_Plot=function(a,ub,figure_title,option,Ylab)
{
  #a=temp_QI;ub=temp_UBQI ;figure_title="A";option="2chip"
  # a = temp_QI; ub = temp_UBQI; figure_title = QI_name[k]; chiptype_choice = option; chipname_list = Ylab
  nind=nrow(a);nchip=ncol(a)
  ub_matrix=rep(ub,nind);ub_matrix=matrix(ub_matrix,nind,nchip,byrow=T)
  ##		z = 1: QI_data >= QI_database
  ##		z = 0: QI_data > QI_database
  z=a-ub_matrix
  z[z>=0]=1;z[z<0]=0
  z=as.matrix(z)
  z=matrix(as.numeric(z),nind,nchip)# Z is a 1/0 matrix
  y=c(1:nchip);x=c(1:nind)
  if(option=="2chip"){nY=3}	# Merge, Chip1, Chip2
  if(option=="1chip"){nY=1}	# Chip1

  ## Color the background to show the chips are good-quality or poor-quality
  if(sum(z)==nind*nchip)	# All chips are bad-quality, so assign the darkest color
  {image(x,y,z,col=c(terrain.colors(18, alpha = 1)[18]),axes=F,ylab="",xlab="",ylim=c(0.5,nY+0.5),main=figure_title,cex.main=3,font.main=4)}
  if(sum(z)==0)	# All chips are good-quality, so assign the lightest color
  {image(x,y,z,col=c(terrain.colors(18, alpha = 1)[1]),axes=F,ylab="",xlab="",ylim=c(0.5,nY+0.5),main=figure_title,cex.main=3,font.main=4)}
  if(sum(z)<nind*nchip & sum(z)>0)	# Assign colors between the lightest and darkest colors
  {image(x,y,z,col=c(terrain.colors(18, alpha = 1)[1],terrain.colors(18, alpha = 1)[18]),axes=F,ylab="",xlab="",ylim=c(0.5,nY+0.5),main=figure_title,cex.main=3,font.main=4)}

  ## This solid skyblue line will be covered by another dotted blue line
  abline(v=x+0.5,col="skyblue")
  abline(h=y+0.5,col="skyblue")

  for(i in 1:3)
  {
    #i=1
    tempQI=a[,i];UB=ub[i]
    b=c(tempQI,UB,0);
	yQI=range(b);	yQI[2]=yQI[2]*1.1
    newY=(i-0.5)+(b-yQI[1])/(yQI[2]-yQI[1])
    tempQI=newY[x]
    newUB=newY[max(x)+1]
    #points(c(x,max(x)+1)-0.5,c(tempQI,tempQI[max(x)]),type="s",col="blue",pch="+")
    points(x,tempQI,col="blue",pch="+")	# "+": the intersected points of the red lines
    polygonX=c(0.5,x,nind+0.5)
    polygonY=c(i-0.5,tempQI,i-0.5)
    polygon(polygonX,polygonY,col="yellow",pch="+",density=45,border="brown")
    UBlab=round(UB*1000)/1000
    abline(h=newUB,col="purple",lwd=2)
    #text(nind*0.1,newUB,labels=UBlab,cex=1.5)
    text(max(nind*0.1, 1),newUB,labels=UBlab,cex=1.5)
    #axis(2,at=c(i-0.5,i+0.5),labels=c("",""),tick=T,line=0.25)
  }
  abline(h=c(0.5,1.5,2.5),lwd=2,col="blue")
  box(col="blue",lwd=2,lty=1)
  cexsize=1.75;ymove=0.5
  ##-------------------------- X axis setting ----------------------------------
  nbroke=10;labsize=1.75
  ntemp=floor(nind/nbroke);
  if(nind>=30)
  {
    xind=c(1:nbroke)*ntemp
    if((nind-max(xind))>=ntemp/2){xind=c(xind,nind)}
    if((nind-max(xind))<ntemp/2){xind=c(xind[-length(xind)],nind)}
  }
  if(nind<30){xind=c(1:nind)}
  axis(1,at=xind,labels=xind,tick=F,line=-0.25,cex.axis=labsize)
  abline(v=xind+0.5,col="blue",lwd=1,lty=2) # to separate each individaul
  abline(v=c(0,nind)+0.5,col="black",lwd=2)	 # The outside box
  abline(h=c(0,1,2,3)+0.5,col="black",lwd=2) # to separate each chip
  ##-------------------------- Y axis setting ----------------------------------
  cexsize=1.75;ymove=0.5
  axis(2,at=c(1:3),labels=Ylab,font.axis=4,tick=F,line=-0.25+ymove,col.axis="blue",cex.axis=cexsize)
  axis(1,at=nind/2+0.5,labels=c("Sample"),tick=F,line=1.5,col.axis="darkgreen",cex.axis=cexsize+0.75,font.axis=4)
  ##--------------------------------------------------
  ##	 z: an 0/1 index matrix, normal_chip_ratio = the average ratio of samples with bad-quality chip (?)
  ##		But this value is not of use currently!!
  normal_chip_ratio=1-apply(z,2,mean)
  normal_chip_ratio=as.numeric(normal_chip_ratio)
  normal_chip_ratio=round(normal_chip_ratio*100)/100
  #for(i in 1:3){text(nind,i+0.45,labels=normal_chip_ratio[i],adj=c(1,1),cex=cexsize)}
  ##--------------------------------------------------
  return(normal_chip_ratio)
}


##==============================================================================
##==============================================================================
##==============================================================================


################################################################################
###                  Numerical Results area                                  ###
################################################################################

#------------------Correlation--------------------------------------------------
QI_corr=function(QI1,QI2)
{
  #QI1=chip_MQI1;QI2=chip_MQI2
  Chip1=corr(cbind(as.numeric(QI1[,3]),as.numeric(QI2[,3])));
  Chip2=corr(cbind(as.numeric(QI1[,4]),as.numeric(QI2[,4])));
  Merge=corr(cbind(as.numeric(QI1[,2]),as.numeric(QI2[,2])));
  Chip1=format(as.numeric(Chip1),nsmall=4)
  Chip2=format(as.numeric(Chip2),nsmall=4)
  Merge=format(as.numeric(Merge),nsmall=4)
  Chip1=sprintf("%1.5s",Chip1)
  Chip2=sprintf("%1.5s",Chip2)
  Merge=sprintf("%1.5s",Merge)
  whole=rbind(Chip1,Chip2,Merge)
  return(whole)
}
##--------------------------------------------------------------------------
spearman_QI_corr=function(QI1,QI2)
{
  #QI1=chip_MQI1;QI2=chip_MQI2
  spearman_X=cor.test(as.numeric(QI1[,4]), as.numeric(QI2[,4]),alternative = "two.sided",method ="spearman",exact = F, conf.level = 0.95)
  spearman_H=cor.test(as.numeric(QI1[,3]), as.numeric(QI2[,3]),alternative = "two.sided",method ="spearman",exact = F, conf.level = 0.95)
  spearman_Combine=cor.test(as.numeric(QI1[,2]), as.numeric(QI2[,2]),alternative = "two.sided",method ="spearman",exact = F, conf.level = 0.95)

  # Chip1=unlist(spearman_X)[3];
  # Chip2=unlist(spearman_H)[3];
  Chip1=unlist(spearman_H)[3]
  Chip2=unlist(spearman_X)[3]
  Merge=unlist(spearman_Combine)[3];
  Chip1=as.matrix(Chip1)
  Chip2=as.matrix(Chip2)
  Merge=as.matrix(Merge)
  Chip1=format(as.numeric(Chip1),nsmall=4)
  Chip2=format(as.numeric(Chip2),nsmall=4)
  Merge=format(as.numeric(Merge),nsmall=4)
  Chip1=sprintf("%1.5s",Chip1)
  Chip2=sprintf("%1.5s",Chip2)
  Merge=sprintf("%1.5s",Merge)

  rho=rbind(Chip1,Chip2,Merge)
  return(rho)
}
##---------------------------------sub-function------------------------------
basic_stat_table=function(chip_QI,UB,QIname,a,chipsize_option)
{
  #chip_QI=chip_WQI1;UB=UBQI
  n=dim(chip_QI)[1]
  QI1_chip1=as.numeric(chip_QI[,3]);
  QI1_chip2=as.numeric(chip_QI[,4]);
  QI1_merge=as.numeric(chip_QI[,2]);
  #QI1=c(QI1_chip1,QI1_chip2,QI1_merge)
  Chip1=summary(QI1_chip1);Chip2=summary(QI1_chip2);Merge=summary(QI1_merge)

  Chip1=sprintf("%.4f",Chip1);Chip2=sprintf("%.4f",Chip2);Merge=sprintf("%.4f",Merge);
  std=cbind(sd(QI1_chip1),sd(QI1_chip2),sd(QI1_merge));std=sprintf("%.4f",std)
  max_chip1=UB[1]; max_chip2=UB[2]; max_merge=UB[3]       #<----from database
  max_value=cbind(max_chip1,max_chip2,max_merge)
  max_numeric=as.numeric(max_value)
  num_over_upper=cbind(sum(as.numeric(QI1_chip1>max_numeric[1])),sum(as.numeric(QI1_chip2>max_numeric[2])),sum(as.numeric(QI1_merge>max_numeric[3])))

  #QI1=cbind("                Chip 1","Chip 2","Merge")
  if(chipsize_option=="Affy 100K"){QI1=cbind("                Hind","Xba","Merge")};
  if(chipsize_option=="Affy 500K"){QI1=cbind("                Nsp","Sty","Merge")};

  # QI1=rbind(QI1,n,cbind(Chip1,Chip2,Merge),std,max_value,num_over_upper)
  QI1=rbind(QI1,n,cbind(Chip1,Chip2,Merge),std,round(max_value, 4),num_over_upper)
  rownames(QI1)[1]=QIname
  rownames(QI1)[2]="Num. of Sample    \t"
  rownames(QI1)[3]="Min.              \t"
  rownames(QI1)[4]="1st Quantile      \t"
  rownames(QI1)[5]="Median            \t"
  rownames(QI1)[6]="Mean              \t"
  rownames(QI1)[7]="3rd Quantile      \t"
  rownames(QI1)[8]="Max               \t"
  rownames(QI1)[9]="Std               \t"
  rownames(QI1)[10]=paste("Upper of ",a*100,"%     \t",sep="")
  rownames(QI1)[11]="Num. of over upper\t"
  # QI1[10,]=sprintf("%.4s",QI1[10,])
  # QI1[10,]=sprintf("%.6s",QI1[10,])
  return(QI1)
}

###---------------------identification poor chip -------------------------------

over_upper_list=function(ub_table,chipQI)
{
  #ub_table=temp_QIUB;chipQI=temp_chip_QI
  ind_name=chipQI[,1]
  nind=length(ind_name)
  QIUB_table=t(matrix(rep(ub_table,nind),3,nind))
  QI=chipQI[,-1];QI=matrix(as.numeric(QI),nrow(QI),ncol(QI))
  a=QI-QIUB_table
  a[a>=0]=1;a[a<0]=0
  idenpoor=apply(a,1,sum)
  idenpoor[idenpoor>0]="  Poor"
  idenpoor[idenpoor==0]=" "
  
  a[a==1]="*";a[a==0]=" "
  QI=matrix(sprintf("%.4f",QI),nrow(a),ncol(a))
  a=matrix(paste(a,QI,sep=""),nrow(a),ncol(a))
  a=cbind(c(1:nind),ind_name,a,idenpoor)
  a=data.frame(a)
  return(a)
}
