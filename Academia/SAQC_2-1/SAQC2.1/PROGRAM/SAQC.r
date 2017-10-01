#rm(list = ls()) # remove whole parameters
library(gtools) ## manage file name list...
library(boot)
cat("---------------------------------------------------------------------  \n")
cat("The main function is starting.  \n")
#main_pathway = paste("J:/working_area/SAQC", sep = "");
main_pathway = strsplit(pathway_program, "PROGRAM");
main_pathway = unlist(main_pathway)[1]
ptm <- proc.time()   ##<----time record
##==============================================================================
## chipsize_choice
##                        Input parameter area
##----Data format------
population_group = c("Asian", "African", "European", "Taiwan", "Combined")    #<---for data description
population_database = c("Asian", "YRI", "CEU", "Taiwan", "Combined")          #<---for data base
chipsize_group = c("Affy 100K", "Affy 500K", "Affy 6.0", "Affy Axiom", "Illumina 550K")
chiptype_group = c("1chip", "2chip") ##<<--------
onechip_group = c("chip1", "chip2")   ##<<-----Hind/Xba ; Nsp/Sty
#chromosome_group = c(1:23)
chromosome_choice = c(1:23)
##----Statistical analysis------
CPA_calculation_option = c(1, 2, 3)
#1: using CPA estimates from SAQC-provided database;	2: based on user-provided data;	3: disable (if input data is AF-based)

AF_calculation_option = c(1, 2, 3);
#1: No;	2: Yes;	3: disable (if input data is AF-based)

AF_ref_calculation_option = c(1, 2, 3)
#1:Database;	2: User_provided;	3: disable (if input data is AF-based)

QI_calculation_option = c(1, 2, 3)
#1: No;	2: Yes;	3: disable (if input data is intensity-based & choose no AF calculation)

Iden_poor_chipchr_option = c(1, 2, 3)
#1: No;	2: Yes;	3: diable

QI_quantile_option = c(0.95, 0.975, 0.99)
##----Ressult visualization-----
AF.plot_option = c(1, 2, 3);#AF.based_option = c(1, 1)      ## Intensity-based, Genotype based,  both
CR.plot_option = c(1, 2, 3)
QI.plot_option = c(1, 2, 3);      ##HeatMap   #QI.based_option = c(1, 1)      ## chip-based, chromosome based,  both
QI_est.based_option = c(1, 1)
QI.scatterplot_option = c(1, 2, 3)##Polygon
##----Numerical output-----------
data_description_option = c(1, 2)
CPA_est_option = c(1, 2, 3);    #3 have been done
AF_est_option = c(1, 2, 3);     #3 have been done
AF_est.based_option = c(1, 1)
#AFbi_est_option = c(1, 2, 3);AFbi_est.based_option = c(1, 1)
QI_est_option = c(1, 2, 3);
#QIbi_est_option = c(1, 2, 3);QIbi_est.based_option = c(1, 1)
Poor_chipchr_option = c(1, 2, 3)                     ## chip-based, chromosome based,  both
#poor_type_option=c(1, 0) ##poor chip: 1/0 ;poor chromosome: 1/0
##==============================================================================
##------------------------------------------------------------------------------

##==============================================================================
##
##                      parameter setting
##------------------------------------------------------------------------------

##=======Input/Output path=======
inputdata_format=as.numeric(inputdata_format) #1:Genotype-based; 2: AF-based

if(data_input == "Test example 1"){
	#if(file.exists(data_output)=F){data_output=paste(pathway,"OUTPUT/Test_Affy100K",sep="")}
	dir.create(data_output,showWarnings=F)
	pathway_input=paste(pathway,"EXAMPLE/Test_Affy500K",sep="");    ##<<---Input pathway
	population_name=c("Test_Affy500K")
} else if(data_input == "Test example 2"){
	#if(file.exists(data_output)=F){data_output=paste(pathway,"OUTPUT/Test_Affy100K",sep="")}
	dir.create(data_output,showWarnings=F)
	pathway_input=paste(pathway,"EXAMPLE/Test_Affy100K",sep="");    ##<<---Input pathway
	population_name=c("Test_Affy100K")
}  else if(data_input == "Test example 3"){
	#if(file.exists(data_output)=F){data_output=paste(pathway,"OUTPUT/Test_Affy100K",sep="")}
	dir.create(data_output,showWarnings=F)
	pathway_input=paste(pathway,"EXAMPLE/Test_Illu550K",sep="");    ##<<---Input pathway
	population_name=c("Test_Illu550K")
} else {
	pathway_input=data_input;
	population_name=strsplit(data_input,"INPUT/")
	population_name=population_name[[1]][2]
}    ##<<---Input pathway

##-----------Input file checking------------------------------------------------
cat("Input file checking is preparing.\n")
infile_checking(infile_path = pathway_input, inputformat = inputdata_format, Ind_Array = chipsize_choice)    #inputdata_format: 1: Genotype/Intensity 2: Allele frequecny
cat("Input file checking is finished.\n")
##------------------------------------------------------------------------------
pathway_output=data_output;  ##<<---Output pathway
dir.create(pathway_output,showWarnings=F)
pathway_output_population=paste(pathway_output,"/",population_name,sep="")
dir.create(pathway_output_population,showWarnings=F)
  
Null_SNP_AFref_info=paste(pathway_output_population,"/Null_SNP_list.txt",sep="")
if(AFref_source_choice[2]==1) #if AFref is user-provided
{
	if( nrow(Null_SNP_AFref_table)>0)
	{write.table(Null_SNP_AFref_table,Null_SNP_AFref_info,col.names=T,row.names=F,quote=F,sep="\t")} 
}
 
logfile=paste(pathway_output_population,"/SAQC-log.txt",sep="")
##------------------------------------------------------------------------------
cat("### ----- SAQC Log file ----- ###\n",file=logfile,append=F)
cat("   ~ The input options are preparing.\n",file=logfile,append=T)
##=======Data format=======
ipopu=as.numeric(group_choice)
population_choice=population_group[ipopu]
population_database_choice=population_database[ipopu]   ## study population : Asia, Africa, Europe, Taiwan, Combine
Ind_Array = chipsize_choice   ## Affy 100K, Affy 500K, Affy 6.0, Affy Axiom, Illumian 550K----
chipsize_choice=chipsize_group[as.numeric(chipsize_choice)]   ## Affy 100K, Affy 500K, Affy 6.0, Affy Axiom, Illumian 550K----
chiptype_choice=chiptype_group[as.numeric(chiptype_choice)]   ##Chip1 Merge ----
if(as.numeric(onechip_choice)==0){
	onechip_choice = chipsize_choice
}else{
	onechip_choice=onechip_group[as.numeric(onechip_choice)] #if user selects one array, this index refers to one of them (e.g. Hind or Xba)
}

# chromosome_choice=as.numeric(chromosome_choice)
# chromosome_choice=chromosome_choice[chromosome_choice>0]
##=======Statistical analysis=======
if(inputdata_format==2) #AF input
{
	CPA_calculation_choice=3
	AF_calculation_choice=3
	CPA_est_choice=3
	AF_est_choice=3
}
CPA_calculation_choice=CPA_calculation_option[as.numeric(CPA_calculation_choice)]  #1:Database ;2 User_provided
AF_calculation_choice=AF_calculation_option[as.numeric(AF_calculation_choice)]

AF_ref_calculation_choice=AF_ref_calculation_option[as.numeric(AF_ref_choice)]
  ##  AFref_resource_choice
  
QI_calculation_choice=QI_calculation_option[as.numeric(QI_calculation_choice)]
Iden_poor_chipchr_choice=Iden_poor_chipchr_option[as.numeric(Iden_poor_chipchr_choice)]

#QI_quantile_choice=QI_quantile_option[1]      ##0.95  0.975 0.99

##=======Result visualization=======
AF.plot_choice<-as.numeric(AF.plot_choice)
AF.based_choice<-as.numeric(AF.based_choice) #intensity-based or genotype-based or both
      ##------ Heatmap     ------##
HM_plot_choice<-as.numeric(HM_plot_choice);#QI.based_choice<-1 #only chip
      ##------ Polygon     ------##
PG_plot_choice=as.numeric(PG_plot_choice)
##=======Numerical output=======
data_description_choice=as.numeric(data_description_choice)
CPA_est_choice=as.numeric(CPA_est_choice)

AF_est_choice<-as.numeric(AF_est_choice);
AF_est.based_choice<-as.numeric(AF_est.based_choice) #intensity-based or genotype-based or both
QI_est_choice<-as.numeric(QI_est_choice);
if(QI_est_choice==2){QI_est.based_choice<-1} #only chip
if(QI_est_choice!=2){QI_est.based_choice<-2} #only chromosome? this option is removed

Poor_chipchr_choice<-as.numeric(poor_chipchr_choice)
if(Poor_chipchr_choice==2){poor_type_choice<-1}        #only chip
if(Poor_chipchr_choice!=2){poor_type_choice<-2}        #only chromosome? this option is removed

##--------------------------input parameters------------------------------------
##==============================================================================

##----------------------------- chip name setting ------------------------------
chipname_list=NULL;
if(chipsize_choice=="Affy 100K")
{
	if(chiptype_choice=="2chip"){chipname_list=c(" Hind"," Xba"," Merge")}
	if(chiptype_choice=="1chip")
	{
		if(onechip_choice=="chip1"){chipname_list=c(" Hind"," Hind"," Hind")}
		if(onechip_choice=="chip2"){chipname_list=c(" Xba"," Xba"," Xba")}
	}
}
if(chipsize_choice=="Affy 500K")
{
	if(chiptype_choice=="2chip"){chipname_list=c(" Nsp"," Sty"," Merge")}
	if(chiptype_choice=="1chip")
	{
		if(onechip_choice=="chip1"){chipname_list=c(" Nsp"," Nsp"," Nsp")}
		if(onechip_choice=="chip2"){chipname_list=c(" Sty"," Sty"," Sty")}
	}
}
if(chipsize_choice=="Affy 6.0"){chipname_list=c(" Affy6.0"," Affy6.0"," Affy6.0")}
if(chipsize_choice=="Affy Axiom"){chipname_list=c(" AffyAxiom"," AffyAxiom"," AffyAxiom")}
if(chipsize_choice=="Illumina 550K"){chipname_list=c(" Illumina 550K"," Illumina 550K"," Illumina 550K")}

##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
##==============================================================================
cat("   ~ The input options are finished.\n",file=logfile,append=T)
cat("\n",file=logfile,append=T)
cat("   ~ The establishment of temp directory is preparing.\n",file=logfile,append=T)
##==================== temp folder for calculation =============================
pathway_temp_population=paste(main_pathway,"TEMP",sep="");
dir.create(pathway_temp_population,showWarnings=F)
pathway_temp=paste(pathway_temp_population,"/",population_name,sep="")##<-----------------
dir.create(pathway_temp,showWarnings=F)
cat("   ~ The establishment of temp directory is finished.\n",file=logfile,append=T)
cat("\n",file=logfile,append=T)
cat("   ~ The databse is preparing.\n",file=logfile,append=T)
##====================database population choice================================
pathway_DB=paste(main_pathway,"DATABASE/",chipsize_choice,"/",population_database_choice,sep="")
pathway_DB_AFref=paste(pathway_DB,"/AF_ref",sep="")
#pathway_DB_AFprop=paste(pathway_DB,"/AF_prop",sep="")
pathway_DB_CPA=paste(pathway_DB,"/CPA",sep="")
pathway_DB_UB=paste(pathway_DB,"/Upper_bound",sep="") #upper bound of QI
cat("   ~ The databse is finished.\n",file=logfile,append=T)
cat("\n",file=logfile,append=T)
##==================== Subfunction =============================================
pathway_program=paste(main_pathway,"PROGRAM",sep="")
pathway_subfunction=paste(pathway_program,"/SAQC_subfunction.r",sep="")
source(pathway_subfunction)

################################################################################
################################################################################

###=============manager the genotype based data to be calculation form==========
pathway_output_population=paste(pathway_output,"/",population_name,sep="")
dir.create(pathway_output_population,showWarnings=F)

pathway_NR=paste(pathway_output_population,"/Numerical results",sep="")
dir.create(pathway_NR,showWarnings=F)
pathway_GR=paste(pathway_output_population,"/Graphical results",sep="")
dir.create(pathway_GR,showWarnings=F)

if(inputdata_format==1) #intensity-based
{
	cat("   ~ The input data format is genotype/intensity-based.\n",file=logfile,append=T)
	cat("     -- The setting of input data format is preparing.\n",file=logfile,append=T)

	# GenoPI=list.files(pathway_input)
	# GenoPI=mixedsort(GenoPI)
	# IndGeno=list.files(path=paste(pathway_input,"/",GenoPI[1],sep=""))
	# IndGeno_file=list.files(path=paste(pathway_input,"/",GenoPI[1],sep=""),full.names=T)
	# IndGeno=mixedsort(IndGeno); ## file name
	# IndGeno_file=mixedsort(IndGeno_file); ## full file name

	# # n=length(IndGeno);
	# #IndChipname=strsplit(IndGeno,"_");IndChipname=unlist(IndChipname)
	# # IndChipname = if(chipsize_choice == "Affy 100K"){ 
		# # # strsplit(IndGeno, "-50khind_mdy.txt|-50kxba_mdy.txt")
		# # sapply(strsplit(IndGeno, "_hind|_xba|_Hind|_Xba"), function(x) x[1]) ## "|" means "or"
	# # } else if(chipsize_choice == "Affy 500K") {
		# # # strsplit(IndGeno, "-250knsp_mdy.txt|-250ksty_mdy.txt")
		# # sapply(strsplit(IndGeno, "_nsp|_sty|_Nsp|_Sty"), function(x) x[1])
	# # }
	# IndChipname = sapply(strsplit(IndGeno, "_"), function(x) x[1])
	# # IndChipname=unlist(IndChipname)
	# #IndChipname=matrix(IndChipname,length(IndChipname)/n,n)
	# #IndChipname=t(IndChipname);temp_chipname_list=table(IndChipname[,2])

	# # if(temp_chipname_list[1]!=temp_chipname_list[2])
	# # {stop("The number of input chips are not equal.")}
	# # above, check the numbers of chips equal or not
	# if(Ind_Array<3){
		# if(!all(table(IndChipname) == 2)) stop("The number of input chips are not equal.")
	# }
	# #derive sample id
	# #n=length(IndGeno)/length(temp_chipname_list); #number of samples
	# n=length(unique(IndChipname)) #number of samples
	nchr=length(chromosome_choice) #number of selected chromosomes
	# # ind_name=list.files(path=paste(pathway_input,"/",GenoPI[1],sep=""))
	# # ind_name=mixedsort(ind_name);ind_name=strsplit(ind_name,"_")
	# # ind_name=unlist(ind_name);
	# ind_name = unique(IndChipname)
	# # ind_name=matrix(ind_name,length(ind_name)/n,n)
	# # ind_name=ind_name[1,]
	# sample_list=cbind(1:n,ind_name)
	# if(chipsize_choice %in% c("Affy 100K", "Affy 500K")){
		# IndPI=list.files(path=paste(pathway_input,"/",GenoPI[2],sep=""))
		# IndPI_file=list.files(path=paste(pathway_input,"/",GenoPI[2],sep=""),full.names=T)
		# IndPI=mixedsort(IndPI);
		# IndPI_file=mixedsort(IndPI_file);
	# }
	cat("     -- The setting of input data format is finished.\n",file=logfile,append=T)
	cat("\n",file=logfile,append=T)
}

##=============manager the AF based data to be calculation form=================
if(inputdata_format==2) #AF-based
{
	cat("   ~ The input data format is Allele frequency-based.\n",file=logfile,append=T)
	cat("     -- The setting of input data format is preparing.\n",file=logfile,append=T)
	ind_name=list.files(pathway_input)
	# pathway_output_population=paste(pathway_output,"/",population_name,sep="")
	# dir.create(pathway_output_population,showWarnings=F)

	# pathway_NR=paste(pathway_output_population,"/Numerical results",sep="")
	# dir.create(pathway_NR,showWarnings=F)
	# pathway_GR=paste(pathway_output_population,"/Graphical results",sep="")
	# dir.create(pathway_GR,showWarnings=F)

	pathway_AFest=pathway_input
	# dir.create(pathway_AFest,showWarnings=F)
	# pathway_input_ind=list.files(pathway_input,full.names=T)
	n=length(list.files(pathway_input,full.names=T))
	nchr=length(chromosome_choice)
	# for(k in 1:n) #for each sample
	# {
	# #k=1
	# pathway_ind_chr=list.files(pathway_input_ind[k],full.names=T)
	# pathway_ind_chr<-mixedsort(pathway_ind_chr)
	# pathway_ind_chr<-pathway_ind_chr[chromosome_choice]
	# pathway_output_ind=paste(pathway_AFest,"/",ind_name[k],sep="")
	# dir.create(pathway_output_ind,showWarning=F)
	# file.copy(pathway_ind_chr,pathway_output_ind,overwrite=T)
	# }
	sample_list=cbind(1:n,ind_name)
	cat("    -- The setting of input data format is finished.\n",file=logfile,append=T)
	cat("\n",file=logfile,append=T)
}

##------------------------------------------------------------------------------
cat("---------------------------------------------------------------------  \n");
cat("\n");
################################################################################
##                              Zone 1                                        ##
##                calculate the AF from intensity data                        ##
################################################################################
##==============================================================================
##111111111111111111111111111111111111111111111111111111111111111111111111111111
if(inputdata_format==1){
	path_ind=paste(pathway_temp,"/",population_name,"_ind",sep="")
	dir.create(path_ind,showWarnings=F)
}else{path_ind<-pathway_AFest}

# gc(TRUE,verbose=F)
      ##--------------Genotype based input format------------------
if(inputdata_format==1)
{
# start.time = Sys.time()
	#if(CPA_calculation_choice==2) #Yes, calculate CPA from input data
	if(CPA_calculation_choice == 2 | AF_calculation_choice == 2) #calculate CPA or estimate AF, we need to derive raw AF first
	{
		Feature_RData_ft(inpath = pathway_input, source_name = population_name, Ind_Array = Ind_Array, NA_str = NA_str, ParaSet = ParaSet, chrlist = chromosome_choice, outpath = main_pathway, chiptype_choice = chiptype_choice, onechip_choice = onechip_choice, logfile = logfile)
		Feature_Stat_Fast_ft(source_name = population_name, Ind_Array = Ind_Array, chrlist = chromosome_choice,  inpath = pathway_input, outpath = main_pathway)
		# ind_name=list.files(path_ind,full.names=F)
		
		
	}
	ind_name=mixedsort(list.files(path_ind,full.names=F)) 
	n = length(ind_name)
	# end.time = Sys.time()
	# cat(end.time-start.time, "\n", sep = "")
}
gc(reset=TRUE,verbose=F)
cat("****************************************************************************************************************\n")

##==============================================================================

##==============================================================================


################################################################################
##                        Data description                                    ##
################################################################################

cat("    ******      ****    ==========               Data description               ==========   ****       ****\n")
cat("   ~ Data description is preparing.\n",file=logfile,append=T)
pathway_output_population=paste(pathway_output,"/",population_name,sep="")
dir.create(pathway_output_population,showWarnings=F)
# file.copy(from = paste(pathway_temp, "/samplelist.txt", sep = ""), to = paste(pathway_output_population, "/Sample list.txt", sep = ""), recursive = T)
sample_list=cbind(1:n,ind_name)
colnames(sample_list)[1]="Obs";
colnames(sample_list)[2]="Ind_ID";
write.table(sample_list,file=paste(pathway_output_population,"/Sample list.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

pathway_NR=paste(pathway_output_population,"/Numerical results",sep="")
dir.create(pathway_NR,showWarnings=F)

if(data_description_choice==2)
{
	# gc(TRUE,verbose=F)
	data_descrption=paste(pathway_NR,"/Data_description.txt",sep="")
	single_line=c("-----------------------------------------------------");
	double_line=c("=========================================================");

	title_part=paste("==========      Parameter/Data Description     ==========",sep="")
	title_part=rbind(double_line,title_part,double_line)

	write.table(title_part,file=data_descrption,row.names =F,quote = F,col.names =F,append=F,sep="\t")
	write.table(c(" "),file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	## 1. ===Directories===
	write.table("1. Directories - ",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	if(inputdata_format==1){inputdata_format_name="Genotype/Intensity-based"}
	if(inputdata_format==2){inputdata_format_name="Allele Frequency-based"}
	input_data_format=paste("   Input data format: ",inputdata_format_name,sep="")
	input_way=paste("   The directory name of input data:  ",pathway_input,sep="")
	output_way=paste("   The directory name of output data: ",pathway_output,sep="")
	write.table(input_way,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(output_way,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(c(" "),file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")


	## 2. ===study group===
	write.table("2. Options -",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	temp_population=paste("   Study population: ",population_choice,sep="");
	write.table(temp_population,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	temp_chiptype=paste("   The type of SNP array: ",chipsize_choice,sep="")
	write.table(temp_chiptype,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	# if(chipsize_choice=="Affy 100K")
	# {
		# if(chiptype_choice=="1chip")
		# {
			# if(onechip_choice=="chip1"){chip_name="Hind"};
			# if(onechip_choice=="chip2"){chip_name="Xba"};
		# }
	# }

	# if(chipsize_choice=="Affy 500K")
	# {
		# if(chiptype_choice=="1chip")
		# {
			# if(onechip_choice=="chip1"){chip_name="Nsp"};
			# if(onechip_choice=="chip2"){chip_name="Sty"};
		# }
	# }
  # if(chipsize_choice=="Affy 6.0"){chip_name="One chip"}
  # if(chipsize_choice=="Illumina 550K"){chip_name="One chip"}
  
	if(chiptype_choice=="1chip")
		temp_chip=paste("   The number of selected array: One array analysis  (",chipname_list[1],")",sep="")
	if(chiptype_choice=="2chip")
		temp_chip=paste("   The number of selected array: Two array analysis",sep="")

	write.table(temp_chip,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	# chr_chromosome=paste(matrix(chromosome_choice,1,length(chromosome_choice)),sep=",")
	# chr_list=chromosome_choice[1];
	# if(length(chromosome_choice)>1)
	# {for(ichr in 2:length(chromosome_choice)){chr_list=paste(chr_list,",",chromosome_choice[ichr],sep="")}}
	# chr_chromosome=paste("   The selected chromosomes: ",chr_list,sep="")
	# write.table(chr_chromosome,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep=",")
	# write.table(c(" "),file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	##===data format===
	write.table("2' Samples and markers -",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	ind_example=list.files(path_ind)

	chr=paste(path_ind,"/",ind_example[1],sep="")
	chr=list.files(chr,full.names=T)
	chr=mixedsort(chr)
	num_chip1=NULL;num_chip2=NULL;
	for(i in 1:length(chr))
	{
		#i=1
		a=read.delim(chr[i],head=T,sep="\t", stringsAsFactors = F)
		SNP_eachchip=table(a[,3])
		num_chip1=c(num_chip1,SNP_eachchip[1])
		num_chip2=c(num_chip2,SNP_eachchip[2])
	}
	num_ind=paste("   The number of individuals: ",n,sep="")
	write.table(num_ind,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	whole_SNP_chip=c(sum(num_chip1),sum(num_chip2))
	if(chipsize_choice=="Affy 100K")
	{
		temp_chip_name=c("Hind","Xba")
		if(chiptype_choice=="1chip")
		{
			if(onechip_choice=="chip1"){whole_SNP_chip[2]=0};
			if(onechip_choice=="chip2"){whole_SNP_chip[1]=0};
		}
	}
	if(chipsize_choice=="Affy 500K")
	{
		temp_chip_name=c("Nsp","Sty")
		if(chiptype_choice=="1chip")
		{
			if(onechip_choice=="chip1"){whole_SNP_chip[2]=0};
			if(onechip_choice=="chip2"){whole_SNP_chip[1]=0};
		}
	}
	if(chipsize_choice %in% c("Affy 6.0" , "Affy Axiom", "Illumina 550K")){temp_chip_name=c("Chip","Chip");whole_SNP_chip[2]=NA}

	sum_chip1=paste("   The number of SNPs in ",temp_chip_name[1],": ",whole_SNP_chip[1],sep="")
	write.table(sum_chip1,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	if(chipsize_choice=="Affy 100K" | chipsize_choice=="Affy 500K")
	{
		sum_chip2=paste("   The number of SNPs in ",temp_chip_name[2],": ",whole_SNP_chip[2],sep="")
		write.table(sum_chip2,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	}
	sum_merge=paste("   The number of SNPs in whole chip: ",sum(whole_SNP_chip,na.rm=T),sep="")
	write.table(sum_merge,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(" ",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

  ##===statistical analysis===
	write.table("3. Statistical analysis -",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	temp_CPA=YesNo_option(CPA_calculation_choice)
	if(CPA_calculation_choice==1)
		temp_CPA=paste(temp_CPA," (using CPA in database)",sep="")
	temp_CPA=paste("   CPA calculation: ",temp_CPA,sep="")
	write.table(temp_CPA,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	temp_AF=YesNo_option(AF_calculation_choice)
	temp_AF=paste("   AF calculation: ",temp_AF,sep="")
	write.table(temp_AF,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	temp_QI=YesNo_option(QI_calculation_choice)
	temp_QI=paste("   QI calculation: ",temp_QI,sep="")
	write.table(temp_QI,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	temp_AF_ref=YesNo_option(AF_ref_calculation_choice)
	if(AF_ref_calculation_choice==1){
		if(AFref_source_choice[1]){
			temp_AF_ref=paste(temp_AF_ref," (using AF reference in SAQC database)",sep="")
		} else temp_AF_ref=paste(temp_AF_ref," (using AF reference in user-provided database)",sep="")
	} else temp_AF_ref=paste(temp_AF_ref,sep="")

    #temp_AF_ref=paste(temp_AF_ref," (using QI reference in database)",sep="")
	temp_QI_ref=paste("   AF reference calculation: ",temp_AF_ref,sep="")
	write.table(temp_QI_ref,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	temp_iden=YesNo_option(Iden_poor_chipchr_choice)
	temp_QI_quantile=paste(as.character(QI_quantile_choice*100),"%",sep="")
	if(QI_quantile_choice!=0){temp_iden=paste("   Identification of poor chip: ",temp_iden," (QI upper quantile: ",temp_QI_quantile,")",sep="")}
	if(QI_quantile_choice==0){temp_iden=paste("   Identification of poor chip: ",temp_iden,sep="")}

	write.table(temp_iden,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(c(" "),file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")


	##===Result visualization===
	write.table("4. Graphical result -",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	if(AF.plot_choice!=3){
		temp = YesNo_option(AF.plot_choice)
		if(temp == "Yes"){
			temp_AFplot=paste("   AF plot: ", temp, " (",based_option(AF.based_choice),")",sep="")
		} else temp_AFplot=paste("   AF plot: ", temp, sep="")	
	}
	if(AF.plot_choice==3){temp_AFplot=paste("   AF plot: ",YesNo_option(AF.plot_choice),sep="")}

	temp_CRplot=paste("   GCR plot: ",YesNo_option(CallRateplot_choice),sep="")
	if(HM_plot_choice!=3){
		HM_col = sum(col.scale_choice)
		HM_label = if(HM_col == 2){
			"SAQC database and user's data"
		} else if(col.scale_choice[1] == 1){
			"SAQC database"
		} else "user's data"
		temp_HM_plot=paste("   QI HeatMap plot: ",YesNo_option(HM_plot_choice)," (", "color scale based on ", HM_label,")",sep="")
	} else {
		temp_HM_plot=paste("   QI HeatMap plot: ",YesNo_option(HM_plot_choice),sep="")
	}
	temp_QISP=paste("   QI Polygon plot : ",YesNo_option(PG_plot_choice),sep="")

	write.table(temp_AFplot,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_CRplot,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_HM_plot,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_QISP,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(c(" "),file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	##===Numerical output===
	write.table("5. Numerical output -",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	temp_data_description=paste("   Parameter/Data description: ",YesNo_option(data_description_choice),sep="")
	temp_CPA_est=paste("   CPA estimate: ",YesNo_option(CPA_est_choice),sep="")
	temp_AF_est=paste("   AF estimate: ",YesNo_option(AF_est_choice),sep="")#" (",based_option(AF_est.based_choice),")",sep="")
	temp_QI_est=paste("   QI estimate: ",YesNo_option(QI_est_choice),sep="")
	temp_poor_chipchr=paste("   Poor SNP array: ",YesNo_option(Poor_chipchr_choice),sep="")
	# temp_poor_type_choice=as.numeric(poor_type_choice)
	#temp_poor_type_choice=temp_poor_type_choice%*%c(1,2)
	# stat_poor_type=NULL;

	# if(sum(temp_poor_type_choice)==1)
	# stat_poor_type="Poor chip";
	# if(sum(temp_poor_type_choice)!=1)
	# stat_poor_type="NA";
	# if(poor_chipchr_choice!=3){temp_poor_chipchr=paste("   Poor chip: ",YesNo_option(Poor_chipchr_choice)," (",stat_poor_type,")",sep="")}
	# if(Poor_chipchr_choice==3){temp_poor_chipchr=paste("   Poor chip: ",YesNo_option(Poor_chipchr_choice),sep="")}


	write.table(temp_data_description,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_CPA_est,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_AF_est,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_QI_est,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(temp_poor_chipchr,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
	write.table(c(" "),file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")

	# gc(reset=TRUE)
}
cat("   ~ Data description is finished.\n",file=logfile,append=T)
cat("\n",file=logfile,append=T)
##===============================================================================================


################################################################################
##                                Part 1                                      ##
##                                  CPA                                       ##
##              CPA_choice==1 meaning the CPA from the datacase               ##
##              CPA_choice==2 meaning the CPA from the original data          ##
##============================================================================##

#CPA_calculation_choice=1        ## As the user choose, the result will be output the specific CPA
if(CPA_calculation_choice==1)
{
    #pathway_CPA=paste(pathway_NR,"/CPA estimate",sep="");##<-----------------
    #dir.create(pathway_CPA,showWarnings=F)
    #file.copy(pathway_DB_AFref, pathway_CPA, overwrite = TRUE)
    #file.copy(list.files(pathway_DB_CPA, full.names = T), pathway_CPA, overwrite = TRUE)
	pathway_CPA = pathway_DB_CPA
}
#CPA_calculation_choice=2
if(CPA_calculation_choice==2)    ## calculate the CPA from the original data
{
    cat("    **	       **  **	==========               CPA calculation                ==========  **  **     ** \n")
    cat("   ~ CPA calculation is preparing.\n",file=logfile,append=T)
    ##==============================================================================
    ##      CPA calculation
    ##  <Part 2> ±N©Ò¦³individualªºGenotype ,Allel frequency ©M individual information
    ##           ¤À§O¦X¨Ö¤@°_¡A¬Goutputªº³¡¥÷¦b¨C¤@­Óchromosome·|¦³¤T­ÓÀÉ®×
    ##
    ##  <Part 3> ­pºâCPA Ku,Kh
    ##==============================================================================
    ##222222222222222222222222222222222222222222222222222222222222222222222222222222
    #ind_name=list.files(path_ind)
    ind_name=mixedsort(list.files(path_ind))
    pathway_GAF=paste(pathway_temp,"/",population_name,"_GAF",sep="");##<-----------------
    dir.create(path=pathway_GAF,showWarnings = F) ##GAF mean genotype and allele frequency
    #chromosome_choice=c(1,2,19,22,23)
    # if(chromosome_choice[length(chromosome_choice)]!=23)
    # {
		# for(k in 1:length(chromosome_choice))
		# {
			# dir.create(path=paste(pathway_GAF,"/Chr_",chromosome_choice[k],sep=""),showWarnings = F)
			# whole_genotype=NULL;
			# whole_AF=NULL;
			# for(i in 1:length(ind_name))
			# {
				# chr_file=list.files(path=paste(path_ind,"/",ind_name[i],sep=""),full.names=T)
				# chr_file=mixedsort(chr_file[1:length(chr_file)])
				# chr=read.delim(chr_file[k],fill=T,header=T,sep="\t")
				# whole_genotype=cbind(whole_genotype,as.character(chr[,4]));
				# whole_AF=cbind(whole_AF,chr[,5]);
			# }
			# write.table(whole_genotype,file=paste(pathway_GAF,"/Chr_",chromosome_choice[k],"/genotype.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
			# write.table(whole_AF,file=paste(pathway_GAF,"/Chr_",chromosome_choice[k],"/allele_frequency.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
			# write.table(chr[,1:3],file=paste(pathway_GAF,"/Chr_",chromosome_choice[k],"/chiptype.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
			# #Probe_Set Phy_position chiptype
		# }
    # }
    # if(chromosome_choice[length(chromosome_choice)]==23)
    # {
		# if(length(chromosome_choice)>1)
		# {
			# for(k in 1:(length(chromosome_choice)-1))
			# {
				# dir.create(path=paste(pathway_GAF,"/Chr_",chromosome_choice[k],sep=""),showWarnings = F)
				# whole_genotype=NULL;
				# whole_AF=NULL;
				# for(i in 1:length(ind_name))
				# {
					# chr_file=list.files(path=paste(path_ind,"/",ind_name[i],sep=""),full.names=T)
					# chr_file=mixedsort(chr_file[1:(length(chr_file))])
					# chr=read.delim(chr_file[k],fill=T,header=T,sep="\t")
					# whole_genotype=cbind(whole_genotype,as.character(chr[,4]));
					# whole_AF=cbind(whole_AF,chr[,5]);
				# }
				# write.table(whole_genotype,file=paste(pathway_GAF,"/Chr_",chromosome_choice[k],"/genotype.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
				# write.table(whole_AF,file=paste(pathway_GAF,"/Chr_",chromosome_choice[k],"/allele_frequency.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
				# write.table(chr[,1:3],file=paste(pathway_GAF,"/Chr_",chromosome_choice[k],"/chiptype.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
			# }
		# }
		# pathway_X_female=paste(pathway_GAF,"/Chr_23",sep="");
		# dir.create(pathway_X_female,showWarnings = F) 
		# pathway_X_male=paste(pathway_GAF,"/Chr_24",sep="");
		# dir.create(pathway_X_male,showWarnings = F) 


		# whole_genotype_female=NULL;  whole_genotype_male=NULL;
		# whole_AF_female=NULL;  whole_AF_male=NULL;
		# for(i in 1:length(ind_name))
		# {
			# chr_file=list.files(path=paste(path_ind,"/",ind_name[i],sep=""),full.names=T)
			# chr_file=mixedsort(chr_file[1:length(chr_file)])
			# chr=read.delim(chr_file[length(chr_file)],fill=T,header=T,sep="\t")

			# AB_code=table(chr[,4])[2]
			# AB_rate=AB_code/length(chr[,4])
			# if(AB_rate>=0.1) #add "="
			# {
				# whole_genotype_female=cbind(whole_genotype_female,as.character(chr[,4]));
				# whole_AF_female=cbind(whole_AF_female,chr[,5]);
			# }
			# if(AB_rate<0.1)
			# {
				# whole_genotype_male=cbind(whole_genotype_male,as.character(chr[,4]));
				# whole_AF_male=cbind(whole_AF_male,chr[,5]);
			# }
		# }

		# #-------------Chromosome X output---------------------------
		# write.table(whole_genotype_female,file=paste(pathway_X_female,"/genotype.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
		# write.table(whole_AF_female,file=paste(pathway_X_female,"/allele_frequency.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
		# write.table(chr[,1:3],file=paste(pathway_X_female,"/chiptype.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

		# #-------------Chromosome Y output---------------------------
		# write.table(whole_genotype_male,file=paste(pathway_X_male,"/genotype.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
		# write.table(whole_AF_male,file=paste(pathway_X_male,"/allele_frequency.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
		# write.table(chr[,1:3],file=paste(pathway_X_male,"/chiptype.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
    # }
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##==============================================================================

    ##===============================================================================
    ##3
    ##   CPA process
    # rm(whole_genotype,whole_AF,chr)
    # rm(whole_genotype_male,whole_AF_male)
    # rm(whole_genotype_female,whole_AF_female)

    chr_name=list.files(pathway_GAF)
    chr_name=mixedsort(chr_name)
#    pathway_CPA=paste(pathway_temp,"/",population_name,"_CPA",sep="");##<-----------------
    pathway_CPA=paste(pathway_NR,"/CPA estimate",sep="");##<-----------------
    dir.create(pathway_CPA,showWarnings=F)
    nchr=length(chromosome_choice)
    for(u in 1:nchr)
    {
		i=chromosome_choice[u]
		chr=list.files(path=paste(pathway_GAF, chr_name[i],sep="/"),full.names=T)
		lapply(chr, load, envir=.GlobalEnv)
		##----modify the approporate CPA-----------------------
		CPA=kh_ku(AF,geno)
		chr_kh=CPA[,1]#Kh
		chr_ku=CPA[,2]#Ku HeterN
		# CPA=cbind(chr_kh,chr_ku) #Kh Ku HeterN
		CPA[CPA[,3]==1,2]=CPA[CPA[,3]==1,1]
		CPA[CPA[,3]==0,2]=1
		##-----------------------------------------------------
		heterN=CPA[,3];
		CPA=CPA[,2] #Ku_modified
		chr_CPA=cbind(chip,heterN,chr_kh,chr_ku,CPA) #Probe_Set Phy_position chiptype HeterN Kh Ku Ku_modified
		colnames(chr_CPA)[5]="K_h" #Probe_Set Phy_position chiptype HeterN K_h K_u Ku_modified
		colnames(chr_CPA)[6]="K_u"
		chr_CPA=chr_CPA[,c(1,2,3,4,7)]#Probe_Set Phy_position chiptype HeterN Ku_modified
		chr_CPA[,5]=sprintf("%1.4f",chr_CPA[,5]) #Ku_modified
		chr_CPA=chr_CPA[,c(1,2,5,3)] #Probe_Set	Phy_position	Ku_modified	chiptype
		# if(i<10){ichr=paste("0",i,sep="")}
		# if(i>9){ichr=paste(i,sep="")}
		ichr = ifelse(i<10, paste("0",i,sep=""), i)
		write.table(chr_CPA,file=paste(pathway_CPA,"/CPA_Chr_",ichr,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
    }
    cat("   ~ CPA calculation is finished.\n",file=logfile,append=T)
    cat("\n",file=logfile,append=T)
    ##---------------remove the GAF file----------------
    ##--------------------------------------------------
}
##==============================================================================


################################################################################
##                            Part 2                                          ##
##                calculate the allele frequency                              ##
##----------------------------------------------------------------------------##
## there are two choices                                                      ##
## one is that the CPA can be calcuated from original data                    ##
## another is that it can be obtained from the database                       ##
##  output file:  population_name_adjAF                                       ##
##  if the input file format is "GI based", the AF have been calculated       ##
##  if the input file format is "AF based", the AF from the input file        ##
##============================================================================##

#CPA_choice=2
#ind_name=list.files(path_ind)
ind_name=mixedsort(list.files(path_ind))

#if(inputdata_format==1){path_ind=paste(pathway_temp,"/AF estimate",sep="")};##<-----------------
#if(inputdata_format==2){path_ind}

pathway_adjAF=paste(pathway_NR,"/AF estimate",sep="");##<-----------------
dir.create(pathway_adjAF,showWarnings=F)


cat("    **         **  **   ==========                AF adjustment       	        ==========  **  **     ** \n")
	##  inputdata_format Genotype/Intensity based
	##  CPA from data base
	##  don't calculate AF
if(inputdata_format==1 & CPA_calculation_choice==1 & AF_calculation_choice==1)
{ #only re-combine genotype data based on raw genotype data, did not calculate AF 
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	cat("Genotype/Intensity based and Do not output AF result!!\n")
	GI_file=paste(pathway_input,"/IndGeno",sep="")
	GIind_file=list.files(GI_file,full.names=T)
	GIind_file=mixedsort(GIind_file)
	#n=length(GIind_file)/2
	n = length(ind_name)
	Genotype_temp=paste(pathway_NR,"/AF estimate",sep="")
	dir.create(Genotype_temp,showWarnings=F)
	for(k in 1:n)
	{
		tempk=k
		temp_ind_genotype=paste(Genotype_temp,"/",ind_name[k],sep="")
		dir.create(temp_ind_genotype,showWarnings=F)
		if(Ind_Array %in% c(1,2)){
			#chip1=read.delim(GIind_file[2*(tempk-1)+1],head=T,sep="\t")
			chip1=read.delim(GIind_file[2*(tempk-1)+1],header=F, skip = ParaSet[1], stringsAsFactors = F);
			chip1 = chip1[, c(2, 3, 4, 6)]
			chip1 = cbind(chip1,"Chiptype" = "Chip1")

			#colnames(chip1)[5]="Chiptype"
			#chip2=read.delim(GIind_file[2*(tempk-1)+2],head=T,sep="\t")
			chip2=read.delim(GIind_file[2*(tempk-1)+2],header=F, skip = ParaSet[1], stringsAsFactors = F);
			chip2 = chip2[, c(2, 3, 4, 6)]
			chip2=cbind(chip2,"Chiptype" = "Chip2")
			#colnames(chip2)[5]="Chiptype"

			ind=rbind(chip1,chip2)
		}else{
			chip1=read.delim(GIind_file[tempk],header=F, skip = ParaSet[1], stringsAsFactors = F);
			chip1 = chip1[, c(2, 3, 4, 6)]
			chip1 = cbind(chip1,"Chiptype" = "Chip1")
			ind=chip1
		}
		geno_name=as.character(names(table(ind[,2])))
		geno_name=geno_name[-grep(" ", geno_name)] 
		geno_name=mixedsort(geno_name)
		nchr=length(geno_name)

		colnames(ind)[1]="Probe_set"
		colnames(ind)[2]="Chr"
		colnames(ind)[3]="Phy_position"
		colnames(ind)[4]="Genotype"
		colnames(ind)[5]="Chiptype"
		for(j in 1:nchr)
		{
			#j=18
			jchr=as.character(geno_name[j])
			chr=ind[ind[,2]==jchr,]
			jchr = if(jchr != "X"){
				as.numeric(jchr)
			} else 23
			chr=chr[,c(1,3,5,4,2)] #Probe_set	Phy_position	Chiptype	Genotype	Chr
			nsnp=nrow(chr)
			probe=as.character(chr[,1]) #<--remove SNP_A-
			# probe=strsplit(probe,"SNP_A-");probe=unlist(probe);probe=probe[(1:nsnp)*2]
			#probe=substr(probe,1,7);
			probe=as.numeric(as.matrix(probe))
			chr[,1]=probe
			chr=sort.append(chr,2)

			if(j<10){jchr=paste("0",j,sep="")}
			if(j>9){jchr=paste(j,sep="")}
			#if(j==23){jchr=paste("23",sep="")}
			#chr_name=paste("Chr_",as.character(jchr),".txt",sep="")
			chr_name=paste("Chr_",jchr,".txt",sep="")
			ind_chr=paste(temp_ind_genotype,"/",chr_name,sep="")
			write.table(chr,ind_chr,col.names=T,row.names=F,quote=F,sep="\t")
		}
	}
}
# cat("break-1", sep = "\n")
	##  inputdata_format Genotype/Intensity based
	##  CPA from data base
	##  don't calculate AF
if(inputdata_format==1 & CPA_calculation_choice==2 & AF_calculation_choice==1)
{	#just copy unadjusted AF to Numerical output
	f=paste(pathway_NR,"/AF estimate",sep="")
	dir.create(NR_AF,showWarnings=F)
	AF_temp=paste(pathway_temp,"/",population_name,"_ind",sep="") #directory with unadjusted AF
	# pathway_temp_population=paste(main_pathway,"/TEMP",sep="");
	# pathway_temp=paste(pathway_temp_population,"/",population_name,sep="")

	ind_file=list.files(AF_temp,full.names=T)
	ind_name=list.files(AF_temp,full.names=F)
	nind=length(ind_name)
	NR_ind=paste(NR_AF,"/",ind_name,sep="")

	for(j in 1:nind)
	{
		dir.create(NR_ind[j],showWarnings=F)
		file.copy(list.files(ind_file[j],full.names=T),NR_ind[j],overwrite=TRUE) #just copy unadjusted AF to Numerical output
	}
}
# cat("break-2", sep = "\n")
if(inputdata_format==1 & AF_est_choice==2) 
{ #AF_est_choice: AF numerical output. If AF_est_choice == 2, it equals AF_calculation_choice == 2
	cat("   ~ Adjust allele frequency calculation is preparing.\n",file=logfile,append=T)
	##--------------------database-------------
	# if(CPA_calculation_choice==1)
	# {
	# CPA_database=list.files(pathway_DB_CPA,full.names=T)
	# chr_file=mixedsort(CPA_database);
	# }
	#--------------------user Provided
	# if(CPA_calculation_choice==2)
	# {
	# chr_file=list.files(pathway_CPA,full.names=T)
	# chr_file=mixedsort(chr_file)
	# }
	chr_file=mixedsort(list.files(pathway_CPA,full.names=T))
	for(k in 1:length(ind_name))
	{
		#k=2
		cat("     ~ The adjust allele frequency calculation of",paste("individual ",k,sep=""),"is preparing.\n",file=logfile,append=T)
		pathway_adjAF_ind=paste(pathway_adjAF,"/",ind_name[k],sep="");
		dir.create(pathway_adjAF_ind,showWarnings=F)

		ind_chr=list.files(paste(path_ind,"/",ind_name[k],sep=""),full.names=T) #read raw AF
		ind_chr=mixedsort(ind_chr)
		for(i in 1:nchr)
		{
			#i=4
			cat("        ~ The adjust allele frequency calculation of",paste("individual ",k," in chromosome ",i,sep=""),"is preparing.\n",file=logfile,append=T)
			ichr=chromosome_choice[i]
			ind=read.delim(ind_chr[i],fill=T,header=T,sep="\t", stringsAsFactors = F)
			##------------- ind modification-------------------------------
			#probe=as.character(ind[,1]);probe=strsplit(probe,"SNP_A-");probe=unlist(probe); probe=probe[2*c(1:(length(probe)/2))];ind[,1]=probe
			if(CPA_calculation_choice==1)     #database
			{CPA=read.delim(chr_file[ichr],fill=T,header=T,sep="\t", stringsAsFactors = F)}
			#ichr is absolute index
			if(CPA_calculation_choice==2)     #user-provided
			{CPA=read.delim(chr_file[i],fill=T,header=T,sep="\t", stringsAsFactors = F)}
			#i is relative index
			colnames(CPA)[1]="Probe_set"

			#ind=merge(ind,CPA[,c(1,3)],by="Probe_set")    ##ind {Probe_set Phy_position Chiptype Genotype AF CPA}
			ind=merge(ind,CPA[,c(1,3)],by="Probe_set", all.x = T) 
			ind$CPA[is.na(ind$CPA)] = 1 #replace CPA of SNPs which are NA by 1

			adjustAF=ind[,5]/(ind[,5]+ind[,6]*(1-ind[,5]))
			ind=cbind(ind,adjustAF)
			#colnames(ind)[6]="CPA"
			ind=ind[,c(1,2,3,4,7)]
			colnames(ind)[2]="Phy_position"
			colnames(ind)[3]="Chiptype"
			colnames(ind)[5]="adjAF"

			ind = sort.append(ind,2)
			  
			# ind[,1]=sprintf("%9s",ind[,1])
			# ind[,2]=sprintf("%9s",ind[,2])
			# ind[,3]=sprintf("%8s",ind[,3])
			# ind[,4]=sprintf("%8s",ind[,4])
			# ind[,5]=sprintf("%1.4f",ind[,5])

			##--- remove the AF estimate numerical result ---
			#temp_AF_est.based_choice=AF_est.based_choice%*%c(1,2)
			#if(sum(temp_AF_est.based_choice)==1){ind=ind[,c(1,2,3,5)]}
			#if(sum(temp_AF_est.based_choice)==2){ind=ind[,c(1,2,3,4)]}
			# ichr=i
			if(ichr<10){Ichr=paste("0",ichr,sep="")}
			if(ichr>9){Ichr=paste(ichr,sep="")}
			write.table(ind,file=paste(pathway_adjAF_ind,"/Chr_",Ichr,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
			# pathway_adjAF=paste(pathway_NR,"/AF estimate",sep=""); under numerical output
			cat("        ~ The adjust allele frequency calculation of",paste("individual ",k," in chromosome ",i,sep=""),"is finished.\n",file=logfile,append=T)
        }
		cat("     ~ The adjust allele frequency calculation of",paste("individual ",k,sep=""),"is finished.\n",file=logfile,append=T)
        cat("\n",file=logfile,append=T)
	}
}
# cat("break-3", sep = "\n")
    ##--------------------- input file is AF based----------------------------------
if(inputdata_format==2 & AF_est_choice==3)
{
	cat("   ~ Adjust allele frequency calculation is preparing.\n",file=logfile,append=T)
	AFref_file=list.files(paste(pathway_DB,"/AF_ref",sep=""),full.names=T)
	##--------------------database-------------
	for(k in 1:length(ind_name))
	{
		#k=1
		cat("     ~ The adjust allele frequency modification of",paste("individual ",k,sep=""),"is preparing.\n",file=logfile,append=T)
		pathway_input_adjAF_ind=paste(pathway_input,"/",ind_name[k],sep="");
		ind_chr=list.files(pathway_input_adjAF_ind,full.names=T)
		ind_chr=mixedsort(ind_chr)

		pathway_adjAF_ind=paste(pathway_adjAF,"/",ind_name[k],sep="");
		# unlink(pathway_adjAF_ind,recursive = T)
		dir.create(pathway_adjAF_ind,showWarnings=F)

		for(i in 1:nchr)
		{
			#i=1
			cat("        ~ The adjust allele frequency calculation of",paste("individual ",k," in chromosome ",i,sep=""),"is preparing.\n",file=logfile,append=T)
			ichr=chromosome_choice[i]
			ind=read.delim(ind_chr[i],fill=T,header=T,sep="\t", stringsAsFactors = F) ##ind {Probe_set Phy_position Chiptype Genotype AF}
			nsnp=nrow(ind)
			colnames(ind)[1]="Probe_set"
			##-------------------------------------------------------------
			##------------- ind modification-------------------------------
			probe=as.character(ind[,1]) # ;probe=strsplit(probe,"SNP_A-");probe=unlist(probe);
			if(length(probe)==2*nsnp){probe=probe[2*c(1:(length(probe)/2))]};
			ind[,1]=probe
			colnames(ind)[2]="Phy_position"
			# colnames(ind)[3]="chiptype"
			colnames(ind)[3]="Chiptype"
			colnames(ind)[5]="adjAF"
			ind[,1]=sprintf("%9s",ind[,1])
			ind[,2]=sprintf("%9s",ind[,2])
			ind[,3]=sprintf("%8s",ind[,3])
			ind[,4]=sprintf("%8s",ind[,4])
			ind[,5]=sprintf("%1.4f",ind[,5])
			#ichr=i
			if(ichr<10){Ichr=paste("0",ichr,sep="")}
			if(ichr>9){Ichr=paste(ichr,sep="")}
			write.table(ind,file=paste(pathway_adjAF_ind,"/Chr_",Ichr,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
			cat("        ~ The adjust allele frequency calculation of",paste("individual ",k," in chromosome ",i,sep=""),"is finished.\n",file=logfile,append=T)
		}
        cat("     ~ The adjust allele frequency calculation of",paste("individual ",k,sep=""),"is finished.\n",file=logfile,append=T)
        cat("\n",file=logfile,append=T)
    }
}
# cat("break-4", sep = "\n")

##--------------------NoCall/ Call Rate calculation-----------------------------
cat("     ~ The Call rate calculation is preparing.\n",file=logfile,append=T)
pathway_Sample_file=paste(pathway_output_population,"/Sample list.txt",sep="")
#Sample_list=read.delim(pathway_Sample_file,head=T,sep="\t")
Sample_list=read.delim(pathway_Sample_file,header=T,sep="\t", colClasses = c("integer", "character"))
#pathway_adjAF=paste(pathway_NR,"/AF estimate",sep="");
CallR=Call_Ratio(pathway_adjAF) ##<----- it have been changed to calculate Call Rate
Sample_list=merge(Sample_list,CallR,by="Ind_ID")
Sample_list=Sample_list[,c(2,1,3,4,5)]
Sample_list=sort.append(Sample_list,1)
colnames(Sample_list)[2]="IndID" 
Sample_list_temp1 = Sample_list
# Obs	IndID           CR_Nsp	CR_Sty	CR_Merge
# 1	01	0.9912	0.9864	0.9889
# 2	02	0.9952	0.9909	0.9932
# 3	14	0.9112	0.9147	0.9129
# 4	16	0.9094	0.8846	0.8976
if(chiptype_choice=="1chip")
{
	if(onechip_choice=="chip2"){
		Sample_list=Sample_list[,-c(3,5)]
	}else{
		Sample_list=Sample_list[,-c(4,5)]
	}
}
write.table(Sample_list,pathway_Sample_file,col.names=T,row.names=F,quote=F,sep="\t")
cat("     ~ The Call rate calculation is finished.\n",file=logfile,append=T)
cat("\n",file=logfile,append=T)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##==============================================================================




################################################################################
##                          Zone 4                                            ##
##                  AF-physcial ploting                                       ##
##----------------------------------------------------------------------------##
pathway_program=paste(main_pathway,"PROGRAM",sep="")
pathway_subfunction=paste(pathway_program,"/SAQC_subfunction.r",sep="")
source(pathway_subfunction)

pathway_GR=paste(pathway_output_population,"/Graphical results",sep="")
dir.create(pathway_GR,showWarnings=F)
# cat("break-5", sep = "\n")
##-----------------------------------Call figure--------------------------------
if(CallRateplot_choice==2)
{

    cat("     ~ The ploting of Call rate calculation is preparing.\n",file=logfile,append=T)
    pathway_GR_CR=paste(pathway_GR,"/CR plot",sep="") #GR stands for graphical output
    dir.create(pathway_GR_CR,showWarnings=F)

    pathway_Sample_file=paste(pathway_output_population,"/Sample list.txt",sep="")
    CR=read.delim(pathway_Sample_file,header=T,sep="\t", colClasses = c("integer", "character", rep("numeric", ncol(Sample_list)-2)))  #Call Rate
    if(chiptype_choice=="1chip") CR=cbind(CR,CR[,3],CR[,3]) #nind x 5

	##--------Call Rate figure
	pathway_Call_figure=paste(pathway_GR_CR,"/Call Ratio.png",sep="")
	png(pathway_Call_figure,width=1440,height=960)
	Call_plot(CR,chiptype_choice,chipname_list)
	dev.off()

	##--------NoCall Rate figure
	NCR=cbind(CR[,c(1,2)],1-matrix(as.matrix(CR[,-c(1,2)]),nrow(CR),3))
	pathway_Call_figure=paste(pathway_GR_CR,"/NoCall Ratio.png",sep="")
	png(pathway_Call_figure,width=1440,height=960)
	NoCall_plot(NCR,chiptype_choice,chipname_list)
	dev.off()
    cat("     ~ The plotting of Call rate is finished.\n",file=logfile,append=T)
    cat("\n",file=logfile,append=T)
}
##------------------------------------------------------------------------------


if(AF.plot_choice==2)
{
	pathway_AFplot=paste(pathway_GR,"/AF plot",sep="")
	dir.create(pathway_AFplot,showWarnings=F)
	pathway_output_figure_intensity=paste(pathway_AFplot,"/Intensity",sep="")
	dir.create(pathway_output_figure_intensity,showWarnings=F)
	pathway_output_figure_genotype=paste(pathway_AFplot,"/Genotype",sep="")
	dir.create(pathway_output_figure_genotype,showWarnings=F)
	cat("    ***        **  **   ==========                AF plotting                   ==========  **  **     **\n")
	cat("   ~ Adjust allele frequency figure is preparing.\n",file=logfile,append=T)
	#pathway_figure_adjAF=pathway_output_population_figure_adjAF
	##----------------------------AF figure-----------------------------------------
	##   k-th chromosome ; i=1 plot for Allele frquency; i=2 plot for genotype
	##
	##------------------------------------------------------------------------------
	##------------------------------------------------------------------------------
	pathway_Geno_AFplot=paste(pathway_output_figure_genotype,"/AFplot",sep="")
	dir.create(pathway_Geno_AFplot,showWarnings=F)
	pathway_Intensity_AFplot=paste(pathway_output_figure_intensity,"/AFplot",sep="")
	dir.create(pathway_Intensity_AFplot,showWarnings=F)

	ind_name=list.files(pathway_adjAF)
	##------------------------------------------------------------------------------
	for(i in 1:length(ind_name))
	{
		#i=1
		cat("     ~ The adjust allele frequency figure of",paste("individual ",i,sep=""),"is preparing.\n",file=logfile,append=T)
		chr_file=list.files(paste(pathway_adjAF,"/",ind_name[i],sep=""),full.names=T)
		chr_file=mixedsort(chr_file)
		temp_ind_name=ind_name[i]
		if(AF.based_choice[1]==1) #intensity based
		{
			cat("       ~ The adjust allele frequency figure based on intensity of",paste("individual ",i,sep=""),"is preparing.\n",file=logfile,append=T)
			AF_plotting(pathway_Intensity_AFplot,chr_file,chromosome_choice,temp_ind_name,1,1)
			AF_plotting(pathway_Intensity_AFplot,chr_file,chromosome_choice,temp_ind_name,2,1)
			AF_plotting(pathway_Intensity_AFplot,chr_file,chromosome_choice,temp_ind_name,3,1)
			cat("       ~ The adjust allele frequency figure based on intensity of",paste("individual ",i,sep=""),"is finished.\n",file=logfile,append=T)
		}
		if(AF.based_choice[2]==1)#genotype based
		{
			cat("       ~ The adjust allele frequency figure based on genotype of",paste("individual ",i,sep=""),"is preparing.\n",file=logfile,append=T)
			AF_plotting(pathway_Geno_AFplot,chr_file,chromosome_choice,temp_ind_name,1,2)
			AF_plotting(pathway_Geno_AFplot,chr_file,chromosome_choice,temp_ind_name,2,2)
			AF_plotting(pathway_Geno_AFplot,chr_file,chromosome_choice,temp_ind_name,3,2)
			cat("       ~ The adjust allele frequency figure based on genotype of",paste("individual ",i,sep=""),"is finished.\n",file=logfile,append=T)
		}
		cat("     ~ The adjust allele frequency figure of",paste("individual ",i,sep=""),"is finished.\n",file=logfile,append=T)
	}
	##-----------remove the empty folder-------------------
	Geno_figure_pathway=list.files(pathway_Geno_AFplot,full.name=T)
	Intensity_figure_pathway=list.files(pathway_Intensity_AFplot,full.name=T)
	if(chiptype_choice=="1chip")
	{
		if(onechip_choice=="chip2"){
			remove_chip=c(1,3)
		}else{
			remove_chip=c(2,3)
		}
		unlink(Geno_figure_pathway[remove_chip],recursive=T)
		unlink(Intensity_figure_pathway[remove_chip],recursive=T)
	}
	cat("   ~ Adjust allele frequency figure is finished.\n",file=logfile,append=T)
	cat("\n",file=logfile,append=T)
}
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##==============================================================================




##==============================================================================

################################################################################
##                                Part 3                                      ##
##                            AF reference                                    ##
##     AF_ref_calculation_choice==1 meaning the AF ref from the datacase      ##
##     AF_ref_calculation_choice==2 meaning the AF ref is obatined the SAQC   ##
##============================================================================##
if(AF_ref_calculation_choice==2)
{
	cat("     ***       **  **   ==========             AF ref. calculation              ==========  **  **     ** \n")
	cat("   ~ The calcaulation of allele frequency reference is preparing.\n",file=logfile,append=T)
	pathway_temp_AF_ref=paste(pathway_temp,"/",population_name,"_AF ref",sep="")
	dir.create(pathway_temp_AF_ref,showWarnings=F)
	##---------------sub function for AF-ref  ------------------------------------
	pathway_NR_AFest=paste(pathway_NR,"/AF estimate",sep="")
	ind_name=list.files(pathway_NR_AFest)
	ind_name=mixedsort(ind_name);

	for(i in 1:nchr) #for each chromosome
	{
		#i=4
		u=chromosome_choice[i] #i: the index from GUI, u: the real choice of chr, e.g. if the user chosen chr 1, 10, 20, then i = 1~3, u = 1, 10, 20
		# gc(TRUE,verbose=F)

		if(u<23)
		{
			AF=NULL;GE=NULL;
			for(k in 1:length(ind_name)) #for each sample
			{
				pathway_chr=paste(pathway_NR_AFest,"/",ind_name[k],sep="")
				pathway_chr=list.files(pathway_chr,full.names=T) #extract paths of chrs
				pathway_chr=mixedsort(pathway_chr)
				pathway_chr=pathway_chr[i] #should be u? A: No, i is correct
				chr=read.delim(pathway_chr,header=T,sep="") #read ith chr of kth sample, #Probe_set	Phy_position	Chiptype	genotype	AF
				AF=cbind(AF,chr[,5]) #AF
				GE=cbind(GE,as.character(chr[,4]))	#genotype
			}
			chip=chr[,1:3] #Probe_set	Phy_position	Chiptype
			##-----------------output first row-----------------------------------------
			if(u<10) u = paste("0",u,sep="")
			name_row=cbind("Probe_set","Phy_Position","Chiptype","num_AA","mean_AA","std_AA","num_AB","mean_AB","std_AB","num_BB","mean_BB","std_BB")
			write.table(name_row,file=paste(pathway_temp_AF_ref,"/Chr_",u,".txt",sep=""),row.names=F,col.names=F,quote=F,append=F,sep="\t")
			ind=dim(GE)[2];nSNP=dim(GE)[1]
			r=2000;remainder=nSNP%%r;tt=ceiling(nSNP/r)
			##----------------------------------------------------------------------------      
			for(yAA in 1:tt) #process by 2000 SNPs
			{
				#yAA=1
				if(yAA<tt)
					R=r
				if(yAA==tt)
					R=remainder
				hAA=(yAA-1)*r+1;eAA=(yAA-1)*r+R;
				GE_temp=GE[hAA:eAA,];AF_temp=AF[hAA:eAA,];chip_temp=chip[hAA:eAA,]

				GE_AA=NULL; GE_AA[GE_temp!="AA"]=0;GE_AA[is.na(GE_AA)==T]=1;
				GE_AA=matrix(GE_AA,dim(GE_temp)[1],dim(GE_temp)[2]);ref_AA=AF_ref(AF_temp,GE_AA)

				GE_AB=NULL; GE_AB[GE_temp!="AB"]=0;GE_AB[is.na(GE_AB)==T]=1;
				GE_AB=matrix(GE_AB,dim(GE_temp)[1],dim(GE_temp)[2]);ref_AB=AF_ref(AF_temp,GE_AB)

				GE_BB=NULL; GE_BB[GE_temp!="BB"]=0;GE_BB[is.na(GE_BB)==T]=1;
				GE_BB=matrix(GE_BB,dim(GE_temp)[1],dim(GE_temp)[2]);ref_BB=AF_ref(AF_temp,GE_BB)

				ref=cbind(chip_temp,ref_AA,ref_AB,ref_BB)
				# ref[,1]=sprintf("%8.0f",ref[,1]);ref[,2]=sprintf("%8.0f",ref[,2])
				ref[,5]=sprintf("%1.5f",ref[,5]);ref[,6]=sprintf("%1.5f",ref[,6])
				ref[,8]=sprintf("%1.5f",ref[,8]);ref[,9]=sprintf("%1.5f",ref[,9])
				ref[,11]=sprintf("%1.5f",ref[,11]);ref[,12]=sprintf("%1.5f",ref[,12])
				write.table(ref,file=paste(pathway_temp_AF_ref,"/Chr_",u,".txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
			}
		}
		if(u==23)
		{
			AF_F=NULL;GE_F=NULL;AF_M=NULL;GE_M=NULL;
			for(k in 1:length(ind_name))
			{
				#k=1
				pathway_chr=paste(pathway_NR_AFest,"/",ind_name[k],sep="")
				pathway_chr=list.files(pathway_chr,full.names=T)
				pathway_chr=mixedsort(pathway_chr)
				pathway_chr=pathway_chr[i] #should be u?
				chr=read.delim(pathway_chr,header=T,sep="")
				temp_genotype=as.character(chr[,4])
				num_AB=length(temp_genotype[temp_genotype=="AB"])
				num_SNP=length(temp_genotype)
				rate_heter=num_AB/num_SNP			## !!!!! 	Change		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if(rate_heter>=0.1)  #female #add "="
				{
					AF_F=cbind(AF_F,chr[,5])
					GE_F=cbind(GE_F,temp_genotype)
				}
				if(rate_heter<0.1)  #male
				{
					AF_M=cbind(AF_M,chr[,5])
					GE_M=cbind(GE_M,temp_genotype)
				}	
			}
			chip=chr[,1:3]
			n_male=ncol(GE_M);if(is.null(n_male)==T){n_male=0}
			n_female=ncol(GE_F);if(is.null(n_female)==T){n_female=0}
			##=======================chromosome X=====================================
			##-----------------output first row-----------------------------------------
			name_row=cbind("Probe_set","Phy_Position","Chiptype","num_AA","mean_AA","std_AA","num_AB","mean_AB","std_AB","num_BB","mean_BB","std_BB")
			write.table(name_row,file=paste(pathway_temp_AF_ref,"/Chr_",u,".txt",sep=""),row.names=F,col.names=F,quote=F,append=F,sep="\t")
			##------intercept the Null article

			# if(n_female<=1) 
			# {
			# ind=dim(GE_M)[2];nSNP=dim(GE_M)[1]
			# nSNP=dim(GE_M)[1]
			# ref_AA=rep(NA,nSNP*3,nSNP,3)
			# ref_AA=matrix(rep(NA,nSNP*3),nSNP,3)
			# ref=cbind(chip,ref_AA,ref_AA,ref_AA)
			# ref[,1]=sprintf("%8.0f",ref[,1]);ref[,2]=sprintf("%8.0f",ref[,2])
			# ref[,5]=sprintf("%1.5f",ref[,5]);ref[,6]=sprintf("%1.5f",ref[,6])
			# ref[,8]=sprintf("%1.5f",ref[,8]);ref[,9]=sprintf("%1.5f",ref[,9])
			# ref[,11]=sprintf("%1.5f",ref[,11]);ref[,12]=sprintf("%1.5f",ref[,12])
			# write.table(ref,file=paste(pathway_temp_AF_ref,"/Chr_",u,".txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
			# }

			if(n_female>1)
			{
				ind=dim(GE_F)[2];nSNP=dim(GE_F)[1]
				r=2000;remainder=nSNP%%r;tt=ceiling(nSNP/r)
				##----------------------------------------------------------------------------
				for(yAA in 1:tt)
				{
					#yAA=1
					if(yAA<tt)
						R=r
					if(yAA==tt)
						R=remainder
					hAA=(yAA-1)*r+1;eAA=(yAA-1)*r+R;
					GE_temp=GE_F[hAA:eAA,];AF_temp=AF_F[hAA:eAA,];chip_temp=chip[hAA:eAA,]

					GE_AA=NULL; GE_AA[GE_temp!="AA"]=0;GE_AA[is.na(GE_AA)==T]=1;
					GE_AA=matrix(GE_AA,dim(GE_temp)[1],dim(GE_temp)[2]);ref_AA=AF_ref(AF_temp,GE_AA)

					GE_AB=NULL; GE_AB[GE_temp!="AB"]=0;GE_AB[is.na(GE_AB)==T]=1;
					GE_AB=matrix(GE_AB,dim(GE_temp)[1],dim(GE_temp)[2]);ref_AB=AF_ref(AF_temp,GE_AB)

					GE_BB=NULL; GE_BB[GE_temp!="BB"]=0;GE_BB[is.na(GE_BB)==T]=1;
					GE_BB=matrix(GE_BB,dim(GE_temp)[1],dim(GE_temp)[2]);ref_BB=AF_ref(AF_temp,GE_BB)

					ref=cbind(chip_temp,ref_AA,ref_AB,ref_BB)
					# ref[,1]=sprintf("%8.0f",ref[,1]);ref[,2]=sprintf("%8.0f",ref[,2])
					ref[,5]=sprintf("%1.5f",ref[,5]);ref[,6]=sprintf("%1.5f",ref[,6])
					ref[,8]=sprintf("%1.5f",ref[,8]);ref[,9]=sprintf("%1.5f",ref[,9])
					ref[,11]=sprintf("%1.5f",ref[,11]);ref[,12]=sprintf("%1.5f",ref[,12])
					write.table(ref,file=paste(pathway_temp_AF_ref,"/Chr_",u,".txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
				}
			}
			##------------------------------------------------------------------------

			##=======================chromosome Y=====================================
			##-----------------output first row-----------------------------------------
			name_row=cbind("Probe_set","Phy_Position","Chiptype","num_AA","mean_AA","std_AA","num_AB","mean_AB","std_AB","num_BB","mean_BB","std_BB")
			write.table(name_row,file=paste(pathway_temp_AF_ref,"/Chr_",as.numeric(u)+1,".txt",sep=""),row.names=F,col.names=F,quote=F,append=F,sep="\t")
			##------intercept the Null article
      
			# if(n_male<=1)
			# {
			# ind=dim(GE_F)[2];nSNP=dim(GE_F)[1]
			# nSNP=dim(GE_F)[1]
			# ref_AA=rep(NA,nSNP*3,nSNP,3)
			# ref_AA=matrix(rep(NA,nSNP*3),nSNP,3) 
			# ref=cbind(chip,ref_AA,ref_AA,ref_AA)
			# ref[,1]=sprintf("%8.0f",ref[,1]);ref[,2]=sprintf("%8.0f",ref[,2])
			# ref[,5]=sprintf("%1.5f",ref[,5]);ref[,6]=sprintf("%1.5f",ref[,6])
			# ref[,8]=sprintf("%1.5f",ref[,8]);ref[,9]=sprintf("%1.5f",ref[,9])
			# ref[,11]=sprintf("%1.5f",ref[,11]);ref[,12]=sprintf("%1.5f",ref[,12])
			# write.table(ref,file=paste(pathway_temp_AF_ref,"/Chr_",u,".txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
			# write.table(ref,file=paste(pathway_temp_AF_ref,"/Chr_",as.numeric(u)+1,".txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
			# }

			if(n_male>1)
			{
				ind=dim(GE_M)[2];nSNP=dim(GE_M)[1]
				r=2000;remainder=nSNP%%r;tt=ceiling(nSNP/r)
				##----------------------------------------------------------------------------
				for(yAA in 1:tt)
				{
					if(yAA<tt)
						R=r
					if(yAA==tt)
						R=remainder
					hAA=(yAA-1)*r+1;eAA=(yAA-1)*r+R;
					GE_temp=GE_M[hAA:eAA,];AF_temp=AF_M[hAA:eAA,];chip_temp=chip[hAA:eAA,]

					GE_AA=NULL; GE_AA[GE_temp!="AA"]=0;GE_AA[is.na(GE_AA)==T]=1;
					GE_AA=matrix(GE_AA,dim(GE_temp)[1],dim(GE_temp)[2]);ref_AA=AF_ref(AF_temp,GE_AA)

					GE_AB=NULL; GE_AB[GE_temp!="AB"]=0;GE_AB[is.na(GE_AB)==T]=1;
					GE_AB=matrix(GE_AB,dim(GE_temp)[1],dim(GE_temp)[2]);ref_AB=AF_ref(AF_temp,GE_AB)

					GE_BB=NULL; GE_BB[GE_temp!="BB"]=0;GE_BB[is.na(GE_BB)==T]=1;
					GE_BB=matrix(GE_BB,dim(GE_temp)[1],dim(GE_temp)[2]);ref_BB=AF_ref(AF_temp,GE_BB)

					ref=cbind(chip_temp,ref_AA,ref_AB,ref_BB)
					# ref[,1]=sprintf("%8.0f",ref[,1]);ref[,2]=sprintf("%8.0f",ref[,2])
					ref[,5]=sprintf("%1.5f",ref[,5]);ref[,6]=sprintf("%1.5f",ref[,6])
					ref[,8]=sprintf("%1.5f",ref[,8]);ref[,9]=sprintf("%1.5f",ref[,9])
					ref[,11]=sprintf("%1.5f",ref[,11]);ref[,12]=sprintf("%1.5f",ref[,12])
					write.table(ref,file=paste(pathway_temp_AF_ref,"/Chr_",as.numeric(u)+1,".txt",sep=""),row.names=F,col.names=F,quote=F,append=T,sep="\t")
				}
			}
		##------------------------------------------------------------------------
		} # end of if(u==23)
		# gc(reset=TRUE,verbose=F)
	}
	cat("   ~ The calcaulation of allele frequency reference is finished.\n",file=logfile,append=T)
	cat("\n",file=logfile,append=T)
	pathway_NR_AFrefest=paste(pathway_NR,"/AF reference estimate",sep="")
	dir.create(pathway_NR_AFrefest,showWarnings = F)

	file.copy(list.files(pathway_temp_AF_ref, full.names = T), pathway_NR_AFrefest, overwrite = T)
}

##============================================================================================================================================================
##============================================================================================================================================================
##============================================================================================================================================================
##============================================================================================================================================================

cat("=========================================================================================================\n",file=logfile,append=T)
cat("===============================                                          ================================\n",file=logfile,append=T)
cat("===============================      The process of QI calculation       --------------------------------\n",file=logfile,append=T)
cat("===============================                                          ================================\n",file=logfile,append=T)
cat("=========================================================================================================\n",file=logfile,append=T)
cat("\n",file=logfile,append=T)

cat("       ***    ********  ==========       The process of QI calculation          ==========  **  **     ** \n")
################################################################################
#######################            Part 2         ##############################
##################     Quality Index Calculation area    #######################
################################################################################

################################################################################
##                                Part 4                                      ##
##                            QI calculation                                  ##
##     QI_calculation_choice==1; the QI calculation doesn't process           ##
##     QI_calculation_choice==2; calculate the QI via the SAQC                ##
##============================================================================##
pathway_program=paste(main_pathway,"PROGRAM",sep="")
pathway_subfunction=paste(pathway_program,"/SAQC_subfunction.r",sep="")
source(pathway_subfunction)
##==============================================================================
##-------------------------  database input ------------------------------------
##                          AF ref selection area
# cat("break-6", sep = "\n")
##	
##	Only "collect" the filename of chr filenames that needs to combine later on  
##		Not including import data
##
if(AF_ref_calculation_choice==1)  ## AF_ref from database
{
	if(AFref_source_choice[1]==1){AFref_database=pathway_DB_AFref} #SAQC database
	if(AFref_source_choice[2]==1){AFref_database=AFref_input} #user-provided database
	chr_name<-list.files(path=AFref_database)
	chr_file<-list.files(path=AFref_database,full.names=T)
	chr_name=mixedsort(chr_name);
	chr_file=mixedsort(chr_file);
	temp_chr_file=chr_file[c(chromosome_choice)]
	temp_chr_name=chr_name[c(chromosome_choice)]
	if(chromosome_choice[nchr]!=23){chr_file=temp_chr_file;chr_name=temp_chr_name}
	## if user chose all chromosomes, then chr_file combines the chr_file of chr 1 to 24 (i.e. chr Y)
	if(chromosome_choice[nchr]==23){chr_file=c(temp_chr_file,chr_file[24]);chr_name=c(temp_chr_name,chr_name[24])}
}
# cat("break-7", sep = "\n")
if(AF_ref_calculation_choice==2)  ## AF_ref calculated from data
{
	AFref_database=pathway_temp_AF_ref ## The AF reference were calculated from the data and output to pathway_temp_AF_ref in the previous step
	chr_name<-list.files(path=AFref_database)
	chr_file<-list.files(path=AFref_database,full.names=T)
	chr_name=mixedsort(chr_name);
	chr_file=mixedsort(chr_file);
	if(chromosome_choice[nchr]==23)
	{
		chr_X=read.delim(chr_file[nchr],head=T,sep="\t", stringsAsFactors = F)
		chr_Y=read.delim(chr_file[nchr+1],head=T,sep="\t", stringsAsFactors = F)
		# AFref_database=pathway_DB_AFref
		# temp_chr_name<-list.files(path=AFref_database)
		# temp_chr_file<-list.files(path=AFref_database,full.names=T)
		# temp_chr_name=mixedsort(temp_chr_name);temp_chr_file=mixedsort(temp_chr_file);
		if(nrow(chr_X)==0) ##  Happens when only one female in the data (Need to pay attention 'cause the column name is exported anyway in the begining)
		{					##	When this happens, use the AF reference from the database
			chr_file[nchr]=mixedsort(dir(pathway_DB_AFref, full.names = T))[23];
			chr_name[nchr]=mixedsort(dir(pathway_DB_AFref))[23];
		}
		if(nrow(chr_Y)==0) ## Happens when only one male in the data (Need to pay attention 'cause the column name is exported anyway in the begining)
		{					##	When this happens, use the AF reference from the database
			chr_file[nchr+1]=mixedsort(dir(pathway_DB_AFref, full.names = T))[24];
			chr_name[nchr+1]=mixedsort(dir(pathway_DB_AFref))[24];
		}
	}
	#chr_file=mixedsort(chr_file)
	#chr_name=mixedsort(chr_name)
}
##==============================================================================
# cat("break-8", sep = "\n")
##==============================================================================
##------------------------  Quality Index of Chromosome ------------------------
##==============================================================================
##------------------------------------------------------------------------------
if(QI_calculation_choice == 2){
	pathway_NR_QI=paste(pathway_NR,"/QI estimate",sep="")
	dir.create(pathway_NR_QI,showWarnings = F)
	pathway_NR_indQI=paste(pathway_NR_QI,"/Ind_QI",sep="")
	dir.create(pathway_NR_indQI,showWarnings = F)
}
# cat("break-9", sep = "\n")
if(QI_calculation_choice==2 & QI_est.based_choice[1]==1)    # chip
{
	cat("   ~ The calcaulation of QI in SNP base is preparing.\n",file=logfile,append=T)
	ind_name<-list.files(path=paste(pathway_NR,"/AF estimate",sep=""))
	ind_name=mixedsort(ind_name);

	##		Import AF data for the following QI calculation
	for(i in 1:length(ind_name))
	{
    #i=1
		cat(paste("     ~  The QI calcaulation of individual ",i," is preparing.\n",sep=""),file=logfile,append=T)
		
		indQI_folder=paste(pathway_NR_indQI,"/",ind_name[i],sep="")
		dir.create(indQI_folder,showWarnings=F)

		pathway_indAF=paste(pathway_NR,"/AF estimate/",ind_name[i],sep="")
		ind_chr_name=list.files(path=pathway_indAF)
		ind_chr_file=list.files(path=pathway_indAF,full.names=T)
		ind_chr_name=mixedsort(ind_chr_name);
		ind_chr_file=mixedsort(ind_chr_file);
		nchr=length(ind_chr_name)
		for(u in 1:nchr) ##---we only calculate the chromosome frome 1 to 23
		{
			#u=2
			k=chromosome_choice[u]
			ind_chr=read.delim(ind_chr_file[u],fill=T,header=T,sep="", stringsAsFactors = F)
			##		chr: AF reference 
			##		ind_chr: individual AF data
			if(k<23)
			{
				chr=read.delim(chr_file[u],fill=T,header=T,sep="", stringsAsFactors = F)
			}
			if(k==23)
			{
				##
				##		Detect gender
				##		
				# AB_num=table(ind_chr$Genotype)[2]
				AB_num=sum(as.character(na.omit(ind_chr$Genotype)) == "AB")
				AB_rate=AB_num/(length(ind_chr[,4]))
				# if(AF_ref_calculation_choice==1){x=23}
				# if(AF_ref_calculation_choice==2){x=u}

				x = u #masaki adds
				if(AB_rate>=0.1)  ## female
				{
					chr=read.delim(chr_file[x],fill=T,header=T,sep="", stringsAsFactors = F)
					gender="female";
				}
				if(AB_rate<0.1)  ## male
				{
					chr=read.delim(chr_file[x+1],fill=T,header=T,sep="", stringsAsFactors = F)
					gender="male";
				}
			}
			#chr=sort.append(chr,2)    ##<<--------------------
			#NoCall=ind_chr[ind_chr$Genotype=="NoCall",];nNoCall=length(NoCall[,1])
			##	ƒÜ	Filter out the SNPs with NoCall!!!
			ind_chr=ind_chr[ind_chr$Genotype!="NoCall",];
			colnames(chr)[1]="Probe_set";
			colnames(ind_chr)[1]="Probe_set";
			##		Only keep the SNPs both existed in AF reference and indiv. AF data	  
			ind_chr=merge(ind_chr,chr,by="Probe_set",all=F)
			chr=ind_chr[,c(1,4,3,2,4,5,9,10,12,13,15,16)]
			#1. Probe_set, 2. Genotype, 3. Chiptype, 4. Phy_position, 5. Genotype.1, 6. adjustAF, 7. mean_AA, 8. std_AA, 9. mean_AB, 10. std_AB, 11. mean_BB, 12.std_BB
			#colnames(chr)[3]="Chiptype"
			geno=as.character(chr[,2]);
			geno[geno=="AA"]=3;geno[geno=="AB"]=2;geno[geno=="BB"]=1;geno[geno=="NoCall"]=0
			## 2. Genotype, 5. Genotype.1    <=====same column?
			chr[,2]=as.numeric(geno);	  chr[,5]=as.numeric(geno)
			##	3. Chiptype
			chip=as.character(chr[,3]);
			chip[chip=="chip1"]=1;chip[chip=="chip2"]=2;
			chr[,3]=as.numeric(chip);
			chr=as.matrix(chr)
			# Nrow=nrow(chr);      Ncol=ncol(chr);
			# chr=matrix(as.numeric(chr),Nrow,Ncol)
			##--------------------QI calculation--------------------------------------
			##	ind_chr_QI1: 1. Probe_set, 2. Genotype, 3. Chiptype, 4. Phy_position, 5. Genotype.1, 6. adjustAF, 7:QI1
			ind_chr_QI1=Quality_Index_1(chr) # genotype-specific 
			##	ind_chr_QI2: 1. Probe_set, 2. Genotype, 3. Chiptype, 4. Phy_position, 5. Genotype.1, 6. adjustAF, 7:QI2
			ind_chr_QI2=Quality_Index_2(chr) # not considering genotype
			##------------------------------------------------------------------------
			QI1=ind_chr_QI1[,7];QI1=as.matrix(QI1);QI1=as.numeric(QI1)
			QI2=ind_chr_QI2[,7];QI2=as.matrix(QI2);QI2=as.numeric(QI2)
			rm(chr)
			ind_chr=cbind(ind_chr_QI1[,c(1:6)],QI1,QI2)
			ind_chr[,7]=round(as.numeric(ind_chr[,7])*100000)/100000;
			ind_chr[,8]=round(as.numeric(ind_chr[,8])*100000)/100000;
			ind_chr=ind_chr[,-2]	  ##		Drop the duplicated column: genotype
			# ind_chr=data.frame(ind_chr)
			colnames(ind_chr)[1]="Probe_set"
			colnames(ind_chr)[2]="Chiptype"
			colnames(ind_chr)[3]="Phy_posi"
			colnames(ind_chr)[4]="genotype"
			colnames(ind_chr)[5]="AF"
			colnames(ind_chr)[6]="QI1"
			colnames(ind_chr)[7]="QI2"
			if(k<10){kchr=paste("0",k,sep="")}
			if(k>9){kchr=paste(k,sep="")}
			QI_table=paste(indQI_folder,"/Chr_",kchr,".txt",sep="")
			write.table(ind_chr,QI_table,row.names=F,col.names=T,quote=F,append=F,sep="\t")
		}
		cat(paste("     ~  The QI calcaulation of individual ",i," is finished.\n",sep=""),file=logfile,append=T)
	}
	cat("   ~ The calcaulation of QI in SNP base is finished.\n",file=logfile,append=T)
	cat("\n",file=logfile,append=T)
}
# cat("break-10", sep = "\n")
##------------------------------------------------------------------------------
##                    QI calculation by chip or chromosome
##------------------------------------------------------------------------------
if(QI_calculation_choice==2)    #chromosome
{
	pathway_QI_table=paste(pathway_NR_QI,"/QI table",sep="")
    cat("   ~ The calculation of QI in chip base is preparing.\n",file=logfile,append=T)
    qtrim=0.95; # For WinsMean: winsorized mean
	qmedian=0.5	# For median
    dir.create(pathway_QI_table,showWarnings=F)
    indQI_file=mixedsort(list.files(pathway_NR_indQI,full.names=T))
    nind=length(indQI_file)
    if(chiptype_choice=="2chip"){first_colname=c("Ind","Merge","Chip1","Chip2",paste("Chr",chromosome_choice,sep=""))}
    if(chiptype_choice=="1chip"){first_colname=c("Ind","Chip1",paste("Chr",chromosome_choice,sep=""))}
    
    first_colname=matrix(first_colname,1,length(first_colname))
    
    mQI1_table=paste(pathway_QI_table,"/WinsMean QI1.txt",sep="")
    mQI2_table=paste(pathway_QI_table,"/WinsMean QI2.txt",sep="")
    MQI1_table=paste(pathway_QI_table,"/Median QI1.txt",sep="")
    MQI2_table=paste(pathway_QI_table,"/Median QI2.txt",sep="")
    QI_table=c(mQI1_table,mQI2_table,MQI1_table,MQI2_table)
    write.table(first_colname,mQI1_table,col.names=F,row.names=F,quote=F,append=F,sep="\t")
    write.table(first_colname,mQI2_table,col.names=F,row.names=F,quote=F,append=F,sep="\t")
    write.table(first_colname,MQI1_table,col.names=F,row.names=F,quote=F,append=F,sep="\t")
    write.table(first_colname,MQI2_table,col.names=F,row.names=F,quote=F,append=F,sep="\t")
        
    for(k in 1:nind)
    {
		#k=1
		cat(paste("     ~  The QI calcaulation of individual ",k," is preparing.\n",sep=""),file=logfile,append=T)
		chr_file=list.files(indQI_file[k],full.names=T)
		nchr=length(chr_file)
		whole_QI=NULL
		QI_chr=NULL
		for(i in 1:nchr)
		{
			ichr=chromosome_choice[i]
			chr=read.delim(chr_file[i],head=T,sep="", stringsAsFactors = F)
			# geno=chr[,4];geno[geno==0]=1;geno[geno!=0]=0
			
			QI=cbind(chr[,c(6,7)],NA)
			QI=QIstat(QI,qtrim,qmedian)
			QI=round(QI[-c(3,6)]*10000)/10000
			if(ichr<10){Ichr=paste("chr0",ichr,sep="")}
			if(ichr>9){Ichr=paste("chr",ichr,sep="")}
			QI=c(Ichr,QI)
			QI_chr=rbind(QI_chr,QI)
			whole_QI=rbind(whole_QI,chr[,c(2,6,7)])
		}
      
		if(chiptype_choice=="2chip")
		{
			## Chip 1
			chip1_QI=whole_QI[whole_QI[,1]==1,]
			chip1_QI=cbind(chip1_QI[,-1],NA)	#	cbind a "NA" 'cause we don't use QI3 in this case, so we fill the third column with NA
			chip1_QI=QIstat(chip1_QI,qtrim,qmedian)
			chip1_QI=round(chip1_QI[-c(3,6)]*10000)/10000
			chip1_QI=c("Chip1",chip1_QI)

			## Chip 2
			chip2_QI=whole_QI[whole_QI[,1]==2,]
			chip2_QI=cbind(chip2_QI[,-1],NA)	#	cbind a "NA" 'cause we don't use QI3 in this case, so we fill the third column with NA
			chip2_QI=QIstat(chip2_QI,qtrim,qmedian)
			chip2_QI=round(chip2_QI[-c(3,6)]*10000)/10000
			chip2_QI=c("Chip2",chip2_QI)

			## Both chips
			whole_QI=cbind(whole_QI[,-1],NA)
			whole_QI=QIstat(whole_QI,qtrim,qmedian)
			whole_QI=round(whole_QI[-c(3,6)]*10000)/10000
			whole_QI=c("Merge",whole_QI)
			whole_QI=rbind(whole_QI,chip1_QI,chip2_QI,QI_chr) # IndID  Merge	 Nsp  	 Sty  	Chr1  	Chr2  	Chr3  	Chr4  	Chr5  	Chr6  	Chr7  	Chr8  	Chr9  	Chr10 	Chr11 	Chr12 	Chr13 	Chr14 	Chr15 	Chr16 	Chr17 	Chr18 	Chr19 	Chr20 	Chr21 	Chr22 	Chr23 

		}
		if(chiptype_choice=="1chip")
		{
			if(onechip_choice=="chip2") Ichip=2 else Ichip=1
			# if(onechip_choice=="chip2"){Ichip=2}
			## Chip 1 or 2 
			chip1_QI=whole_QI[whole_QI[,1]==Ichip,]
			chip1_QI=cbind(chip1_QI[,-1],NA)
			chip1_QI=QIstat(chip1_QI,qtrim,qmedian)
			chip1_QI=round(chip1_QI[-c(3,6)]*10000)/10000
			chip1_QI=c("Chip1",chip1_QI)
		  
			## Both chips
			whole_QI=chip1_QI
			whole_QI=as.matrix(whole_QI)
			whole_QI=matrix(whole_QI,ncol(whole_QI),nrow(whole_QI))
			#whole_QI=t(whole_QI)
			whole_QI=rbind(whole_QI,QI_chr) # Copy chip1(or chip 2) to the column "merge" 
			# e.g. IndID  Chip1	 Chip1  	 Chip1  	Chr1  	Chr2  	Chr3  	Chr4  	Chr5  	Chr6  	Chr7  	Chr8  	Chr9  	Chr10 	Chr11 	Chr12 	Chr13 	Chr14 	Chr15 	Chr16 	Chr17 	Chr18 	Chr19 	Chr20 	Chr21 	Chr22 	Chr23 
		}

		whole_QI=as.matrix(whole_QI)
		whole_QI=matrix(whole_QI,nrow(whole_QI),ncol(whole_QI))
		whole_QI=data.frame(whole_QI)

		colnames(whole_QI)[1]="Chiptype"
		colnames(whole_QI)[2]="WinsMeanQI1"
		colnames(whole_QI)[3]="WinsMeanQI2"
		colnames(whole_QI)[4]="MedianQI1"
		colnames(whole_QI)[5]="MedianQI2"
		QIname=c("WinsMeanQI1","WinsMeanQI2","MedianQI1","MedianQI2")
		for(j in 1:4)
		{
			cat(paste("       ~  The calcaulation of ",QIname[j]," is preparing.\n",sep=""),file=logfile,append=T)
			ind_QI=as.matrix(whole_QI[,j+1])
			ind_QI=as.numeric(ind_QI)
			ind_QI=c(ind_name[k],ind_QI)
			ind_QI=matrix(ind_QI,1,length(ind_QI))
			write.table(ind_QI,QI_table[j],col.names=F,row.names=F,quote=F,append=T,sep="\t")	# QI_table: outpath
			cat(paste("       ~  The calcaulation of ",QIname[j]," is finished.\n",sep=""),file=logfile,append=T)
		}
		cat(paste("     ~  The QI calcaulation of individual ",k," is finished.\n",sep=""),file=logfile,append=T)
    }
    cat("   ~ The calculation of QI in chip base is finished.\n",file=logfile,append=T)
    cat("\n",file=logfile,append=T)
}
############################################################################################### YTLin
# pathway_DBUB=paste(pathway_DB,"/Upper_bound",sep="")
# pathway_DBMS = paste(pathway_DB, "/Mean_Std", sep = "")
# ind_winQ1 = read.table(QI_table[1], sep = "\t", header = T)
# MS_winQ1 = ind_winQ1[,c(1,3:25,2,2,2)]
# colnames(MS_winQ1)[c(2:10, 25:27)] = c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chip1", "Chip2", "Merge")
# MS_winQ1 = data.frame(Chip = colnames(MS_winQ1)[2:27], Mean = round(colMeans(MS_winQ1[,2:27]),4), Std = round(sapply(MS_winQ1[2:27], sd),4))
# write.table(MS_winQ1, file = paste(pathway_DBMS, "/WinsMean_QI1.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)

# ind_winQ2 = read.table(QI_table[2], sep = "\t", header = T)
# MS_winQ2 = ind_winQ2[,c(1,3:25,2,2,2)]
# colnames(MS_winQ2)[c(2:10, 25:27)] = c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chip1", "Chip2", "Merge")
# MS_winQ2 = data.frame(Chip = colnames(MS_winQ2)[2:27], Mean = round(colMeans(MS_winQ2[,2:27]),4), Std = round(sapply(MS_winQ2[2:27], sd),4))
# write.table(MS_winQ2, file = paste(pathway_DBMS, "/WinsMean_QI2.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)

# ind_medQ1 = read.table(QI_table[3], sep = "\t", header = T)
# MS_medQ1 = ind_medQ1[,c(1,3:25,2,2,2)]
# colnames(MS_medQ1)[c(2:10, 25:27)] = c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chip1", "Chip2", "Merge")
# MS_medQ1 = data.frame(Chip = colnames(MS_medQ1)[2:27], Mean = round(colMeans(MS_medQ1[,2:27]),4), Std = round(sapply(MS_medQ1[2:27], sd),4))
# write.table(MS_medQ1, file = paste(pathway_DBMS, "/Median_QI1.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)

# ind_medQ2 = read.table(QI_table[4], sep = "\t", header = T)
# MS_medQ2 = ind_medQ2[,c(1,3:25,2,2,2)]
# colnames(MS_medQ2)[c(2:10, 25:27)] = c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chip1", "Chip2", "Merge")
# MS_medQ2 = data.frame(Chip = colnames(MS_medQ2)[2:27], Mean = round(colMeans(MS_medQ2[,2:27]),4), Std = round(sapply(MS_medQ2[2:27], sd),4))
# write.table(MS_medQ2, file = paste(pathway_DBMS, "/Median_QI2.txt", sep = ""), row.names = F, col.names = T, sep = "\t", quote = F)

# winQ1 = t(apply(ind_winQ1[,2:ncol(ind_winQ1)], 2, function(a) quantile(x = a, probs = c(0.95, 0.975, 0.99))))
# winQ2 = t(apply(ind_winQ2[,2:ncol(ind_winQ2)], 2, function(a) quantile(x = a, probs = c(0.95, 0.975, 0.99))))
# medQ1 = t(apply(ind_medQ1[,2:ncol(ind_medQ1)], 2, function(a) quantile(x = a, probs = c(0.95, 0.975, 0.99))))
# medQ2 = t(apply(ind_medQ2[,2:ncol(ind_medQ2)], 2, function(a) quantile(x = a, probs = c(0.95, 0.975, 0.99))))
# Chip = c("Mean QI1", "Mean QI2", "Mean QI1", "Median QI1", "Median QI2", "Median QI1")
# upper = c("UB_95", "UB_975", "UB_99")
# temp = cbind(winQ1, winQ2, winQ1, medQ1, medQ2, medQ1)
# for(chr in 1:23){
	# chrname = ifelse(chr<10, paste("0", chr, sep = ""), chr)
	# Chr_upper = data.frame(Chip = Chip, t(matrix(temp[chr+1,], 3)))
	# colnames(Chr_upper)[2:4] = upper
	# write.table(Chr_upper, file = paste(pathway_DBUB, "/Chr", chrname, ".txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
# }
# chip_upper = data.frame(Chip = Chip, t(matrix(temp[1,], 3)))
# colnames(chip_upper)[2:4] = upper
# write.table(chip_upper, file = paste(pathway_DBUB, "/Chip1.txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
# write.table(chip_upper, file = paste(pathway_DBUB, "/Chip2.txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
# write.table(chip_upper, file = paste(pathway_DBUB, "/Merge.txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")

# cat("break-11", sep = "\n")
################################################################################
#===============================================================================
#===============================================================================
################################################################################

################################################################################
###                  Graphical Results area                                  ###
################################################################################

cat("        **     **  **   ==========       The process of QI graphical output     ==========  **  ***    **\n")
pathway_program=paste(main_pathway,"PROGRAM",sep="")
pathway_subfunction=paste(pathway_program,"/SAQC_subfunction.r",sep="")
source(pathway_subfunction)
##-----------------------------------------------------------------------------

####--------------------------Graphical Part 1------------------------------####
####--------------------------Heat Map--------------------------------------####
if(QI_calculation_choice == 2){
	pathway_output_GR_QIplot=paste(pathway_GR,"/QI plot",sep="")
	dir.create(pathway_output_GR_QIplot,showWarnings=F)

	pathway_output_GR_HeatMap=paste(pathway_output_GR_QIplot,"/QI-HeatMap plot",sep="")
	dir.create(pathway_output_GR_HeatMap,showWarnings=F)

	pathway_output_CRplot=paste(pathway_output_GR_QIplot,"/QI-Polygon plot",sep="")
	dir.create(pathway_output_CRplot,showWarnings=F)
}

if(HM_plot_choice==2)      #whether plotting or not!!
{
	cat("   ~ The HeatMap plotting of QI in chip base is preparing.\n",file=logfile,append=T)
	QI_file=list.files(pathway_QI_table,full.names=T)	# 	pathway_QI_table=paste(pathway_NR_QI,"/QI table",sep="")
	QI_name=list.files(pathway_QI_table,full.names=F)
	QI_name=strsplit(QI_name,".txt")
	QI_name=unlist(QI_name)
	nstat=length(QI_file)
    ##--------------------------------------------------------------------------
    if(col.scale_choice[1]==1)      # scale from database
    {
		pathway_output_GR_HeatMap_Ref=paste(pathway_output_GR_HeatMap,"/Reference scale",sep="")
		dir.create(pathway_output_GR_HeatMap_Ref,showWarnings=F)
		DB_QIinfo=paste(pathway_DB,"/Mean_Std",sep="")
		DB_QIinfo_file=list.files(DB_QIinfo,full.names=T)
		for(k in 1:nstat)
		{
			#k=1
			QI=read.delim(QI_file[k],head=T,sep="\t")	# QI of sample data
			DB_QI=read.delim(DB_QIinfo_file[k],head=T,sep="\t")
			DB_QI=DB_QI[c(24:26),] # Mean and std of 24. Chip1, 25. Chip2 and 26. Merge
								 # Currently, SAQC doesn't support drawing plot by chromosome
			if(chiptype_choice=="2chip"){QI=QI[,c(3,4,2)];	#chip1 chip2 merge, row: an individual
			DB_QI=as.numeric(DB_QI[3,-1])} 	# Only keep result of "merge", col 1: Chip, 2. Mean, 3. Std 
			if(chiptype_choice=="1chip")
			{
				QI=QI[,c(2,2,2)] # In the previous step, column 2 is the one-chip QI value
				if(onechip_choice=="chip2") DB_QI=as.numeric(DB_QI[2,-1]) else DB_QI=as.numeric(DB_QI[1,-1])
			}
        
			QI=as.matrix(QI);nind=nrow(QI);nchip=ncol(QI)
			pathway_QI_figure=paste(pathway_output_GR_HeatMap_Ref,"/",QI_name[k],".png",sep="")
			png(pathway_QI_figure,width=1440,height=960)
			HeatMap_Plot(QI,QI_name[k],chiptype_choice,chipname_list,DB_QI)
			dev.off()
		}
	}
    ##--------------------------------------------------------------------------
    if(col.scale_choice[2]==1)    # scale from user provided
    {
		pathway_output_GR_HeatMap_User=paste(pathway_output_GR_HeatMap,"/Provided scale",sep="")
		dir.create(pathway_output_GR_HeatMap_User,showWarnings=F)
		for(k in 1:nstat)
		{
			#k=1
			QI=read.delim(QI_file[k],head=T,sep="\t")
			if(chiptype_choice=="2chip")
			{
				QI=QI[,c(3,4,2)];QI=as.matrix(QI);nind=nrow(QI);nchip=ncol(QI)
				QIinfo=c(mean(as.numeric(QI)),sd(as.numeric(QI)))
				#QIinfo=c(mean(as.numeric(QI[,3])),sd(as.numeric(QI[,3])))

			}
			if(chiptype_choice=="1chip")
			{
				QI=QI[,c(2,2,2)]
				QI=as.matrix(QI);nind=nrow(QI);nchip=ncol(QI)
				QIinfo=c(mean(as.numeric(QI[,1])),sd(as.numeric(QI[,1])))
			}
			pathway_QI_figure=paste(pathway_output_GR_HeatMap_User,"/",QI_name[k],".png",sep="")
			png(pathway_QI_figure,width=1440,height=960)
			HeatMap_Plot(QI,QI_name[k],chiptype_choice,chipname_list,QIinfo)
			dev.off()
		}
    }
    cat("   ~ The HeatMap of QI in chip base is finished.\n",file=logfile,append=T)
    cat("\n",file=logfile,append=T)
}

##------------------------------------------------------------------------------
####---------------------------Polygon Plot--------------------------------####

# if(Iden_poor_chipchr_choice==2)
if(PG_plot_choice==2)
{
    cat("   ~ The Polygon plotting of QI in chip base is preparing.\n",file=logfile,append=T)
    ##--------------------------------------------------------------------------
    QI_file=list.files(pathway_QI_table,full.names=T)	  #	pathway_NR_QI=paste(pathway_NR,"/QI estimate",sep="")
														  # pathway_QI_table=paste(pathway_NR_QI,"/QI table",sep="")
    QI_name=list.files(pathway_QI_table,full.names=F)
    QI_name=strsplit(QI_name,".txt")
    QI_name=unlist(QI_name)
    nstat=length(QI_file)
    ##--------------------------------------------------------------------------
    
    pathway_DBUB=paste(pathway_DB,"/Upper_bound",sep="")

	##	Column 1. Chip, 2. UB_95, 3. UB_975, 4. UB_99	
    UB_Merge=read.delim(paste(pathway_DBUB,"/Merge.txt",sep=""),head=T,sep="\t")
    UB_Chip1=read.delim(paste(pathway_DBUB,"/Chip1.txt",sep=""),head=T,sep="\t")
    UB_Chip2=read.delim(paste(pathway_DBUB,"/Chip2.txt",sep=""),head=T,sep="\t")
    
	# if(QI_quantile_choice==0.95){UB_QI=cbind(UB_Chip1[,1],UB_Chip2[,1],UB_Merge[,1])}
    # if(QI_quantile_choice==0.975){UB_QI=cbind(UB_Chip1[,2],UB_Chip2[,2],UB_Merge[,2])}
    # if(QI_quantile_choice==0.99){UB_QI=cbind(UB_Chip1[,3],UB_Chip2[,3],UB_Merge[,3])}
	
	if(QI_quantile_choice==0.95){UB_QI=cbind(UB_Chip1[,1 + 1],UB_Chip2[,1 + 1],UB_Merge[,1 + 1])}
    if(QI_quantile_choice==0.975){UB_QI=cbind(UB_Chip1[,2 + 1],UB_Chip2[,2 + 1],UB_Merge[,2 + 1])}
    if(QI_quantile_choice==0.99){UB_QI=cbind(UB_Chip1[,3 + 1],UB_Chip2[,3 + 1],UB_Merge[,3 + 1])}
	##	Original idea has 3 QI values, but QI3 is removed from analyses
	##	Hsin-Chi copied the values of QI1 to replace empty QI3
	UB_QI = UB_QI[-c(3,6), ];	##	Removed the row of QI3
	UB_QI = UB_QI[c(3,4,1,2), ] ##  row:  MedianQI1, MedianQI2, MeanQI1, MeanQI2;  col: QI_chip1  QI_chip2  QI_merge

    for(k in 1:nstat)	# MedianQI1, MedianQI2, MeanQI1, MeanQI2
    {
		#k=1
		QI=read.delim(QI_file[k],head=T,sep="\t")
		if(chiptype_choice=="2chip")
		{
			QI=QI[,c(1,3,4,2)] # IndID, Chip1, Chip2, Merge
			temp_QI=QI[,-1] 	 # Chip1, Chip2, Merge
			temp_UBQI=UB_QI[k,]
		}
		if(chiptype_choice=="1chip")
		{
			QI=QI[,c(1,2,2,2)] # IndID, Chip1, Chip1, Chip1
			temp_QI=QI[,-1]
			if(onechip_choice=="chip2") ichip=2 else ichip=1
			temp_UBQI=rep(UB_QI[k,ichip],3)
		}
		pathway_QI_idenplot=paste(pathway_output_CRplot,"/",QI_name[k],".png",sep="")
		png(pathway_QI_idenplot,width=1440,height=960)
		Polygon_Plot(temp_QI,temp_UBQI,QI_name[k],chiptype_choice,chipname_list)
		dev.off()
    }
    cat("   ~ The Polygon plotting of QI in chip base is finished.\n",file=logfile,append=T)
    cat("\n",file=logfile,append=T)
}


cat("    ******     **  **   ==========       The process of QI numerical output     ==========   **** **    ****\n")
################################################################################
###                  Numerical Results area                                  ###
################################################################################

##------------------numerical result -------------------------
##  calculate the QI_poor_chip.txt Table_QI_summary
##  upper_rate=0.95,0.975,0.99  there are three choices
##------------------------------------------------------------
pathway_program=paste(main_pathway,"PROGRAM",sep="")
pathway_subfunction=paste(pathway_program,"/SAQC_subfunction.r",sep="")
source(pathway_subfunction)
##==============================================================================
n=length(ind_name)

if(Iden_poor_chipchr_choice==2)
{
	cat("   ~ The summary table of QI is preparing.\n",file=logfile,append=T)
	##--------------------------------------------------------------------------
    QI_file=list.files(pathway_QI_table,full.names=T)
    QI_name=list.files(pathway_QI_table,full.names=F)
    QI_name=strsplit(QI_name,".txt")
    QI_name=unlist(QI_name)
    nstat=length(QI_file)
    ##--------------------------------------------------------------------------
    pathway_DBUB=paste(pathway_DB,"/Upper_bound",sep="")
    UB_Merge=read.delim(paste(pathway_DBUB,"/Merge.txt",sep=""),head=T,sep="\t")
    UB_Chip1=read.delim(paste(pathway_DBUB,"/Chip1.txt",sep=""),head=T,sep="\t")
    UB_Chip2=read.delim(paste(pathway_DBUB,"/Chip2.txt",sep=""),head=T,sep="\t")
    
	# if(QI_quantile_choice==0.95){UB_QI=cbind(UB_Chip1[,1],UB_Chip2[,1],UB_Merge[,1])}
    # if(QI_quantile_choice==0.975){UB_QI=cbind(UB_Chip1[,2],UB_Chip2[,2],UB_Merge[,2])}
    # if(QI_quantile_choice==0.99){UB_QI=cbind(UB_Chip1[,3],UB_Chip2[,3],UB_Merge[,3])}
	
	if(QI_quantile_choice==0.95){UB_QI=cbind(UB_Chip1[,1 + 1],UB_Chip2[,1 + 1],UB_Merge[,1 + 1])}
    if(QI_quantile_choice==0.975){UB_QI=cbind(UB_Chip1[,2 + 1],UB_Chip2[,2 + 1],UB_Merge[,2 + 1])}
    if(QI_quantile_choice==0.99){UB_QI=cbind(UB_Chip1[,3 + 1],UB_Chip2[,3 + 1],UB_Merge[,3 + 1])}
    UB_QI=UB_QI[-c(3,6),];UB_QI=UB_QI[c(3,4,1,2),] ##  row:  MedianQI1, MedianQI2,MeanQI1, MeanQI2;  col: QI_chip1  QI_chip2  QI_merge

	if(chiptype_choice=="1chip")
    {
        if(onechip_choice=="chip2")ichip=2 else ichip=1
        temp_UBQI=rep(UB_QI[,ichip],3)
        UB_QI=matrix(temp_UBQI,nrow(UB_QI),ncol(UB_QI))
    }

	chip_MQI1=read.delim(QI_file[1],header=T,sep="") #Median QI1
	chip_MQI2=read.delim(QI_file[2],header=T,sep="") #Median QI2
	chip_WQI1=read.delim(QI_file[3],header=T,sep="") #WinsMean QI1
	chip_WQI2=read.delim(QI_file[4],header=T,sep="") #WinsMean QI2

	if(chiptype_choice=="2chip"){specific_chip=c(1:4)}
	if(chiptype_choice=="1chip"){specific_chip=c(1,2,2,2)}

	chip_MQI1=as.matrix(chip_MQI1);chip_MQI1=chip_MQI1[,specific_chip];
	chip_MQI2=as.matrix(chip_MQI2);chip_MQI2=chip_MQI2[,specific_chip];
	chip_WQI1=as.matrix(chip_WQI1);chip_WQI1=chip_WQI1[,specific_chip];
	chip_WQI2=as.matrix(chip_WQI2);chip_WQI2=chip_WQI2[,specific_chip];

	chip_QI=list(chip_WQI1,chip_WQI2,chip_MQI1,chip_MQI2)
	##==============================================================================
	##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	##                                Table_QI_summary.txt
	#################### Chapter 2 Correlation between whole QI index ##############
	corrM1M2=QI_corr(chip_MQI1,chip_MQI2)
	corrM1W2=QI_corr(chip_MQI1,chip_WQI2)
	corrW1M2=QI_corr(chip_WQI1,chip_MQI2)
	corrW1W2=QI_corr(chip_WQI1,chip_WQI2)
	corrM1W1=QI_corr(chip_MQI1,chip_WQI1)
	corrM2W2=QI_corr(chip_MQI2,chip_WQI2)
	corr_whole=cbind(corrM1M2,corrM1W2,corrW1M2,corrW1W2,corrM1W1,corrM2W2)
	corr_whole=t(corr_whole)
	rownames(corr_whole)=c(" Pearson_corr(MQI1,MQI2)"," Pearson_corr(MQI1,WQI2)"," Pearson_corr(WQI1,MQI2)"," Pearson_corr(WQI1,WQI2)"," Pearson_corr(MQI1,WQI1)"," Pearson_corr(MQI2,WQI2)")
	##-------------------spearman correlation---------------------------------------
	##==========================================================================
	##--------------------------------------------------------------------------
	ScorrM1M2=spearman_QI_corr(chip_MQI1,chip_MQI2)
	ScorrM1W2=spearman_QI_corr(chip_MQI1,chip_WQI2)
	ScorrW1M2=spearman_QI_corr(chip_WQI1,chip_MQI2)
	ScorrW1W2=spearman_QI_corr(chip_WQI1,chip_WQI2)
	ScorrM1W1=spearman_QI_corr(chip_MQI1,chip_WQI1)
	ScorrM2W2=spearman_QI_corr(chip_MQI2,chip_WQI2)
	Scorr_whole=cbind(ScorrM1M2,ScorrM1W2,ScorrW1M2,ScorrW1W2,ScorrM1W1,ScorrM2W2)
	Scorr_whole=t(Scorr_whole)
	rownames(Scorr_whole)=c("Spearman_corr(MQI1,MQI2)","Spearman_corr(MQI1,WQI2)","Spearman_corr(WQI1,MQI2)","Spearman_corr(WQI1,WQI2)","Spearman_corr(MQI1,WQI1)","Spearman_corr(MQI2,WQI2)")
	#chip_name=cbind("                        Chip 1","Chip 2","Merge")
	if(chipsize_choice=="Affy 100K"){
		chip_name=cbind("                        Hind","Xba","Merge")
	}else if(chipsize_choice=="Affy 500K"){
		chip_name=cbind("                        Nsp","Sty","Merge")
	}else{
		chip_name=cbind(paste("                        ", chipname_list[1], sep = ""), chipname_list[1], chipname_list[1])
	}
	corr_whole=rbind(chip_name,corr_whole,Scorr_whole)
	##--------------------------------------------------
	pathway_NR_QIsummary=paste(pathway_NR,"/QI summary",sep="")
	dir.create(pathway_NR_QIsummary,showWarnings=F)
	TQS=paste(pathway_NR_QIsummary,"/Table_QI_Summary.txt",sep="") ##TQS=Tabele_QI_Summary

	write.table(corr_whole,file=TQS,row.names =T,quote = F,col.names =F,append=F,sep="\t")
	split_line=c("=========================================================")
	write.table(split_line,file=TQS,row.names =F,quote = F,col.names = F,append=T,sep="\t")
	##==============================================================================

	## WinsMean QI1   WinsMean QI2   Median QI1 Median QI2
	QI_name=c("WinsMean QI1","WinsMean QI2","Median QI1","Median QI2")
	UB_QI=UB_QI[c(3,4,1,2),]
	for(k in 1:4)
	{
		UBQI=UB_QI[k,]
		QI1=basic_stat_table(chip_QI[[k]],UBQI,QI_name[k],QI_quantile_choice,chipsize_choice)
		write.table(QI1,file=TQS,row.names =T,quote = F,col.names =F,append=T,sep="\t")
		write.table(split_line,file=TQS,row.names =F,quote = F,col.names = F,append=T,sep="\t")
	}
	##---------------------------------------------------------------------------

	##==========================================================================
	##------------  result modified as one chip --------------------------------
	if(chiptype_choice=="1chip")  ##<---------------
	{
		temp_summary_table=read.delim(TQS,fill=T,sep="\t")#temp_summary_table[,1]=sprintf("%s",temp_summary_table[,1])
		temp_summary_table=temp_summary_table[,c(1,3)]
		dashline=paste("=====================================",sep="")
		output_pathway=TQS
		#temp_summary_table_A=temp_summary_table[c(2:9),]
		temp_summary_table_A=temp_summary_table[1:12,]
		title_name=paste(" "," ",chipname_list[1],sep="\t\t")
		write.table(title_name,file=output_pathway,quote=F,col.names=F,row.names=F,append=F,sep="\t\t")
		write.table(temp_summary_table_A,file=output_pathway,quote=F,col.names=F,row.names=F,append=T,sep="\t")
		write.table(dashline,file=output_pathway,quote=F,col.names=F,row.names=F,append=T,sep=" ")
		for(k in 1:4)
		{
			Head=15+(k-1)*22;End=Head+19
			B=temp_summary_table[c(Head:End),]
			title_name=paste(QI_name[k],chipname_list[1],sep="\t\t\t")
			write.table(title_name,file=output_pathway,quote=F,col.names=F,row.names=F,append=T,sep="\t\t")
			B=B[-c(2*c(1:10)),]
			write.table(B,file=output_pathway,quote=F,col.names=F,row.names=F,append=T,sep="\t\t")
			write.table(dashline,file=output_pathway,quote=F,col.names=F,row.names=F,append=T,sep=" ")
		}
	}
	cat("   ~ The summary table of QI is finished.\n",file=logfile,append=T)
	cat("\n",file=logfile,append=T)

      
	##==========================================================================
	##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	################ Chapter 3 List of over upper of chip ######################
	#---------------------------------------------------------------------------
	# if(Poor_chipchr_option[1]==1)
	if(poor_chipchr_choice==2)
	{
        cat("   ~ The identification of poor chip is preparing.\n",file=logfile,append=T)
        ##--------------------------------------------------------------------------
        QI_file=list.files(pathway_QI_table,full.names=T)
        QI_name=list.files(pathway_QI_table,full.names=F)
        QI_name=strsplit(QI_name,".txt")
        QI_name=unlist(QI_name)
        nstat=length(QI_file)
        chip_MQI1=read.delim(QI_file[1],header=T,sep="") #Median QI1
        chip_MQI2=read.delim(QI_file[2],header=T,sep="") #Median QI2
        chip_WQI1=read.delim(QI_file[3],header=T,sep="") #WinsMean QI1
        chip_WQI2=read.delim(QI_file[4],header=T,sep="") #WinsMean QI2
        
        if(chiptype_choice=="2chip"){specific_chip=c(1:4)}
        if(chiptype_choice=="1chip"){specific_chip=c(1,2,2,2)}
        chip_MQI1=as.matrix(chip_MQI1);chip_MQI1=chip_MQI1[,specific_chip];
        chip_MQI2=as.matrix(chip_MQI2);chip_MQI2=chip_MQI2[,specific_chip];
        chip_WQI1=as.matrix(chip_WQI1);chip_WQI1=chip_WQI1[,specific_chip];
        chip_WQI2=as.matrix(chip_WQI2);chip_WQI2=chip_WQI2[,specific_chip];
        chip_QI=list(chip_MQI1,chip_MQI2,chip_WQI1,chip_WQI2)
        
          ##--------------------------------------------------------------------------
        pathway_DBUB=paste(pathway_DB,"/Upper_bound",sep="")
        UB_Merge=read.delim(paste(pathway_DBUB,"/Merge.txt",sep=""),head=T,sep="\t")
        UB_Chip1=read.delim(paste(pathway_DBUB,"/Chip1.txt",sep=""),head=T,sep="\t")
        UB_Chip2=read.delim(paste(pathway_DBUB,"/Chip2.txt",sep=""),head=T,sep="\t")
        
		# if(QI_quantile_choice==0.95){UB_QI=cbind(UB_Chip1[,1],UB_Chip2[,1],UB_Merge[,1])}
		# if(QI_quantile_choice==0.975){UB_QI=cbind(UB_Chip1[,2],UB_Chip2[,2],UB_Merge[,2])}
		# if(QI_quantile_choice==0.99){UB_QI=cbind(UB_Chip1[,3],UB_Chip2[,3],UB_Merge[,3])}
		
		if(QI_quantile_choice==0.95){UB_QI=cbind(UB_Chip1[,1 + 1],UB_Chip2[,1 + 1],UB_Merge[,1 + 1])}
		if(QI_quantile_choice==0.975){UB_QI=cbind(UB_Chip1[,2 + 1],UB_Chip2[,2 + 1],UB_Merge[,2 + 1])}
		if(QI_quantile_choice==0.99){UB_QI=cbind(UB_Chip1[,3 + 1],UB_Chip2[,3 + 1],UB_Merge[,3 + 1])}
		UB_QI=UB_QI[-c(3,6),];UB_QI=UB_QI[c(3,4,1,2),] ##  row:  MedianQI1, MedianQI2,MeanQI1, MeanQI2;  col: QI_chip1  QI_chip2  QI_merge

        if(chiptype_choice=="1chip")
        {
			if(onechip_choice=="chip2") ichip=2 else ichip=1
			temp_UBQI=rep(UB_QI[,ichip],3)
			UB_QI=matrix(temp_UBQI,nrow(UB_QI),ncol(UB_QI))
        }

      ##----------------------------------------------------------
		pathway_NR_QIpoor=paste(pathway_NR,"/QI poor_chip",sep="")
		dir.create(pathway_NR_QIpoor,showWarnings=F)
        if(chiptype_choice=="2chip")
        {
			split_line=c("========================================================="); dash_line=c("---------------------------------------------------------")
        }
        if(chiptype_choice=="1chip")
        {
			split_line=c("==============================");dash_line=c("------------------------------")
        }
      ##----------------------------------------------------------
		for(k in 1:4)
		{
			#k=1
			IdenQItable=paste(pathway_NR_QIpoor,"/Poor_chip_list_",QI_name[k],".txt",sep="") ##TQS=Tabele_QI_Summary
			if(chiptype_choice=="2chip")
			{
				if(chipsize_choice=="Affy 100K"){head_title=paste("Obs","IndID"," Hind"," Xba"," Merge",sep="\t")}
				if(chipsize_choice=="Affy 500K"){head_title=paste("Obs","IndID"," Nsp"," Sty"," Merge",sep="\t")}
			}
			if(chiptype_choice=="1chip"){head_title=paste("Obs","IndID",chipname_list[1],sep="\t")}
			write.table(head_title,file=IdenQItable,row.names =F,quote = F,col.names = F,append=F,sep="\t")
			write.table(split_line,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
			temp_QIUB=as.numeric(UB_QI[k,])
			temp_chip_QI=chip_QI[[k]]
			temp_chip_QI=temp_chip_QI[,c(1,3,4,2)]  #Ind Chip1 Chip2 Merge
			list_ind_QI=over_upper_list(temp_QIUB,temp_chip_QI)
			if(chiptype_choice=="2chip")
			{
				UB_table=paste("","",paste(" ",sprintf("%.4f",temp_QIUB[1]),sep=""),paste(" ",sprintf("%.4f",temp_QIUB[2]),sep=""),paste(" ",sprintf("%.4f",temp_QIUB[3]),sep=""),sep="\t")
				write.table(UB_table,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
				write.table(dash_line,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
			}
			if(chiptype_choice=="1chip")
			{
				list_ind_QI=list_ind_QI[,-c(4,5)]
				UB_table=paste("","",paste(" ",sprintf("%.4f",temp_QIUB[1]),sep=""),sep="\t")
				write.table(UB_table,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
			  write.table(dash_line,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
			}
        
			write.table(list_ind_QI,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
			write.table(dash_line,file=IdenQItable,row.names =F,quote = F,col.names = F,append=T,sep="\t")
			cat(" * - The marker of poor chip.\n",file=IdenQItable,append=T)
		}
		cat("   ~ The identification of poor chip is finished.\n",file=logfile,append=T)
		cat("\n",file=logfile,append=T)
    }
}
##==============================================================================
cat("****************************************************************************************************************\n")



###-------------------Numerical result modification region-----------------------
### Due to the format of result could not management, i develope the sub-fucntion
### which can manage the result format and make the result clear.
#
library(base)
NR_name=list.files(pathway_NR)

pathway_NR_QI=paste(pathway_NR,"/QI estimate",sep="")
if(QI_calculation_choice == 2 & QI_est_choice == 2)
{
	pathway_QI_table=paste(pathway_NR_QI,"/QI table",sep="")
	if(file.exists(pathway_QI_table)==T){
		QI_file=list.files(pathway_QI_table,full.names=T)
		QI_name=list.files(pathway_QI_table,full.names=F)
		nQIstat=length(QI_name)
		for(k in 1:nQIstat){
			chr=read.delim(QI_file[k],head=T,sep="\t", stringsAsFactors = F)
			ind_name=as.character(chr[,1]);max_width_indname=max(nchar(ind_name,type="width"))
			ind_name=format(ind_name,width=max_width_indname)
			QI=chr[,-1];
			QIname=names(QI)
			if(chiptype_choice=="2chip"){QIname[1]=chipname_list[3];QIname[2]=chipname_list[1];QIname[3]=chipname_list[2];}
			if(chiptype_choice=="1chip"){QIname[1]=chipname_list[1];}
			IndIDtitle=format("IndID",width=max(5,max_width_indname))
			QIname=format(QIname,width=6)
			QIname=c(IndIDtitle,QIname)
			QI=as.matrix(QI);Nrow=nrow(QI);Ncol=ncol(QI)
			QI=matrix(sprintf("%.4f",QI),Nrow,Ncol)
			QI=format(QI,width=6)
			QI=cbind(ind_name,QI)
			QI=data.frame(QI);colnames(QI)=QIname
			#write.table(QI,paste("J:/test.txt"),col.names=T,row.names=F,quote=F,sep="\t")
			write.table(QI,QI_file[k],col.names=T,row.names=F,quote=F,sep="\t")
		}
	}
	##----------------------------------------------------------------------------
	pathway_IndQI=paste(pathway_NR_QI,"/Ind_QI",sep="")
	if(file.exists(pathway_IndQI)==T){
		IndQI_file=list.files(pathway_IndQI,full.names=T) 
		IndQI_name=list.files(pathway_IndQI,full.names=F)     
		nind=length(IndQI_name)
		for(k in 1:nind){
			#k=1
			chr_file=list.files(IndQI_file[k],full.names=T)
			nchr=length(chr_file) 
			for(i in 1:nchr){
				#i=1
				chr=read.delim(chr_file[i],head=T,sep="\t", stringsAsFactors = F) 
				chiptype=chr[,2]
				chiptype[chiptype==1]="Chip1" 
				chiptype[chiptype==2]="Chip2"
				chr[,2]=chiptype         
				Genotype=chr[,4] 
				Genotype[Genotype==1]="BB" 
				Genotype[Genotype==2]="AB" 
				Genotype[Genotype==3]="AA" 
				chr[,4]=Genotype 

				AF=as.numeric(chr[,5]);AF=round(AF*1000)/1000;AF=sprintf("%1.3f",AF)
				chr[,5]=AF                        

				# QI1=as.numeric(chr[,6]);QI1=round(QI1*1000)/1000;QI1=sprintf("%3.4f",QI1)
				QI1=as.numeric(chr[,6]);QI1=round(QI1*10000)/10000;QI1=sprintf("%3.4f",QI1)
				chr[,6]=QI1                        
				# QI2=as.numeric(chr[,7]);QI2=round(QI2*1000)/1000;QI2=sprintf("%3.4f",QI2)                          
				QI2=as.numeric(chr[,7]);QI2=round(QI2*10000)/10000;QI2=sprintf("%3.4f",QI2)                                
				chr[,7]=QI2

				chr=sort.append(chr,3) 
				chr=data.frame(chr)
				colnames(chr)[2]="Chiptype" 
				colnames(chr)[3]="Phy_position" 
				colnames(chr)[4]="Genotype" 
				colnames(chr)[5]="   AdjAF" 
				colnames(chr)[6]="     QI1"                 
				colnames(chr)[7]="     QI2"                         
				pathway_result=chr_file[i]#paste("F:/text.txt",sep="")
				chr=format(chr,width=8,justify="right") 
				write.table(chr,pathway_result,col.names=T,row.names=F,quote=F,sep="\t") 
			}
		} 
	} 
}
chipname=unlist(strsplit(chipname_list," "))[c(1:3)*2]

pathway_sample_list=paste(pathway_output_population,"/Sample list.txt",sep="")
if(file.exists(pathway_sample_list)==T)
{
	A=read.delim(pathway_sample_list,head=T,sep="\t")
	Colnames=names(A)
	if(chiptype_choice=="2chip"){Colnames[3:5]=paste("CR_",chipname[1:3],sep="")}
	if(chiptype_choice=="1chip"){Colnames[3]=paste("CR_",chipname[1],sep="")}
	ind_name=as.character(A[,2]);max_width_indname=max(nchar(ind_name,type="width"))
	ind_name=format(ind_name,width=max_width_indname)
	Colnames[2]=format(Colnames[2],width=max_width_indname)
	max_width_NCR=max(nchar(Colnames[4],type="width"))
	if(chiptype_choice=="2chip"){for(k in 3:5){A[,k]=format(A[,k],width=max_width_NCR)}}
	if(chiptype_choice=="1chip"){for(k in 3){A[,k]=format(A[,k],width=max_width_NCR)}}
	colnames(A)=Colnames
	sample_list_table=pathway_sample_list
	write.table(A,sample_list_table,col.names=T,row.names=F,quote=F,sep="\t")
}
##------------------------------------------------------------------------------
if(CPA_est_choice == 1) unlink(pathway_CPA,recursive = T)

#if(AF_est_choice != 3){
#if(AF_est_choice == 1){
if(AF_est_choice != 2){
	unlink(pathway_adjAF,recursive = T)
} else {
	AF_file = mixedsort(list.files(pathway_adjAF,full.names=T))
	AF_name = mixedsort(list.files(pathway_adjAF,full.names=F))
	nind = length(AF_name)
	for(i in 1:nind){
		chr_file = list.files(AF_file[i], full.names = T)
		chr_name = list.files(AF_file[i], full.names = F)
		for(k in 1:length(chr_name)){ #for each chromosome
			chr=read.delim(chr_file[k], stringsAsFactors = F)
			if(AF_est.based_choice[1]){
				colnames(chr)[ncol(chr)] = "AF(Intensity-based)"
				chr$"AF(Intensity-based)" = sprintf("%1.4f", chr$"AF(Intensity-based)")
			} else {
				chr$adjAF = NULL
			}
			if(AF_est.based_choice[2]){
				AF.geno = chr$Genotype
				AF.geno = sapply(strsplit(AF.geno, " "), function(x) x[length(x)])
				AF.geno = as.character(factor(AF.geno, labels = c("0", "0.5", "1"), levels = c("BB", "AB", "AA")))
				chr$"AF(Genotype-based)" = sprintf("%.1f", as.numeric(AF.geno))
				#chr$"AF(Genotype-based)" = AF.geno
			}
			chr[,1]=sprintf("%9s",chr[,1])
			chr[,2]=sprintf("%9s",chr[,2])
			chr[,3]=sprintf("%8s",chr[,3])
			chr[,4]=sprintf("%8s",chr[,4])
			write.table(chr, chr_file[k], col.names = T, row.names = F, quote = F, sep = "\t")
		}
	}
}
#}

if(QI_est_choice == 1) unlink(pathway_NR_QI,recursive = T)

##------------------------------------------------------------------------------
u=proc.time() - ptm
sec=u[3]
hour=sec/3600
minute=(hour-floor(hour))*60
second=(minute-floor(minute))*60
time_stance=paste("Elapsed time: ",floor(hour),"-H, ",floor(minute),"-M, ",floor(second),"-S ",sep="")
write.table(" ",file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
write.table(time_stance,file=data_descrption,row.names =F,quote = F,col.names =F,append=T,sep="\t")
cat("\n");cat("The calculation has been finished!!!\n")
cat("   ~ Whole calculation have been done.\n",file=logfile,append=T)
#rm(list=ls())

