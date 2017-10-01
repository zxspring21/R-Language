##------------------parameter setting---------------------------

################################################################################
##-------------------------- 1 Input/output path:-----------------------------##
################################################################################
inputdata_format=tclvalue(inputdata_format.statsValue)
  ##inputdata_format=1 : Genotype/Intensity based
  ##inputdata_format=2 : Allele Frequency based
data_input=tclvalue(pathway_input)
if(nchar(data_input)==0)
{stop("SAQC error message: The directory of data input is incorrect.")}

data_output=tclvalue(pathway_output)
if(data_input == "Test example 2"){
	temp_data_output=data_output
	if(file.exists(temp_data_output)==T){data_output=paste(temp_data_output,"/Test_Affy100K",sep="")}
	if(file.exists(temp_data_output)==F){data_output=paste(pathway,"OUTPUT/Test_Example_Output",sep="")}
} else if(data_input == "Test example 1"){
	temp_data_output=data_output
	if(file.exists(temp_data_output)==T){data_output=paste(temp_data_output,"/Test_Affy500K",sep="")}
	if(file.exists(temp_data_output)==F){data_output=paste(pathway,"OUTPUT/Test_Example_Output",sep="")}
} else if(data_input == "Test example 3"){
	temp_data_output=data_output
	if(file.exists(temp_data_output)==T){data_output=paste(temp_data_output,"/Test_Illu550K",sep="")}
	if(file.exists(temp_data_output)==F){data_output=paste(pathway,"OUTPUT/Test_Example_Output",sep="")}
} else if(file.exists(data_output)==F) stop("SAQC error message: The directory of result output is incorrect.")
##==============================================================================

################################################################################
##-------------------------- 2 Data format:-----------------------------------##
################################################################################
group_choice=tclvalue(groupValue)
  ##  group_choice=1  :Asia (CHB+JPT)
  ##  group_choice=2  :Africa (YRI)
  ##  group_choice=3  :Europe (CEU)
  ##  group_choice=4  :Taiwan (TWN)
  ##  group_choice=5  :Combine
chipsize_choice=c(tclvalue(chipValue))
  ##  chipsize_choice=1 : Affy 100K
  ##  chipsize_choice=2 : Affy 500K
  ##  chipsize_choice=3 : Affy 6.0
  ##  chipsize_choice=4 : Affy Axiom
  ##  chipsize_choice=5 : Illumina 550K
chiptype_choice=c(tclvalue(chip1.val),tclvalue(chip2.val))
chiptype_choice=sum(as.numeric(chiptype_choice)*c(1:2))
onechip_choice=c(tclvalue(hind.val),tclvalue(xba.val),tclvalue(nsp.val),tclvalue(sty.val))
onechip_choice=sum(as.numeric(onechip_choice)*c(1,2,1,2))

# NA_str = c(tclvalue(NA_str))
# index.probe = c(tclvalue(index.probe))
##  onechip_choice: chip1= "Hind" "Nsp";  chip2=" Xbar" "Sty"
  ##  chipsize=c("chip1","chip2","merge")
# chromosome_choice=c(tclvalue(chr1.val),tclvalue(chr2.val),tclvalue(chr3.val),
# tclvalue(chr4.val),tclvalue(chr5.val),tclvalue(chr6.val),tclvalue(chr7.val),
# tclvalue(chr8.val),tclvalue(chr9.val),tclvalue(chr10.val),tclvalue(chr11.val),
# tclvalue(chr12.val),tclvalue(chr13.val),tclvalue(chr14.val),tclvalue(chr15.val),
# tclvalue(chr16.val),tclvalue(chr17.val),tclvalue(chr18.val),tclvalue(chr19.val),
# tclvalue(chr20.val),tclvalue(chr21.val),tclvalue(chr22.val),tclvalue(chr23.val),
# tclvalue(chrall.val))
# chromosome_choice=as.numeric(chromosome_choice[-24])
# chromosome_choice=chromosome_choice*(1:23)

##==============================================================================
##==============================================================================


################################################################################
##-------------------------- 3 Statistical analysis:--------------------------##
################################################################################

    ##++++++++ 3-1 CPA calculation +++++++
CPA_calculation_choice=tclvalue(CPA.statsValue)
  ##  CPA_calculataion_choice=1 NO,from database; 2 Yes,User-provided data
CPA_calculation_choice=as.numeric(CPA_calculation_choice)
if(CPA_calculation_choice==0){CPA_calculation_choice=3}
    ##++++++++ 3-2 AF calculation +++++++
AF_calculation_choice=tclvalue(AF.statsValue)
  ##  AF_calculation_choice=1 NO; 2 Yes;
AF_calculation_choice=as.numeric(AF_calculation_choice)
if(AF_calculation_choice==0){AF_calculation_choice=3}
    ##++++++++ 3-3 AF reference  +++++++
AF_ref_choice=tclvalue(AF_ref.Value)
AF_ref_choice=as.numeric(AF_ref_choice)
if(AF_ref_choice==0){AF_ref_choice=3}
  ##  AF_ref_choice=1 NO; 2 Yes;
   AFref_source1_choice=tclvalue(AFref_source1.val)
   AFref_source2_choice=tclvalue(AFref_source2.val)
   AFref_source_choice=c(AFref_source1_choice,AFref_source2_choice)
   AFref_source_choice=as.numeric(AFref_source_choice)
   AFref_input=tclvalue(pathway_AFrefinput)
   if(AFref_source_choice[2]==1)
   {
      Null_SNP_AFref_table=AFref_checking(AFref_input)
      cat("SAQC message: The AF reference of user provided is suitable for SAQC calculation.\n")
   }
   
    ##++++++++ 3-4 QI calculation ++++++
QI_calculation_choice=tclvalue(QIValue)
  ##  QI_calculation_choice=1 NO; 2 Yes;
QI_calculation_choice=as.numeric(QI_calculation_choice)
if(QI_calculation_choice==0){QI_calculation_choice=3}
    ##++++++++ 3-5 QI identification +++++
Iden_poor_chipchr_choice=as.numeric(tclvalue(QI_iden.statsValue))
if(Iden_poor_chipchr_choice==0){Iden_poor_chipchr_choice=3}
  ##  Iden_poor_chipchr_choice=1 NO; 2 Yes;
QI_quantile_choice=c(tclvalue(QI_quantile1.val),tclvalue(QI_quantile2.val),tclvalue(QI_quantile3.val))
QI_quantile_choice=as.numeric(QI_quantile_choice);QI_quantile_choice=QI_quantile_choice%*%c(1,2,3)

if(Iden_poor_chipchr_choice==2)
{
  if(QI_quantile_choice==0)
  {
    {stop("SAQC error message: Please enter the QI quantile of upper bound.")}
  }
  if(!QI_quantile_choice==0)
  {
    QI_quantile=c(0.95,0.975,0.99)
    QI_quantile_choice=QI_quantile[QI_quantile_choice]
  }
}
  ##  QI_quantile_choice=c("0.95","0.975","0.99")
##==============================================================================
##==============================================================================


################################################################################
##--------------------------- 4 Graphical output------------------------------##
################################################################################
##
AF.plot_choice=as.numeric(tclvalue(AF_plot.statsValue))   ## No:1 ,Yes:2
if(AF.plot_choice==0){AF.plot_choice=3}
AF.based_choice=as.numeric(c(tclvalue(AF_type1.val),tclvalue(AF_type2.val)))  ##  Intensity_based,Genotype_based
CallRateplot_choice=as.numeric(tclvalue(CallRate_plot.statsValue))    ## No:1 ,Yes:2
    ## HeatMap
#QI.plot_choice=as.numeric(tclvalue(QI_plot.statsValue))   ## No:1 ,Yes:2
#QI.plot_choice is replaced by HM_plot_choice
HM_plot_choice=as.numeric(tclvalue(HM_plot.statsValue))   ## No:1 ,Yes:2
if(HM_plot_choice==0){HM_plot_choice=3}
col.scale_choice=c(as.numeric(tclvalue(col_scale_type1.val)),as.numeric(tclvalue(col_scale_type2.val)))

    ## Polygon
#CR_plot.statsValue is replaced by PG_plot.statsValue
#QI.scatterplot_choice is replaced by PG_plot_choice
PG_plot_choice=as.numeric(tclvalue(PG_plot.statsValue))
if(PG_plot_choice==0){PG_plot_choice=3}
##==============================================================================
##==============================================================================

################################################################################
##--------------------------- 5 Numerical output------------------------------##
################################################################################
data_description_choice=as.numeric(tclvalue(DD_plot.statsValue))
CPA_est_choice=as.numeric(tclvalue(CPA_est_plot.statsValue))
AF_est_choice=as.numeric(tclvalue(AF_numeric.statsValue))
if(AF_est_choice==0){AF_est_choice=3}
AF_est.based_choice=as.numeric(c(tclvalue(AF_numeric_type1.val),tclvalue(AF_numeric_type2.val)))
QI_est_choice=as.numeric(tclvalue(QI_numeric.statsValue))
if(QI_est_choice==0){QI_est_choice=3}
poor_chipchr_choice=as.numeric(tclvalue(poor_numeric.statsValue))
if(poor_chipchr_choice==0){poor_chipchr_choice=3}
##==============================================================================

##----------------------chromosome debuging area--------------------------------
#if(sum(chromosome_choice)==0){stop("SAQC error message: Please select the chromosome which you want!!.")}
#if(sum(chromosome_choice)!=0){source(SAQC.file)}
if(inputdata_format == "1"){
	if(data_input == "Test example 1"){
		NA_str = "null"
		ParaSet = c(2, 2, 3, 4, 6, NA, NA, NA, NA)
	}else if(data_input == "Test example 3"){
		NA_str = c("NaN", "-")
		ParaSet = c(11, 1, 2, 3, 10, 11, 6, 7, 4)
	}else{
		NA_str = unlist(strsplit(tclvalue(NA_str), split = ","))
		ParaSet = as.numeric(c(tclvalue(Skiprow), tclvalue(index.probe), tclvalue(index.chr), tclvalue(index.posi), tclvalue(index.CallA), tclvalue(index.CallB), tclvalue(index.PI_A), tclvalue(index.PI_B), tclvalue(index.BAF))) 
	}
	# cat(ParaSet)
}
source(SAQC.file)