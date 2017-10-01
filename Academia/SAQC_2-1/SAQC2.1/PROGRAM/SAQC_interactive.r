#-----------------------setting the libaray path-------------------------------
# download the packages, as the R-gui don't have
libname = .packages(all.available = TRUE)
if(!"tcltk" %in% libname) install.packages("tcltk", repos = "http://cran.csie.ntu.edu.tw")
if(!"tcltk2" %in% libname) install.packages("tcltk2", repos = "http://cran.csie.ntu.edu.tw")
if(!"gtools" %in% libname) install.packages("gtools", repos = "http://cran.csie.ntu.edu.tw")
if(!"boot" %in% libname) install.packages("boot", repos="http://cran.csie.ntu.edu.tw")
##==============================================================================
##R-GUI for Quality Index
library(tcltk); 
library(tcltk2) # for multi-page GUI

pathway_program = strsplit(SAQCgui, "SAQC_interactive.r")
pathway_program = unlist(pathway_program)

SAQC.QIpara = paste(pathway_program, "QI_parameter_setting.r", sep = "")
#addTclPath("C:\\Tcl\\lib")
#tclRequire("BWidget")
####################################################################################
							#Define Subfunction by Masaki
###################################################################################
button.disable = function(x) for(i in x) tkconfigure(i, state = "disable") #按鈕設成無效

button.enable = function(x) for(i in x) tkconfigure(i, state = "normal") #按鈕設成無效

interactive_polygon = function(index.decide.ref = F, button.stat = "disable", button.stat1 = "disable", button.stat2 = "disable", 
								default.val1.1 = 0,	default.val1.2 = 0, default.val1.3 = 0, default.val1.4 = 0,
								default.val2.1 = "", default.val2.2 = "", default.val2.3 = ""){
	#---------------------------------- SAQC database -----------------------------------------------#
	if(!index.decide.ref){
		button.stat = button.stat1
		tkconfigure(database.but, variable = QI_resouce.val, stat = button.stat)
		tkconfigure(userprovided.but, variable = QI_resouce.val, stat = button.stat)
		tclvalue(QI_resouce.val) = default.val1.1
	}
	##Study population
	tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, stat = button.stat1)
	tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, stat = button.stat1)
	tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, stat = button.stat1)
	tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, stat = button.stat1)
	tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, stat = button.stat1)
	tclvalue(QIfigure_groupValue) = default.val1.2
	##Genome-wide SNP array
	tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, stat = button.stat1)
	tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, stat = button.stat1)
	tclvalue(QIfigure_chipValue) = default.val1.3
	##QI upper quantile
	tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = button.stat1)
	tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = button.stat1)
	tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = button.stat1)
	tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = button.stat1)
	tclvalue(QIfigure_quantile.Value1) = default.val1.4; tclvalue(QIfigure_quantile.Value2) = default.val1.4
	tclvalue(QIfigure_quantile.Value3) = default.val1.4; tclvalue(QIfigure_quantile.Value4) = default.val1.4
	#--------------------------------- User provided -----------------------------------------------#
	tkconfigure(MQI1_chip1.q95.box, state = button.stat2);tclvalue(MQI1_chip1.q95.val) = default.val2.1
	tkconfigure(MQI1_chip2.q95.box, state = button.stat2);tclvalue(MQI1_chip2.q95.val) = default.val2.1
	tkconfigure(MQI1_merge.q95.box, state = button.stat2);tclvalue(MQI1_merge.q95.val) = default.val2.1
	tkconfigure(MQI2_chip1.q95.box, state = button.stat2);tclvalue(MQI2_chip1.q95.val) = default.val2.1
	tkconfigure(MQI2_chip2.q95.box, state = button.stat2);tclvalue(MQI2_chip2.q95.val) = default.val2.1
	tkconfigure(MQI2_merge.q95.box, state = button.stat2);tclvalue(MQI2_merge.q95.val) = default.val2.1
	tkconfigure(WQI1_chip1.q95.box, state = button.stat2);tclvalue(WQI1_chip1.q95.val) = default.val2.1
	tkconfigure(WQI1_chip2.q95.box, state = button.stat2);tclvalue(WQI1_chip2.q95.val) = default.val2.1
	tkconfigure(WQI1_merge.q95.box, state = button.stat2);tclvalue(WQI1_merge.q95.val) = default.val2.1
	tkconfigure(WQI2_chip1.q95.box, state = button.stat2);tclvalue(WQI2_chip1.q95.val) = default.val2.1
	tkconfigure(WQI2_chip2.q95.box, state = button.stat2);tclvalue(WQI2_chip2.q95.val) = default.val2.1
	tkconfigure(WQI2_merge.q95.box, state = button.stat2);tclvalue(WQI2_merge.q95.val) = default.val2.1
	tkconfigure(MQI1_chip1.q975.box, state = button.stat2);tclvalue(MQI1_chip1.q975.val) = default.val2.2
	tkconfigure(MQI1_chip2.q975.box, state = button.stat2);tclvalue(MQI1_chip2.q975.val) = default.val2.2
	tkconfigure(MQI1_merge.q975.box, state = button.stat2);tclvalue(MQI1_merge.q975.val) = default.val2.2
	tkconfigure(MQI2_chip1.q975.box, state = button.stat2);tclvalue(MQI2_chip1.q975.val) = default.val2.2
	tkconfigure(MQI2_chip2.q975.box, state = button.stat2);tclvalue(MQI2_chip2.q975.val) = default.val2.2
	tkconfigure(MQI2_merge.q975.box, state = button.stat2);tclvalue(MQI2_merge.q975.val) = default.val2.2
	tkconfigure(WQI1_chip1.q975.box, state = button.stat2);tclvalue(WQI1_chip1.q975.val) = default.val2.2
	tkconfigure(WQI1_chip2.q975.box, state = button.stat2);tclvalue(WQI1_chip2.q975.val) = default.val2.2
	tkconfigure(WQI1_merge.q975.box, state = button.stat2);tclvalue(WQI1_merge.q975.val) = default.val2.2
	tkconfigure(WQI2_chip1.q975.box, state = button.stat2);tclvalue(WQI2_chip1.q975.val) = default.val2.2
	tkconfigure(WQI2_chip2.q975.box, state = button.stat2);tclvalue(WQI2_chip2.q975.val) = default.val2.2
	tkconfigure(WQI2_merge.q975.box, state = button.stat2);tclvalue(WQI2_merge.q975.val) = default.val2.2
	tkconfigure(MQI1_chip1.q99.box, state = button.stat2);tclvalue(MQI1_chip1.q99.val) = default.val2.3
	tkconfigure(MQI1_chip2.q99.box, state = button.stat2);tclvalue(MQI1_chip2.q99.val) = default.val2.3
	tkconfigure(MQI1_merge.q99.box, state = button.stat2);tclvalue(MQI1_merge.q99.val) = default.val2.3
	tkconfigure(MQI2_chip1.q99.box, state = button.stat2);tclvalue(MQI2_chip1.q99.val) = default.val2.3
	tkconfigure(MQI2_chip2.q99.box, state = button.stat2);tclvalue(MQI2_chip2.q99.val) = default.val2.3
	tkconfigure(MQI2_merge.q99.box, state = button.stat2);tclvalue(MQI2_merge.q99.val) = default.val2.3
	tkconfigure(WQI1_chip1.q99.box, state = button.stat2);tclvalue(WQI1_chip1.q99.val) = default.val2.3
	tkconfigure(WQI1_chip2.q99.box, state = button.stat2);tclvalue(WQI1_chip2.q99.val) = default.val2.3
	tkconfigure(WQI1_merge.q99.box, state = button.stat2);tclvalue(WQI1_merge.q99.val) = default.val2.3
	tkconfigure(WQI2_chip1.q99.box, state = button.stat2);tclvalue(WQI2_chip1.q99.val) = default.val2.3
	tkconfigure(WQI2_chip2.q99.box, state = button.stat2);tclvalue(WQI2_chip2.q99.val) = default.val2.3
	tkconfigure(WQI2_merge.q99.box, state = button.stat2);tclvalue(WQI2_merge.q99.val) = default.val2.3
}
##################################################################################
##---------------color changing start part------------------------
color = "aliceblue"
#color = "azure"#"mistyrose"#"aquamarine2"#"#80ff00"
##---------------color changing start part------------------------
##------------
QI_interface = tktoplevel()#(bg = color)

fontHeading = tkfont.create(family = "times", size = 14, weight = "bold", slant = "italic")
fontIntro = tkfont.create(family = "arial", size = 10)
fontIntro_para = tkfont.create(family = "arial", size = 1)
fontsubtitle = tkfont.create(family = "times", size = 14, weight = "bold")
fontTextLabel = tkfont.create(family = "times", size = 10)
fontTextLabel_para = tkfont.create(family = "times", size = 11)
fontFixedWidth = tkfont.create(family = "courier", size = 10)
fontRun = tkfont.create(family = "courier", size = 12, weight = "bold", slant = "italic", underline = F, overstrike = F)
color.run = "lavender"; button.space = 14; blankspcae.run = 45
##---------------------------------------------------------------------------

##--------------------------
##--------------------------
tktitle(QI_interface) = "SAQC v2.1"

#####Ex. 1: Add a text labels in a window ####
##------------------------------------------------------------------------
QI_interface_A = tkframe(QI_interface)  #, background = color
QI_interface_B = tkframe(QI_interface)  #, background = color
a = tklabel(QI_interface_A, text = "     SNP Array Quality Control (SAQC)", width = 29, font = fontHeading)
b = tklabel(QI_interface_A, text = "is developed for quality examination of SNP arrays. SAQC incorporates", font = fontIntro)
tkgrid(a, b)
tkgrid.configure(a, sticky = "es")
tkgrid.configure(b, sticky = "w")
##------------------
tkgrid(tklabel(QI_interface_B, text = "       allele frequency (AF), quality indices (QI) and graphical visualization to identify abnormal samples/arrays/chromosomes for ", font = fontIntro), sticky = "w")#, bg = color
tkgrid(tklabel(QI_interface_B, text = "       genome-wide SNP array data.", font = fontIntro), sticky = "w")
tkgrid(QI_interface_A, sticky = "w")
tkgrid(QI_interface_B, sticky = "w")


##==============================================================================
##---------------------------   Notebook 1  ------------------------------------
noteframe = tkframe(QI_interface)
nb <- tk2notebook(noteframe, tabs = c("Main functions", "Interactive visualization"), width = 800, height = 740) #height default: 670
tkpack(nb, fill = "both", expand = 0) # expand = 0 maybe not work

##
##	Modify the content of the first tab: "Main functions"
##
tb1 <- tk2notetab(nb, "Main functions")

frame <- tkframe(tb1, bg = color)
frame.dataformat <- tkframe(frame, bg = color)
frame.para1 <- tkframe(frame, bg = color)
frame.para1A <- tkframe(frame, bg = color)
frame.inputformat <- tkframe(frame, bg = color)
frame.path <- tkframe(frame, bg = color)

frame.para2 <- tkframe(frame, bg = color)
frame.chip <- tkframe(frame, bg = color)
frame.chip2 <- tkframe(frame, bg = color)
frame.group <- tkframe(frame, bg = color)
frame.para <- tkframe(frame, bg = color)
frame.paratitle3 <- tkframe(frame, bg = color)

frame.SA <- tkframe(frame, bg = color)
frame.para3 <- tkframe(frame, bg = color)
	frame.para3A <- tkframe(frame.para3, bg = color)
	frame.para3C <- tkframe(frame.para3, bg = color)
	frame.para3B <- tkframe(frame.para3, bg = color)
	frame.para3D <- tkframe(frame.para3, bg = color)
frame.para3E <- tkframe(frame, bg = color)
frame.para4 <- tkframe(frame, bg = color)
	frame.para4A <- tkframe(frame, bg = color)
	frame.para4B <- tkframe(frame, bg = color)
	frame.para4C <- tkframe(frame, bg = color)
frame.para5 <- tkframe(frame, bg = color)
	frame.para5A <- tkframe(frame, bg = color)
##==============================================================================
##==============================================================================
##-----------------------------------------------------------------------------
							##1. Input/output path
##-----------------------------------------------------------------------------

tkgrid(tklabel(frame.para1, text = "   1. Input/output path:", font = fontsubtitle, bg = color), sticky = "w") # title
tkgrid.configure(frame.para1, sticky = "w")  ##???????????????
##-----------------------------------------------------------------------------------------

#-------------------------------# 1-1 Input Data format selection #-----------------------#
inputdata_format.port1 <- tkradiobutton(frame.para1A, bg = color, command = function(){
	if((tclvalue(pathway_input) == "Test example 2")|(tclvalue(pathway_input) == "Test example 3")){
		tclvalue(pathway_input) = "Test example 1"
		tclvalue(chipValue) = 2 # set to 500K
		tclvalue(groupValue) = 5 #set to combine
		tclvalue(chip1.val) = 0
		tclvalue(chip2.val) = 1
		tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(chip2.port, state = "normal")
		tkconfigure(hind.port, variable = hind.val, state = "disable")
		tkconfigure(xba.port, variable = xba.val, state = "disable")
		tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		tkconfigure(sty.port, variable = sty.val, state = "disable")
	}else if((tclvalue(pathway_input) == "Test example 2")|(tclvalue(pathway_input) == "Test example 1")){
		tclvalue(pathway_input) = "Test example 3"
		tclvalue(chipValue) = 5 # set to 550K
		tclvalue(groupValue) = 4 #set to TWN
		tclvalue(chip1.val) = 1
		tclvalue(chip2.val) = 0
		tkconfigure(chip2.port, state = "disable")
		tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(hind.port, variable = hind.val, state = "disable")
		tkconfigure(xba.port, variable = xba.val, state = "disable")
		tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		tkconfigure(sty.port, variable = sty.val, state = "disable")
	}
	tkconfigure(group.pp1, variable = groupValue, state = "normal")
	tkconfigure(group.pp2, variable = groupValue, state = "normal")
	tkconfigure(group.pp3, variable = groupValue, state = "normal")
	tkconfigure(group.pp5, variable = groupValue, state = "normal")
	#Statistical Analysis
	##CPA calculation
	tkconfigure(CPA.stats.port1, variable = CPA.statsValue, state = "normal")
	tkconfigure(CPA.stats.port2, variable = CPA.statsValue, state = "normal")
	tclvalue(CPA.statsValue) = 1
	##AF calculation
	tkconfigure(AF.stats.port1, variable = AF.statsValue, state = "normal")
	tkconfigure(AF.stats.port2, variable = AF.statsValue, state = "normal")
	tclvalue(AF.statsValue) = 2
	##AF reference
	tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "normal")
    tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "normal")
    tclvalue(AF_ref.Value) = 1
    tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
    tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
    tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
    ##QI calculation
	tkconfigure(QI.port1, variable = QIValue, state = "normal")
    tkconfigure(QI.port2, variable = QIValue, state = "normal")
    tclvalue(QIValue) = 2
    ##Identification of poor array
    tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
    tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
    tclvalue(QI_iden.statsValue) = 2
    ##-----------------------------
    tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
    tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
    tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
    tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 1
    
	#Graphical result
	##AF plot
	tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "normal")
    tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "normal")
    tclvalue(AF_plot.statsValue) = 2
    tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
    tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
    tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
    ##GCR plot
		##	No matter which input data format, GCR plot can draw or not 
		##	So the radiobutton would not reset
    ##HeatMap plot
    tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
    tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
    tclvalue(HM_plot.statsValue) = 2
    tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "normal")
    tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "normal")
    tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
    ##Polygon plot
    tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
    tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
    tclvalue(PG_plot.statsValue) = 2
   	
	#Numerical output
     ##Data description
		##	No matter which input data format, Data description can output or not 
		##	So the radiobutton would not reset
 	##CPA estimate
	# tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "normal")
	# tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "normal")
	# tclvalue(CPA_est_plot.statsValue) = 2
	tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "disable")
	tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "disable")
	tclvalue(CPA_est_plot.statsValue) = 0 # could be any number rather than 1/2
	##AF estimate
	tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue, state = "normal")
	tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue, state = "normal")
	tclvalue(AF_numeric.statsValue) = 2
	###Intensity-based & Genotype-based
	tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "normal")
	tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "normal")
	tclvalue(AF_numeric_type1.val) = 1; tclvalue(AF_numeric_type2.val) = 1
	
	##QI estimate
	tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
    tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
    tclvalue(QI_numeric.statsValue) = 2
    ##PoorSNP array
	tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
    tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
    tclvalue(poor_numeric.statsValue) = 2
    
	# if(0){
		# #Affy6.0
		# tkconfigure(chip.port3, variable = chipValue, state = "disable")
		# #Illumina 550k
		# tkconfigure(chip.port4, variable = chipValue, state = "disable")
		# tclvalue(chipValue) = 1
	# }
})



inputdata_format.port2 <- tkradiobutton(frame.para1A, bg = color, command = function(){
	if((tclvalue(pathway_input) == "Test example 1")|(tclvalue(pathway_input) == "Test example 3")){
		tclvalue(pathway_input) = "Test example 2"
		tclvalue(chipValue) = 1 #set to 100K
		tclvalue(groupValue) = 5 #set to combine
		tclvalue(chip1.val) = 0
		tclvalue(chip2.val) = 1
		tkconfigure(chip2.port, state = "normal")
		tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(hind.port, variable = hind.val, state = "disable")
		tkconfigure(xba.port, variable = xba.val, state = "disable")
		tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		tkconfigure(sty.port, variable = sty.val, state = "disable")
	}
	tkconfigure(group.pp1, variable = groupValue, state = "normal")
	tkconfigure(group.pp2, variable = groupValue, state = "normal")
	tkconfigure(group.pp3, variable = groupValue, state = "normal")
	tkconfigure(group.pp5, variable = groupValue, state = "normal")
	#Statistical Analysis
	##CPA calculation
	tkconfigure(CPA.stats.port1, variable = CPA.statsValue, state = "disable")
	tkconfigure(CPA.stats.port2, variable = CPA.statsValue, state = "disable")
	tclvalue(CPA.statsValue) = 3
	##AF calculation
	tkconfigure(AF.stats.port1, variable = CPA.statsValue, state = "disable")
	tkconfigure(AF.stats.port2, variable = CPA.statsValue, state = "disable")
	tclvalue(AF.statsValue) = 3
	##AF reference
	tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "normal")
    tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "normal")
    tclvalue(AF_ref.Value) = 1
    tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
    tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
    tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
    ##QI calculation
	tkconfigure(QI.port1, variable = QIValue, state = "normal")
    tkconfigure(QI.port2, variable = QIValue, state = "normal")
    tclvalue(QIValue) = 2
    ##Identification of poor array
    tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
    tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
    tclvalue(QI_iden.statsValue) = 2
    ##-----------------------------
    tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
    tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
    tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
    tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 1
    
	#Graphical result
	##AF plot
	tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "normal")
    tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "normal")
    tclvalue(AF_plot.statsValue) = 2
    tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
    tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
    tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
    ##HeatMap plot
    tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
    tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
    tclvalue(HM_plot.statsValue) = 2
    tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "normal")
    tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "normal")
    tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
    ##Polygon plot
    tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
    tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
    tclvalue(PG_plot.statsValue) = 2
   
	#Numerical Result
	##CPA estimate
	tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "disable")
	tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "disable")
	tclvalue(CPA_est_plot.statsValue) = 3
	#AF estimate
	tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue, state = "disable")
	tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue, state = "disable")
	tclvalue(AF_numeric.statsValue) = 3
	tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "disable")
	tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "disable")
	tclvalue(AF_numeric_type1.val) = 3; tclvalue(AF_numeric_type2.val) = 3
	##QI estimate
	tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
    tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
    tclvalue(QI_numeric.statsValue) = 2
    ##PoorSNP array
	tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
    tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
    tclvalue(poor_numeric.statsValue) = 2
    
	# if(0){
		# #Affy6.0
		# tkconfigure(chip.port3, variable = chipValue, state = "normal")
		# #Illumina 550k
		# tkconfigure(chip.port4, variable = chipValue, state = "normal")
		# tclvalue(chipValue) = 1
	# }
})

#Input data format
inputdata_format.statsValue <- tclVar("2") #default is AF-based
tkconfigure(inputdata_format.port1, variable = inputdata_format.statsValue, value = "1") 
tkconfigure(inputdata_format.port2, variable = inputdata_format.statsValue, value = "2")
inputdata_format.stats.Label <- tklabel(frame.para1A, bg = color, text = "           Input data format:         ", font = fontTextLabel, width = 23)
inputdata_format.LeftLabel <- tklabel(frame.para1A, bg = color, text = "Genotype/Intensity-based", font = fontTextLabel)
inputdata_format.RightLabel <- tklabel(frame.para1A, bg = color, text = "AF-based", font = fontTextLabel)
blank.space <- tklabel(frame.para1A, bg = color, text = "  ", font = fontTextLabel, width = 50)

tkgrid(inputdata_format.stats.Label, inputdata_format.port1, inputdata_format.LeftLabel, inputdata_format.port2, inputdata_format.RightLabel, blank.space, sticky = "w")
tkgrid.configure(frame.para1A, sticky = "w")

#Directory of data input
pathlab <- tklabel(frame.path, text = "     Directory of data input:", font = fontTextLabel, bg = color, width = 23)
pathway_input <- tclVar("Test example 2")
entry.Path <- tkentry(frame.path, width = "85", textvariable = pathway_input, font = fontTextLabel, xscrollcommand = function(...){
	if(tclvalue(pathway_input) == "Test example 1"){
		tclvalue(chipValue) = 2 # set to 500K
		tclvalue(groupValue) = 5 #set to combine
		tclvalue(chip1.val) = 0
		tclvalue(chip2.val) = 1
		tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(chip2.port, state = "normal")
		tkconfigure(hind.port, variable = hind.val, state = "disable")
		tkconfigure(xba.port, variable = xba.val, state = "disable")
		tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		tkconfigure(sty.port, variable = sty.val, state = "disable")
		tkconfigure(group.pp1, variable = groupValue, state = "normal")
		tkconfigure(group.pp2, variable = groupValue, state = "normal")
		tkconfigure(group.pp3, variable = groupValue, state = "normal")
		tkconfigure(group.pp5, variable = groupValue, state = "normal")
		#Statistical Analysis
		##CPA calculation
		tkconfigure(CPA.stats.port1, variable = CPA.statsValue, state = "normal")
		tkconfigure(CPA.stats.port2, variable = CPA.statsValue, state = "normal")
		tclvalue(CPA.statsValue) = 1
		##AF calculation
		tkconfigure(AF.stats.port1, variable = AF.statsValue, state = "normal")
		tkconfigure(AF.stats.port2, variable = AF.statsValue, state = "normal")
		tclvalue(AF.statsValue) = 2
		##AF reference
		tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "normal")
		tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "normal")
		tclvalue(AF_ref.Value) = 1
		tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
		tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
		tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
		##QI calculation
		tkconfigure(QI.port1, variable = QIValue, state = "normal")
		tkconfigure(QI.port2, variable = QIValue, state = "normal")
		tclvalue(QIValue) = 2
		##Identification of poor array
		tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
		tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
		tclvalue(QI_iden.statsValue) = 2
		##-----------------------------
		tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
		tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
		tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
		tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 1
		
		#Graphical result
		##AF plot
		tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "normal")
		tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "normal")
		tclvalue(AF_plot.statsValue) = 2
		tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
		tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
		tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
		##GCR plot
			##	No matter which input data format, GCR plot can draw or not 
			##	So the radiobutton would not reset
		##HeatMap plot
		tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
		tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
		tclvalue(HM_plot.statsValue) = 2
		tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "normal")
		tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "normal")
		tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
		##Polygon plot
		tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
		tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
		tclvalue(PG_plot.statsValue) = 2
		
		#Numerical output
		 ##Data description
			##	No matter which input data format, Data description can output or not 
			##	So the radiobutton would not reset
		##CPA estimate
		# tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "normal")
		# tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "normal")
		# tclvalue(CPA_est_plot.statsValue) = 2
		tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "disable")
		tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "disable")
		tclvalue(CPA_est_plot.statsValue) = 0 # could be any number rather than 1/2
		##AF estimate
		tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue, state = "normal")
		tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue, state = "normal")
		tclvalue(AF_numeric.statsValue) = 2
		###Intensity-based & Genotype-based
		tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "normal")
		tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "normal")
		tclvalue(AF_numeric_type1.val) = 1; tclvalue(AF_numeric_type2.val) = 1
		
		##QI estimate
		tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
		tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
		tclvalue(QI_numeric.statsValue) = 2
		##PoorSNP array
		tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
		tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
		tclvalue(poor_numeric.statsValue) = 2

	}else if(tclvalue(pathway_input) == "Test example 2"){
		tclvalue(chipValue) = 1 #set to 100K
		tclvalue(groupValue) = 5 #set to combine
		tclvalue(chip1.val) = 0
		tclvalue(chip2.val) = 1
		tkconfigure(chip2.port, state = "normal")
		tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(hind.port, variable = hind.val, state = "disable")
		tkconfigure(xba.port, variable = xba.val, state = "disable")
		tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		tkconfigure(sty.port, variable = sty.val, state = "disable")
		tkconfigure(group.pp1, variable = groupValue, state = "normal")
		tkconfigure(group.pp2, variable = groupValue, state = "normal")
		tkconfigure(group.pp3, variable = groupValue, state = "normal")
		tkconfigure(group.pp5, variable = groupValue, state = "normal")
		#Statistical Analysis
		##CPA calculation
		tkconfigure(CPA.stats.port1, variable = CPA.statsValue, state = "disable")
		tkconfigure(CPA.stats.port2, variable = CPA.statsValue, state = "disable")
		tclvalue(CPA.statsValue) = 3
		##AF calculation
		tkconfigure(AF.stats.port1, variable = CPA.statsValue, state = "disable")
		tkconfigure(AF.stats.port2, variable = CPA.statsValue, state = "disable")
		tclvalue(AF.statsValue) = 3
		##AF reference
		tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "normal")
		tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "normal")
		tclvalue(AF_ref.Value) = 1
		tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
		tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
		tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
		##QI calculation
		tkconfigure(QI.port1, variable = QIValue, state = "normal")
		tkconfigure(QI.port2, variable = QIValue, state = "normal")
		tclvalue(QIValue) = 2
		##Identification of poor array
		tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
		tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
		tclvalue(QI_iden.statsValue) = 2
		##-----------------------------
		tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
		tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
		tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
		tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 1
		
		#Graphical result
		##AF plot
		tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "normal")
		tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "normal")
		tclvalue(AF_plot.statsValue) = 2
		tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
		tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
		tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
		##HeatMap plot
		tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
		tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
		tclvalue(HM_plot.statsValue) = 2
		tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "normal")
		tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "normal")
		tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
		##Polygon plot
		tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
		tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
		tclvalue(PG_plot.statsValue) = 2
	   
		#Numerical Result
		##CPA estimate
		tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "disable")
		tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "disable")
		tclvalue(CPA_est_plot.statsValue) = 3
		#AF estimate
		tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue, state = "disable")
		tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue, state = "disable")
		tclvalue(AF_numeric.statsValue) = 3
		tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "disable")
		tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "disable")
		tclvalue(AF_numeric_type1.val) = 3; tclvalue(AF_numeric_type2.val) = 3
		##QI estimate
		tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
		tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
		tclvalue(QI_numeric.statsValue) = 2
		##PoorSNP array
		tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
		tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
		tclvalue(poor_numeric.statsValue) = 2
	}else if(tclvalue(pathway_input) == "Test example 3"){
		tclvalue(chipValue) = 5 # set to 550K
		tclvalue(groupValue) = 4 #set to TWN
		tclvalue(chip1.val) = 1
		tclvalue(chip2.val) = 0
		tkconfigure(chip2.port, state = "disable")
		tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(hind.port, variable = hind.val, state = "disable")
		tkconfigure(xba.port, variable = xba.val, state = "disable")
		tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		tkconfigure(sty.port, variable = sty.val, state = "disable")
		tkconfigure(group.pp1, variable = groupValue, state = "normal")
		tkconfigure(group.pp2, variable = groupValue, state = "normal")
		tkconfigure(group.pp3, variable = groupValue, state = "normal")
		tkconfigure(group.pp5, variable = groupValue, state = "normal")
		#Statistical Analysis
		##CPA calculation
		tkconfigure(CPA.stats.port1, variable = CPA.statsValue, state = "normal")
		tkconfigure(CPA.stats.port2, variable = CPA.statsValue, state = "normal")
		tclvalue(CPA.statsValue) = 1
		##AF calculation
		tkconfigure(AF.stats.port1, variable = AF.statsValue, state = "normal")
		tkconfigure(AF.stats.port2, variable = AF.statsValue, state = "normal")
		tclvalue(AF.statsValue) = 2
		##AF reference
		tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "normal")
		tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "normal")
		tclvalue(AF_ref.Value) = 1
		tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
		tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
		tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
		##QI calculation
		tkconfigure(QI.port1, variable = QIValue, state = "normal")
		tkconfigure(QI.port2, variable = QIValue, state = "normal")
		tclvalue(QIValue) = 2
		##Identification of poor array
		tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
		tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
		tclvalue(QI_iden.statsValue) = 2
		##-----------------------------
		tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
		tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
		tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
		tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 1
		
		#Graphical result
		##AF plot
		tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "normal")
		tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "normal")
		tclvalue(AF_plot.statsValue) = 2
		tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
		tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
		tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
		##GCR plot
			##	No matter which input data format, GCR plot can draw or not 
			##	So the radiobutton would not reset
		##HeatMap plot
		tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
		tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
		tclvalue(HM_plot.statsValue) = 2
		tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "normal")
		tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "normal")
		tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
		##Polygon plot
		tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
		tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
		tclvalue(PG_plot.statsValue) = 2
		
		#Numerical output
		 ##Data description
			##	No matter which input data format, Data description can output or not 
			##	So the radiobutton would not reset
		##CPA estimate
		# tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "normal")
		# tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "normal")
		# tclvalue(CPA_est_plot.statsValue) = 2
		tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "disable")
		tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "disable")
		tclvalue(CPA_est_plot.statsValue) = 0 # could be any number rather than 1/2
		##AF estimate
		tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue, state = "normal")
		tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue, state = "normal")
		tclvalue(AF_numeric.statsValue) = 2
		###Intensity-based & Genotype-based
		tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "normal")
		tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "normal")
		tclvalue(AF_numeric_type1.val) = 1; tclvalue(AF_numeric_type2.val) = 1
		
		##QI estimate
		tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
		tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
		tclvalue(QI_numeric.statsValue) = 2
		##PoorSNP array
		tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
		tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
		tclvalue(poor_numeric.statsValue) = 2
	}
})
box.input <- tkbutton(frame.path, text = "...",  command = function() tclvalue(pathway_input) <- tkchooseDirectory())
tkgrid(pathlab, entry.Path, box.input)
tkgrid.configure(frame.path, sticky = "w")

#Directory of result output
pathlab <- tklabel(frame.path, text = "         Directory of result output:", font = fontTextLabel, bg = color, width = 26)
pathway_output <- tclVar("")
entry.Path <- tkentry(frame.path, width = "85", textvariable = pathway_output, font = fontTextLabel)
box.output <- tkbutton(frame.path, text = "...",  command = function() tclvalue(pathway_output) <- tkchooseDirectory())
tkgrid(pathlab, entry.Path, box.output)
tkgrid.configure(frame.path, sticky = "w")

##-----------------------------------------------------------------------------
								##2. Data format
##-----------------------------------------------------------------------------
tkgrid(tklabel(frame.dataformat, text = "   2. Data format:", font = fontsubtitle, bg = color), sticky = "w") # title
#-------------------------------# 2.1  study population  #-----------------------#
group.Left <- tklabel(frame.dataformat, text = "            Study population:  ", font = fontTextLabel, bg = color)

group.pp1 <- tkradiobutton(frame.dataformat, bg = color)
group.pp2 <- tkradiobutton(frame.dataformat, bg = color)
group.pp3 <- tkradiobutton(frame.dataformat, bg = color)
group.pp4 <- tkradiobutton(frame.dataformat, bg = color)
group.pp5 <- tkradiobutton(frame.dataformat, bg = color)
groupValue <- tclVar("5") #default study population is Combine
tkconfigure(group.pp1, variable = groupValue, value = "1")
tkconfigure(group.pp2, variable = groupValue, value = "2")
tkconfigure(group.pp3, variable = groupValue, value = "3")
tkconfigure(group.pp4, variable = groupValue, value = "4")
tkconfigure(group.pp5, variable = groupValue, value = "5")

group1 <- tklabel(frame.dataformat, text = "Asia (CHB+JPT) ", width = 13, font = fontTextLabel, bg = color)
group2 <- tklabel(frame.dataformat, text = "Africa (YRI)", width = 10, font = fontTextLabel, bg = color)
group3 <- tklabel(frame.dataformat, text = "Europe (CEU)", width = 10, font = fontTextLabel, bg = color)
group4 <- tklabel(frame.dataformat, text = "Taiwan (TWN)", width = 12, font = fontTextLabel, bg = color)
group5 <- tklabel(frame.dataformat, text = "Combine", width = 6, font = fontTextLabel, bg = color)

tkgrid(group.Left, group.pp1, group1, group.pp2, group2, group.pp3, group3, group.pp4, group4, group.pp5, group5, sticky = "w")

#-------------------------------# 2.2 Genome-wide SNP array  #-----------------------#
chip.lab <- tklabel(frame.dataformat, text = "            Genome-wide SNP array:", font = fontTextLabel, bg = color)

##Affy 100K
chip.port1 <- tkradiobutton(frame.dataformat, bg = color, command = function(){
	tclvalue(chip1.val) = 0; tclvalue(chip2.val) = 1 #two array analysis
	tkconfigure(group.pp1, variable = groupValue, state = "normal")
	tkconfigure(group.pp2, variable = groupValue, state = "normal")
	tkconfigure(group.pp3, variable = groupValue, state = "normal")
	tkconfigure(group.pp5, variable = groupValue, state = "normal")
	#tkconfigure(chip1.port, variable = chip1.val, state = "normal")
	#tkconfigure(chip2.port, variable = chip2.val, state = "normal")
	tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
	tkconfigure(hind.port, variable = hind.val, state = "disable")
	tkconfigure(xba.port, variable = xba.val, state = "disable")
	tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
	tkconfigure(nsp.port, variable = nsp.val, state = "disable")
	tkconfigure(sty.port, variable = sty.val, state = "disable")
	if(c(tclvalue(inputdata_format.statsValue))=="1") source(SAQC.QIpara)
})
##Affy 500K
chip.port2 <- tkradiobutton(frame.dataformat, bg = color, command = function(){
	tclvalue(chip1.val) = 0; tclvalue(chip2.val) = 1
	tkconfigure(group.pp1, variable = groupValue, state = "normal")
	tkconfigure(group.pp2, variable = groupValue, state = "normal")
	tkconfigure(group.pp3, variable = groupValue, state = "normal")
	tkconfigure(group.pp5, variable = groupValue, state = "normal")
	#tkconfigure(chip1.port, variable = chip1.val, state = "normal")
	#tkconfigure(chip2.port, variable = chip2.val, state = "normal")
	tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
	tkconfigure(hind.port, variable = hind.val, state = "disable")
	tkconfigure(xba.port, variable = xba.val, state = "disable")
	tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
	tkconfigure(nsp.port, variable = nsp.val, state = "disable")
	tkconfigure(sty.port, variable = sty.val, state = "disable")
	if(c(tclvalue(inputdata_format.statsValue))=="1") source(SAQC.QIpara)
})

##Affy 6.0
chip.port3 <- tkradiobutton(frame.dataformat, bg = color, command = function(){
	tclvalue(chip1.val) = 1;tclvalue(chip2.val) = 0
	tkconfigure(group.pp1, variable = groupValue, state = "normal")
	tkconfigure(group.pp2, variable = groupValue, state = "normal")
	tkconfigure(group.pp3, variable = groupValue, state = "normal")
	tkconfigure(group.pp5, variable = groupValue, state = "normal")
	tkconfigure(chip1.port, variable = chip1.val, state = "normal")
	tkconfigure(chip2.port, variable = chip2.val, state = "disable")
	tclvalue(hind.val) = 0;tclvalue(xba.val) = 0
	tkconfigure(hind.port, variable = hind.val, state = "disable")
	tkconfigure(xba.port, variable = xba.val, state = "disable")
	tclvalue(nsp.val) = 0;tclvalue(sty.val) = 0
	tkconfigure(nsp.port, variable = nsp.val, state = "disable")
	tkconfigure(sty.port, variable = sty.val, state = "disable")
	if(c(tclvalue(inputdata_format.statsValue))=="1") source(SAQC.QIpara)
})

##Affy Axiom
chip.port4 <- tkradiobutton(frame.dataformat, bg = color, command = function(){
	tclvalue(chip1.val) = 1;tclvalue(chip2.val) = 0
	tkconfigure(group.pp1, variable = groupValue, state = "disable")
	tkconfigure(group.pp2, variable = groupValue, state = "disable")
	tkconfigure(group.pp3, variable = groupValue, state = "disable")
	tkconfigure(group.pp5, variable = groupValue, state = "disable")
	tclvalue(groupValue) = 4
	tkconfigure(chip1.port, variable = chip1.val, state = "normal")
	tkconfigure(chip2.port, variable = chip2.val, state = "disable")
	tclvalue(hind.val) = 0;tclvalue(xba.val) = 0
	tkconfigure(hind.port, variable = hind.val, state = "disable")
	tkconfigure(xba.port, variable = xba.val, state = "disable")
	tclvalue(nsp.val) = 0;tclvalue(sty.val) = 0
	tkconfigure(nsp.port, variable = nsp.val, state = "disable")
	tkconfigure(sty.port, variable = sty.val, state = "disable")
	if(c(tclvalue(inputdata_format.statsValue))=="1") source(SAQC.QIpara)
})

##Illumina
chip.port5 <- tkradiobutton(frame.dataformat, bg = color, command = function(){
	tclvalue(chip1.val) = 1;tclvalue(chip2.val) = 0
	tkconfigure(group.pp1, variable = groupValue, state = "disable")
	tkconfigure(group.pp2, variable = groupValue, state = "disable")
	tkconfigure(group.pp3, variable = groupValue, state = "disable")
	tkconfigure(group.pp5, variable = groupValue, state = "disable")
	tclvalue(groupValue) = 4
	tkconfigure(chip1.port, variable = chip1.val, state = "normal")
	tkconfigure(chip2.port, variable = chip2.val, state = "disable")
	tclvalue(hind.val) = 0;tclvalue(xba.val) = 0
	tkconfigure(hind.port, variable = hind.val, state = "disable")
	tkconfigure(xba.port, variable = xba.val, state = "disable")
	tclvalue(nsp.val) = 0;tclvalue(sty.val) = 0
	tkconfigure(nsp.port, variable = nsp.val, state = "disable")
	tkconfigure(sty.port, variable = sty.val, state = "disable")
	if(c(tclvalue(inputdata_format.statsValue))=="1") source(SAQC.QIpara)
})


chipValue <- tclVar("1") #default 100k  
tkconfigure(chip.port1, variable = chipValue, value = "1") #Affy 100k
tkconfigure(chip.port2, variable = chipValue, value = "2") #Affy 500k
tkconfigure(chip.port3, variable = chipValue, value = "3") #Affy 6.0
tkconfigure(chip.port4, variable = chipValue, value = "4") #Affy axiom
tkconfigure(chip.port5, variable = chipValue, value = "5") #illumina 550k


chip.affy100 <- tklabel(frame.dataformat, text = " Affymetrix 100K", width = 12, font = fontTextLabel, bg = color)
chip.affy500 <- tklabel(frame.dataformat, text = " Affymetrix 500K", width = 13, font = fontTextLabel, bg = color)
chip.affy6 <- tklabel(frame.dataformat, text = "Affymetrix 6.0", width = 10, font = fontTextLabel, bg = color)
chip.affyaxiom = tklabel(frame.dataformat, text = "Affymetrix Axiom", width = 14, font = fontTextLabel, bg = color)
chip.illumina550 <- tklabel(frame.dataformat, text = "Illumina 550K", width = 10, font = fontTextLabel, bg = color)


#tkgrid.configure(frame.dataformat, sticky = "w")
tkgrid(chip.lab, chip.port1, chip.affy100, chip.port2, chip.affy500, chip.port3, chip.affy6, chip.port4, chip.affyaxiom, chip.port5, chip.illumina550, sticky = "w")

#-------------------------------# 2.3 Array selection  #-----------------------#
#one array analysis
chip1.port <- tkcheckbutton(frame.dataformat, bg = color, command = function(){
	tclvalue(chip1.val) = 1; tclvalue(chip2.val) = 0
	if(tclvalue(chipValue) == 1){ # 100k
		#tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		tkconfigure(hind.port, variable = hind.val, state = "normal")
		tkconfigure(xba.port, variable = xba.val, state = "normal")
		#tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		#tkconfigure(nsp.port, variable = nsp.val, state = "disable")
		#tkconfigure(sty.port, variable = sty.val, state = "disable")
	}
	if(tclvalue(chipValue) == 2){ # 500k
		#tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
		#tkconfigure(hind.port, variable = hind.val, state = "disable")
		#tkconfigure(xba.port, variable = xba.val, state = "disable")
		#tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
		tkconfigure(nsp.port, variable = nsp.val, state = "normal")
		tkconfigure(sty.port, variable = sty.val, state = "normal")
	}
})

#two array analysis
chip2.port <- tkcheckbutton(frame.dataformat, bg = color, command = function()
{
	tclvalue(chip1.val) = 0; tclvalue(chip2.val) = 1
	tclvalue(hind.val) = 0; tclvalue(xba.val) = 0
	tkconfigure(hind.port, variable = hind.val, state = "disable")
	tkconfigure(xba.port, variable = xba.val, state = "disable")
	tclvalue(nsp.val) = 0; tclvalue(sty.val) = 0
	tkconfigure(nsp.port, variable = nsp.val, state = "disable")
	tkconfigure(sty.port, variable = sty.val, state = "disable")
})

##
##	default two array analysis
##
chip1.val <- tclVar("0"); chip2.val <- tclVar("1") 
tkconfigure(chip1.port, variable = chip1.val) #one array analysis
tkconfigure(chip2.port, variable = chip2.val) #two array analysis

hind.port <- tkradiobutton(frame.dataformat, bg = color, state = "disable", command = function() {tclvalue(hind.val) = 1;tclvalue(xba.val) = 0})
xba.port <- tkradiobutton(frame.dataformat, bg = color, state = "disable", command = function() {tclvalue(xba.val) = 1;tclvalue(hind.val) = 0})
nsp.port <- tkradiobutton(frame.dataformat, bg = color, state = "disable", command = function() {tclvalue(nsp.val) = 1;tclvalue(sty.val) = 0})
sty.port <- tkradiobutton(frame.dataformat, bg = color, state = "disable", command = function() {tclvalue(sty.val) = 1;tclvalue(nsp.val) = 0})

hind.val <- tclVar("0"); xba.val <- tclVar("0"); nsp.val <- tclVar("0"); sty.val <- tclVar("0")

tkconfigure(hind.port, variable = hind.val, value = "1")
tkconfigure(xba.port, variable = xba.val, value = "1")
tkconfigure(nsp.port, variable = nsp.val, value = "1")
tkconfigure(sty.port, variable = sty.val, value = "1")

chip.Label <- tklabel(frame.dataformat, text = "            Array selection:", font = fontTextLabel, bg = color)
chip1.Label <- tklabel(frame.dataformat, bg = color, text = "One array analysis  (", width = 16, font = fontTextLabel)
hind.Label <- tklabel(frame.dataformat, bg = color, text = "Hind", width = 5, font = fontTextLabel)
xba.Label <- tklabel(frame.dataformat, bg = color, text = "Xba", font = fontTextLabel)
nsp.Label <- tklabel(frame.dataformat, bg = color, text = "Nsp", font = fontTextLabel)
sty.Label <- tklabel(frame.dataformat, bg = color, text = "Sty )", font = fontTextLabel)
tkgrid(chip.Label, chip1.port, chip1.Label, hind.port, hind.Label, xba.port, xba.Label, nsp.port, nsp.Label, sty.port, sty.Label, sticky = "w")

chip.Label.second <- tklabel(frame.dataformat, text = "                    ", bg = color)
chip2.Label <- tklabel(frame.dataformat, bg = color, text = "Two array analysis", font = fontTextLabel)
tkgrid(chip.Label.second, chip2.port, chip2.Label, sticky = "w")

#tkgrid.configure(chip1.Label, sticky = "w")
#tkgrid.configure(chip2.Label, sticky = "w")

##==============================================================================

##-------------------------  parameter setting  --------------------------------
##      it have been removed
##----- Chromosome -----
if(0){
	chr1.lab <- tklabel(frame.para,bg=color,text="1")
	chr2.lab <- tklabel(frame.para,bg=color,text="2")
	chr3.lab <- tklabel(frame.para,bg=color,text="3")
	chr4.lab <- tklabel(frame.para,bg=color,text="4")
	chr5.lab <- tklabel(frame.para,bg=color,text="5")
	chr6.lab <- tklabel(frame.para,bg=color,text="6")
	chr7.lab <- tklabel(frame.para,bg=color,text="7")
	chr8.lab <- tklabel(frame.para,bg=color,text="8")
	chr9.lab <- tklabel(frame.para,bg=color,text="9")
	chr10.lab <- tklabel(frame.para,bg=color,text="10")
	chr11.lab <- tklabel(frame.para,bg=color,text="11")
	chr12.lab <- tklabel(frame.para,bg=color,text="12")
	chr13.lab <- tklabel(frame.para,bg=color,text="13")
	chr14.lab <- tklabel(frame.para,bg=color,text="14")
	chr15.lab <- tklabel(frame.para,bg=color,text="15")
	chr16.lab <- tklabel(frame.para,bg=color,text="16")
	chr17.lab <- tklabel(frame.para,bg=color,text="17")
	chr18.lab <- tklabel(frame.para,bg=color,text="18")
	chr19.lab <- tklabel(frame.para,bg=color,text="19")
	chr20.lab <- tklabel(frame.para,bg=color,text="20")
	chr21.lab <- tklabel(frame.para,bg=color,text="21")
	chr22.lab <- tklabel(frame.para,bg=color,text="22")
	chr23.lab <- tklabel(frame.para,bg=color,text="23")
	chrall.lab <- tklabel(frame.para,bg=color,text="All")

	chr1.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr1.val)=="0"){tclvalue(chrall.val)=0}})
	chr2.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr2.val)=="0"){tclvalue(chrall.val)=0}})
	chr3.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr3.val)=="0"){tclvalue(chrall.val)=0}})
	chr4.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr4.val)=="0"){tclvalue(chrall.val)=0}})
	chr5.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr5.val)=="0"){tclvalue(chrall.val)=0}})
	chr6.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr6.val)=="0"){tclvalue(chrall.val)=0}})
	chr7.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr7.val)=="0"){tclvalue(chrall.val)=0}})
	chr8.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr8.val)=="0"){tclvalue(chrall.val)=0}})
	chr9.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr9.val)=="0"){tclvalue(chrall.val)=0}})
	chr10.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr10.val)=="0"){tclvalue(chrall.val)=0}})
	chr11.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr11.val)=="0"){tclvalue(chrall.val)=0}})
	chr12.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr12.val)=="0"){tclvalue(chrall.val)=0}})
	chr13.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr13.val)=="0"){tclvalue(chrall.val)=0}})
	chr14.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr14.val)=="0"){tclvalue(chrall.val)=0}})
	chr15.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr15.val)=="0"){tclvalue(chrall.val)=0}})
	chr16.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr16.val)=="0"){tclvalue(chrall.val)=0}})
	chr17.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr17.val)=="0"){tclvalue(chrall.val)=0}})
	chr18.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr18.val)=="0"){tclvalue(chrall.val)=0}})
	chr19.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr19.val)=="0"){tclvalue(chrall.val)=0}})
	chr20.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr20.val)=="0"){tclvalue(chrall.val)=0}})
	chr21.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr21.val)=="0"){tclvalue(chrall.val)=0}})
	chr22.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr22.val)=="0"){tclvalue(chrall.val)=0}})
	chr23.port <- tkcheckbutton(frame.para,bg=color,command=function(){if(tclvalue(chr23.val)=="0"){tclvalue(chrall.val)=0}})
	chrall.port <- tkcheckbutton(frame.para,bg=color,command=function(){
		if(tclvalue(chrall.val)=="1"){
			tclvalue(chr1.val)=1;tclvalue(chr2.val)=1;tclvalue(chr3.val)=1;
			tclvalue(chr4.val)=1;tclvalue(chr5.val)=1;tclvalue(chr6.val)=1;
			tclvalue(chr7.val)=1;tclvalue(chr8.val)=1;tclvalue(chr9.val)=1;
			tclvalue(chr10.val)=1;tclvalue(chr11.val)=1;tclvalue(chr12.val)=1;
			tclvalue(chr13.val)=1;tclvalue(chr14.val)=1;tclvalue(chr15.val)=1;
			tclvalue(chr16.val)=1;tclvalue(chr17.val)=1;tclvalue(chr18.val)=1;
			tclvalue(chr19.val)=1;tclvalue(chr20.val)=1;tclvalue(chr21.val)=1;
			tclvalue(chr22.val)=1;tclvalue(chr23.val)=1;
		}else{
			tclvalue(chr1.val)=0;tclvalue(chr2.val)=0;tclvalue(chr3.val)=0;
			tclvalue(chr4.val)=0;tclvalue(chr5.val)=0;tclvalue(chr6.val)=0;
			tclvalue(chr7.val)=0;tclvalue(chr8.val)=0;tclvalue(chr9.val)=0;
			tclvalue(chr10.val)=0;tclvalue(chr11.val)=0;tclvalue(chr12.val)=0;
			tclvalue(chr13.val)=0;tclvalue(chr14.val)=0;tclvalue(chr15.val)=0;
			tclvalue(chr16.val)=0;tclvalue(chr17.val)=0;tclvalue(chr18.val)=0;
			tclvalue(chr19.val)=0;tclvalue(chr20.val)=0;tclvalue(chr21.val)=0;
			tclvalue(chr22.val)=0;tclvalue(chr23.val)=0;
		}
	})

	chr1.val <- tclVar("1")
	chr2.val <- tclVar("1")
	chr3.val <- tclVar("1")
	chr4.val <- tclVar("1")
	chr5.val <- tclVar("1")
	chr6.val <- tclVar("1")
	chr7.val <- tclVar("1")
	chr8.val <- tclVar("1")
	chr9.val <- tclVar("1")
	chr10.val <- tclVar("1")
	chr11.val <- tclVar("1")
	chr12.val <- tclVar("1")
	chr13.val <- tclVar("1")
	chr14.val <- tclVar("1")
	chr15.val <- tclVar("1")
	chr16.val <- tclVar("1")
	chr17.val <- tclVar("1")
	chr18.val <- tclVar("1")
	chr19.val <- tclVar("1")
	chr20.val <- tclVar("1")
	chr21.val <- tclVar("1")
	chr22.val <- tclVar("1")
	chr23.val <- tclVar("1")
	chrall.val <- tclVar("1")

	tkconfigure(chr1.port,variable=chr1.val)
	tkconfigure(chr2.port,variable=chr2.val)
	tkconfigure(chr3.port,variable=chr3.val)
	tkconfigure(chr4.port,variable=chr4.val)
	tkconfigure(chr5.port,variable=chr5.val)
	tkconfigure(chr6.port,variable=chr6.val)
	tkconfigure(chr7.port,variable=chr7.val)
	tkconfigure(chr8.port,variable=chr8.val)
	tkconfigure(chr9.port,variable=chr9.val)
	tkconfigure(chr10.port,variable=chr10.val)
	tkconfigure(chr11.port,variable=chr11.val)
	tkconfigure(chr12.port,variable=chr12.val)
	tkconfigure(chr13.port,variable=chr13.val)
	tkconfigure(chr14.port,variable=chr14.val)
	tkconfigure(chr15.port,variable=chr15.val)
	tkconfigure(chr16.port,variable=chr16.val)
	tkconfigure(chr17.port,variable=chr17.val)
	tkconfigure(chr18.port,variable=chr18.val)
	tkconfigure(chr19.port,variable=chr19.val)
	tkconfigure(chr20.port,variable=chr20.val)
	tkconfigure(chr21.port,variable=chr21.val)
	tkconfigure(chr22.port,variable=chr22.val)
	tkconfigure(chr23.port,variable=chr23.val)
	tkconfigure(chrall.port,variable=chrall.val)
}
#tkgrid(chromosome.lab)
#tkgrid(tklabel(frame.para,bg=color,text="            Chromosome: ",font=fontTextLabel),chr1.port,chr1.lab,chr2.port,chr2.lab,chr3.port,chr3.lab,chr4.port,chr4.lab,chr5.port,chr5.lab,chr6.port,chr6.lab,chr7.port,chr7.lab,chr8.port,chr8.lab,chr9.port,chr9.lab,chr10.port,chr10.lab,chr11.port,chr11.lab,chr12.port,chr12.lab)
#tkgrid(tklabel(frame.para,bg=color,text=""),chr13.port,chr13.lab,chr14.port,chr14.lab,chr15.port,chr15.lab,chr16.port,chr16.lab,chr17.port,chr17.lab,chr18.port,chr18.lab,chr19.port,chr19.lab,chr20.port,chr20.lab,chr21.port,chr21.lab,chr22.port,chr22.lab,chr23.port,chr23.lab,chrall.port,chrall.lab)
##==============================================================================

##-----------------------------------------------------------------------------
								##3. Statistical Analysis
##-----------------------------------------------------------------------------
tkgrid(tklabel(frame.SA, text = "   3. Statistical analysis:", font = fontsubtitle, bg = color), sticky = "w") # title
#-------------------------------# 3-1 CPA calculation  #-----------------------#
CPA.stats.port1 <- tkradiobutton(frame.para3A, bg = color, state = "disable", command = function(){ #default disable
	tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "disable")
	tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "disable")
	tclvalue(CPA_est_plot.statsValue) = 0
}) 
CPA.stats.port2 <- tkradiobutton(frame.para3A, bg = color, state = "disable", command = function(){
	tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, state = "normal")
	tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, state = "normal")
	tclvalue(CPA_est_plot.statsValue) = 2
}) #default disable

CPA.statsValue <- tclVar("0") #default disable
tkconfigure(CPA.stats.port1, variable = CPA.statsValue, value = "1")
tkconfigure(CPA.stats.port2, variable = CPA.statsValue, value = "2")

CPA.stats.Label <- tklabel(frame.para3A, bg = color, text = "            CPA calculation:                  ", font = fontTextLabel)
CPA.stats.LeftLabel <- tklabel(frame.para3A, bg = color, text = "No (using CPA in SAQC database)      ", width = 29, font = fontTextLabel)
CPA.stats.RightLabel <- tklabel(frame.para3A, bg = color, text = "Yes", font = fontTextLabel)

tkgrid(CPA.stats.Label, CPA.stats.port1, CPA.stats.LeftLabel, CPA.stats.port2, CPA.stats.RightLabel, sticky = "w")

#------------------------------#  3-2 AF calculation  #----------------------------#
AF.stats.port1 <- tkradiobutton(frame.para3A, bg = color, state = "disable", command = function(){ #default disable
	#Options in 3. Statistical analysis
	#disable AF reference calculation
	tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "disable")
    tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "disable")
    tclvalue(AF_ref.Value) = 0
    ##AF reference source
    tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "disable"); tclvalue(AFref_source1.val) = 0
    tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "disable"); tclvalue(AFref_source2.val) = 0
    ##AF reference source directory path
	tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
    
	#disable QI calculation
	tkconfigure(QI.port1, variable = QIValue, state = "disable")
    tkconfigure(QI.port2, variable = QIValue, state = "disable")
    tclvalue(QIValue) = 0
    
	#disable Identification of poor array
    tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "disable")
    tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "disable")
    tclvalue(QI_iden.statsValue) = 0
    ##QI upper quantile
    tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "disable")
    tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "disable")
    tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "disable")
    tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 0
    
	#Options in 4. Graphical output
	#disable AF plot
    tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "disable")
    tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "disable")
    tclvalue(AF_plot.statsValue) = 0
    tkconfigure(AF_type1.port, variable = AF_type1.val, state = "disable")
    tkconfigure(AF_type2.port, variable = AF_type2.val, state = "disable")
    tclvalue(AF_type1.val) = 0; tclvalue(AF_type2.val) = 0
    #disable QI HeatMap plot
    tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "disable")
    tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "disable")
    tclvalue(HM_plot.statsValue) = 0
    tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "disable")
    tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "disable")
    tclvalue(col_scale_type1.val) = 0; tclvalue(col_scale_type2.val) = 0
    #disable QI Polygon plot
    tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "disable")
    tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "disable")
    tclvalue(PG_plot.statsValue) = 0
    
	#Options in 5. Numerical output
	#disable AF estimate
    tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue , state = "disable")
    tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue , state = "disable")
    tclvalue(AF_numeric.statsValue) = 0
    tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "disable")
    tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "disable")
    tclvalue(AF_numeric_type1.val) = 0; tclvalue(AF_numeric_type2.val) = 0
    #disable QI estimate
    tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue , state = "disable")
    tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue , state = "disable")
    tclvalue(QI_numeric.statsValue) = 0
    #disable Poor SNP array
    tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "disable")
    tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "disable")
    tclvalue(poor_numeric.statsValue) = 0
})
AF.stats.port2 <- tkradiobutton(frame.para3A, bg = color, state = "disable", command = function(){
    #Options in 3. Statistical analysis
	#enable AF reference calculation
	tkconfigure(AF_ref.port1, variable = AF_ref.Value, state = "normal")
    tkconfigure(AF_ref.port2, variable = AF_ref.Value, state = "normal")
    tclvalue(AF_ref.Value) = 1
    ##AF reference source
    tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
    tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
    ##AF reference source directory path
	tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
    
	#enable QI calculation
	tkconfigure(QI.port1, variable = QIValue, state = "normal")
    tkconfigure(QI.port2, variable = QIValue, state = "normal")
    tclvalue(QIValue) = 2
    
	#enable Identification of poor array
    tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
    tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
    tclvalue(QI_iden.statsValue) = 2
    ##QI upper quantile
    tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
    tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
    tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
    tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 1
    
	#Options in 4. Graphical output
	#enable AF plot
    tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, state = "normal")
    tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, state = "normal")
    tclvalue(AF_plot.statsValue) = 2
    tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
    tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
    tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
    
	#enable QI HeatMap plot
    tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
    tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
    tclvalue(HM_plot.statsValue) = 2
    tkconfigure(col_scale_type1.port, variable = col_scale_type1.val , state = "normal")
    tkconfigure(col_scale_type2.port, variable = col_scale_type2.val , state = "normal")
    tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
    #enable QI Polygon plot
    tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
    tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
    tclvalue(PG_plot.statsValue) = 2
    
	#Options in 5. Numerical output
	#enable AF estimate
    tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue , state = "normal")
    tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue , state = "normal")
    tclvalue(AF_numeric.statsValue) = 2
    tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "normal") ##<< it could be removed
    tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "normal") ##<<
    tclvalue(AF_numeric_type1.val) = 1; tclvalue(AF_numeric_type2.val) = 1
    #enable QI estimate
    tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
    tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
    tclvalue(QI_numeric.statsValue) = 2
	#enable Poor SNP array
    tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
    tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
    tclvalue(poor_numeric.statsValue) = 2
})

AF.statsValue <- tclVar("0") #default disable #???????????????????????????????????????????????????????
tkconfigure(AF.stats.port1, variable = AF.statsValue, value = "1")
tkconfigure(AF.stats.port2, variable = AF.statsValue, value = "2")

AF.stats.Label <- tklabel(frame.para3A, bg = color, text = "            AF calculation:", font = fontTextLabel)
AF.stats.LeftLabel <- tklabel(frame.para3A, bg = color, text = "No", width = 2, font = fontTextLabel)
AF.stats.RightLabel <- tklabel(frame.para3A, bg = color, text = "Yes", font = fontTextLabel)
tkgrid(AF.stats.Label, AF.stats.port1, AF.stats.LeftLabel, AF.stats.port2, AF.stats.RightLabel, sticky = "w")

#-------------------------------# 3-3 AF reference calculation#---------------------------#
AF_ref.port1 <- tkradiobutton(frame.para3D, bg = color, command = function(){
	#SAQC database
    tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "normal"); tclvalue(AFref_source1.val) = 1
    #User-provided
	tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "normal"); tclvalue(AFref_source2.val) = 0
	#Directory of AF ref
    tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
})
AF_ref.port2 <- tkradiobutton(frame.para3D, bg = color, command = function(){
    #SAQC database
	tkconfigure(AFref_source1.port, variable = AFref_source1.val, state = "disable"); tclvalue(AFref_source1.val) = 0
    #User-provided
	tkconfigure(AFref_source2.port, variable = AFref_source2.val, state = "disable"); tclvalue(AFref_source2.val) = 0
    #Directory of AF ref
	tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
})

AF_ref.Value <- tclVar("1") #default No
tkconfigure(AF_ref.port1, variable = AF_ref.Value, value = "1")
tkconfigure(AF_ref.port2, variable = AF_ref.Value, value = "2")

AFref_source1.port <- tkcheckbutton(frame.para3D, bg = color, command = function(){
	tclvalue(AFref_source1.val) = 1; tclvalue(AFref_source2.val) = 0;
	#Directory of AF ref
	tkconfigure(entry.AFrefPath, state = "disable"); tclvalue(pathway_AFrefinput) = ""
})

AFref_source2.port <- tkcheckbutton(frame.para3D, bg = color, command = function(){
	tclvalue(AFref_source1.val) = 0; tclvalue(AFref_source2.val) = 1
	#Directory of AF ref
	tkconfigure(entry.AFrefPath, state = "normal"); tclvalue(pathway_AFrefinput) = ""
})

AFref_source1.val <- tclVar("1") #default SAQC database
AFref_source2.val <- tclVar("0")

tkconfigure(AFref_source1.port, variable = AFref_source1.val)
tkconfigure(AFref_source2.port, variable = AFref_source2.val)

AF_ref.Label <- tklabel(frame.para3D, bg = color, text = "            AF reference calculation:   ", font = fontTextLabel)
AF_ref.LeftLabel <- tklabel(frame.para3D, bg = color, text = "No ( Source: ", font = fontTextLabel)
AFref_source1.lab <- tklabel(frame.para3D, bg = color, text = "SAQC database", font = fontTextLabel)
AFref_source2.lab <- tklabel(frame.para3D, bg = color, text = "User-provided )", font = fontTextLabel)
AF_ref.RightLabel <- tklabel(frame.para3D, bg = color, text = "Yes", font = fontTextLabel)

tkgrid(AF_ref.Label, AF_ref.port1, AF_ref.LeftLabel, AFref_source1.port, AFref_source1.lab, AFref_source2.port, AFref_source2.lab, AF_ref.port2, AF_ref.RightLabel, sticky = "w")

#-------------------------------# 3-3-1 User-provided AF reference#---------------------------#
frame.AFrefpath <- tkframe(frame.para3C, bg = color)

AFrefpathlab <- tklabel(frame.AFrefpath, text = "              --- Directory of AF ref. :  ", font = fontTextLabel, bg = color)
pathway_AFrefinput <- tclVar("")
entry.AFrefPath <- tkentry(frame.AFrefpath, width = "85", textvariable = pathway_AFrefinput, font = fontTextLabel, state = "disable")
box.AFrefinput <- tkbutton(frame.AFrefpath, text = "...",  command = function() tclvalue(pathway_AFrefinput) <- tkchooseDirectory())
tkgrid(AFrefpathlab, entry.AFrefPath, box.AFrefinput)
#tkgrid.configure(frame.AFrefpath, sticky = "e")
tkgrid(frame.AFrefpath, sticky = "w")

#-------------------------------# 3-4 QI calculation  #-----------------------#
QI.port1 <- tkradiobutton(frame.para3B, bg = color, command = function(){
	#Options in 3. Statistical analysis
	#disable Identification of poor array
	tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "disable")
	tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "disable")
	tclvalue(QI_iden.statsValue) = 0
	##QI upper quantile
	tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "disable")
	tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "disable")
	tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "disable")
	tclvalue(QI_quantile1.val) = 0; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 0

	#Options in 4. Graphical output
	#disable QI HeatMap plot
	tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "disable")
	tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "disable")
	tclvalue(HM_plot.statsValue) = 0
	tkconfigure(col_scale_type1.port, variable = col_scale_type1.val, state = "disable")
	tkconfigure(col_scale_type2.port, variable = col_scale_type2.val, state = "disable")
	tclvalue(col_scale_type1.val) = 0; tclvalue(col_scale_type2.val) = 0

	#disable QI Polygon plot
	tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "disable")
	tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "disable")
	tclvalue(PG_plot.statsValue) = 0

	#Options in 5. Numerical output
	#disable QI estimate
	tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "disable")
	tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "disable")
	tclvalue(QI_numeric.statsValue) = 0

	#disable Poor SNP array estimate
	tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "disable")
	tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "disable")
	tclvalue(poor_numeric.statsValue) = 0
})
QI.port2 <- tkradiobutton(frame.para3B, bg = color, command = function(){
	#Options in 3. Statistical analysis
	#enable Identification of poor array
	tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, state = "normal")
	tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, state = "normal")
	tclvalue(QI_iden.statsValue) = 2
	##QI upper quantile
	tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal")
	tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal")
	tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal")
	tclvalue(QI_quantile1.val) = 1; tclvalue(QI_quantile2.val) = 0; tclvalue(QI_quantile3.val) = 0

	#Options in 4. Graphical output
	#enable QI HeatMap plot
	tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, state = "normal")
	tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, state = "normal")
	tclvalue(HM_plot.statsValue) = 2

	tkconfigure(col_scale_type1.port, variable = col_scale_type1.val, state = "normal")
	tkconfigure(col_scale_type2.port, variable = col_scale_type2.val, state = "normal")
	tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1

	#enable QI Polygon plot
	tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
	tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
	tclvalue(PG_plot.statsValue) = 2

	#Options in 5. Numerical output
	#enable QI estimate
	tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, state = "normal")
	tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, state = "normal")
	tclvalue(QI_numeric.statsValue) = 2

	#enable Poor SNP array estimate
	tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
	tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
	tclvalue(poor_numeric.statsValue) = 2
})

QIValue <- tclVar("2") #default Yes
tkconfigure(QI.port1, variable = QIValue, value = "1")
tkconfigure(QI.port2, variable = QIValue, value = "2")

QI.Label <- tklabel(frame.para3B, bg = color, text = "            QI calculation: ", font = fontTextLabel)
QI.LeftLabel <- tklabel(frame.para3B, bg = color, text = "No", font = fontTextLabel)
QI.RightLabel <- tklabel(frame.para3B, bg = color, text = "Yes", font = fontTextLabel)
tkgrid(QI.Label, QI.port1, QI.LeftLabel, QI.port2, QI.RightLabel, sticky = "w")

#-------------------------------# 3-5 Identification of poor array #------------------------#
QI_iden.stats.port1 <- tkradiobutton(frame.para3B, bg = color, command = function(){
	#disable QI upper quantile
	tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "disable"); tclvalue(QI_quantile1.val) = 0
	tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "disable"); tclvalue(QI_quantile2.val) = 0
	tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "disable"); tclvalue(QI_quantile3.val) = 0
	
	#disable QI Polygon plot
	tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "disable")
	tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "disable")
	tclvalue(PG_plot.statsValue) = 0
	
	#disable Poor SNP array	
	tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "disable")
	tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "disable")
	tclvalue(poor_numeric.statsValue) = 0
})
QI_iden.stats.port2 <- tkradiobutton(frame.para3B, bg = color, command = function(){
	#enable QI upper quantile
	tkconfigure(QI_quantile1.port, variable = QI_quantile1.val, state = "normal"); tclvalue(QI_quantile1.val) = 1
	tkconfigure(QI_quantile2.port, variable = QI_quantile2.val, state = "normal"); tclvalue(QI_quantile2.val) = 0
	tkconfigure(QI_quantile3.port, variable = QI_quantile3.val, state = "normal"); tclvalue(QI_quantile3.val) = 0
	
	#enable QI Polygon plot
	tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, state = "normal")
	tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, state = "normal")
	tclvalue(PG_plot.statsValue) = 2
	
	#enable Poor SNP array
	tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, state = "normal")
	tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, state = "normal")
	tclvalue(poor_numeric.statsValue) = 2
})

QI_iden.statsValue <- tclVar("2") #default Yes
tkconfigure(QI_iden.stats.port1, variable = QI_iden.statsValue, value = "1")
tkconfigure(QI_iden.stats.port2, variable = QI_iden.statsValue, value = "2")

QI_iden.stats.Label <- tklabel(frame.para3B, bg = color, text = "            Identification of poor array:", font = fontTextLabel)
QI_iden.stats.LeftLabel <- tklabel(frame.para3B, bg = color, text = "No", font = fontTextLabel)
QI_iden.stats.RightLabel <- tklabel(frame.para3B, bg = color, text = "Yes", font = fontTextLabel)

##-------------------------quantile of QI-----------------------------------##
QI_quantile.lab <- tklabel(frame.para3B, bg = color, text = "  (QI upper quantile:", font = fontTextLabel)
QI_quantile1.lab <- tklabel(frame.para3B, bg = color, text = "95%", font = fontTextLabel)
QI_quantile2.lab <- tklabel(frame.para3B, bg = color, text = "97.5%", font = fontTextLabel)
QI_quantile3.lab <- tklabel(frame.para3B, bg = color, text = "99%)", font = fontTextLabel)

QI_quantile1.port <- tkcheckbutton(frame.para3B, bg = color, command = function(){
	tclvalue(QI_quantile1.val) = "1"
	tclvalue(QI_quantile2.val) = 0
	tclvalue(QI_quantile3.val) = 0
})
QI_quantile2.port <- tkcheckbutton(frame.para3B, bg = color, command = function(){
	tclvalue(QI_quantile1.val) = 0
	tclvalue(QI_quantile2.val) = "1"
	tclvalue(QI_quantile3.val) = 0
})
QI_quantile3.port <- tkcheckbutton(frame.para3B, bg = color, command = function(){
	tclvalue(QI_quantile1.val) = 0
	tclvalue(QI_quantile2.val) = 0
	tclvalue(QI_quantile3.val) = "1"
})

QI_quantile1.val <- tclVar("0")
QI_quantile2.val <- tclVar("0")
QI_quantile3.val <- tclVar("1") #default 99%

tkconfigure(QI_quantile1.port, variable = QI_quantile1.val)
tkconfigure(QI_quantile2.port, variable = QI_quantile2.val)
tkconfigure(QI_quantile3.port, variable = QI_quantile3.val)

tkgrid(QI_iden.stats.Label, QI_iden.stats.port1, QI_iden.stats.LeftLabel, QI_iden.stats.port2, QI_iden.stats.RightLabel, 
		QI_quantile.lab, QI_quantile1.port, QI_quantile1.lab, QI_quantile2.port, QI_quantile2.lab, 
		QI_quantile3.port, QI_quantile3.lab, sticky = "w")
#tkfocus(frame.para3B)

##-----------------------------------------------------------------------------
								##4. Graphical output
##-----------------------------------------------------------------------------
tkgrid(tklabel(frame.para4, text = "   4. Graphical output:", font = fontsubtitle, bg = color), sticky = "w") # title
#-------------------------------# 4-1 AF plot  #-----------------------#
AF_plot.port1 <- tkradiobutton(frame.para4B, bg = color, command = function(){
	tkconfigure(AF_type1.port, variable = AF_type1.val, state = "disable")
	tkconfigure(AF_type2.port, variable = AF_type2.val, state = "disable")
	tclvalue(AF_type1.val) = 0; tclvalue(AF_type2.val) = 0
})
AF_plot.port2 <- tkradiobutton(frame.para4B, bg = color, command = function(){
	tkconfigure(AF_type1.port, variable = AF_type1.val, state = "normal")
	tkconfigure(AF_type2.port, variable = AF_type2.val, state = "normal")
	tclvalue(AF_type1.val) = 1; tclvalue(AF_type2.val) = 1
})

AF_plot.statsValue <- tclVar("2") #default Yes
tkconfigure(AF_plot.port1, variable = AF_plot.statsValue, value = "1")
tkconfigure(AF_plot.port2, variable = AF_plot.statsValue, value = "2")

##-------------------------Intensity-based and genotype-based----------------------------##
AF_type1.port <- tkcheckbutton(frame.para4B, bg = color)
AF_type2.port <- tkcheckbutton(frame.para4B, bg = color)

AF_type1.val <- tclVar("1") #default both plot intensity-based & genotype-based AF plots
AF_type2.val <- tclVar("1")

tkconfigure(AF_type1.port, variable = AF_type1.val)
tkconfigure(AF_type2.port, variable = AF_type2.val)

AF_plot.stats.Label <- tklabel(frame.para4B, bg = color, text = "            AF plot:                ", font = fontTextLabel)
AF_plot.stats.LeftLabel <- tklabel(frame.para4B, bg = color, text = "No", font = fontTextLabel)
AF_plot.stats.RightLabel <- tklabel(frame.para4B, bg = color, text = "Yes   ( AF estimate:", font = fontTextLabel)

AF_type1.lab <- tklabel(frame.para4B, bg = color, text = "Intensity-based", font = fontTextLabel)
AF_type2.lab <- tklabel(frame.para4B, bg = color, text = "Genotype-based )", font = fontTextLabel)

tkgrid(AF_plot.stats.Label, AF_plot.port1, AF_plot.stats.LeftLabel, AF_plot.port2, AF_plot.stats.RightLabel, 
		AF_type1.port, AF_type1.lab, AF_type2.port, AF_type2.lab, sticky = "w")

#-------------------------------# 4-2 GCR plot #-----------------------#
CallRate_plot.stats.port1 <- tkradiobutton(frame.para4A, bg = color)
CallRate_plot.stats.port2 <- tkradiobutton(frame.para4A, bg = color)

CallRate_plot.statsValue <- tclVar("2") #default Yes
tkconfigure(CallRate_plot.stats.port1, variable = CallRate_plot.statsValue, value = "1")
tkconfigure(CallRate_plot.stats.port2, variable = CallRate_plot.statsValue, value = "2")

CallRate_plot.stats.Label <- tklabel(frame.para4A, bg = color, text = "            GCR plot:", font = fontTextLabel)
CallRate_plot.stats.LeftLabel <- tklabel(frame.para4A, bg = color, text = "No", font = fontTextLabel)
CallRate_plot.stats.RightLabel <- tklabel(frame.para4A, bg = color, text = "Yes", font = fontTextLabel)

tkgrid(CallRate_plot.stats.Label, CallRate_plot.stats.port1, CallRate_plot.stats.LeftLabel, CallRate_plot.stats.port2, CallRate_plot.stats.RightLabel, sticky = "w")

#-------------------------------# 4-3 HeatMap plot  #-----------------------#
HM_plot.port1 <- tkradiobutton(frame.para4A, bg = color, command = function(){
	tkconfigure(col_scale_type1.port, variable = col_scale_type1.val, state = "disable")
	tkconfigure(col_scale_type2.port, variable = col_scale_type2.val, state = "disable")
	tclvalue(col_scale_type1.val) = 0; tclvalue(col_scale_type2.val) = 0
})
HM_plot.port2 <- tkradiobutton(frame.para4A, bg = color, command = function(){
	tkconfigure(col_scale_type1.port, variable = col_scale_type1.val, state = "normal")
	tkconfigure(col_scale_type2.port, variable = col_scale_type2.val, state = "normal")
	tclvalue(col_scale_type1.val) = 1; tclvalue(col_scale_type2.val) = 1
})

HM_plot.statsValue <- tclVar("2") #default Yes
tkconfigure(HM_plot.port1, variable = HM_plot.statsValue, value = "1")
tkconfigure(HM_plot.port2, variable = HM_plot.statsValue, value = "2")

HM_plot.stats.Label <- tklabel(frame.para4A, bg = color, text = "            QI HeatMap plot:", font = fontTextLabel)
HM_plot.stats.LeftLabel <- tklabel(frame.para4A, bg = color, text = "No", font = fontTextLabel)
HM_plot.stats.RightLabel <- tklabel(frame.para4A, bg = color, text = "Yes   ( Base of color scale:", font = fontTextLabel)

##--------------------------  color scale  -----------------------------------##
col_scale_type1.lab <- tklabel(frame.para4A, bg = color, text = "Reference in SAQC database", font = fontTextLabel)
col_scale_type2.lab <- tklabel(frame.para4A, bg = color, text = "User's data )", font = fontTextLabel)

col_scale_type1.port <- tkcheckbutton(frame.para4A, bg = color)
col_scale_type2.port <- tkcheckbutton(frame.para4A, bg = color)

col_scale_type1.val <- tclVar("1") #default both scales are used
col_scale_type2.val <- tclVar("1")

tkconfigure(col_scale_type1.port, variable = col_scale_type1.val)
tkconfigure(col_scale_type2.port, variable = col_scale_type2.val)

tkgrid(HM_plot.stats.Label, HM_plot.port1, HM_plot.stats.LeftLabel, HM_plot.port2, HM_plot.stats.RightLabel, 
		col_scale_type1.port, col_scale_type1.lab, col_scale_type2.port, col_scale_type2.lab, sticky = "w")

#-------------------------------# 4-4 Polygon plot  #-----------------------#
PG_plot.port1 <- tkradiobutton(frame.para4A, bg = color)
PG_plot.port2 <- tkradiobutton(frame.para4A, bg = color)

PG_plot.statsValue <- tclVar("2") #default Yes
tkconfigure(PG_plot.port1, variable = PG_plot.statsValue, value = "1")
tkconfigure(PG_plot.port2, variable = PG_plot.statsValue, value = "2")

PG_plot.stats.Label <- tklabel(frame.para4A, bg = color, text = "            QI Polygon plot:", font = fontTextLabel)
PG_plot.stats.LeftLabel <- tklabel(frame.para4A, bg = color, text = "No", font = fontTextLabel)
PG_plot.stats.RightLabel <- tklabel(frame.para4A, bg = color, text = "Yes", font = fontTextLabel)

tkgrid(PG_plot.stats.Label, PG_plot.port1, PG_plot.stats.LeftLabel, PG_plot.port2, PG_plot.stats.RightLabel, sticky = "w")

##-----------------------------------------------------------------------------
								##5. Numerical output
##-----------------------------------------------------------------------------
tkgrid(tklabel(frame.para5, text = "   5. Numerical output:", font = fontsubtitle, bg = color), sticky = "w") # title
#-------------------------------# 5-1 Data description #-----------------------#
DD_plot.stats.port1 <- tkradiobutton(frame.para5A, bg = color)
DD_plot.stats.port2 <- tkradiobutton(frame.para5A, bg = color)

DD_plot.statsValue <- tclVar("2") #default Yes
tkconfigure(DD_plot.stats.port1, variable = DD_plot.statsValue, value = "1")
tkconfigure(DD_plot.stats.port2, variable = DD_plot.statsValue, value = "2")

DD_plot.stats.Label <- tklabel(frame.para5A, bg = color, text = "            Data description:", font = fontTextLabel)
DD_plot.stats.LeftLabel <- tklabel(frame.para5A, bg = color, text = "No", font = fontTextLabel)
DD_plot.stats.RightLabel <- tklabel(frame.para5A, bg = color, text = "Yes", font = fontTextLabel)

tkgrid(DD_plot.stats.Label, DD_plot.stats.port1, DD_plot.stats.LeftLabel, DD_plot.stats.port2, DD_plot.stats.RightLabel, sticky = "w")

#-------------------------------# 5-2 CPA estimate #-----------------------#
CPA_est_plot.port1 <- tkradiobutton(frame.para5A, bg = color, state = "disable")
CPA_est_plot.port2 <- tkradiobutton(frame.para5A, bg = color, state = "disable")

CPA_est_plot.statsValue <- tclVar("0") #default disable
tkconfigure(CPA_est_plot.port1, variable = CPA_est_plot.statsValue, value = "1")
tkconfigure(CPA_est_plot.port2, variable = CPA_est_plot.statsValue, value = "2")

CPA_est_plot.stats.Label <- tklabel(frame.para5A, bg = color, text = "            CPA estimate:", font = fontTextLabel)
CPA_est_plot.stats.LeftLabel <- tklabel(frame.para5A, bg = color, text = "No", font = fontTextLabel)
CPA_est_plot.stats.RightLabel <- tklabel(frame.para5A, bg = color, text = "Yes", font = fontTextLabel)

tkgrid(CPA_est_plot.stats.Label, CPA_est_plot.port1, CPA_est_plot.stats.LeftLabel, 
		CPA_est_plot.port2, CPA_est_plot.stats.RightLabel, sticky = "w")

#-------------------------------# 5-3 AF estimate #----------------------------#
AF_numeric.port1 <- tkradiobutton(frame.para5A, bg = color, state = "disable", command = function(){
	tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "disable")
	tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "disable")
	tclvalue(AF_numeric_type1.val) = 0; tclvalue(AF_numeric_type2.val) = 0
})
AF_numeric.port2 <- tkradiobutton(frame.para5A, bg = color, state = "disable", command = function(){
	tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val, state = "normal")
	tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val, state = "normal")
	tclvalue(AF_numeric_type1.val) = 1; tclvalue(AF_numeric_type2.val) = 1

})

AF_numeric.statsValue <- tclVar("0") #default disable
tkconfigure(AF_numeric.port1, variable = AF_numeric.statsValue, value = "1")
tkconfigure(AF_numeric.port2, variable = AF_numeric.statsValue, value = "2")

#-------------------------Intensity and genotype----------------------------
AF_numeric_type1.port <- tkcheckbutton(frame.para5A, bg = color, state = "disable", command = function(){
	if(tclvalue(AF_numeric_type1.val) == 0 & tclvalue(AF_numeric_type2.val) == 0) tclvalue(AF_numeric_type1.val) = 1
})
AF_numeric_type2.port <- tkcheckbutton(frame.para5A, bg = color, state = "disable", command = function(){
	if(tclvalue(AF_numeric_type1.val) == 0 & tclvalue(AF_numeric_type2.val) == 0) tclvalue(AF_numeric_type2.val) = 1
})

AF_numeric_type1.val <- tclVar("0") #default disable
AF_numeric_type2.val <- tclVar("0")

tkconfigure(AF_numeric_type1.port, variable = AF_numeric_type1.val)
tkconfigure(AF_numeric_type2.port, variable = AF_numeric_type2.val)

AF_numeric.stats.Label <- tklabel(frame.para5A, bg = color, text = "            AF estimate:", font = fontTextLabel)
AF_numeric.stats.LeftLabel <- tklabel(frame.para5A, bg = color, text = "No", font = fontTextLabel)
AF_numeric.stats.RightLabel <- tklabel(frame.para5A, bg = color, text = "Yes   ( AF estimate:", font = fontTextLabel)
AF_numeric_type1.lab <- tklabel(frame.para5A, bg = color, text = "Intensity-based", font = fontTextLabel)
AF_numeric_type2.lab <- tklabel(frame.para5A, bg = color, text = "Genotype-based ) ", font = fontTextLabel)

tkgrid(AF_numeric.stats.Label, AF_numeric.port1, AF_numeric.stats.LeftLabel, AF_numeric.port2, AF_numeric.stats.RightLabel,
		AF_numeric_type1.port, AF_numeric_type1.lab, AF_numeric_type2.port, AF_numeric_type2.lab, sticky = "w")

#-------------------------------# 5-4 QI estimate #----------------------------#
QI_numeric.port1 <- tkradiobutton(frame.para5A, bg = color)
QI_numeric.port2 <- tkradiobutton(frame.para5A, bg = color)

QI_numeric.statsValue <- tclVar("2") #default Yes
tkconfigure(QI_numeric.port1, variable = QI_numeric.statsValue, value = "1")
tkconfigure(QI_numeric.port2, variable = QI_numeric.statsValue, value = "2")

QI_numeric.stats.Label <- tklabel(frame.para5A, bg = color, text = "            QI estimate:", font = fontTextLabel)
QI_numeric.stats.LeftLabel <- tklabel(frame.para5A, bg = color, text = "No", font = fontTextLabel)
QI_numeric.stats.RightLabel <- tklabel(frame.para5A, bg = color, text = "Yes", font = fontTextLabel)
tkgrid(QI_numeric.stats.Label, QI_numeric.port1, QI_numeric.stats.LeftLabel, QI_numeric.port2, QI_numeric.stats.RightLabel, sticky = "w")

#-------------------------------# 5-5 Poor SNP array #----------------------------#
poor_numeric.port1 <- tkradiobutton(frame.para5A, bg = color)
poor_numeric.port2 <- tkradiobutton(frame.para5A, bg = color)

poor_numeric.statsValue <- tclVar("2") #default Yes
tkconfigure(poor_numeric.port1, variable = poor_numeric.statsValue, value = "1")
tkconfigure(poor_numeric.port2, variable = poor_numeric.statsValue, value = "2")

poor_numeric.stats.Label <- tklabel(frame.para5A, bg = color, text = "            Poor SNP array:", font = fontTextLabel)
poor_numeric.stats.LeftLabel <- tklabel(frame.para5A, bg = color, text = "No", font = fontTextLabel)
poor_numeric.stats.RightLabel <- tklabel(frame.para5A, bg = color, text = "Yes", font = fontTextLabel)
tkgrid(poor_numeric.stats.Label, poor_numeric.port1, poor_numeric.stats.LeftLabel, 
		poor_numeric.port2, poor_numeric.stats.RightLabel, sticky = "w")
##-----------------------------------------------------------------------------
								##tkgrid all frames
##-----------------------------------------------------------------------------
#tkgrid(frame.dataformat)
tkgrid.configure(frame.dataformat, sticky = "w")  

#tkgrid(frame.para1) #have been tkgrid before
#tkgrid.configure(frame.para1, sticky = "w")
#tkgrid(frame.para1A) #have been tkgrid before
#tkgrid.configure(frame.para1A, sticky = "w")

#tkgrid.configure(frame.dataformat, sticky = "w")
##-----------------------------------------------
#tkgrid(frame.group)
#tkgrid.configure(frame.group, sticky = "w")
#tkgrid(frame.chip)
#tkgrid.configure(frame.chip, sticky = "w")
#tkgrid(frame.chip2)
#tkgrid.configure(frame.chip2, sticky = "w")
#tkgrid(frame.para)
#tkgrid.configure(frame.para, sticky = "w")
##-----------------------------------------------
#tkgrid(frame.paratitle3)
#tkgrid.configure(frame.paratitle3, sticky = "w")
#tkgrid(frame.SA)
tkgrid.configure(frame.SA, sticky = "w")
#tkgrid(frame.para2)
tkgrid.configure(frame.para2, sticky = "w")
##----------------------------
#tkgrid(frame.para3)
tkgrid.configure(frame.para3, sticky = "w")
#tkgrid(frame.para3A)
tkgrid.configure(frame.para3A, sticky = "w")
#tkgrid(frame.para3D)
tkgrid.configure(frame.para3D, sticky = "w")
#tkgrid(frame.para3C)
tkgrid.configure(frame.para3C, sticky = "w")
#tkgrid(frame.para3B)
tkgrid.configure(frame.para3B, sticky = "w")
##----------------------------
#tkgrid(frame.para4)
tkgrid.configure(frame.para4, sticky = "w")
#tkgrid(frame.para4B)
tkgrid.configure(frame.para4B, sticky = "w")
#tkgrid(frame.para4A)
tkgrid.configure(frame.para4A, sticky = "w")
#tkgrid(frame.para4C)
#tkgrid.configure(frame.para4C, sticky = "w")
##----------------------------
#tkgrid(frame.para5)
tkgrid.configure(frame.para5, sticky = "w")
#tkgrid(frame.para5A)
tkgrid.configure(frame.para5A, sticky = "w")
tkgrid(tklabel(frame, text = "", bg = color))

##-----------------------------------------------------------------------------
								##Run botton
##-----------------------------------------------------------------------------
pathway_subfunction = paste(pathway_program, "SAQC_subfunction.r", sep = "")
source(pathway_subfunction)

pathway = strsplit(pathway_program, "PROGRAM")
pathway = unlist(pathway)[1]

SAQC.file <- paste(pathway_program, "SAQC.r", sep = "")  ## Used in SAQC.r
SAQC.para <- paste(pathway_program, "parameter_setting.r", sep = "")

##------------------------------------------------------------------------------
Runframe <- tkframe(frame, bg = color)
Run.button <- tkbutton(Runframe, text = "   Run   ", width = button.space, bg = color.run, font = fontRun, command = function(){
	cat("Please wait a moment...\n")
	source(SAQC.para)
})
blank1.space <- tklabel(Runframe, bg = color, text = "  ", font = fontTextLabel, width = blankspcae.run)
tkgrid(blank1.space, Run.button, sticky = "w")
tkgrid(tklabel(Runframe, text = "", bg = color))
tkgrid.configure(Runframe, sticky = "w")


##-----------------------------------------------------------------------------
						##tkgrid the whole frame!
##-----------------------------------------------------------------------------
tkgrid(frame)

##==============================================================================
##---------------------------   Notebook 2  ------------------------------------

##---------------color changing end part------------------------
#noteframe2 <- tkframe(QI_interface, background = color)
tb2 <- tk2notetab(nb, "Interactive visualization")

Frame <- tkframe(tb2, bg = color)

##-----------------------------------------------------------------------------
							##1. Input/output path
##-----------------------------------------------------------------------------
Frame.para1 <- tkframe(Frame, bg = color)
tkgrid(tklabel(Frame.para1, text = "   1. Input/Output path:", font = fontsubtitle, bg = color), sticky = "w") # title
tkgrid.configure(Frame.para1, sticky = "w")
#tkgrid(Frame.para1)

#-------------------------------# 1-1 Directory of QI data input #----------------------------#
Frame.path <- tkframe(Frame, bg = color)

QIfigure_pathlab <- tklabel(Frame.path, text = "            Directory of QI data input:", font = fontTextLabel, bg = color)
QIfigure_inputpath <- tclVar("Example")
entry.QI_Path <- tkentry(Frame.path, width = "85", textvariable = QIfigure_inputpath, font = fontTextLabel)
QIfigure_box.input <- tkbutton(Frame.path, text = "...", command = function() tclvalue(QIfigure_inputpath) <- tkchooseDirectory())
QIfigure_blank.space <- tklabel(Frame.path, bg = color, text = "  ", font = fontsubtitle, width = 50)

tkgrid(QIfigure_pathlab, entry.QI_Path, QIfigure_box.input, QIfigure_blank.space)
tkgrid.configure(Frame.path, sticky = "w")
#tkgrid(Frame.path)

#-------------------------------# 1-2 Directory of figure output #----------------------------#
Frame.outputpath <- tkframe(Frame, bg = color)

QIfigure_outputpathlab <- tklabel(Frame.outputpath, text = "            Directory of figure output:", font = fontTextLabel, bg = color)
QIfigure_outputpath <- tclVar(" ")
output.QI_Path <- tkentry(Frame.outputpath, width = "85", textvariable = QIfigure_outputpath, font = fontTextLabel, state = "normal")
QIfigure_box.output <- tkbutton(Frame.outputpath, text = "...",  command = function() tclvalue(QIfigure_outputpath) <- tkchooseDirectory())
QIfigure_outputblank.space <- tklabel(Frame.outputpath, bg = color, text = "  ", font = fontsubtitle, width = 50)

tkgrid(QIfigure_outputpathlab, output.QI_Path, QIfigure_box.output, QIfigure_outputblank.space)
tkgrid.configure(Frame.outputpath, sticky = "w")
tkgrid(Frame.outputpath)

##-----------------------------------------------------------------------------
							##2. Plot parameters
##-----------------------------------------------------------------------------
Frame.para2 <- tkframe(Frame, bg = color)

tkgrid(tklabel(Frame.para2, text = "   2. Plot parameters:", font = fontsubtitle, bg = color), sticky = "w") # title
tkgrid.configure(Frame.para2, sticky = "w")
#tkgrid(Frame.para2)

#-------------------------------# 2-1 Plot selection #----------------------------#
Frame.paraFigureSelection <- tkframe(Frame, bg = color)

QIfigure.port1 <- tkcheckbutton(Frame.paraFigureSelection, bg = color, command = function(){
	if(tclvalue(QIfigure.Value1) == 0){
		if(tclvalue(QIfigure.Value3) == 1) tclvalue(QIfigure.Value3) = 0
	} else if(tclvalue(QIfigure.Value2) == 1) tclvalue(QIfigure.Value3) = 1
	#if(tclvalue(QIfigure.Value1) == "1" & tclvalue(QIfigure.Value2) == "1") tclvalue(QIfigure.Value3) = 1	
	#if(tclvalue(QIfigure.Value1) == "0") {
	#	if(tclvalue(QIfigure.Value2) == "0"){
	#		tclvalue(QIfigure.Value1) = 1
	#	} else tclvalue(QIfigure.Value3) = 0
	#}
	if(0){
		#SAQC database
		tkconfigure(database.but, variable = QI_resouce.val, stat = "normal"); tclvalue(QI_resouce.val) = 1
		#User provided
		tkconfigure(userprovided.but, variable = QI_resouce.val, stat = "normal");
		tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, value = "1", stat = "normal")
		tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, value = "2", stat = "normal")
		tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, value = "3", stat = "normal")
		tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, value = "4", stat = "normal")
		tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, value = "5", stat = "normal")
		tclvalue(QIfigure_groupValue) = 5
		tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, value = "1", stat = "normal")
		tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, value = "2", stat = "normal")
		tclvalue(QIfigure_chipValue) = 1
		tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "normal")
		tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "normal")
		tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "normal")
		tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "normal")
		tclvalue(QIfigure_quantile.Value1) = 1;tclvalue(QIfigure_quantile.Value2) = 1
		tclvalue(QIfigure_quantile.Value3) = 1;tclvalue(QIfigure_quantile.Value4) = 1
		##----------------------------------------------------------------------------------
		tkconfigure(MQI1_chip1.q95.box, state = "disable");tclvalue(MQI1_chip1.q95.val) = ""
		tkconfigure(MQI1_chip2.q95.box, state = "disable");tclvalue(MQI1_chip2.q95.val) = ""
		tkconfigure(MQI1_merge.q95.box, state = "disable");tclvalue(MQI1_merge.q95.val) = ""
		tkconfigure(MQI2_chip1.q95.box, state = "disable");tclvalue(MQI2_chip1.q95.val) = ""
		tkconfigure(MQI2_chip2.q95.box, state = "disable");tclvalue(MQI2_chip2.q95.val) = ""
		tkconfigure(MQI2_merge.q95.box, state = "disable");tclvalue(MQI2_merge.q95.val) = ""
		tkconfigure(WQI1_chip1.q95.box, state = "disable");tclvalue(WQI1_chip1.q95.val) = ""
		tkconfigure(WQI1_chip2.q95.box, state = "disable");tclvalue(WQI1_chip2.q95.val) = ""
		tkconfigure(WQI1_merge.q95.box, state = "disable");tclvalue(WQI1_merge.q95.val) = ""
		tkconfigure(WQI2_chip1.q95.box, state = "disable");tclvalue(WQI2_chip1.q95.val) = ""
		tkconfigure(WQI2_chip2.q95.box, state = "disable");tclvalue(WQI2_chip2.q95.val) = ""
		tkconfigure(WQI2_merge.q95.box, state = "disable");tclvalue(WQI2_merge.q95.val) = ""
		tkconfigure(MQI1_chip1.q975.box, state = "disable");tclvalue(MQI1_chip1.q975.val) = ""
		tkconfigure(MQI1_chip2.q975.box, state = "disable");tclvalue(MQI1_chip2.q975.val) = ""
		tkconfigure(MQI1_merge.q975.box, state = "disable");tclvalue(MQI1_merge.q975.val) = ""
		tkconfigure(MQI2_chip1.q975.box, state = "disable");tclvalue(MQI2_chip1.q975.val) = ""
		tkconfigure(MQI2_chip2.q975.box, state = "disable");tclvalue(MQI2_chip2.q975.val) = ""
		tkconfigure(MQI2_merge.q975.box, state = "disable");tclvalue(MQI2_merge.q975.val) = ""
		tkconfigure(WQI1_chip1.q975.box, state = "disable");tclvalue(WQI1_chip1.q975.val) = ""
		tkconfigure(WQI1_chip2.q975.box, state = "disable");tclvalue(WQI1_chip2.q975.val) = ""
		tkconfigure(WQI1_merge.q975.box, state = "disable");tclvalue(WQI1_merge.q975.val) = ""
		tkconfigure(WQI2_chip1.q975.box, state = "disable");tclvalue(WQI2_chip1.q975.val) = ""
		tkconfigure(WQI2_chip2.q975.box, state = "disable");tclvalue(WQI2_chip2.q975.val) = ""
		tkconfigure(WQI2_merge.q975.box, state = "disable");tclvalue(WQI2_merge.q975.val) = ""
		tkconfigure(MQI1_chip1.q99.box, state = "disable");tclvalue(MQI1_chip1.q99.val) = ""
		tkconfigure(MQI1_chip2.q99.box, state = "disable");tclvalue(MQI1_chip2.q99.val) = ""
		tkconfigure(MQI1_merge.q99.box, state = "disable");tclvalue(MQI1_merge.q99.val) = ""
		tkconfigure(MQI2_chip1.q99.box, state = "disable");tclvalue(MQI2_chip1.q99.val) = ""
		tkconfigure(MQI2_chip2.q99.box, state = "disable");tclvalue(MQI2_chip2.q99.val) = ""
		tkconfigure(MQI2_merge.q99.box, state = "disable");tclvalue(MQI2_merge.q99.val) = ""
		tkconfigure(WQI1_chip1.q99.box, state = "disable");tclvalue(WQI1_chip1.q99.val) = ""
		tkconfigure(WQI1_chip2.q99.box, state = "disable");tclvalue(WQI1_chip2.q99.val) = ""
		tkconfigure(WQI1_merge.q99.box, state = "disable");tclvalue(WQI1_merge.q99.val) = ""
		tkconfigure(WQI2_chip1.q99.box, state = "disable");tclvalue(WQI2_chip1.q99.val) = ""
		tkconfigure(WQI2_chip2.q99.box, state = "disable");tclvalue(WQI2_chip2.q99.val) = ""
		tkconfigure(WQI2_merge.q99.box, state = "disable");tclvalue(WQI2_merge.q99.val) = ""
		#---------------------------------------------------------------------------------------------------#
        tkconfigure(database.but, variable = QI_resouce.val, stat = "disable");tclvalue(QI_resouce.val) = 0
        tkconfigure(userprovided.but, variable = QI_resouce.val, stat = "disable");
        tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, stat = "disable")
        tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, stat = "disable")
        tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, stat = "disable")
        tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, stat = "disable")
        tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, stat = "disable")
        tclvalue(QIfigure_groupValue) = 0
        tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, stat = "disable")
        tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, stat = "disable")
        tclvalue(QIfigure_chipValue) = 0
        tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "disable")
        tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "disable")
        tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "disable")
        tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "disable")
        tclvalue(QIfigure_quantile.Value1) = 0;tclvalue(QIfigure_quantile.Value2) = 0
        tclvalue(QIfigure_quantile.Value3) = 0;tclvalue(QIfigure_quantile.Value4) = 0
        ##----------------------------------------------------------------------------------
        tkconfigure(MQI1_chip1.q95.box, state = "disable");tclvalue(MQI1_chip1.q95.val) = ""
        tkconfigure(MQI1_chip2.q95.box, state = "disable");tclvalue(MQI1_chip2.q95.val) = ""
        tkconfigure(MQI1_merge.q95.box, state = "disable");tclvalue(MQI1_merge.q95.val) = ""
        tkconfigure(MQI2_chip1.q95.box, state = "disable");tclvalue(MQI2_chip1.q95.val) = ""
        tkconfigure(MQI2_chip2.q95.box, state = "disable");tclvalue(MQI2_chip2.q95.val) = ""
        tkconfigure(MQI2_merge.q95.box, state = "disable");tclvalue(MQI2_merge.q95.val) = ""
        tkconfigure(WQI1_chip1.q95.box, state = "disable");tclvalue(WQI1_chip1.q95.val) = ""
        tkconfigure(WQI1_chip2.q95.box, state = "disable");tclvalue(WQI1_chip2.q95.val) = ""
        tkconfigure(WQI1_merge.q95.box, state = "disable");tclvalue(WQI1_merge.q95.val) = ""
        tkconfigure(WQI2_chip1.q95.box, state = "disable");tclvalue(WQI2_chip1.q95.val) = ""
        tkconfigure(WQI2_chip2.q95.box, state = "disable");tclvalue(WQI2_chip2.q95.val) = ""
        tkconfigure(WQI2_merge.q95.box, state = "disable");tclvalue(WQI2_merge.q95.val) = ""
        tkconfigure(MQI1_chip1.q975.box, state = "disable");tclvalue(MQI1_chip1.q975.val) = ""
        tkconfigure(MQI1_chip2.q975.box, state = "disable");tclvalue(MQI1_chip2.q975.val) = ""
        tkconfigure(MQI1_merge.q975.box, state = "disable");tclvalue(MQI1_merge.q975.val) = ""
        tkconfigure(MQI2_chip1.q975.box, state = "disable");tclvalue(MQI2_chip1.q975.val) = ""
        tkconfigure(MQI2_chip2.q975.box, state = "disable");tclvalue(MQI2_chip2.q975.val) = ""
        tkconfigure(MQI2_merge.q975.box, state = "disable");tclvalue(MQI2_merge.q975.val) = ""
        tkconfigure(WQI1_chip1.q975.box, state = "disable");tclvalue(WQI1_chip1.q975.val) = ""
        tkconfigure(WQI1_chip2.q975.box, state = "disable");tclvalue(WQI1_chip2.q975.val) = ""
        tkconfigure(WQI1_merge.q975.box, state = "disable");tclvalue(WQI1_merge.q975.val) = ""
        tkconfigure(WQI2_chip1.q975.box, state = "disable");tclvalue(WQI2_chip1.q975.val) = ""
        tkconfigure(WQI2_chip2.q975.box, state = "disable");tclvalue(WQI2_chip2.q975.val) = ""
        tkconfigure(WQI2_merge.q975.box, state = "disable");tclvalue(WQI2_merge.q975.val) = ""
        tkconfigure(MQI1_chip1.q99.box, state = "disable");tclvalue(MQI1_chip1.q99.val) = ""
        tkconfigure(MQI1_chip2.q99.box, state = "disable");tclvalue(MQI1_chip2.q99.val) = ""
        tkconfigure(MQI1_merge.q99.box, state = "disable");tclvalue(MQI1_merge.q99.val) = ""
        tkconfigure(MQI2_chip1.q99.box, state = "disable");tclvalue(MQI2_chip1.q99.val) = ""
        tkconfigure(MQI2_chip2.q99.box, state = "disable");tclvalue(MQI2_chip2.q99.val) = ""
        tkconfigure(MQI2_merge.q99.box, state = "disable");tclvalue(MQI2_merge.q99.val) = ""
        tkconfigure(WQI1_chip1.q99.box, state = "disable");tclvalue(WQI1_chip1.q99.val) = ""
        tkconfigure(WQI1_chip2.q99.box, state = "disable");tclvalue(WQI1_chip2.q99.val) = ""
        tkconfigure(WQI1_merge.q99.box, state = "disable");tclvalue(WQI1_merge.q99.val) = ""
        tkconfigure(WQI2_chip1.q99.box, state = "disable");tclvalue(WQI2_chip1.q99.val) = ""
        tkconfigure(WQI2_chip2.q99.box, state = "disable");tclvalue(WQI2_chip2.q99.val) = ""
        tkconfigure(WQI2_merge.q99.box, state = "disable");tclvalue(WQI2_merge.q99.val) = ""
	}
})
QIfigure.port2 <- tkcheckbutton(Frame.paraFigureSelection, bg = color, command = function(){
	if(tclvalue(QIfigure.Value2) == "0"){
		if(tclvalue(QIfigure.Value3) == 1) tclvalue(QIfigure.Value3) = 0
		
		#disable SAQC database & User provided
		interactive_polygon(button.stat1 = "disable", button.stat2 = "disable", default.val1.1 = 0, 
								default.val1.2 = 0, default.val1.3 = 0, default.val1.4 = 0)	
		#if(tclvalue(QIfigure.Value1) == "0"){
		#	tclvalue(QIfigure.Value2) = 1
		#} else {
		#	tclvalue(QIfigure.Value3) = 0			
		#	#disable SAQC database & User provided
		#	interactive_polygon(button.stat1 = "disable", button.stat2 = "disable", default.val1.1 = 0, 
		#						default.val1.2 = 0, default.val1.3 = 0, default.val1.4 = 0, default.val2 = "")
		#}
	} else { #tclvalue(QIfigure.Value2) == "1"
		if(tclvalue(QIfigure.Value1) == "1") tclvalue(QIfigure.Value3) = 1
		#enable SAQC database & disable User provided
		interactive_polygon(button.stat1 = "normal", button.stat2 = "disable", default.val1.1 = 1, 
								default.val1.2 = 5, default.val1.3 = 1, default.val1.4 = 1)
		if(0){		
			tkconfigure(database.but, variable = QI_resouce.val, stat = "normal");tclvalue(QI_resouce.val) = 1
			tkconfigure(userprovided.but, variable = QI_resouce.val, stat = "normal");
			tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, value = "1", stat = "normal")
			tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, value = "2", stat = "normal")
			tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, value = "3", stat = "normal")
			tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, value = "4", stat = "normal")
			tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, value = "5", stat = "normal")
			tclvalue(QIfigure_groupValue) = 5
			tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, value = "1", stat = "normal")
			tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, value = "2", stat = "normal")
			tclvalue(QIfigure_chipValue) = 1
			tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "normal")
			tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "normal")
			tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "normal")
			tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "normal")
			tclvalue(QIfigure_quantile.Value1) = 1;tclvalue(QIfigure_quantile.Value2) = 1
			tclvalue(QIfigure_quantile.Value3) = 1;tclvalue(QIfigure_quantile.Value4) = 1
			##----------------------------------------------------------------------------------
			tkconfigure(MQI1_chip1.q95.box, state = "disable");tclvalue(MQI1_chip1.q95.val) = ""
			tkconfigure(MQI1_chip2.q95.box, state = "disable");tclvalue(MQI1_chip2.q95.val) = ""
			tkconfigure(MQI1_merge.q95.box, state = "disable");tclvalue(MQI1_merge.q95.val) = ""
			tkconfigure(MQI2_chip1.q95.box, state = "disable");tclvalue(MQI2_chip1.q95.val) = ""
			tkconfigure(MQI2_chip2.q95.box, state = "disable");tclvalue(MQI2_chip2.q95.val) = ""
			tkconfigure(MQI2_merge.q95.box, state = "disable");tclvalue(MQI2_merge.q95.val) = ""
			tkconfigure(WQI1_chip1.q95.box, state = "disable");tclvalue(WQI1_chip1.q95.val) = ""
			tkconfigure(WQI1_chip2.q95.box, state = "disable");tclvalue(WQI1_chip2.q95.val) = ""
			tkconfigure(WQI1_merge.q95.box, state = "disable");tclvalue(WQI1_merge.q95.val) = ""
			tkconfigure(WQI2_chip1.q95.box, state = "disable");tclvalue(WQI2_chip1.q95.val) = ""
			tkconfigure(WQI2_chip2.q95.box, state = "disable");tclvalue(WQI2_chip2.q95.val) = ""
			tkconfigure(WQI2_merge.q95.box, state = "disable");tclvalue(WQI2_merge.q95.val) = ""
			tkconfigure(MQI1_chip1.q975.box, state = "disable");tclvalue(MQI1_chip1.q975.val) = ""
			tkconfigure(MQI1_chip2.q975.box, state = "disable");tclvalue(MQI1_chip2.q975.val) = ""
			tkconfigure(MQI1_merge.q975.box, state = "disable");tclvalue(MQI1_merge.q975.val) = ""
			tkconfigure(MQI2_chip1.q975.box, state = "disable");tclvalue(MQI2_chip1.q975.val) = ""
			tkconfigure(MQI2_chip2.q975.box, state = "disable");tclvalue(MQI2_chip2.q975.val) = ""
			tkconfigure(MQI2_merge.q975.box, state = "disable");tclvalue(MQI2_merge.q975.val) = ""
			tkconfigure(WQI1_chip1.q975.box, state = "disable");tclvalue(WQI1_chip1.q975.val) = ""
			tkconfigure(WQI1_chip2.q975.box, state = "disable");tclvalue(WQI1_chip2.q975.val) = ""
			tkconfigure(WQI1_merge.q975.box, state = "disable");tclvalue(WQI1_merge.q975.val) = ""
			tkconfigure(WQI2_chip1.q975.box, state = "disable");tclvalue(WQI2_chip1.q975.val) = ""
			tkconfigure(WQI2_chip2.q975.box, state = "disable");tclvalue(WQI2_chip2.q975.val) = ""
			tkconfigure(WQI2_merge.q975.box, state = "disable");tclvalue(WQI2_merge.q975.val) = ""
			tkconfigure(MQI1_chip1.q99.box, state = "disable");tclvalue(MQI1_chip1.q99.val) = ""
			tkconfigure(MQI1_chip2.q99.box, state = "disable");tclvalue(MQI1_chip2.q99.val) = ""
			tkconfigure(MQI1_merge.q99.box, state = "disable");tclvalue(MQI1_merge.q99.val) = ""
			tkconfigure(MQI2_chip1.q99.box, state = "disable");tclvalue(MQI2_chip1.q99.val) = ""
			tkconfigure(MQI2_chip2.q99.box, state = "disable");tclvalue(MQI2_chip2.q99.val) = ""
			tkconfigure(MQI2_merge.q99.box, state = "disable");tclvalue(MQI2_merge.q99.val) = ""
			tkconfigure(WQI1_chip1.q99.box, state = "disable");tclvalue(WQI1_chip1.q99.val) = ""
			tkconfigure(WQI1_chip2.q99.box, state = "disable");tclvalue(WQI1_chip2.q99.val) = ""
			tkconfigure(WQI1_merge.q99.box, state = "disable");tclvalue(WQI1_merge.q99.val) = ""
			tkconfigure(WQI2_chip1.q99.box, state = "disable");tclvalue(WQI2_chip1.q99.val) = ""
			tkconfigure(WQI2_chip2.q99.box, state = "disable");tclvalue(WQI2_chip2.q99.val) = ""
			tkconfigure(WQI2_merge.q99.box, state = "disable");tclvalue(WQI2_merge.q99.val) = ""
		}
	}
})
      
QIfigure.port3 <- tkcheckbutton(Frame.paraFigureSelection, bg = color, command = function(){
	#tclvalue(QIfigure.Value3) = 1
	#if(tclvalue(QIfigure.Value1) == 0) tclvalue(QIfigure.Value1) = 1
	#if(tclvalue(QIfigure.Value2) == 0){
	#	tclvalue(QIfigure.Value2) = 1
	#	interactive_polygon(button.stat1 = "normal", button.stat2 = "disable", default.val1.1 = 1, 
	#							default.val1.2 = 5, default.val1.3 = 1, default.val1.4 = 1, default.val2 = "")
	#} else tclvalue(QIfigure.Value1) = 1
	
	if(tclvalue(QIfigure.Value3) == "1"){
		#tclvalue(QIfigure.Value1) = 1; tclvalue(QIfigure.Value2) = 1
		if(tclvalue(QIfigure.Value1) == 0) tclvalue(QIfigure.Value1) = 1
		if(tclvalue(QIfigure.Value2) == 0){
			tclvalue(QIfigure.Value2) = 1
			interactive_polygon(button.stat1 = "normal", button.stat2 = "disable", default.val1.1 = 1, 
									default.val1.2 = 5, default.val1.3 = 1, default.val1.4 = 1)		
		}
		if(0){
			tkconfigure(database.but, variable = QI_resouce.val, stat = "normal");tclvalue(QI_resouce.val) = 0
			tkconfigure(userprovided.but, variable = QI_resouce.val, stat = "normal");tclvalue(QI_resouce.val) = 1
			tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, value = "1", stat = "normal")
			tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, value = "2", stat = "normal")
			tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, value = "3", stat = "normal")
			tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, value = "4", stat = "normal")
			tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, value = "5", stat = "normal")
			tclvalue(QIfigure_groupValue) = 5
			tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, value = "1", stat = "normal")
			tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, value = "2", stat = "normal")
			tclvalue(QIfigure_chipValue) = 1
			tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "normal")
			tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "normal")
			tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "normal")
			tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "normal")
			tclvalue(QIfigure_quantile.Value1) = 1;tclvalue(QIfigure_quantile.Value2) = 1
			tclvalue(QIfigure_quantile.Value3) = 1;tclvalue(QIfigure_quantile.Value4) = 1
			##----------------------------------------------------------------------------------
			tkconfigure(MQI1_chip1.q95.box, state = "disable");tclvalue(MQI1_chip1.q95.val) = ""
			tkconfigure(MQI1_chip2.q95.box, state = "disable");tclvalue(MQI1_chip2.q95.val) = ""
			tkconfigure(MQI1_merge.q95.box, state = "disable");tclvalue(MQI1_merge.q95.val) = ""
			tkconfigure(MQI2_chip1.q95.box, state = "disable");tclvalue(MQI2_chip1.q95.val) = ""
			tkconfigure(MQI2_chip2.q95.box, state = "disable");tclvalue(MQI2_chip2.q95.val) = ""
			tkconfigure(MQI2_merge.q95.box, state = "disable");tclvalue(MQI2_merge.q95.val) = ""
			tkconfigure(WQI1_chip1.q95.box, state = "disable");tclvalue(WQI1_chip1.q95.val) = ""
			tkconfigure(WQI1_chip2.q95.box, state = "disable");tclvalue(WQI1_chip2.q95.val) = ""
			tkconfigure(WQI1_merge.q95.box, state = "disable");tclvalue(WQI1_merge.q95.val) = ""
			tkconfigure(WQI2_chip1.q95.box, state = "disable");tclvalue(WQI2_chip1.q95.val) = ""
			tkconfigure(WQI2_chip2.q95.box, state = "disable");tclvalue(WQI2_chip2.q95.val) = ""
			tkconfigure(WQI2_merge.q95.box, state = "disable");tclvalue(WQI2_merge.q95.val) = ""
			tkconfigure(MQI1_chip1.q975.box, state = "disable");tclvalue(MQI1_chip1.q975.val) = ""
			tkconfigure(MQI1_chip2.q975.box, state = "disable");tclvalue(MQI1_chip2.q975.val) = ""
			tkconfigure(MQI1_merge.q975.box, state = "disable");tclvalue(MQI1_merge.q975.val) = ""
			tkconfigure(MQI2_chip1.q975.box, state = "disable");tclvalue(MQI2_chip1.q975.val) = ""
			tkconfigure(MQI2_chip2.q975.box, state = "disable");tclvalue(MQI2_chip2.q975.val) = ""
			tkconfigure(MQI2_merge.q975.box, state = "disable");tclvalue(MQI2_merge.q975.val) = ""
			tkconfigure(WQI1_chip1.q975.box, state = "disable");tclvalue(WQI1_chip1.q975.val) = ""
			tkconfigure(WQI1_chip2.q975.box, state = "disable");tclvalue(WQI1_chip2.q975.val) = ""
			tkconfigure(WQI1_merge.q975.box, state = "disable");tclvalue(WQI1_merge.q975.val) = ""
			tkconfigure(WQI2_chip1.q975.box, state = "disable");tclvalue(WQI2_chip1.q975.val) = ""
			tkconfigure(WQI2_chip2.q975.box, state = "disable");tclvalue(WQI2_chip2.q975.val) = ""
			tkconfigure(WQI2_merge.q975.box, state = "disable");tclvalue(WQI2_merge.q975.val) = ""
			tkconfigure(MQI1_chip1.q99.box, state = "disable");tclvalue(MQI1_chip1.q99.val) = ""
			tkconfigure(MQI1_chip2.q99.box, state = "disable");tclvalue(MQI1_chip2.q99.val) = ""
			tkconfigure(MQI1_merge.q99.box, state = "disable");tclvalue(MQI1_merge.q99.val) = ""
			tkconfigure(MQI2_chip1.q99.box, state = "disable");tclvalue(MQI2_chip1.q99.val) = ""
			tkconfigure(MQI2_chip2.q99.box, state = "disable");tclvalue(MQI2_chip2.q99.val) = ""
			tkconfigure(MQI2_merge.q99.box, state = "disable");tclvalue(MQI2_merge.q99.val) = ""
			tkconfigure(WQI1_chip1.q99.box, state = "disable");tclvalue(WQI1_chip1.q99.val) = ""
			tkconfigure(WQI1_chip2.q99.box, state = "disable");tclvalue(WQI1_chip2.q99.val) = ""
			tkconfigure(WQI1_merge.q99.box, state = "disable");tclvalue(WQI1_merge.q99.val) = ""
			tkconfigure(WQI2_chip1.q99.box, state = "disable");tclvalue(WQI2_chip1.q99.val) = ""
			tkconfigure(WQI2_chip2.q99.box, state = "disable");tclvalue(WQI2_chip2.q99.val) = ""
			tkconfigure(WQI2_merge.q99.box, state = "disable");tclvalue(WQI2_merge.q99.val) = ""
		}
	} else {
		if(tclvalue(QIfigure.Value1) == 1) tclvalue(QIfigure.Value1) = 0
		if(tclvalue(QIfigure.Value2) == 1) {
			tclvalue(QIfigure.Value2) = 0
			interactive_polygon(button.stat1 = "disable", button.stat2 = "disable", default.val1.1 = 0, 
									default.val1.2 = 0, default.val1.3 = 0, default.val1.4 = 0)
		}
		if(0){
			tkconfigure(database.but, variable = QI_resouce.val, stat = "disable");tclvalue(QI_resouce.val) = 0
			tkconfigure(userprovided.but, variable = QI_resouce.val, stat = "disable");tclvalue(QI_resouce.val) = 0
			tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, stat = "disable")
			tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, stat = "disable")
			tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, stat = "disable")
			tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, stat = "disable")
			tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, stat = "disable")
			tclvalue(QIfigure_groupValue) = 0
			tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, stat = "disable")
			tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, stat = "disable")
			tclvalue(QIfigure_chipValue) = 0
			tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "disable")
			tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "disable")
			tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "disable")
			tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "disable")
			tclvalue(QIfigure_quantile.Value1) = 0;tclvalue(QIfigure_quantile.Value2) = 0
			tclvalue(QIfigure_quantile.Value3) = 0;tclvalue(QIfigure_quantile.Value4) = 0
			##----------------------------------------------------------------------------------
			tkconfigure(MQI1_chip1.q95.box, state = "disable");tclvalue(MQI1_chip1.q95.val) = ""
			tkconfigure(MQI1_chip2.q95.box, state = "disable");tclvalue(MQI1_chip2.q95.val) = ""
			tkconfigure(MQI1_merge.q95.box, state = "disable");tclvalue(MQI1_merge.q95.val) = ""
			tkconfigure(MQI2_chip1.q95.box, state = "disable");tclvalue(MQI2_chip1.q95.val) = ""
			tkconfigure(MQI2_chip2.q95.box, state = "disable");tclvalue(MQI2_chip2.q95.val) = ""
			tkconfigure(MQI2_merge.q95.box, state = "disable");tclvalue(MQI2_merge.q95.val) = ""
			tkconfigure(WQI1_chip1.q95.box, state = "disable");tclvalue(WQI1_chip1.q95.val) = ""
			tkconfigure(WQI1_chip2.q95.box, state = "disable");tclvalue(WQI1_chip2.q95.val) = ""
			tkconfigure(WQI1_merge.q95.box, state = "disable");tclvalue(WQI1_merge.q95.val) = ""
			tkconfigure(WQI2_chip1.q95.box, state = "disable");tclvalue(WQI2_chip1.q95.val) = ""
			tkconfigure(WQI2_chip2.q95.box, state = "disable");tclvalue(WQI2_chip2.q95.val) = ""
			tkconfigure(WQI2_merge.q95.box, state = "disable");tclvalue(WQI2_merge.q95.val) = ""
			tkconfigure(MQI1_chip1.q975.box, state = "disable");tclvalue(MQI1_chip1.q975.val) = ""
			tkconfigure(MQI1_chip2.q975.box, state = "disable");tclvalue(MQI1_chip2.q975.val) = ""
			tkconfigure(MQI1_merge.q975.box, state = "disable");tclvalue(MQI1_merge.q975.val) = ""
			tkconfigure(MQI2_chip1.q975.box, state = "disable");tclvalue(MQI2_chip1.q975.val) = ""
			tkconfigure(MQI2_chip2.q975.box, state = "disable");tclvalue(MQI2_chip2.q975.val) = ""
			tkconfigure(MQI2_merge.q975.box, state = "disable");tclvalue(MQI2_merge.q975.val) = ""
			tkconfigure(WQI1_chip1.q975.box, state = "disable");tclvalue(WQI1_chip1.q975.val) = ""
			tkconfigure(WQI1_chip2.q975.box, state = "disable");tclvalue(WQI1_chip2.q975.val) = ""
			tkconfigure(WQI1_merge.q975.box, state = "disable");tclvalue(WQI1_merge.q975.val) = ""
			tkconfigure(WQI2_chip1.q975.box, state = "disable");tclvalue(WQI2_chip1.q975.val) = ""
			tkconfigure(WQI2_chip2.q975.box, state = "disable");tclvalue(WQI2_chip2.q975.val) = ""
			tkconfigure(WQI2_merge.q975.box, state = "disable");tclvalue(WQI2_merge.q975.val) = ""
			tkconfigure(MQI1_chip1.q99.box, state = "disable");tclvalue(MQI1_chip1.q99.val) = ""
			tkconfigure(MQI1_chip2.q99.box, state = "disable");tclvalue(MQI1_chip2.q99.val) = ""
			tkconfigure(MQI1_merge.q99.box, state = "disable");tclvalue(MQI1_merge.q99.val) = ""
			tkconfigure(MQI2_chip1.q99.box, state = "disable");tclvalue(MQI2_chip1.q99.val) = ""
			tkconfigure(MQI2_chip2.q99.box, state = "disable");tclvalue(MQI2_chip2.q99.val) = ""
			tkconfigure(MQI2_merge.q99.box, state = "disable");tclvalue(MQI2_merge.q99.val) = ""
			tkconfigure(WQI1_chip1.q99.box, state = "disable");tclvalue(WQI1_chip1.q99.val) = ""
			tkconfigure(WQI1_chip2.q99.box, state = "disable");tclvalue(WQI1_chip2.q99.val) = ""
			tkconfigure(WQI1_merge.q99.box, state = "disable");tclvalue(WQI1_merge.q99.val) = ""
			tkconfigure(WQI2_chip1.q99.box, state = "disable");tclvalue(WQI2_chip1.q99.val) = ""
			tkconfigure(WQI2_chip2.q99.box, state = "disable");tclvalue(WQI2_chip2.q99.val) = ""
			tkconfigure(WQI2_merge.q99.box, state = "disable");tclvalue(WQI2_merge.q99.val) = ""
		}
	}

})

QIfigure.Value1 <- tclVar("1") #default both HeatMap & Polygon
QIfigure.Value2 <- tclVar("1")
QIfigure.Value3 <- tclVar("1")

tkconfigure(QIfigure.port1, variable = QIfigure.Value1)
tkconfigure(QIfigure.port2, variable = QIfigure.Value2)
tkconfigure(QIfigure.port3, variable = QIfigure.Value3)

QIfigure.lab<-tklabel(Frame.paraFigureSelection, text = "            Plot selection:", font = fontTextLabel, bg = color)
HeatMap <- tklabel(Frame.paraFigureSelection, text = "HeatMap", font = fontTextLabel, bg = color)
Polygon <- tklabel(Frame.paraFigureSelection, text = "Polygon", font = fontTextLabel, bg = color)
Allplot <- tklabel(Frame.paraFigureSelection, text = "Both", font = fontTextLabel, bg = color)

tkgrid(QIfigure.lab, QIfigure.port1, HeatMap, QIfigure.port2, Polygon, QIfigure.port3, Allplot, sticky = "w")
#chip.port3, chip.affy6, chip.port4, chip.illumina550, sticky = "w")

tkgrid.configure(Frame.paraFigureSelection, sticky = "w")

#-------------------------------# 2-2 QI reference #----------------------------#
Frame.para2A <- tkframe(Frame, bg = color)
tkgrid(tklabel(Frame.para2A, text = "            QI reference:", font = fontTextLabel, bg = color), sticky = "w") # title

tkgrid.configure(Frame.para2A, sticky = "w")
#tkgrid(Frame.para2A)

#-------------------------------# 2-2-1 SAQC database #----------------------------#
Frame.title1.para2A <- tkframe(Frame.para2A, bg = color)
database.but <- tkradiobutton(Frame.title1.para2A, bg = color, command = function(){
	#tclvalue(QI_resouce1.val) = 1;tclvalue(QI_resouce2.val) = 0
	interactive_polygon(index.decide.ref = T, button.stat1 = "normal", button.stat2 = "disable", 
								default.val1.2 = 5, default.val1.3 = 1, default.val1.4 = 1)
	if(0){
		tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, value = "1", stat = "normal")
		tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, value = "2", stat = "normal")
		tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, value = "3", stat = "normal")
		tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, value = "4", stat = "normal")
		tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, value = "5", stat = "normal")
		tclvalue(QIfigure_groupValue) = 5
		tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, value = "1", stat = "normal")
		tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, value = "2", stat = "normal")
		tclvalue(QIfigure_chipValue) = 1
		tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "normal")
		tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "normal")
		tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "normal")
		tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "normal")
		tclvalue(QIfigure_quantile.Value1) = 1;tclvalue(QIfigure_quantile.Value2) = 1
		tclvalue(QIfigure_quantile.Value3) = 1;tclvalue(QIfigure_quantile.Value4) = 1
		##----------------------------------------------------------------------------------
		tkconfigure(MQI1_chip1.q95.box, state = "disable");tclvalue(MQI1_chip1.q95.val) = ""
		tkconfigure(MQI1_chip2.q95.box, state = "disable");tclvalue(MQI1_chip2.q95.val) = ""
		tkconfigure(MQI1_merge.q95.box, state = "disable");tclvalue(MQI1_merge.q95.val) = ""
		tkconfigure(MQI2_chip1.q95.box, state = "disable");tclvalue(MQI2_chip1.q95.val) = ""
		tkconfigure(MQI2_chip2.q95.box, state = "disable");tclvalue(MQI2_chip2.q95.val) = ""
		tkconfigure(MQI2_merge.q95.box, state = "disable");tclvalue(MQI2_merge.q95.val) = ""
		tkconfigure(WQI1_chip1.q95.box, state = "disable");tclvalue(WQI1_chip1.q95.val) = ""
		tkconfigure(WQI1_chip2.q95.box, state = "disable");tclvalue(WQI1_chip2.q95.val) = ""
		tkconfigure(WQI1_merge.q95.box, state = "disable");tclvalue(WQI1_merge.q95.val) = ""
		tkconfigure(WQI2_chip1.q95.box, state = "disable");tclvalue(WQI2_chip1.q95.val) = ""
		tkconfigure(WQI2_chip2.q95.box, state = "disable");tclvalue(WQI2_chip2.q95.val) = ""
		tkconfigure(WQI2_merge.q95.box, state = "disable");tclvalue(WQI2_merge.q95.val) = ""
		tkconfigure(MQI1_chip1.q975.box, state = "disable");tclvalue(MQI1_chip1.q975.val) = ""
		tkconfigure(MQI1_chip2.q975.box, state = "disable");tclvalue(MQI1_chip2.q975.val) = ""
		tkconfigure(MQI1_merge.q975.box, state = "disable");tclvalue(MQI1_merge.q975.val) = ""
		tkconfigure(MQI2_chip1.q975.box, state = "disable");tclvalue(MQI2_chip1.q975.val) = ""
		tkconfigure(MQI2_chip2.q975.box, state = "disable");tclvalue(MQI2_chip2.q975.val) = ""
		tkconfigure(MQI2_merge.q975.box, state = "disable");tclvalue(MQI2_merge.q975.val) = ""
		tkconfigure(WQI1_chip1.q975.box, state = "disable");tclvalue(WQI1_chip1.q975.val) = ""
		tkconfigure(WQI1_chip2.q975.box, state = "disable");tclvalue(WQI1_chip2.q975.val) = ""
		tkconfigure(WQI1_merge.q975.box, state = "disable");tclvalue(WQI1_merge.q975.val) = ""
		tkconfigure(WQI2_chip1.q975.box, state = "disable");tclvalue(WQI2_chip1.q975.val) = ""
		tkconfigure(WQI2_chip2.q975.box, state = "disable");tclvalue(WQI2_chip2.q975.val) = ""
		tkconfigure(WQI2_merge.q975.box, state = "disable");tclvalue(WQI2_merge.q975.val) = ""
		tkconfigure(MQI1_chip1.q99.box, state = "disable");tclvalue(MQI1_chip1.q99.val) = ""
		tkconfigure(MQI1_chip2.q99.box, state = "disable");tclvalue(MQI1_chip2.q99.val) = ""
		tkconfigure(MQI1_merge.q99.box, state = "disable");tclvalue(MQI1_merge.q99.val) = ""
		tkconfigure(MQI2_chip1.q99.box, state = "disable");tclvalue(MQI2_chip1.q99.val) = ""
		tkconfigure(MQI2_chip2.q99.box, state = "disable");tclvalue(MQI2_chip2.q99.val) = ""
		tkconfigure(MQI2_merge.q99.box, state = "disable");tclvalue(MQI2_merge.q99.val) = ""
		tkconfigure(WQI1_chip1.q99.box, state = "disable");tclvalue(WQI1_chip1.q99.val) = ""
		tkconfigure(WQI1_chip2.q99.box, state = "disable");tclvalue(WQI1_chip2.q99.val) = ""
		tkconfigure(WQI1_merge.q99.box, state = "disable");tclvalue(WQI1_merge.q99.val) = ""
		tkconfigure(WQI2_chip1.q99.box, state = "disable");tclvalue(WQI2_chip1.q99.val) = ""
		tkconfigure(WQI2_chip2.q99.box, state = "disable");tclvalue(WQI2_chip2.q99.val) = ""
		tkconfigure(WQI2_merge.q99.box, state = "disable");tclvalue(WQI2_merge.q99.val) = ""
	}
})
#QI_resouce1.val <- tclVar("1")
#tkconfigure(database.but, variable = QI_resouce.val)  #
QI_resouce.val <- tclVar("1") #default SAQC database
tkconfigure(database.but, variable = QI_resouce.val, value = "1")

nblank1 = 3
blank1.space <- tklabel(Frame.title1.para2A, bg = color, text = "  ", font = fontsubtitle, width = nblank1)
tkgrid(blank1.space, database.but, tklabel(Frame.title1.para2A, text = "SAQC database:", font = fontTextLabel, bg = color), sticky = "w") # title

tkgrid.configure(Frame.title1.para2A, sticky = "w")
#tkgrid(Frame.title1.para2A)

#-------------------------------# 2-2-1-1 Study population #----------------------------#
Frame.para2Aa <- tkframe(Frame.para2A, bg = color)
nblank2 = 6
blank2.space <- tklabel(Frame.para2Aa, bg = color, text = "  ", font = fontsubtitle, width = nblank2)
QIfigure_group.pp1 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_group.pp2 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_group.pp3 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_group.pp4 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_group.pp5 <- tkradiobutton(Frame.para2Aa, bg = color)

QIfigure_groupValue <- tclVar("5") #default Combine
tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, value = "1")
tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, value = "2")
tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, value = "3")
tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, value = "4")
tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, value = "5")

QIfigure_group.Left <- tklabel(Frame.para2Aa, text = "                      Study population:  ", font = fontTextLabel, bg = color)
QIfigure_group1 <- tklabel(Frame.para2Aa, text = "Asia (CHB+JPT) ", width = 13, font = fontTextLabel, bg = color)
QIfigure_group2 <- tklabel(Frame.para2Aa, text = "Africa (YRI)", width = 9, font = fontTextLabel, bg = color)
QIfigure_group3 <- tklabel(Frame.para2Aa, text = "Europe (CEU)", width = 10, font = fontTextLabel, bg = color)
QIfigure_group4 <- tklabel(Frame.para2Aa, text = "Taiwan (TWN)", width = 12, font = fontTextLabel, bg = color)
QIfigure_group5 <- tklabel(Frame.para2Aa, text = "Combine", width = 6, font = fontTextLabel, bg = color)

QIfigure_group.space <- tklabel(Frame.para2Aa, text = "  ", font = fontTextLabel, bg = color)
tkgrid(QIfigure_group.Left, QIfigure_group.pp1, QIfigure_group1, QIfigure_group.pp2, QIfigure_group2, 
	QIfigure_group.pp3, QIfigure_group3, QIfigure_group.pp4, QIfigure_group4, QIfigure_group.pp5, QIfigure_group5, sticky = "w")

#-------------------------------# 2-2-1-2 Genome-wide SNP array #----------------------------#
QIfigure_chip.port1 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_chip.port2 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_chip.port3 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_chip.port4 <- tkradiobutton(Frame.para2Aa, bg = color)
QIfigure_chip.port5 <- tkradiobutton(Frame.para2Aa, bg = color)

QIfigure_chipValue <- tclVar("1") #default Affymetrix 100K
tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, value = "1")
tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, value = "2")
tkconfigure(QIfigure_chip.port3, variable = QIfigure_chipValue, value = "3")
tkconfigure(QIfigure_chip.port4, variable = QIfigure_chipValue, value = "4")
tkconfigure(QIfigure_chip.port5, variable = QIfigure_chipValue, value = "5")

QIfigure_chip.lab <- tklabel(Frame.para2Aa, text = "                      Genome-wide SNP array:", font = fontTextLabel, bg = color)
QIfigure_chip.affy100 <- tklabel(Frame.para2Aa, text = "Affymetrix 100K", font = fontTextLabel, bg = color)
QIfigure_chip.affy500 <- tklabel(Frame.para2Aa, text = "Affymetrix 500K", font = fontTextLabel, bg = color)
QIfigure_chip.array6 <- tklabel(Frame.para2Aa, text = "Affymetrix 6.0", font = fontTextLabel, bg = color)
QIfigure_chip.axiom <- tklabel(Frame.para2Aa, text = "Affymetrix Axiom", font = fontTextLabel, bg = color)
QIfigure_chip.illumina <- tklabel(Frame.para2Aa, text = "Illumina 550K", font = fontTextLabel, bg = color)
tkgrid.configure(Frame.para2Aa, sticky = "w")
tkgrid(QIfigure_chip.lab, QIfigure_chip.port1, QIfigure_chip.affy100, QIfigure_chip.port2, QIfigure_chip.affy500, QIfigure_chip.port3, QIfigure_chip.array6, QIfigure_chip.port4, QIfigure_chip.axiom, QIfigure_chip.port5, QIfigure_chip.illumina, sticky = "w")

#-------------------------------# 2-2-1-3 QI upper quantile #----------------------------#
QIfigure_quantile.port1 <- tkcheckbutton(Frame.para2Aa, bg = color, command = function(){
	if(tclvalue(QIfigure_quantile.Value1) == "0" & tclvalue(QIfigure_quantile.Value4) == "1") tclvalue(QIfigure_quantile.Value4) = 0
	if(tclvalue(QIfigure_quantile.Value1) == "1" & tclvalue(QIfigure_quantile.Value2) == "1" & tclvalue(QIfigure_quantile.Value3) == "1") tclvalue(QIfigure_quantile.Value4) = 1
})

QIfigure_quantile.port2 <- tkcheckbutton(Frame.para2Aa, bg = color, command = function(){
	if(tclvalue(QIfigure_quantile.Value2) == "0" & tclvalue(QIfigure_quantile.Value4) == "1") tclvalue(QIfigure_quantile.Value4) = 0
	if(tclvalue(QIfigure_quantile.Value1) == "1" & tclvalue(QIfigure_quantile.Value2) == "1" & tclvalue(QIfigure_quantile.Value3) == "1") tclvalue(QIfigure_quantile.Value4) = 1
})

QIfigure_quantile.port3 <- tkcheckbutton(Frame.para2Aa, bg = color, command = function(){ 
	if(tclvalue(QIfigure_quantile.Value3) == "0" & tclvalue(QIfigure_quantile.Value4) == "1") tclvalue(QIfigure_quantile.Value4) = 0
	if(tclvalue(QIfigure_quantile.Value1) == "1" & tclvalue(QIfigure_quantile.Value2) == "1" & tclvalue(QIfigure_quantile.Value3) == "1") tclvalue(QIfigure_quantile.Value4) = 1
})

QIfigure_quantile.port4 <- tkcheckbutton(Frame.para2Aa, bg = color, command = function()
	if(tclvalue(QIfigure_quantile.Value4) == "1"){
		tclvalue(QIfigure_quantile.Value1) = 1; tclvalue(QIfigure_quantile.Value2) = 1; tclvalue(QIfigure_quantile.Value3) = 1
	} else {
		tclvalue(QIfigure_quantile.Value1) = 0; tclvalue(QIfigure_quantile.Value2) = 0; tclvalue(QIfigure_quantile.Value3) = 0
	}
)

QIfigure_quantile.Value1 <- tclVar("1") #default all selected
QIfigure_quantile.Value2 <- tclVar("1")
QIfigure_quantile.Value3 <- tclVar("1")
QIfigure_quantile.Value4 <- tclVar("1")

tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1)
tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2)
tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3)
tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4)

QIfigure_quantile.lab <- tklabel(Frame.para2Aa, text = "                      QI upper quantile:", font = fontTextLabel, bg = color)
QIfigure_q95 <- tklabel(Frame.para2Aa, text = "95%", font = fontTextLabel, bg = color)
QIfigure_q975 <- tklabel(Frame.para2Aa, text = "97.5%", font = fontTextLabel, bg = color)
QIfigure_q99 <- tklabel(Frame.para2Aa, text = "99%", font = fontTextLabel, bg = color)
QIfigure_qAll <- tklabel(Frame.para2Aa, text = "All", font = fontTextLabel, bg = color)

#tkgrid.configure(Frame.para2Aa, sticky = "w")
tkgrid(QIfigure_quantile.lab, QIfigure_quantile.port1, QIfigure_q95, QIfigure_quantile.port2, QIfigure_q975, 
	QIfigure_quantile.port3, QIfigure_q99, QIfigure_quantile.port4, QIfigure_qAll, sticky = "w")#chip.port3, chip.affy6, chip.port4, chip.illumina550, sticky = "w")
tkgrid.configure(Frame.para2Aa, sticky = "w")
#tkgrid(Frame.para2Aa)

#-------------------------------# 2-2-2 User provided #----------------------------#
Frame.title2.para2A <- tkframe(Frame.para2A, bg = color)
blank1.space <- tklabel(Frame.title2.para2A, bg = color, text = "  ", font = fontsubtitle, width = nblank1)
userprovided.but <- tkradiobutton(Frame.title2.para2A, bg = color, command = function(){
	#tclvalue(QI_resouce1.val) = 0;tclvalue(QI_resouce2.val) = 1
	interactive_polygon(index.decide.ref = T, button.stat1 = "disable", button.stat2 = "normal",  
								default.val2.1 = 1, default.val2.2 = 2, default.val2.3 = 3)
	if(0){
		tkconfigure(QIfigure_group.pp1, variable = QIfigure_groupValue, stat = "disable")
		tkconfigure(QIfigure_group.pp2, variable = QIfigure_groupValue, stat = "disable")
		tkconfigure(QIfigure_group.pp3, variable = QIfigure_groupValue, stat = "disable")
		tkconfigure(QIfigure_group.pp4, variable = QIfigure_groupValue, stat = "disable")
		tkconfigure(QIfigure_group.pp5, variable = QIfigure_groupValue, stat = "disable")
		tclvalue(QIfigure_groupValue) = 0
		tkconfigure(QIfigure_chip.port1, variable = QIfigure_chipValue, stat = "disable")
		tkconfigure(QIfigure_chip.port2, variable = QIfigure_chipValue, stat = "disable")
		tclvalue(QIfigure_chipValue) = 0
		tkconfigure(QIfigure_quantile.port1, variable = QIfigure_quantile.Value1, stat = "disable")
		tkconfigure(QIfigure_quantile.port2, variable = QIfigure_quantile.Value2, stat = "disable")
		tkconfigure(QIfigure_quantile.port3, variable = QIfigure_quantile.Value3, stat = "disable")
		tkconfigure(QIfigure_quantile.port4, variable = QIfigure_quantile.Value4, stat = "disable")
		tclvalue(QIfigure_quantile.Value1) = 0;tclvalue(QIfigure_quantile.Value2) = 0
		tclvalue(QIfigure_quantile.Value3) = 0;tclvalue(QIfigure_quantile.Value4) = 0
		##------------------------------------------------------------------------
		tkconfigure(MQI1_chip1.q95.box, state = "normal");tclvalue(MQI1_chip1.q95.val) = "1"
		tkconfigure(MQI1_chip2.q95.box, state = "normal");tclvalue(MQI1_chip2.q95.val) = "1"
		tkconfigure(MQI1_merge.q95.box, state = "normal");tclvalue(MQI1_merge.q95.val) = "1"
		tkconfigure(MQI2_chip1.q95.box, state = "normal");tclvalue(MQI2_chip1.q95.val) = "1"
		tkconfigure(MQI2_chip2.q95.box, state = "normal");tclvalue(MQI2_chip2.q95.val) = "1"
		tkconfigure(MQI2_merge.q95.box, state = "normal");tclvalue(MQI2_merge.q95.val) = "1"
		tkconfigure(WQI1_chip1.q95.box, state = "normal");tclvalue(WQI1_chip1.q95.val) = "1"
		tkconfigure(WQI1_chip2.q95.box, state = "normal");tclvalue(WQI1_chip2.q95.val) = "1"
		tkconfigure(WQI1_merge.q95.box, state = "normal");tclvalue(WQI1_merge.q95.val) = "1"
		tkconfigure(WQI2_chip1.q95.box, state = "normal");tclvalue(WQI2_chip1.q95.val) = "1"
		tkconfigure(WQI2_chip2.q95.box, state = "normal");tclvalue(WQI2_chip2.q95.val) = "1"
		tkconfigure(WQI2_merge.q95.box, state = "normal");tclvalue(WQI2_merge.q95.val) = "1"
		tkconfigure(MQI1_chip1.q975.box, state = "normal");tclvalue(MQI1_chip1.q975.val) = "2"
		tkconfigure(MQI1_chip2.q975.box, state = "normal");tclvalue(MQI1_chip2.q975.val) = "2"
		tkconfigure(MQI1_merge.q975.box, state = "normal");tclvalue(MQI1_merge.q975.val) = "2"
		tkconfigure(MQI2_chip1.q975.box, state = "normal");tclvalue(MQI2_chip1.q975.val) = "2"
		tkconfigure(MQI2_chip2.q975.box, state = "normal");tclvalue(MQI2_chip2.q975.val) = "2"
		tkconfigure(MQI2_merge.q975.box, state = "normal");tclvalue(MQI2_merge.q975.val) = "2"
		tkconfigure(WQI1_chip1.q975.box, state = "normal");tclvalue(WQI1_chip1.q975.val) = "2"
		tkconfigure(WQI1_chip2.q975.box, state = "normal");tclvalue(WQI1_chip2.q975.val) = "2"
		tkconfigure(WQI1_merge.q975.box, state = "normal");tclvalue(WQI1_merge.q975.val) = "2"
		tkconfigure(WQI2_chip1.q975.box, state = "normal");tclvalue(WQI2_chip1.q975.val) = "2"
		tkconfigure(WQI2_chip2.q975.box, state = "normal");tclvalue(WQI2_chip2.q975.val) = "2"
		tkconfigure(WQI2_merge.q975.box, state = "normal");tclvalue(WQI2_merge.q975.val) = "2"
		tkconfigure(MQI1_chip1.q99.box, state = "normal");tclvalue(MQI1_chip1.q99.val) = "3"
		tkconfigure(MQI1_chip2.q99.box, state = "normal");tclvalue(MQI1_chip2.q99.val) = "3"
		tkconfigure(MQI1_merge.q99.box, state = "normal");tclvalue(MQI1_merge.q99.val) = "3"
		tkconfigure(MQI2_chip1.q99.box, state = "normal");tclvalue(MQI2_chip1.q99.val) = "3"
		tkconfigure(MQI2_chip2.q99.box, state = "normal");tclvalue(MQI2_chip2.q99.val) = "3"
		tkconfigure(MQI2_merge.q99.box, state = "normal");tclvalue(MQI2_merge.q99.val) = "3"
		tkconfigure(WQI1_chip1.q99.box, state = "normal");tclvalue(WQI1_chip1.q99.val) = "3"
		tkconfigure(WQI1_chip2.q99.box, state = "normal");tclvalue(WQI1_chip2.q99.val) = "3"
		tkconfigure(WQI1_merge.q99.box, state = "normal");tclvalue(WQI1_merge.q99.val) = "3"
		tkconfigure(WQI2_chip1.q99.box, state = "normal");tclvalue(WQI2_chip1.q99.val) = "3"
		tkconfigure(WQI2_chip2.q99.box, state = "normal");tclvalue(WQI2_chip2.q99.val) = "3"
		tkconfigure(WQI2_merge.q99.box, state = "normal");tclvalue(WQI2_merge.q99.val) = "3"
	}
})
#QI_resouce2.val <- tclVar("0")
tkconfigure(userprovided.but, variable = QI_resouce.val, value = "2") 

tkgrid(blank1.space, userprovided.but, 
		tklabel(Frame.title2.para2A, text = "User provided:", font = fontTextLabel, bg = color), sticky = "w") # title
tkgrid.configure(Frame.title2.para2A, sticky = "w")
#tkgrid(Frame.title2.para2A)
#------------------------------------------------------------------------------------------------------------------#
QIspace = 9
Frame.MQI1 <- tkframe(Frame.para2A, bg = color)
nblank_table = 4
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
MQI1.lab <- tklabel(Frame.MQI1, text = "   Median QI1", width = 14, font = fontTextLabel, bg = color)
MQI2.lab <- tklabel(Frame.MQI1, text = "Median QI2", width = 12, font = fontTextLabel, bg = color)
MQI1_q95.lab <- tklabel(Frame.MQI1, text = "95%", width = QIspace, font = fontTextLabel, bg = color)
MQI1_q975.lab <- tklabel(Frame.MQI1, text = "97.5%", width = QIspace, font = fontTextLabel, bg = color)
MQI1_q99.lab <- tklabel(Frame.MQI1, text = "99%", width = QIspace, font = fontTextLabel, bg = color)
MQI2_q95.lab <- tklabel(Frame.MQI1, text = "95%", width = QIspace, font = fontTextLabel, bg = color)
MQI2_q975.lab <- tklabel(Frame.MQI1, text = "97.5%", width = QIspace, font = fontTextLabel, bg = color)
MQI2_q99.lab <- tklabel(Frame.MQI1, text = "99%", width = QIspace, font = fontTextLabel, bg = color)

tkgrid(blank3.space, MQI1.lab, MQI1_q95.lab, MQI1_q975.lab, MQI1_q99.lab, MQI2.lab, MQI2_q95.lab, MQI2_q975.lab, MQI2_q99.lab, sticky = "w")

#-------------------------------# 2-2-2-1 Median QI1 & QI2 #----------------------------#
#Chip1
#Median QI1
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
MQI1_chip1.lab <- tklabel(Frame.MQI1, text = "    Chip1", width = 12, font = fontTextLabel, bg = color)
MQI1_chip1.q95.val <- tclVar(" "); MQI1_chip1.q975.val <- tclVar(" "); MQI1_chip1.q99.val <- tclVar(" ")
MQI1_chip1.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_chip1.q95.val, state = "disable") #default disable
MQI1_chip1.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_chip1.q975.val, state = "disable")
MQI1_chip1.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_chip1.q99.val, state = "disable")
#Median QI2
MQI2_chip1.lab <- tklabel(Frame.MQI1, text = "Chip1", width = 12, font = fontTextLabel, bg = color)
MQI2_chip1.q95.val <- tclVar(" "); MQI2_chip1.q975.val <- tclVar(" "); MQI2_chip1.q99.val <- tclVar(" ")
MQI2_chip1.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_chip1.q95.val, state = "disable")
MQI2_chip1.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_chip1.q975.val, state = "disable")
MQI2_chip1.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_chip1.q99.val, state = "disable")

tkgrid(blank3.space, MQI1_chip1.lab, MQI1_chip1.q95.box, MQI1_chip1.q975.box, MQI1_chip1.q99.box, 
		MQI2_chip1.lab, MQI2_chip1.q95.box, MQI2_chip1.q975.box, MQI2_chip1.q99.box, sticky = "w")

#Chip2
#Median QI1
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
MQI1_chip2.lab <- tklabel(Frame.MQI1, text = "    Chip2", width = 12, font = fontTextLabel, bg = color)
MQI1_chip2.q95.val <- tclVar(" "); MQI1_chip2.q975.val <- tclVar(" "); MQI1_chip2.q99.val <- tclVar(" ")
MQI1_chip2.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_chip2.q95.val, state = "disable")
MQI1_chip2.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_chip2.q975.val, state = "disable")
MQI1_chip2.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_chip2.q99.val, state = "disable")
#Median QI2
MQI2_chip2.lab <- tklabel(Frame.MQI1, text = "Chip2", width = 12, font = fontTextLabel, bg = color)
MQI2_chip2.q95.val <- tclVar(" "); MQI2_chip2.q975.val <- tclVar(" "); MQI2_chip2.q99.val <- tclVar(" ")
MQI2_chip2.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_chip2.q95.val, state = "disable")
MQI2_chip2.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_chip2.q975.val, state = "disable")
MQI2_chip2.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_chip2.q99.val, state = "disable")

tkgrid(blank3.space, MQI1_chip2.lab, MQI1_chip2.q95.box, MQI1_chip2.q975.box, MQI1_chip2.q99.box, 
		MQI2_chip2.lab, MQI2_chip2.q95.box, MQI2_chip2.q975.box, MQI2_chip2.q99.box, sticky = "w")

#Merge
#Median QI1
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
MQI1_merge.lab <- tklabel(Frame.MQI1, text = "    Merge", width = 12, font = fontTextLabel, bg = color)
MQI1_merge.q95.val <- tclVar(" "); MQI1_merge.q975.val <- tclVar(" "); MQI1_merge.q99.val <- tclVar(" ")
MQI1_merge.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_merge.q95.val, state = "disable")
MQI1_merge.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_merge.q975.val, state = "disable")
MQI1_merge.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI1_merge.q99.val, state = "disable")
#Median QI2
MQI2_merge.lab <- tklabel(Frame.MQI1, text = "Merge", width = 12, font = fontTextLabel, bg = color)
MQI2_merge.q95.val<- tclVar(" ");MQI2_merge.q975.val<- tclVar(" ");MQI2_merge.q99.val<- tclVar(" ")
MQI2_merge.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_merge.q95.val, state = "disable")
MQI2_merge.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_merge.q975.val, state = "disable")
MQI2_merge.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = MQI2_merge.q99.val, state = "disable")

tkgrid(blank3.space, MQI1_merge.lab, MQI1_merge.q95.box, MQI1_merge.q975.box, MQI1_merge.q99.box, 
		MQI2_merge.lab, MQI2_merge.q95.box, MQI2_merge.q975.box, MQI2_merge.q99.box, sticky = "w")

#------------------------------------------------------------------------------------------------------------------#
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
WQI1.lab <- tklabel(Frame.MQI1, text = "     WinsMean QI1", width = 14, font = fontTextLabel, bg = color)
WQI2.lab <- tklabel(Frame.MQI1, text = "WinsMean QI2", width = 12, font = fontTextLabel, bg = color)
WQI1_q95.lab <- tklabel(Frame.MQI1, text = "95%", width = QIspace, font = fontTextLabel, bg = color)
WQI1_q975.lab <- tklabel(Frame.MQI1, text = "97.5%", width = QIspace, font = fontTextLabel, bg = color)
WQI1_q99.lab <- tklabel(Frame.MQI1, text = "99%", width = QIspace, font = fontTextLabel, bg = color)
WQI2_q95.lab <- tklabel(Frame.MQI1, text = "95%", width = QIspace, font = fontTextLabel, bg = color)
WQI2_q975.lab <- tklabel(Frame.MQI1, text = "97.5%", width = QIspace, font = fontTextLabel, bg = color)
WQI2_q99.lab <- tklabel(Frame.MQI1, text = "99%", width = QIspace, font = fontTextLabel, bg = color)
tkgrid(blank3.space, WQI1.lab, WQI1_q95.lab, WQI1_q975.lab, WQI1_q99.lab, 
		WQI2.lab, WQI2_q95.lab, WQI2_q975.lab, WQI2_q99.lab, sticky = "w")

#-------------------------------# 2-2-2-2 WinsMean QI1 & QI2 #----------------------------#
#Chip1
#WinsMean QI1
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
WQI1_chip1.lab <- tklabel(Frame.MQI1, text = "    Chip1", width = 12, font = fontTextLabel, bg = color)
WQI1_chip1.q95.val <- tclVar(" "); WQI1_chip1.q975.val <- tclVar(" "); WQI1_chip1.q99.val <- tclVar(" ")
WQI1_chip1.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_chip1.q95.val, state = "disable")
WQI1_chip1.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_chip1.q975.val, state = "disable")
WQI1_chip1.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_chip1.q99.val, state = "disable")
#WinsMean QI2
WQI2_chip1.lab <- tklabel(Frame.MQI1, text = "Chip1", width = 12, font = fontTextLabel, bg = color)
WQI2_chip1.q95.val <- tclVar(" "); WQI2_chip1.q975.val <- tclVar(" "); WQI2_chip1.q99.val <- tclVar(" ")
WQI2_chip1.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_chip1.q95.val, state = "disable")
WQI2_chip1.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_chip1.q975.val, state = "disable")
WQI2_chip1.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_chip1.q99.val, state = "disable")

tkgrid(blank3.space, WQI1_chip1.lab, WQI1_chip1.q95.box, WQI1_chip1.q975.box, WQI1_chip1.q99.box, 
		WQI2_chip1.lab, WQI2_chip1.q95.box, WQI2_chip1.q975.box, WQI2_chip1.q99.box, sticky = "w")

#Chip2
#WinsMean QI1
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
WQI1_chip2.lab <- tklabel(Frame.MQI1, text = "    Chip2", width = 12, font = fontTextLabel, bg = color)
WQI1_chip2.q95.val <- tclVar(" "); WQI1_chip2.q975.val <- tclVar(" "); WQI1_chip2.q99.val <- tclVar(" ")
WQI1_chip2.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_chip2.q95.val, state = "disable")
WQI1_chip2.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_chip2.q975.val, state = "disable")
WQI1_chip2.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_chip2.q99.val, state = "disable")
#WinsMean QI2
WQI2_chip2.lab <- tklabel(Frame.MQI1, text = "Chip2", width = 12, font = fontTextLabel, bg = color)
WQI2_chip2.q95.val <- tclVar(" "); WQI2_chip2.q975.val <- tclVar(" "); WQI2_chip2.q99.val <- tclVar(" ")
WQI2_chip2.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_chip2.q95.val, state = "disable")
WQI2_chip2.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_chip2.q975.val, state = "disable")
WQI2_chip2.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_chip2.q99.val, state = "disable")

tkgrid(blank3.space, WQI1_chip2.lab, WQI1_chip2.q95.box, WQI1_chip2.q975.box, WQI1_chip2.q99.box, 
		WQI2_chip2.lab, WQI2_chip2.q95.box, WQI2_chip2.q975.box, WQI2_chip2.q99.box, sticky = "w")

#Merge
#WinsMean QI1
blank3.space <- tklabel(Frame.MQI1, bg = color, text = "  ", font = fontsubtitle, width = nblank_table)
WQI1_merge.lab <- tklabel(Frame.MQI1, text = "    Merge", width = 12, font = fontTextLabel, bg = color)
WQI1_merge.q95.val <- tclVar(" "); WQI1_merge.q975.val <- tclVar(" "); WQI1_merge.q99.val <- tclVar(" ")
WQI1_merge.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_merge.q95.val, state = "disable")
WQI1_merge.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_merge.q975.val, state = "disable")
WQI1_merge.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI1_merge.q99.val, state = "disable")
#WinsMean QI2
WQI2_merge.lab <- tklabel(Frame.MQI1, text = "Merge", width = 12, font = fontTextLabel, bg = color)
WQI2_merge.q95.val <- tclVar(" "); WQI2_merge.q975.val <- tclVar(" "); WQI2_merge.q99.val <- tclVar(" ")
WQI2_merge.q95.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_merge.q95.val, state = "disable")
WQI2_merge.q975.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_merge.q975.val, state = "disable")
WQI2_merge.q99.box <- tkentry(Frame.MQI1, width = QIspace, textvariable = WQI2_merge.q99.val, state = "disable")

tkgrid(blank3.space, WQI1_merge.lab, WQI1_merge.q95.box, WQI1_merge.q975.box, WQI1_merge.q99.box, 
		WQI2_merge.lab, WQI2_merge.q95.box, WQI2_merge.q975.box, WQI2_merge.q99.box, sticky = "w")
tkgrid.configure(Frame.MQI1, sticky = "w")
#tkgrid(Frame.MQI1)

Space.frame <- tkframe(Frame, bg = color)
tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontsubtitle, width = 50), sticky = "w")
# tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontsubtitle, width = 50), sticky = "w")
#tkgrid.configure(Space.frame, sticky = "w")

##-----------------------------------------------------------------------------
					##Functions of interactive gui
##-----------------------------------------------------------------------------
SAQC.QIinteractive = paste(pathway_program, "/HeatMap.r", sep = "")

##-----------------------------------------------------------------------------
								##Run botton
##-----------------------------------------------------------------------------
Runframe <- tkframe(Frame, bg = color)
tkgrid(blank1.space, sticky = "w")
Run.button <- tkbutton(Runframe, text = "   Run   ", width = button.space, bg = color.run, font = fontRun, command = function(){
	cat("Please wait a moment...\n")
	source(SAQC.QIinteractive)
})

blank1.space <- tklabel(Runframe, bg = color, text = "  ", font = fontTextLabel, width = blankspcae.run)
tkgrid(blank1.space, Run.button, sticky = "w")
tkgrid(tklabel(Runframe, text = "", bg = color), sticky = "w")
tkgrid(tklabel(Runframe, text = "", bg = color), sticky = "w")
tkgrid.configure(Runframe, sticky = "w")

Space.frame <- tkframe(Frame, bg = color)
tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontsubtitle, width = 50), sticky = "w")
tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontsubtitle, width = 50), sticky = "w")
tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontsubtitle, width = 50), sticky = "w")
tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontTextLabel, width = 50), sticky = "w")
tkgrid(tklabel(Space.frame, bg = color, text = "  ", font = fontTextLabel, width = 50), sticky = "w")
tkgrid.configure(Space.frame, sticky = "w")
tkgrid(Frame)
##------------------------------------------------------------------------------
tk2notetab.select(nb, "Main function")# default "Main functions"
tk2notetab.text(nb) # Text of the currently selected tab
tkgrid(noteframe, sticky = "w")

