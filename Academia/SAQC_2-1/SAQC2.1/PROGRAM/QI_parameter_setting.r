# pathway_subfunction = paste(pathway_program, "SAQC_subfunction.r", sep = "")
# source(pathway_subfunction)

# pathway = strsplit(pathway_program, "PROGRAM")
# pathway = unlist(pathway)[1]

# SAQC.file <- paste(pathway_program, "SAQC.r", sep = "")  ## Used in SAQC.r
# SAQC.para <- paste(pathway_program, "parameter_setting.r", sep = "")

if(exists("QI_para")) tkdestroy(QI_para)
QI_para = tktoplevel()
tktitle(QI_para) = "Import DATA Setting"
tkgrid(tklabel(QI_para, text = "   Parameter Setting", width = 14, font = fontHeading), sticky = "w")
tkgrid(tklabel(QI_para, text = "", font = fontIntro_para, height = 0), sticky = "w")
###### Skiprow ######
skip.lab = tklabel(QI_para, text = "   Number of header lines to skip: ", font = fontTextLabel_para)
Skiprow = tclVar("")
skip_entry = tkentry(QI_para, width = "5", textvariable = Skiprow, font = fontTextLabel_para)	

###### NA string ######
na.lab <- tklabel(QI_para, text = "   Missing value: ", font = fontTextLabel_para)
NA_str = tclVar("")
NA_entry = tkentry(QI_para, width = "5", textvariable = NA_str, font = fontTextLabel_para)

###### Probe set ######
probe.lab = tklabel(QI_para, text = "   Column index of SNP ID: ", font = fontTextLabel_para)
index.probe = tclVar("")
probe_entry = tkentry(QI_para, width = "5", textvariable = index.probe, font = fontTextLabel_para)

###### Phy posi ######
posi.lab = tklabel(QI_para, text = "   Column index of physical position: ", font = fontTextLabel_para)
index.posi = tclVar("")
posi_entry = tkentry(QI_para, width = "5" ,textvariable = index.posi, font = fontTextLabel_para)

###### Chromosome ######
chr.lab = tklabel(QI_para, text = "   Column index of chromosome: ", font = fontTextLabel_para)
index.chr = tclVar("")
chr_entry = tkentry(QI_para, width = "5", textvariable = index.chr, font = fontTextLabel_para)

###### Genotype of allele A ######
CallA.lab = tklabel(QI_para, text = "   Column index of allele A: ", font = fontTextLabel_para)
index.CallA = tclVar("")
CallA_entry = tkentry(QI_para, width = "5", textvariable = index.CallA, font = fontTextLabel_para)

###### Genotype of allele B ######
CallB.lab = tklabel(QI_para, text = "   Column index of allele B: ", font = fontTextLabel_para)
index.CallB = tclVar("")
CallB_entry = tkentry(QI_para, width = "5", textvariable = index.CallB, font = fontTextLabel_para)

###### AF of allele B ######
BAF.lab = tklabel(QI_para, text = "   Column index of AF of allele B: ", font = fontTextLabel_para)
index.BAF = tclVar("")
BAF_entry = tkentry(QI_para, width = "5", textvariable = index.BAF, font = fontTextLabel_para)	

###### PI_A ######
PI_A.lab = tklabel(QI_para, text = "   Column index of raw intensity of X: ", font = fontTextLabel_para)
index.PI_A = tclVar("")
PI_A_entry = tkentry(QI_para, width = "5", textvariable = index.PI_A, font = fontTextLabel_para)

###### PI_B ######
PI_B.lab = tklabel(QI_para, text = "   Column index of raw intensity of Y: ", font = fontTextLabel_para)
index.PI_B = tclVar("")
PI_B_entry = tkentry(QI_para, width = "5", textvariable = index.PI_B, font = fontTextLabel_para)

###### default setting for Test example 1 ######

if(tclvalue(pathway_input)=="Test example 1"){
	tclvalue(Skiprow) = 2
	tclvalue(NA_str) = "null"
	tclvalue(index.probe) = 2
	tclvalue(index.chr) = 3
	tclvalue(index.posi) = 4
	tclvalue(index.CallA) = 6
}else if(tclvalue(pathway_input)=="Test example 3"){
	tclvalue(Skiprow) = 11
	tclvalue(NA_str) = "NaN,-"
	tclvalue(index.probe) = 1
	tclvalue(index.chr) = 2
	tclvalue(index.posi) = 3
	tclvalue(index.CallA) = 10
	tclvalue(index.CallB) = 11
	tclvalue(index.BAF) = 4
	tclvalue(index.PI_A) = 6
	tclvalue(index.PI_B) = 7
	
}

####### Turn off no use entry respect to different chips #######

if(c(tclvalue(chipValue)) %in% c("1", "2")){
	tkconfigure(CallA.lab, text = "   Column index of allele A: ")
	tkconfigure(CallB_entry, textvariable = index.CallB, state = "disable")
	tkconfigure(BAF_entry, textvariable = index.BAF, state = "disable")
	tkconfigure(PI_A_entry, textvariable = index.PI_A, state = "disable")
	tkconfigure(PI_B_entry, textvariable = index.PI_B, state = "disable")
}else if(c(tclvalue(chipValue)) %in% c("3")){
	tkconfigure(CallA.lab, text = "   Column index of allele A: ")
	tkconfigure(CallB_entry, textvariable = index.CallB, state = "disable")
	tkconfigure(BAF_entry, textvariable = index.BAF, state = "disable")
}else if(c(tclvalue(chipValue))=="4"){
	tkconfigure(posi_entry, textvariable = index.posi, state = "disable")
	tkconfigure(chr_entry, textvariable = index.chr, state = "disable")
	tkconfigure(CallA.lab, text = "   Column index of allele A: ")
	tkconfigure(CallB_entry, textvariable = index.CallB, state = "disable")
	tkconfigure(BAF_entry, textvariable = index.BAF, state = "disable")
}



tkgrid(skip.lab, skip_entry, tklabel(QI_para, text = "   "), sticky = "w")


tkgrid(na.lab, NA_entry, sticky = "w")


tkgrid(probe.lab, probe_entry, sticky = "w")


tkgrid(posi.lab, posi_entry, sticky = "w")


tkgrid(chr.lab, chr_entry, sticky = "w")


tkgrid(CallA.lab, CallA_entry, sticky = "w")


tkgrid(CallB.lab, CallB_entry, sticky = "w")


tkgrid(BAF.lab, BAF_entry, sticky = "w")


tkgrid(PI_A.lab, PI_A_entry, sticky = "w")


tkgrid(PI_B.lab, PI_B_entry, sticky = "w")
tkgrid(tklabel(QI_para, text = ""), sticky = "w")

##-----------------------------------------------------------------------------
								##Run botton
##-----------------------------------------------------------------------------

saveframe = tkframe(QI_para)

save.button <- tkbutton(saveframe, text = "   SAVE   ", width = 8, bg = color.run, font = fontRun, command = function(){
	cat("The input parameters are saved.\n")
	tkdestroy(QI_para)
})
blank1.space <- tklabel(saveframe, text = "  ", font = fontTextLabel_para, width = "12")
tkgrid(blank1.space, save.button, sticky = "w")
tkgrid(tklabel(saveframe, text = ""), sticky = "w")
tkgrid.configure(saveframe, sticky = "w")

