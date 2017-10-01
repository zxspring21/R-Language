pathway_program = strsplit(ALOHA.gui, "ALOHA_interface.r"); pathway_program = unlist(pathway_program)
pathway = strsplit(pathway_program, "PROGRAM"); pathway = unlist(pathway)[1]
setwd(pathway)
library(tcltk)
button.disable = function(x) for(i in x) tkconfigure(i, state = "disable"); 
button.enable = function(x) for(i in x) tkconfigure(i, state = "normal"); 
set.default1 = function(){
	tclvalue(database.val) = 0; 
	tclvalue(userprovided.val) = 0; 
	tclvalue(pathway_AFrefinput) = ""; 
	tclvalue(groupValue) = 0; 
	tclvalue(chipValue) = 0; 
	tclvalue(Alphascaling.val) = ""; 
	tclvalue(windowsize.val) = ""; 
	tclvalue(qRef.val) = ""; 
}
set.default2 = function(){
	tclvalue(database.val) = 1;
	tclvalue(userprovided.val) = 0;
	tclvalue(pathway_AFrefinput) = "";
	tclvalue(groupValue) = 5;
	tclvalue(chipValue) = 1;
	tclvalue(Alphascaling.val) = 0.95;
	tclvalue(windowsize.val) = 100;
	tclvalue(qRef.val) = 0.95;
}
color="palegreen"
QI_interface = tktoplevel(bg = color)
fontHeading = tkfont.create(family = "arial", size = 14, weight = "bold")
fontIntro = tkfont.create(family = "arial", size = 10, weight = "bold")
fontIntro_ALOHA = tkfont.create(family = "times", size = 12, weight = "bold", slant = "italic")
fontsubtitle = tkfont.create(family = "times", size = 14, weight = "bold")
fontTextLabel = tkfont.create(family = "times", size = 10)
tktitle(QI_interface) = "ALOHA v1.0"
QI_interface_A = tkframe(QI_interface, background = color) 
QI_interface_B = tkframe(QI_interface, background = color)
a = tklabel(QI_interface_A, text = "       ALOHA  ", font = fontHeading, bg = color)
b = tklabel(QI_interface_A, 
	text = " (Allele-frequency/Loss-of-heterozygosity/Allele-imbalance; AF/LOH/AI),  written in R and R GUI,            ",
	font = fontIntro, bg = color)
tkgrid(a, b)
tkgrid(tklabel(QI_interface_B, 
	text = "        provides for genome-wide analysis of allele frequency and detection of both loss of heterozygosity (LOH) and            ",
	font = fontIntro, bg = color), sticky = "w")
tkgrid(tklabel(QI_interface_B, 
	text = "        allelic imbalance (AI). An allele frequency biplot is also provided for sample classification,  outlier detection,        ",
	font = fontIntro, bg = color), sticky = "w")
tkgrid(tklabel(QI_interface_B, text = "        and SNP clustering.", font = fontIntro, bg = color), sticky = "w")
tkgrid(tklabel(QI_interface_B, text = " ", font = fontIntro, bg = color))
tkgrid(QI_interface_A, sticky = "w")
tkgrid(QI_interface_B)
Frame = tkframe(QI_interface, bg = color)     
frame.para1 = tkframe(Frame, bg = color)
frame.path = tkframe(Frame, bg = color)
frame.dataformat = tkframe(Frame, bg = color)
frame.dataformat_database = tkframe(Frame, bg = color)
frame.dataformat_A = tkframe(Frame, bg = color)
frame.dataformat_userprovided = tkframe(Frame, bg = color)
frame.AFrefpath = tkframe(Frame, bg = color)
frame.AILOH = tkframe(Frame, bg = color)
frame.Biplot = tkframe(Frame, bg = color)
tkgrid(tklabel(frame.para1, text = "     1. Input/Output path:", font = fontsubtitle, bg = color)) 
inputdata_format.port1 = tkradiobutton(frame.para1, bg = color, command = function(){
	tclvalue(pathway_input) = "Test2"; 
	button.disable(list(database.port, userprovided.port, entry.AFrefPath, group.pp1, group.pp2, group.pp3, group.pp4, group.pp5, chip.port1,
		chip.port2, chip.port3, Alphascaling.box, windowsize.box, qRef.box))
	set.default1()
})
inputdata_format.port2 = tkradiobutton(frame.para1, bg = color, command = function(){
	tclvalue(pathway_input) = "Test1";
	button.enable(list(database.port, userprovided.port, entry.AFrefPath, group.pp1, group.pp2, group.pp3, group.pp4, group.pp5, chip.port1,
		chip.port2, chip.port3, Alphascaling.box, windowsize.box, qRef.box))
	set.default2()
})
inputdata_format.statsValue = tclVar("2")
tkconfigure(inputdata_format.port1, variable = inputdata_format.statsValue, value = "1") 
tkconfigure(inputdata_format.port2, variable = inputdata_format.statsValue, value = "2")
inputdata_format.stats.Label = tklabel(frame.para1, bg = color, text = "               Study group:         ", font = fontTextLabel)
inputdata_format.LeftLabel = tklabel(frame.para1, bg = color, text = "One group (One or multiple populations)", font = fontTextLabel)
inputdata_format.RightLabel = tklabel(frame.para1, bg = color, text = "Two groups (Case/Control)", font = fontTextLabel)
tkgrid(inputdata_format.stats.Label, inputdata_format.port1, inputdata_format.LeftLabel,
		inputdata_format.port2, inputdata_format.RightLabel, sticky = "w")
pathway_input = tclVar("Test1") 
pathlab = tklabel(frame.path, text = "          Directory of data input:", font = fontTextLabel, bg = color)
entry.Path = tkentry(frame.path, width = 85, textvariable = pathway_input, font = fontTextLabel) 
box.input = tkbutton(frame.path, text = "...",  command = function() tclvalue(pathway_input) = tkchooseDirectory()) 
tkgrid(pathlab, entry.Path, box.input)
pathway_output = tclVar("")
pathlab = tklabel(frame.path, text = "               Directory of result output:", font = fontTextLabel, bg = color)
entry.Path = tkentry(frame.path, width = 85, textvariable = pathway_output, font = fontTextLabel)
box.output = tkbutton(frame.path, text = "...",  command = function() tclvalue(pathway_output) = tkchooseDirectory())
tkgrid(pathlab, entry.Path, box.output)
tkgrid(tklabel(frame.dataformat, text = "     2. AF reference:", font = fontsubtitle, bg = color)) 
database.space = tklabel(frame.dataformat_database, bg = color, text = "             ")
database.port = tkcheckbutton(frame.dataformat_database, bg = color, command = function()
{
	tclvalue(database.val) = 1
	tclvalue(userprovided.val) = 0
	button.enable(list(group.pp1, group.pp2, group.pp3, group.pp4, group.pp5, chip.port1, chip.port2, chip.port3))
	tclvalue(groupValue) = 5;
	tclvalue(chipValue) = 1;
	tkconfigure(entry.AFrefPath, state = "disable");
	tclvalue(pathway_AFrefinput) = "";
})
database.val = tclVar("1")
database.lab = tklabel(frame.dataformat_database, bg = color, text = " ALOHA database -")
tkconfigure(database.port, variable = database.val)
tkgrid(database.space, database.port, database.lab, sticky = "w")
group.pp1 = tkradiobutton(frame.dataformat_A, bg = color)
group.pp2 = tkradiobutton(frame.dataformat_A, bg = color)
group.pp3 = tkradiobutton(frame.dataformat_A, bg = color)
group.pp4 = tkradiobutton(frame.dataformat_A, bg = color)
group.pp5 = tkradiobutton(frame.dataformat_A, bg = color)
groupValue = tclVar("5")
tkconfigure(group.pp1, variable = groupValue, value = "1")
tkconfigure(group.pp2, variable = groupValue, value = "2")
tkconfigure(group.pp3, variable = groupValue, value = "3")
tkconfigure(group.pp4, variable = groupValue, value = "4")
tkconfigure(group.pp5, variable = groupValue, value = "5")
group1 = tklabel(frame.dataformat_A, text = "Taiwan (TWN)", font = fontTextLabel, bg = color)
group2 = tklabel(frame.dataformat_A, text = "Asia (CHB+JPT) ", font = fontTextLabel, bg = color)
group3 = tklabel(frame.dataformat_A, text = "Africa (YRI)", font = fontTextLabel, bg = color)
group4 = tklabel(frame.dataformat_A, text = "Europe (CEU)", font = fontTextLabel, bg = color)
group5 = tklabel(frame.dataformat_A, text = "Combined", font = fontTextLabel, bg = color)
group.Left = tklabel(frame.dataformat_A, text = "                         Study population:  ", font = fontTextLabel, bg = color)
tkgrid(group.Left, group.pp2, group2, group.pp3, group3, group.pp4, group4, group.pp1, group1, group.pp5, group5, sticky = "w")
chip.port1 = tkradiobutton(frame.dataformat_A, bg = color)
chip.port2 = tkradiobutton(frame.dataformat_A, bg = color)
chip.port3 = tkradiobutton(frame.dataformat_A, bg = color)
chip.port4 = tkradiobutton(frame.dataformat_A, bg = color)
chipValue = tclVar("1")
tkconfigure(chip.port1, variable = chipValue, value = "1")
tkconfigure(chip.port2, variable = chipValue, value = "2")
tkconfigure(chip.port3, variable = chipValue, value = "3")
tkconfigure(chip.port4, variable = chipValue, value = "4")
chip.lab = tklabel(frame.dataformat_A, text = "                         Genome-wide SNP chip:", font = fontTextLabel, bg = color)
chip.affy100 = tklabel(frame.dataformat_A, text = "Affymetrix 100K", font = fontTextLabel, bg = color)
chip.affy500 = tklabel(frame.dataformat_A, text = "Affymetrix 500K", font = fontTextLabel, bg = color)
tkgrid(chip.lab, chip.port1, chip.affy100, chip.port2, chip.affy500, sticky = "w")
userprovided.space = tklabel(frame.dataformat_userprovided, bg = color, text = "             ")
userprovided.port = tkcheckbutton(frame.dataformat_userprovided, bg = color, command = function()
{
	tclvalue(userprovided.val) = 1 
	tclvalue(database.val) = 0
	button.disable(list(group.pp1, group.pp2, group.pp3, group.pp4, group.pp5, chip.port1, chip.port2, chip.port3))
	tclvalue(groupValue) = 0;
	tclvalue(chipValue) = 0;
	tkconfigure(entry.AFrefPath, state = "normal");
	tclvalue(pathway_AFrefinput) = "";
})
userprovided.val = tclVar("0")
userprovided.lab = tklabel(frame.dataformat_userprovided, bg = color, text = " User provided -")
tkconfigure(userprovided.port, variable = userprovided.val)
tkgrid(userprovided.space, userprovided.port, userprovided.lab, sticky = "w")
pathway_AFrefinput=tclVar("")
box.AFrefinput = tkbutton(frame.AFrefpath, text = "...",  command = function() tclvalue(pathway_AFrefinput) = tkchooseDirectory())
entry.AFrefPath=tkentry(frame.AFrefpath, width = "95", textvariable = pathway_AFrefinput, font = fontTextLabel)
AFrefpathlab=tklabel(frame.AFrefpath, text = "                         Directory:", font = fontTextLabel, bg = color)
tkgrid(AFrefpathlab, entry.AFrefPath, box.AFrefinput)
tkgrid(tklabel(frame.AILOH, text = "     3. AI/LOH calculation:", font = fontsubtitle, bg = color), sticky = "w") 
Alphascaling.val= tclVar("0.95")
Alphascaling.Label = tklabel(frame.AILOH, bg = color, text = "               Confidence level :  ")
Alphascaling.box = tkentry(frame.AILOH, width = 6, textvariable = Alphascaling.val)
Alphascaling.ps = tklabel(frame.AILOH, bg = color, text = "    ( between 0.9 and 1 )")
tkgrid(Alphascaling.Label, Alphascaling.box, Alphascaling.ps, sticky = "w")
windowsize.val= tclVar("100")
windowsize.Label = tklabel(frame.AILOH, bg = color, text = "               Window size :  ")
windowsize.box = tkentry(frame.AILOH, width = 6, textvariable = windowsize.val)
windowsize.ps = tklabel(frame.AILOH, bg = color, text = "    ( at least 10 )")
tkgrid(windowsize.Label, windowsize.box, windowsize.ps, sticky = "w")
qRef.val= tclVar("0.95")
qRef.Label = tklabel(frame.AILOH, bg = color, text = "               Upper bound (Quantile) of reference :  ")
qRef.box = tkentry(frame.AILOH, width = 6, textvariable = qRef.val)
qRef.ps = tklabel(frame.AILOH, bg = color, text = "    ( between 0.9 and 1 )")
tkgrid(qRef.Label, qRef.box, qRef.ps, sticky = "w")
tkgrid(tklabel(frame.Biplot, text = "     4. AF biplot:", font = fontsubtitle, bg = color),sticky = "w") 
biplot.Alphascaling.val= tclVar("1")
biplot.Alphascaling.Label = tklabel(frame.Biplot, bg = color, text = "               Alpha scaling :  ")
biplot.Alphascaling.box = tkentry(frame.Biplot, width = 6, textvariable = biplot.Alphascaling.val)
biplot.Alphascaling.ps = tklabel(frame.Biplot, bg = color, text = "    (  0 or 1 )")
tkgrid(biplot.Alphascaling.Label, biplot.Alphascaling.box, biplot.Alphascaling.ps, sticky = "w")
tkgrid(frame.para1, sticky = "w")
tkgrid(frame.path, sticky = "w")
tkgrid(frame.dataformat, sticky = "w")
tkgrid(frame.dataformat_database, sticky = "w")
tkgrid(frame.dataformat_A, sticky = "w")
tkgrid(frame.dataformat_userprovided, sticky = "w")
tkgrid(frame.AFrefpath, sticky = "w")
tkgrid.configure(frame.AILOH, sticky = "w")
tkgrid.configure(frame.Biplot, sticky = "w")
tkgrid(Frame, sticky = "w")
tkgrid(tklabel(QI_interface, text = "", bg = color))
ALOHA.para = paste(pathway_program, "ALOHA_parameter_setting.r", sep = "")
Run.button = tkbutton(QI_interface, text = "  Run  ", command = function(){cat("Please wait a while,  ALOHA is running...\n");source(ALOHA.para)})
tkgrid(Run.button)
tkgrid(tklabel(QI_interface, text = "", bg = color))
