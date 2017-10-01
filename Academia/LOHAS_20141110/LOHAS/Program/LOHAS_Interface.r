###########################
#####      Modified version      #####
##### Independent sample page only (20121030)  #####
##### Partition genotype file before coding, save/rm/load mAFchr to free memory, renew coding method (20130611) #####
##### Implement debug part into "RUN" button (20130730) #####
###################################################################################
##### Input .RData folder to draw plots without run LOHAS again               #####
##### Add interface for variable selection of Regression/GEE model (20130815) #####
##### Add example5 plot example(.RData of example3)                           #####
###################################################################################
##### Remove the "Info" button for variable selection and change              #####
##### into showing the variable selection window as soon as users (20130902)  #####
##### selecting suitable covariate format text file                           #####
###################################################################################
###################################################################################
#####                Combine GUI and Batch version (20140310)                 #####
###################################################################################
###########################
### Running flow: GUI ---> debug ---> LOHAS.main ---> subfunction
### Running flow: GUI ---> LOHAS.main ---> subfunction 20130730
### Load the package ####

##### Code to detect whether batch mode is executed ##### 20130223

libname = .packages(all.available = TRUE)


library(tcltk)
if(!"tcltk2" %in% libname) install.packages("tcltk2", repos = "http://cran.csie.ntu.edu.tw")
library(tcltk2)


# install geeM package for GEE (available for 2.15 or later)
if((!"geeM" %in% libname)) install.packages("geeM", repos = "http://cran.csie.ntu.edu.tw")
library(geeM)


  # addTclPath("C:/Tcl/lib")
  # tclRequire("BWidget")
  # tclRequire("Iwidgets")
  
  ##### background color #####
  color <- colors()[531]
  
  ##### Create a new toplevel window #####
  
  LOH <- tktoplevel(background=color)
  
  ##### Fontsize setting #####
  ##### font selection: "arial", "times", "courier", "helvetica"
  ##### slant selection: "roman", "italic"
  ##### weight selection: "normal", "bold"
  ##### underline selection: TRUE, FALSE
  ##### overstrike selection: TRUE, FALSE
  fonttitle <- tkfont.create(family="arial",size=9)
  fontHeading <- tkfont.create(family="times",size=16,weight="bold")
  
  ##### parameter setting #####
  ##### relief selections:"raised", "sunken", "flat", "ridge", "solid", and "groove" #####
  frame.up <- tkframe(LOH,background=color)
  
  frame.title <- tkframe(frame.up,background=color,relief="ridge",borderwidth=5)
  
  ####################
  ##### function #####
  ####################
  Tkgrid <- function(Tklabel) tkgrid(Tklabel, sticky="w")
  
  Disabled <- function(TEXTNAME){
    for(i in TEXTNAME){
      tkconfigure(get(i), state="disabled")
    }
  }
  Normal <- function(TEXTNAME){
    for(i in TEXTNAME){
      tkconfigure(get(i), state="normal")
    }
  }
  
  
  ##### Title  #####
  
  tktitle(LOH)<- "LOHAS Analysis"
  tkgrid(tklabel(frame.title,text="Welcome to use LOHAS",background=color,font=fontHeading),sticky="w")
  #tkgrid(tklabel(frame.title,text="LOHAS (Loss-Of-Heterozygosity Analysis Suite) is a package designed to detect loss of (genotype) heterozygosity in human genome.",font=fonttitle,background=color),sticky="w")
  #tkgrid(tklabel(frame.title,background=color,text="The phenomenon of LOH can be found in cancer patients (due to DNA deletion, uniparental disomy, mitotic recombination, gene conversionor, etc.) ",font=fonttitle),sticky="w")
  #tkgrid(tklabel(frame.title,background=color,text="or normal individuals (due to consanguineous marriage, etc.). LOHAS provides non-parametric estimates of LOH intensity, chromosome-wise ",font=fonttitle),sticky="w")
  #tkgrid(tklabel(frame.title,background=color,text="LOH tests, and biplot visualization tools for detections of aberrant chromosomes, unusual samples and LOH genomic regions.",font=fonttitle),sticky="w")
  tkgrid(tklabel(frame.title, text = "LOHAS (Loss-of-heterozygosity analysis suite) was developed for studies of homozygosity disequilibrium. LOHAS can be used to identify                    ", font = fonttitle, background = color), sticky = "w")
  tkgrid(tklabel(frame.title, text = "run of homozygosity associated with complex disorders, detect loss of heterozygosity in cancer research, and characterize long contiguous", font = fonttitle, background = color), sticky = "w")
  tkgrid(tklabel(frame.title, text = "stretches of homozygosity in population genetics studies using whole-genome SNP data and DNA sequencing data.", font = fonttitle, background = color), sticky = "w")
  ##tkgrid(tklabel(frame.title,background=color,text=""),sticky="w")
  ##tkgrid(tklabel(frame.title,background=color,text="Reference: Richard M Huggins, Ling-Hui Li, You-Chin Lin, Alice L Yu & Hsin-Chou Yang* (2008). Nonparametric estimation of LOH using Affymetrix",font=fonttitle),sticky="w")
  ##tkgrid(tklabel(frame.title,background=color,text="SNP genotyping arrays for unpaired samples. Journal of Human Genetics 53: 983-990. (*Correspondence to hsinchou@stat.sinica.edu.tw)",font=fonttitle),sticky="w")
  
  ########################################################################################
  ##### Multi pages: page1 ###############################################################
  ########################################################################################
  tn = tkframe(LOH)
  #nb <- tk2notebook(tn, tabs = c("  Independent samples  ", "   Matched samples   "), width = 820, height = 635) #height default: 580, width = 820
  #nb <- tk2notebook(tn, tabs = c("  Independent samples  "), width = 820, height = 635) #height default: 580, width = 820 # 20121030
  #tkpack(nb, fill = "both", expand = 0)
  #tb1 <- tk2notetab(nb, "  Independent samples  ") # no tab 20121030
  frame <- tkframe(tn, bg = color)
  ### change the background color of tabs (http://r.789695.n4.nabble.com/One-application-of-Tcltk2-td2233546.html)
  ### 20111229
  tcl("ttk::style", "configure", "TNotebook", background=color)
  
  # tn <- tkwidget(LOH, "iwidgets::tabnotebook")
  # tkconfigure(tn,tabpos="n",width=820,height=580,
   		 # angle=0,bevelamount=8,gap=2,margin=1,tabborders=2,
  		 # tabbackground=color,background=color,backdrop=color)
  
  # tb1 <- tclvalue(tkadd(tn,label="  Independent samples  "))
  # frame <- .Tk.newwin(tb1)
  
  #################################
  frame.path <- tkframe(frame,background=color)  
  frame.group <- tkframe(frame,background=color)
  frame.group.1 <-tkframe(frame.group,background=color)
  frame.group.1.1 <-tkframe(frame.group.1,background=color)
  #frame.group.1.2<-tkframe(frame.group.1,background=color)
  frame.chip <- tkframe(frame,background=color) 
  frame.chip.1 <- tkframe(frame.chip,background=color)
  frame.chip.2 <- tkframe(frame.chip,background=color)
  
  frame.para <- tkframe(frame,background=color)
  frame.para2 <- tkframe(frame,background=color)
  frame.para2.1 <-tkframe(frame,background=color)
  frame.para3 <- tkframe(frame,background=color)
  frame.para4 <- tkframe(frame,background=color)
  frame.para4.1 <- tkframe(frame.para4,background=color)
  frame.para4.1.1 <- tkframe(frame.para4,background=color)
  frame.para4.1.2 <- tkframe(frame.para4,background=color)
  frame.para4.1.3 <- tkframe(frame.para4,background=color)
  frame.para4.2 <- tkframe(frame.para4,background=color)
  #frame.para4.3 <- tkframe(frame.para4, background=color)
  
  frame.para8 <- tkframe(frame,background=color)
  frame.para5 <- tkframe(frame.para8,background=color)
  frame.para6 <- tkframe(frame.para8,background=color)
  frame.para7 <- tkframe(frame.para8,background=color)
  
  ##### input and output path title #####
  
  font.path<- tkfont.create(family="times",size=12,weight="bold") # text type and size of title
  #tkgrid(tklabel(frame.path,background=color,text=""),sticky="w")
  tkgrid(tklabel(frame.path,background=color,text="1. Input/output path:",font=font.path),sticky="w") # title
  
  ##### input and output box #####
  
  dir.input <- tclVar("")
  box.input <- tkbutton(frame.path,text="...", command=function() tclvalue(dir.input) <- tkchooseDirectory())
  #entry.Path <- tkentry(frame.path, width = "94", textvariable = dir.input)
  # gui parameter setting for examples 20130109
  entry.Path <- tkentry(frame.path, width = "94", textvariable = dir.input, xscrollcommand = function(...){
    plot.only <<- FALSE
    ### example5 is the plot example
    example.groupValue <- c("1", "1", "2", "2", "2")
    example.ethnicity.val <- c("1", "2", "1", "1", "1")
    example.patient1.val <- c("1", "1", "1", "1", "1")
    example.chipValue <- c("3", "3", "1", "4", "4")
    example.chip1.val <- c("0", "0", "1", "0", "0")
    example.chip2.val <- c("0", "0", "1", "0", "0")
    example.chip1_2.val <- c("0", "0", "1", "0", "0")
    #example.chr.val <- c("1", "1", "1", "0", "1", "1", "1", "0")
    example.chr.val <- list(example1 = rep(1, 23), example2 = rep(1, 23), example3 = rep(1, 23), example4 = rep(0, 23),
      example5 = rep(0, 23))
    example.statsValue <- c("1", "1", "1", "1", "1")
    example.LOHintensityValue <- c("1", "1", "1", "1", "1")
    example.SNPthinningvalue <- c("2", "2", "2", "2", "2")
    example.SNPthinningsize.val <- c("2", "2", "2", "2", "2")
    example.SNPselectionValue <- c("2", "2", "2", "2", "2")
    example.fixSNPsize.val <- c("50", "50", "50", "50", "50")
    example.fixSNPprop.val <- c("0.05", "0.05", "0.05", "0.05", "0.05")
    example.chromotestValue <- c("2", "2", "1", "1", "1")
    example.associationcheck <- c("2", "2", "1", "1", "1")
    example.associationValue <- c("0", "0", "1", "1", "1")
    #example.ftestValue <- c("0", "0", "0", "1")
    example.ftestValue <- c("0", "0", "0", "0", "0")
    example.geeValue <- c("0", "0", "0", "0", "0")
    #example.dir.cov <- c("", "", "", paste(gsub("/Program/LOH_gui_new110307.r$", "", gsub("\\\\", "/", LOHgui)), "/Examples/Example4/covariates.txt", sep = ""))
    example.dir.cov <- c("", "", "", "", "")
    example.DataVisualValue <- c("1", "1", "1", "1", "1")
    example.Plotselection1.val <- c("1", "1", "1", "1", "1")
    example.Plotselection2.val <- c("1", "1", "1", "1", "1")
    example.Plotselection3.val <- c("0", "1", "1", "1", "1")
    example.QLOHintensity.val <- c("0.9", "0.9", "0.9", "0.9", "0.9")
    example.Alphascaling.val <- c("0", "0", "0", "0", "0")
    #example.dir.output <- c(paste(gsub(paste("/Program/LOH_gui_new", tempdate, ".r$", sep = ""), "", gsub("\\\\", "/", LOHgui)), "/Output/", rep(c("", "Plot"), c(4, 4)), "Example", rep(1:4, 2), sep = ""))
    example.dir.output <- c(paste(gsub(paste("/Program/LOHAS_Interface.r", sep = ""), "", gsub("\\\\", "/", LOHAS.gui)), "/Output/", "Example", 1:5, sep = ""))
    example.disabled <- c("output.lab", "entry.Path2", "box.output",
      "group.Left", "group.port1", "group.LeftLabel", "ethnicity.Label", "ethnicity.port1", "ethnicity1.Label", "ethnicity.port2", "ethnicity2.Label",
      "group.port2", "group.RightLabel", "patienttitle.Label", "patient.port1", "patient1.Label", "patient.port2", "patient2.Label",
      "user.Label", "User.chip",
      "chip.lab", "chip.port1", "chip.LeftLabel", "chip.port2", "chip.MiddleLabel", "chip.port3", "chip.RightLabel", "chip.Label", "chip1.port", "chip1.Label", "chip2.port", "chip2.Label", "chip1_2.port", "chip1_2.Label",
      "chrselect.name", "chr1.port", "chr1.lab", "chr2.port", "chr2.lab", "chr3.port", "chr3.lab", "chr4.port", "chr4.lab", "chr5.port", "chr5.lab", "chr6.port", "chr6.lab", "chr7.port", "chr7.lab", "chr8.port", "chr8.lab", "chr9.port", "chr9.lab", "chr10.port", "chr10.lab", "chr11.port", "chr11.lab", "chr12.port", "chr12.lab",
      "chrselect.space1", "chr13.port", "chr13.lab", "chr14.port", "chr14.lab", "chr15.port", "chr15.lab", "chr16.port", "chr16.lab", "chr17.port", "chr17.lab", "chr18.port", "chr18.lab", "chr19.port", "chr19.lab", "chr20.port", "chr20.lab", "chr21.port", "chr21.lab", "chr22.port", "chr22.lab", "chrall.port", "chrall.lab",
      "stats.Label", "stats.port1", "stats.LeftLabel", "stats.port2", "stats.RightLabel",
      "LOHintensity.Label", "LOHintensity.port1", "LOHintensity.LeftLabel", "LOHintensity.port2", "LOHintensity.RightLabel",
      "SNP.thinning.space", "SNP.thinning.label", "SNP.thinning.port1", "SNPthinning.LeftLabel", "SNPthinningsize.Label", "SNPthinningsize.box", "SNPthinningsize.Label2", "SNP.thinning.port2", "SNPthinning.RightLabel",
      "SNPselection.Label", "SNPselection.port1", "fixSNPsize.Label", "fixSNPsize.box", "fixSNPsize.ps",
      "fixsnp.space1", "SNPselection.port2", "fixSNPprop.Label", "fixSNPprop.box", "fixSNPprop.ps",
      "chromotest.Label", "chromotest.port1", "chromotest.LeftLabel", "chromotest.port2", "chromotest.RightLabel",
      "associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.port2", "associationtest.RightLabel",
      "associationtest.Label1", "associationtest.port1.1", "associationtest.Wilcoxon", "associationtest.port1.2", "associationtest.Ftest", "associationtest.port1.3", "associationtest.GEE", "input.cov", "entry.Path3", "box.input3", "input.cov_end",
      "DataVisual.Label", "DataVisual.port1", "DataVisual.LeftLabel", "DataVisual.port2", "DataVisual.RightLabel",
      "Plotselection.space1", "Plotselection1.port", "Plotselection1.lab",
      "Plotselection.space2", "Plotselection2.port", "Plotselection2.lab", "QLOHintensity.Label", "QLOHintensity.box", "QLOHintensity.ps",
      "Plotselection.space3", "Plotselection3.port", "Plotselection3.lab", "Alphascaling.Label", "Alphascaling.box", "Alphascaling.ps")
    example.normal <- c("output.lab", "entry.Path2", "box.output",
      "group.Left", "group.port1", "group.LeftLabel",
      "group.port2", "group.RightLabel", "patienttitle.Label", "patient.port1", "patient1.Label", "patient.port2", "patient2.Label",
      "user.Label", "User.chip",
      "chip.lab", "chip.port1", "chip.LeftLabel", "chip.port2", "chip.MiddleLabel", "chip.port3", "chip.RightLabel", "chip.Label", "chip1.port", "chip1.Label", "chip2.port", "chip2.Label", "chip1_2.port", "chip1_2.Label",
      "chrselect.name", "chr1.port", "chr1.lab", "chr2.port", "chr2.lab", "chr3.port", "chr3.lab", "chr4.port", "chr4.lab", "chr5.port", "chr5.lab", "chr6.port", "chr6.lab", "chr7.port", "chr7.lab", "chr8.port", "chr8.lab", "chr9.port", "chr9.lab", "chr10.port", "chr10.lab", "chr11.port", "chr11.lab", "chr12.port", "chr12.lab",
      "chrselect.space1", "chr13.port", "chr13.lab", "chr14.port", "chr14.lab", "chr15.port", "chr15.lab", "chr16.port", "chr16.lab", "chr17.port", "chr17.lab", "chr18.port", "chr18.lab", "chr19.port", "chr19.lab", "chr20.port", "chr20.lab", "chr21.port", "chr21.lab", "chr22.port", "chr22.lab", "chrall.port", "chrall.lab",
      "stats.Label", "stats.port1", "stats.LeftLabel", "stats.port2", "stats.RightLabel",
      "LOHintensity.Label", "LOHintensity.port1", "LOHintensity.LeftLabel", "LOHintensity.port2", "LOHintensity.RightLabel",
      "SNP.thinning.space", "SNP.thinning.label", "SNP.thinning.port1", "SNPthinning.LeftLabel", "SNPthinningsize.Label", "SNPthinningsize.box", "SNPthinningsize.Label2", "SNP.thinning.port2", "SNPthinning.RightLabel",
      "SNPselection.Label", "SNPselection.port1", "fixSNPsize.Label", "fixSNPsize.box", "fixSNPsize.ps",
      "fixsnp.space1", "SNPselection.port2", "fixSNPprop.Label", "fixSNPprop.box", "fixSNPprop.ps",
      "chromotest.Label", "chromotest.port1", "chromotest.LeftLabel", "chromotest.port2", "chromotest.RightLabel",
      "associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.port2", "associationtest.RightLabel",
      "associationtest.Label1", "associationtest.port1.1", "associationtest.Wilcoxon", "associationtest.port1.2", "associationtest.Ftest", "associationtest.port1.3", "associationtest.GEE",
      "DataVisual.Label", "DataVisual.port1", "DataVisual.LeftLabel", "DataVisual.port2", "DataVisual.RightLabel",
      "Plotselection.space1", "Plotselection1.port", "Plotselection1.lab",
      "Plotselection.space2", "Plotselection2.port", "Plotselection2.lab", "QLOHintensity.Label", "QLOHintensity.box", "QLOHintensity.ps",
      "Plotselection.space3", "Plotselection3.port", "Plotselection3.lab", "Alphascaling.Label", "Alphascaling.box", "Alphascaling.ps")
    plotonly.disabled <- c("group.Left", "group.port1", "group.LeftLabel", "ethnicity.Label", "ethnicity.port1", "ethnicity1.Label", "ethnicity.port2", "ethnicity2.Label",
      "group.port2", "group.RightLabel", "patienttitle.Label", "patient.port1", "patient1.Label", "patient.port2", "patient2.Label",
      "user.Label", "User.chip",
      "chip.lab", "chip.port1", "chip.LeftLabel", "chip.port2", "chip.MiddleLabel", "chip.port3", "chip.RightLabel", "chip.Label", "chip1.port", "chip1.Label", "chip2.port", "chip2.Label", "chip1_2.port", "chip1_2.Label",
      "chrselect.name", "chr1.port", "chr1.lab", "chr2.port", "chr2.lab", "chr3.port", "chr3.lab", "chr4.port", "chr4.lab", "chr5.port", "chr5.lab", "chr6.port", "chr6.lab", "chr7.port", "chr7.lab", "chr8.port", "chr8.lab", "chr9.port", "chr9.lab", "chr10.port", "chr10.lab", "chr11.port", "chr11.lab", "chr12.port", "chr12.lab",
      "chrselect.space1", "chr13.port", "chr13.lab", "chr14.port", "chr14.lab", "chr15.port", "chr15.lab", "chr16.port", "chr16.lab", "chr17.port", "chr17.lab", "chr18.port", "chr18.lab", "chr19.port", "chr19.lab", "chr20.port", "chr20.lab", "chr21.port", "chr21.lab", "chr22.port", "chr22.lab", "chrall.port", "chrall.lab",
      "stats.Label", "stats.port1", "stats.LeftLabel", "stats.port2", "stats.RightLabel",
      "LOHintensity.Label", "LOHintensity.port1", "LOHintensity.LeftLabel", "LOHintensity.port2", "LOHintensity.RightLabel",
      "SNP.thinning.space", "SNP.thinning.label", "SNP.thinning.port1", "SNPthinning.LeftLabel", "SNPthinningsize.Label", "SNPthinningsize.box", "SNPthinningsize.Label2", "SNP.thinning.port2", "SNPthinning.RightLabel",
      "SNPselection.Label", "SNPselection.port1", "fixSNPsize.Label", "fixSNPsize.box", "fixSNPsize.ps",
      "fixsnp.space1", "SNPselection.port2", "fixSNPprop.Label", "fixSNPprop.box", "fixSNPprop.ps",
      "chromotest.Label", "chromotest.port1", "chromotest.LeftLabel", "chromotest.port2", "chromotest.RightLabel",
      "associationtest.Label1", "associationtest.port1.1", "associationtest.Wilcoxon", "associationtest.port1.2", "associationtest.Ftest", "associationtest.port1.3", "associationtest.GEE")
    if (any(toupper(tclvalue(dir.input)) %in% paste("EXAMPLE", 1:5, sep = ""))){
    #if (any(toupper(tclvalue(dir.input)) %in% c("EXAMPLE1", "EXAMPLE2", "EXAMPLE3", "EXAMPLE4"))){
       #examplei <- which(c("EXAMPLE1", "EXAMPLE2", "EXAMPLE3", "EXAMPLE4") == toupper(tclvalue(dir.input)))
       #examplei <- which(paste(rep(c("", "PLOT"), c(4, 4)), "EXAMPLE", rep(1:4, 2), sep = "") == toupper(tclvalue(dir.input)))
       examplei <- which(paste("EXAMPLE", 1:5, sep = "") == toupper(tclvalue(dir.input)))
       tclvalue(dir.output) <- example.dir.output[examplei]
       tclvalue(groupValue) <- example.groupValue[examplei]
       tclvalue(ethnicity.val) <- example.ethnicity.val[examplei]
       tclvalue(chipValue) <- example.chipValue[examplei]
       tclvalue(chip1.val) <- example.chip1.val[examplei]
       tclvalue(chip2.val) <- example.chip2.val[examplei]
       tclvalue(chip1_2.val) <- example.chip1_2.val[examplei]
       for (i in c(1:22, "all")){
         #eval(parse(text = paste("tclvalue(chr", i, ".val) <- example.chr.val[examplei]", sep = "")))
         eval(parse(text = paste("tclvalue(chr", i, ".val) <- example.chr.val[[examplei]][", which(c(1:22, "all") == i), "]", sep = "")))
       }
       tclvalue(statsValue) <- example.statsValue[examplei]
       tclvalue(LOHintensityValue) <- example.LOHintensityValue[examplei]
       tclvalue(SNPthinningvalue) <- example.SNPthinningvalue[examplei]
       tclvalue(SNPthinningsize.val) <- example.SNPthinningsize.val[examplei]
       tclvalue(SNPselectionValue) <- example.SNPselectionValue[examplei]
       tclvalue(fixSNPsize.val) <- example.fixSNPsize.val[examplei]
       tclvalue(fixSNPprop.val) <- example.fixSNPprop.val[examplei]
       tclvalue(chromotestValue) <- example.chromotestValue[examplei]
       tclvalue(associationcheck) <- example.associationcheck[examplei]
       tclvalue(associationValue) <- example.associationValue[examplei]
       tclvalue(ftestValue) <- example.ftestValue[examplei]
       tclvalue(geeValue) <- example.geeValue[examplei]
       tclvalue(dir.cov) <- example.dir.cov[examplei]
       tclvalue(DataVisualValue) <- example.DataVisualValue[examplei]
       tclvalue(Plotselection1.val) <- example.Plotselection1.val[examplei]
       tclvalue(Plotselection2.val) <- example.Plotselection2.val[examplei]
       tclvalue(Plotselection3.val) <- example.Plotselection3.val[examplei]
       tclvalue(QLOHintensity.val) <- example.QLOHintensity.val[examplei]
       tclvalue(Alphascaling.val) <- example.Alphascaling.val[examplei]
       Disabled(example.disabled)
    }else if(any(c(paste("chr", 1:22, ".part1.RData", sep = ""),
                   paste("chr", 1:22, "_Chip1.part1.RData", sep = ""),
                   paste("chr", 1:22, "_Chip2.part1.RData", sep = ""),
                   paste("chr", 1:22, "_Chip12.part1.RData", sep = "")) %in% dir(tclvalue(dir.input), pattern = ".RData"))){
       plot.only <<- TRUE # stand for indicating using .RData draw plots
       ### drop needed .RData
       plot.rdfile <<- dir(tclvalue(dir.input), pattern = ".RData")
       plot.rdfile <<- plot.rdfile[substr(plot.rdfile, 1, 3) %in% "chr"]
       plot.rdfile <<- plot.rdfile[unlist(lapply(1:length(plot.rdfile), function(k) substr(strsplit(plot.rdfile[k], ".", fixed = TRUE)[[1]][2], 1, 4))) %in% "part"]
  
       tclvalue(dir.output) <- ""
       tclvalue(groupValue) <- "2"
       tclvalue(ethnicity.val) <- "1"
       ### find out whether there are more than one chip for a person(Aff100k/Affy500k or Affy6.0/Customized)
       plot.chip <<- unique(unlist(lapply(1:length(plot.rdfile), function(k) gsub("[.]$", "", substr(strsplit(plot.rdfile[k], "_Chip", fixed = TRUE)[[1]][2], 1, 2)))))
       tclvalue(chipValue) <- if(all(is.na(plot.chip))){"4"}else{"1"}
       tclvalue(chip1.val) <- if("1" %in% plot.chip){"1"}else{"0"}
       tclvalue(chip2.val) <- if("2" %in% plot.chip){"1"}else{"0"}
       tclvalue(chip1_2.val) <- if("12" %in% plot.chip){"1"}else{"0"}
       ### only turn on chips which have .RData
       plotonly.disabled <- plotonly.disabled[!(plotonly.disabled %in% c("chip1.port", "chip1.Label", "chip2.port", "chip2.Label", "chip1_2.port", "chip1_2.Label")[rep(c("1", "2", "12") %in% plot.chip, 2)])]
       if (tclvalue(chipValue) == "1"){plotonly.disabled <- plotonly.disabled[!(plotonly.disabled %in% "chip.Label")]}else{}
  
       ### find which chromosome could draw biplot
       plot.chr <<- unique(unlist(lapply(1:length(plot.rdfile), function(k) strsplit(gsub("chr", "", gsub("_Chip[1-2]{1,2}", "", plot.rdfile[k])), ".", fixed = TRUE)[[1]][1])))
       for (i in 1:22){
         eval(parse(text = paste("tclvalue(chr", i, ".val) <- ", if(i %in% plot.chr){1}else{0}, sep = "")))
       }
       eval(parse(text = paste("tclvalue(chrall.val) <- ", if(all(1:22 %in% plot.chr)){1}else{0}, sep = "")))
       ### only turn on chromosome selected buttons which have .RData
       plotonly.disabled <- plotonly.disabled[!(plotonly.disabled %in% "chrselect.name")]
       plotonly.disabled <- plotonly.disabled[!(plotonly.disabled %in% c(paste("chr", plot.chr, ".lab", sep = ""), paste("chr", plot.chr, ".port", sep = "")))]
       if (all(1:22 %in% plot.chr)){plotonly.disabled <- plotonly.disabled[!(plotonly.disabled %in% c("chrall.lab", "chrall.port"))]}else{}
  
       tclvalue(statsValue) <- "1"
       tclvalue(LOHintensityValue) <- "2" # turn off this because there's .RData
       tclvalue(SNPthinningvalue) <- "2"
       tclvalue(SNPthinningsize.val) <- "2"
       tclvalue(SNPselectionValue) <- "1"
       tclvalue(fixSNPsize.val) <- "50"
       tclvalue(fixSNPprop.val) <- "0.05"
       tclvalue(chromotestValue) <- "2"
       tclvalue(associationcheck) <- "2"
       tclvalue(associationValue) <- "0"
       tclvalue(ftestValue) <- "0"
       tclvalue(geeValue) <- "0"
       tclvalue(dir.cov) <- ""
       tclvalue(DataVisualValue) <- "1"
       tclvalue(Plotselection1.val) <- "0"
       tclvalue(Plotselection2.val) <- "0"
       tclvalue(Plotselection3.val) <- "1"
       tclvalue(QLOHintensity.val) <- "0.9"
       tclvalue(Alphascaling.val) <- "0"
       Disabled(plotonly.disabled)
    }else{
       tclvalue(dir.output) <- ""
       tclvalue(groupValue) <- "2"
       tclvalue(ethnicity.val) <- "1"
       tclvalue(chipValue) <- "1"
       tclvalue(chip1.val) <- "1"
       tclvalue(chip2.val) <- "1"
       tclvalue(chip1_2.val) <- "1"
       tclvalue(chr1.val) <- "1"
       for (i in c(2:22, "all")){
         eval(parse(text = paste("tclvalue(chr", i, ".val) <- 0", sep = "")))
       }
       tclvalue(statsValue) <- "1"
       tclvalue(LOHintensityValue) <- "1"
       tclvalue(SNPthinningvalue) <- "2"
       tclvalue(SNPthinningsize.val) <- "2"
       tclvalue(SNPselectionValue) <- "1"
       tclvalue(fixSNPsize.val) <- "50"
       tclvalue(fixSNPprop.val) <- "0.05"
       tclvalue(chromotestValue) <- "1"
       tclvalue(associationcheck) <- "1"
       tclvalue(associationValue) <- "1"
       tclvalue(ftestValue) <- "0"
       tclvalue(geeValue) <- "0"
       tclvalue(dir.cov) <- ""
       tclvalue(DataVisualValue) <- "1"
       tclvalue(Plotselection1.val) <- "1"
       tclvalue(Plotselection2.val) <- "1"
       tclvalue(Plotselection3.val) <- "1"
       tclvalue(QLOHintensity.val) <- "0.9"
       tclvalue(Alphascaling.val) <- "0"
       Normal(example.normal)
    }
  })
  input.lab<-tklabel(frame.path,background=color,font=fonttitle,text=" Input directory:        ")
  tkgrid(input.lab,entry.Path,box.input)
  
  dir.output<- tclVar("")
  box.output <- tkbutton(frame.path,text="...", command=function() tclvalue(dir.output) <- tkchooseDirectory())
  entry.Path2<-tkentry(frame.path,width="94",textvariable=dir.output)
  output.lab<-tklabel(frame.path,background=color,font=fonttitle,text=" Output directory:       ")
  tkgrid(output.lab,entry.Path2,box.output)
  
  ##### group selection title #####
  
  font.group<- tkfont.create(family="times",size=12,weight="bold") # text type and size of title
  group.Left <- tklabel(frame.group.1,background=color,text="The number of study groups  ")
  tkgrid(tklabel(frame.group.1,background=color,text="2. Study group:",font=font.group),sticky="w") 
  
  ##### group selection box #####
  ### turn off the buttons of group 2 if 1 is chosen, and vise versa 20111223
  group.port1 <- tkradiobutton(frame.group.1.1,background=color,command=function(){
    tclvalue(chromotestValue)="2";
    tclvalue(associationValue)="2";  # Allow Regression test and GEE in one group
    tclvalue(associationcheck)="2";
    Normal(c("ethnicity.Label", "ethnicity1.Label", "ethnicity2.Label", "ethnicity.port1", "ethnicity.port2"))
    Normal(c("associationtest.port2"))
    Disabled(c("chromotest.Label", "chromotest.port1", "chromotest.LeftLabel", "chromotest.port2", "chromotest.RightLabel"))
    #Disabled(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
    #"associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
    #"associationtest.port2", "input.cov", "input.cov", "input.cov_end", "box.input3", "entry.Path3", "associationtest.port1.3", "associationtest.GEE"))
    Disabled(c("associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
    "input.cov", "input.cov_end", "box.input3", "entry.Path3", "associationtest.port1.3", "associationtest.GEE"))
    Disabled(c("patienttitle.Label", "patient1.Label", "patient2.Label", "patient.port1", "patient.port2"))
  })
  
  group.port2 <- tkradiobutton(frame.group.1.1,background=color,command=function(){
    tclvalue(associationValue)="1";
    tclvalue(chromotestValue)="1";
    tclvalue(associationcheck)="1";
    Normal(c("chromotest.Label", "chromotest.port1", "chromotest.LeftLabel", "chromotest.port2", "chromotest.RightLabel"))
    Normal(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
    "associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
    "associationtest.port2", "associationtest.port1.3", "associationtest.GEE"))
    Normal(c("patienttitle.Label", "patient1.Label", "patient2.Label", "patient.port1", "patient.port2"))
    Disabled(c("ethnicity.Label", "ethnicity1.Label", "ethnicity2.Label", "ethnicity.port1", "ethnicity.port2"))
  })
  
  patient.port1 <- tkradiobutton(frame.group.1.1,background=color,command=function(){})
  patient.port2 <- tkradiobutton(frame.group.1.1,background=color,command=function(){})
  patient1.val <- tclVar("1")
  tkconfigure(patient.port1,variable=patient1.val,value="1")
  tkconfigure(patient.port2,variable=patient1.val,value="2")
  
  patienttitle.Label <- tklabel(frame.group.1.1,background=color,font=fonttitle,text="   ( Patient group: ")
  patient1.Label <- tklabel(frame.group.1.1,background=color,font=fonttitle,text="1st group")
  patient2.Label <- tklabel(frame.group.1.1,background=color,font=fonttitle,text="2nd group )")
  
  groupValue <- tclVar("2")
  tkconfigure(group.port1,variable=groupValue,value="1")
  tkconfigure(group.port2,variable=groupValue,value="2")
  group.LeftLabel <- tklabel(frame.group.1.1,text="One group",font=fonttitle,background=color)
  group.RightLabel <- tklabel(frame.group.1.1,text=" Two groups",font=fonttitle,background=color)
  group.Left <- tklabel(frame.group.1.1,text="     The number of study groups:",font=fonttitle,background=color)
  
  #tkgrid(frame.group.1.1,frame.group.1.2)
  
  ##### ethnities selection #####
  	 ethnicity.port1 <-tkradiobutton(frame.group.1.1,background=color,command=function(){
  	   tclvalue(associationValue)="2"; 
  	   tclvalue(associationcheck)="2";
  	   Normal(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel", "associationtest.port2"))
  	 }, state = "disabled")
  	 ethnicity.port2 <-tkradiobutton(frame.group.1.1,background=color,command=function(){
  	   tclvalue(associationValue)="2"; 
  	   tclvalue(associationcheck)="2";
  	   Disabled(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
  	   "associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
  	   "associationtest.port2", "input.cov", "input.cov", "input.cov_end", "box.input3", "entry.Path3", "associationtest.port1.3", "associationtest.GEE"))
           }, state = "disabled")
  	 ethnicity.val <- tclVar("1")
  	 tkconfigure(ethnicity.port1,variable=ethnicity.val,value="1")
  	 tkconfigure(ethnicity.port2,variable=ethnicity.val,value="2")
  	 
  
  	 ethnicity.Label <- tklabel(frame.group.1.1,background=color,font=fonttitle,text=" ( Populations: ", state = "disabled")
  	 ethnicity1.Label <- tklabel(frame.group.1.1,background=color,font=fonttitle,text="Single     ", state = "disabled")
  	 ethnicity2.Label <- tklabel(frame.group.1.1,background=color,font=fonttitle,text="Mutiple )     ", state = "disabled")
  
  tkgrid(group.Left,group.port1,group.LeftLabel,ethnicity.Label,ethnicity.port1,ethnicity1.Label,ethnicity.port2,ethnicity2.Label)
  tkgrid(tklabel(frame.group.1.1,background=color,text="           "),group.port2,group.RightLabel,patienttitle.Label,patient.port1,patient1.Label,patient.port2,patient2.Label)
  
  
  
  
  ##### chip selection title #####
  
  font.chip<- tkfont.create(family="times",size=12,weight="bold") # text type and size od title
  tkgrid(tklabel(frame.chip,background=color,text="3. Data format:",font=font.chip),sticky="w") # title
  
  #####  Gene chip selection box #####
  
  User.chip<- tkradiobutton(frame.chip.1,background=color,command=function(){tclvalue(chr1.val)=0;tclvalue(chr2.val)=0;tclvalue(chr3.val)=0;
  		tclvalue(chr4.val)=0;tclvalue(chr5.val)=0;tclvalue(chr6.val)=0;
  		tclvalue(chr7.val)=0;tclvalue(chr8.val)=0;tclvalue(chr9.val)=0;
  		tclvalue(chr10.val)=0;tclvalue(chr11.val)=0;tclvalue(chr12.val)=0;
  		tclvalue(chr13.val)=0;tclvalue(chr14.val)=0;tclvalue(chr15.val)=0;
  		tclvalue(chr16.val)=0;tclvalue(chr17.val)=0;tclvalue(chr18.val)=0;
  		tclvalue(chr19.val)=0;tclvalue(chr20.val)=0;tclvalue(chr21.val)=0;
  		tclvalue(chr22.val)=0;tclvalue(chip1.val)=0;
  		#tclvalue(chr23.val)=0; # chr23 are not analyzed 20130109
  		tclvalue(chrall.val)=0;tclvalue(chip2.val)=0;tclvalue(chip1_2.val)=0;
      Disabled(c("chip.Label", "chip1.port", "chip2.port", "chip1_2.port", "chip1.Label", "chip2.Label", "chip1_2.Label"))
      #Disabled(paste("chr", c(1:23, "all"), ".lab", sep=""))
      Disabled(paste("chr", c(1:22, "all"), ".lab", sep=""))
      #Disabled(paste("chr", c(1:23, "all"), ".port", sep=""))
      Disabled(paste("chr", c(1:22, "all"), ".port", sep=""))
  })
  chip.port1 <- tkradiobutton(frame.chip.1,background=color,command=function(){
      tclvalue(chip1.val) <- "1";tclvalue(chip2.val) <- "1";tclvalue(chip1_2.val) <- "1";
      Normal(c("chip.Label", "chip1.port", "chip2.port", "chip1_2.port", "chip1.Label", "chip2.Label", "chip1_2.Label"))
      #Normal(paste("chr", c(1:23, "all"), ".lab", sep="")) # chr23 are not analyzed 20130109
      Normal(paste("chr", c(1:22, "all"), ".lab", sep=""))
      #Normal(paste("chr", c(1:23, "all"), ".port", sep=""))
      Normal(paste("chr", c(1:22, "all"), ".port", sep=""))
  })
  chip.port2 <- tkradiobutton(frame.chip.1,background=color,command=function(){
      tclvalue(chip1.val) <- "1";tclvalue(chip2.val) <- "1";tclvalue(chip1_2.val) <- "1";
      Normal(c("chip.Label", "chip1.port", "chip2.port", "chip1_2.port", "chip1.Label", "chip2.Label", "chip1_2.Label"))
      #Normal(paste("chr", c(1:23, "all"), ".lab", sep="")) # chr23 are not analyzed 20130109
      Normal(paste("chr", c(1:22, "all"), ".lab", sep=""))
      #Normal(paste("chr", c(1:23, "all"), ".port", sep=""))
      Normal(paste("chr", c(1:22, "all"), ".port", sep=""))
  })
  chip.port3 <- tkradiobutton(frame.chip.1,background=color,command=function(){
      tclvalue(chip1.val)=0;tclvalue(chip2.val)=0;tclvalue(chip1_2.val)=0;
      Disabled(c("chip.Label", "chip1.port", "chip2.port", "chip1_2.port", "chip1.Label", "chip2.Label", "chip1_2.Label"))
      #Normal(paste("chr", c(1:23, "all"), ".lab", sep="")) # chr23 are not analyzed 20130109
      Normal(paste("chr", c(1:22, "all"), ".lab", sep=""))
      #Normal(paste("chr", c(1:23, "all"), ".port", sep=""))
      Normal(paste("chr", c(1:22, "all"), ".port", sep=""))
  })
  
  chipValue <- tclVar("1")
  tkconfigure(chip.port1,variable=chipValue,value="1")
  tkconfigure(chip.port2,variable=chipValue,value="2")
  tkconfigure(chip.port3,variable=chipValue,value="3")
  tkconfigure(User.chip,variable=chipValue,value="4")
  user.Label <- tklabel(frame.chip.1,background=color,font=fonttitle,text="  Customized SNP panel:")
  chip.LeftLabel <- tklabel(frame.chip.1,background=color,font=fonttitle,text="Affy 100K")
  chip.MiddleLabel <- tklabel(frame.chip.1,background=color,font=fonttitle,text="Affy 500K")
  chip.RightLabel <- tklabel(frame.chip.1,background=color,font=fonttitle,text="Affy 6.0")
  chip.lab<-tklabel(frame.chip.1,background=color,font=fonttitle,text="    Genome-wide gene chip:")
  
  tkgrid(user.Label,User.chip)
  chip1.lab <- tklabel(frame.chip.1,background=color,font=fonttitle,text="1")
  chip2.lab <- tklabel(frame.chip.1,background=color,font=fonttitle,text="2")
  chrp1_2.lab <- tklabel(frame.chip.1,background=color,font=fonttitle,text="3")
  
  chip.Label <- tklabel(frame.chip.1,background=color,font=fonttitle,text="(")
  
  chip1.port <- tkcheckbutton(frame.chip.1,background=color)
  chip2.port <- tkcheckbutton(frame.chip.1,background=color)
  chip1_2.port <- tkcheckbutton(frame.chip.1,background=color)
  
  chip1.val <- tclVar("1")
  chip2.val <- tclVar("1")
  chip1_2.val <- tclVar("1")
  
  tkconfigure(chip1.port,variable=chip1.val)
  tkconfigure(chip2.port,variable=chip2.val)
  tkconfigure(chip1_2.port,variable=chip1_2.val)
  
  ### example ###
  
  chip1.Label <- tklabel(frame.chip.1,background=color,font=fonttitle,text="Chip1")
  chip2.Label <- tklabel(frame.chip.1,background=color,font=fonttitle,text="Chip2")
  chip1_2.Label<-tklabel(frame.chip.1,background=color,font=fonttitle,text="Merge 2 chips )")
  
  tkgrid(chip.lab,chip.port1,chip.LeftLabel,chip.port2,chip.MiddleLabel,chip.port3,chip.RightLabel,chip.Label,chip1.port,chip1.Label,chip2.port,chip2.Label,chip1_2.port,chip1_2.Label)
  #tkgrid(chip.Label,chip1.port,chip1.Label,chip2.port,chip2.Label,chip1_2.port,chip1_2.Label)
  
  tkgrid.configure(chip1.Label,sticky="w")
  tkgrid.configure(chip2.Label,sticky="w")
  tkgrid.configure(chip1_2.Label,sticky="w")
  
  ##### chromosome and rank test selection title #####
  
  ##### chromosome and rank test selection box #####
  
  #for(i in 1:23){ # chr23 are not analyzed 20130109
  for(i in 1:22){
    assign(paste("chr", i, ".lab", sep=""), tklabel(frame.para, background = color, font = fonttitle, text = paste(i)))
  }
  chrall.lab <- tklabel(frame.para,background=color,font=fonttitle, text="All")
  
  chr1.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr1.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr2.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr2.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr3.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr3.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr4.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr4.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr5.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr5.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr6.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr6.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr7.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr7.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr8.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr8.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr9.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr9.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr10.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr10.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr11.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr11.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr12.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr12.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr13.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr13.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr14.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr14.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr15.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr15.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr16.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr16.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr17.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr17.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr18.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr18.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr19.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr19.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr20.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr20.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr21.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr21.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  chr22.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr22.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  #chr23.port <- tkcheckbutton(frame.para,background=color,command=function(){if(tclvalue(chr23.val)=="0"){tclvalue(chrall.val)=0}else{tclvalue(LOHintensityValue)=1}})
  # chr23 are not analyzed 20130109
  chrall.port <- tkcheckbutton(frame.para,background=color,command=function(){
    if(tclvalue(chrall.val)=="1"){
      tclvalue(chr1.val)=1;tclvalue(chr2.val)=1;tclvalue(chr3.val)=1;
      tclvalue(chr4.val)=1;tclvalue(chr5.val)=1;tclvalue(chr6.val)=1;
      tclvalue(chr7.val)=1;tclvalue(chr8.val)=1;tclvalue(chr9.val)=1;
      tclvalue(chr10.val)=1;tclvalue(chr11.val)=1;tclvalue(chr12.val)=1;
      tclvalue(chr13.val)=1;tclvalue(chr14.val)=1;tclvalue(chr15.val)=1;
      tclvalue(chr16.val)=1;tclvalue(chr17.val)=1;tclvalue(chr18.val)=1;
      tclvalue(chr19.val)=1;tclvalue(chr20.val)=1;tclvalue(chr21.val)=1;
      tclvalue(chr22.val)=1;tclvalue(LOHintensityValue)=1;
      #tclvalue(chr23.val)=1; # chr23 are not analyzed 20130109
    }else{
      tclvalue(chr1.val)=0;tclvalue(chr2.val)=0;tclvalue(chr3.val)=0;
      tclvalue(chr4.val)=0;tclvalue(chr5.val)=0;tclvalue(chr6.val)=0;
      tclvalue(chr7.val)=0;tclvalue(chr8.val)=0;tclvalue(chr9.val)=0;
      tclvalue(chr10.val)=0;tclvalue(chr11.val)=0;tclvalue(chr12.val)=0;
      tclvalue(chr13.val)=0;tclvalue(chr14.val)=0;tclvalue(chr15.val)=0;
      tclvalue(chr16.val)=0;tclvalue(chr17.val)=0;tclvalue(chr18.val)=0;
      tclvalue(chr19.val)=0;tclvalue(chr20.val)=0;tclvalue(chr21.val)=0;
      tclvalue(chr22.val)=0;
      #tclvalue(chr23.val)=0; # chr23 are not analyzed 20130109
    }
  })
  
  chr1.val <- tclVar("1")
  #for(i in 2:23){
  for(i in 2:22){
    assign(paste("chr", i, ".val", sep=""), tclVar("0"))
  }
  chrall.val <- tclVar("0")
  
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
  #tkconfigure(chr23.port,variable=chr23.val) # chr23 are not analyzed 20130109
  tkconfigure(chrall.port,variable=chrall.val)
  
  chrselect.name <- tklabel(frame.para,background=color,font=fonttitle,text="    Chromosome: ")
  chrselect.space1 <- tklabel(frame.para,background=color,text="")
  #chrselect.space2<-tklabel(frame.para,background=color,text="")
  
  tkgrid(chrselect.name,chr1.port,chr1.lab,chr2.port,chr2.lab,chr3.port,chr3.lab,chr4.port,chr4.lab,chr5.port,chr5.lab,chr6.port,chr6.lab,chr7.port,chr7.lab,chr8.port,chr8.lab,chr9.port,chr9.lab,chr10.port,chr10.lab,chr11.port,chr11.lab,chr12.port,chr12.lab)
  #tkgrid(chrselect.space1,chr13.port,chr13.lab,chr14.port,chr14.lab,chr15.port,chr15.lab,chr16.port,chr16.lab,chr17.port,chr17.lab,chr18.port,chr18.lab,chr19.port,chr19.lab,chr20.port,chr20.lab,chr21.port,chr21.lab,chr22.port,chr22.lab,chr23.port,chr23.lab,chrall.port,chrall.lab) # chr23 are not analyzed 20130109
  tkgrid(chrselect.space1,chr13.port,chr13.lab,chr14.port,chr14.lab,chr15.port,chr15.lab,chr16.port,chr16.lab,chr17.port,chr17.lab,chr18.port,chr18.lab,chr19.port,chr19.lab,chr20.port,chr20.lab,chr21.port,chr21.lab,chr22.port,chr22.lab,chrall.port,chrall.lab)
  #tkgrid(chrselect.space2,chr21.port,chr21.lab,chr22.port,chr22.lab,chr23.port,chr23.lab,chrall.port,chrall.lab)
  ##### Statistical analysis title #####
  
  font.para2<- tkfont.create(family="times",size=12,weight="bold") # text type and size od title
  tkgrid(tklabel(frame.para2,background=color,text="4. Statistical analysis:",font=font.para2),sticky="w") # title
  
  ##### Statistical analysis #####
  # Descriptive statistics #
  stats.port1 <- tkradiobutton(frame.para2,background=color)
  stats.port2 <- tkradiobutton(frame.para2,background=color)
  statsValue <- tclVar("1")
  tkconfigure(stats.port1,variable=statsValue,value="1")
  tkconfigure(stats.port2,variable=statsValue,value="2")
  stats.LeftLabel <- tklabel(frame.para2,background=color,font=fonttitle,text="Yes")
  stats.RightLabel <- tklabel(frame.para2,background=color,font=fonttitle,text="No")
  stats.Label <- tklabel(frame.para2,background=color,font=fonttitle,text="Data description:           ")
  tkgrid(stats.Label,stats.port1,stats.LeftLabel,stats.port2,stats.RightLabel)
  
  # Estimate of LOH intensity #
  LOHintensity.port1 <- tkradiobutton(frame.para2,background=color,command=function(){tclvalue(DataVisualValue)=1;
      tclvalue(associationcheck)="1";
      tclvalue(chromotestValue)="1";
      tclvalue(Plotselection1.val)=1;tclvalue(Plotselection2.val)=1;tclvalue(Plotselection3.val)=1;
      Normal(c("SNPselection.Label", "SNPselection.port1", "SNPselection.port2"))
      Normal(c("fixSNPsize.Label", "fixSNPsize.box", "fixSNPsize.ps"))
      Normal(c("fixSNPprop.Label", "fixSNPprop.box", "fixSNPprop.ps"))
      Normal(c("chromotest.port1", "chromotest.port2", "chromotest.LeftLabel", "chromotest.RightLabel", "chromotest.Label"))
      Normal(c("DataVisual.port1", "DataVisual.port2", "DataVisual.LeftLabel", "DataVisual.Label", "DataVisual.RightLabel"))
      Normal(c("Plotselection1.lab", "Plotselection2.lab", "Plotselection3.lab", "Plotselection1.port", "Plotselection2.port", "Plotselection3.port"))
      Normal(c("QLOHintensity.Label", "QLOHintensity.box", "QLOHintensity.ps"))
      Normal(c("Alphascaling.Label", "Alphascaling.box", "Alphascaling.ps"))
      Normal(c("SNP.thinning.space", "SNP.thinning.label", "SNP.thinning.port1", "SNPthinning.LeftLabel", "SNPthinningsize.Label",
        "SNPthinningsize.box", "SNPthinningsize.Label2", "SNP.thinning.port2", "SNPthinning.RightLabel"))
      if (tclvalue(groupValue) == "2" & tclvalue(ethnicity.val) == "1"){
        tclvalue(associationValue)="1";
        Normal(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
        "associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
        "associationtest.port2", "associationtest.port1.3", "associationtest.GEE"))
      }else{
        tclvalue(associationValue)="2";
        Normal(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
        "associationtest.Ftest", "associationtest.port1.2",
        "associationtest.port2", "associationtest.port1.3", "associationtest.GEE"))
      }
      #Normal(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
      #  "associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
      #  "associationtest.port2", "associationtest.port1.3", "associationtest.GEE"))
  })
  LOHintensity.port2 <- tkradiobutton(frame.para2,background=color,command=function()
  		{tclvalue(chromotestValue)=2;tclvalue(DataVisualValue)=2;
  		tclvalue(Plotselection1.val)=0;tclvalue(Plotselection2.val )=0;
  		tclvalue(Plotselection3.val)=0;tclvalue(chrall.val)=0;
  		tclvalue(chr1.val)=0;tclvalue(chr2.val)=0;tclvalue(chr3.val)=0;
  		tclvalue(chr4.val)=0;tclvalue(chr5.val)=0;tclvalue(chr6.val)=0;
  		tclvalue(chr7.val)=0;tclvalue(chr8.val)=0;tclvalue(chr9.val)=0;
  		tclvalue(chr10.val)=0;tclvalue(chr11.val)=0;tclvalue(chr12.val)=0;
  		tclvalue(chr13.val)=0;tclvalue(chr14.val)=0;tclvalue(chr15.val)=0;
  		tclvalue(chr16.val)=0;tclvalue(chr17.val)=0;tclvalue(chr18.val)=0;
  		tclvalue(chr19.val)=0;tclvalue(chr20.val)=0;tclvalue(chr21.val)=0;
  		tclvalue(chr22.val)=0#;tclvalue(chr23.val)=0;
  		tclvalue(associationValue)="2";
  		tclvalue(associationcheck)="2";
  		tclvalue(Plotselection1.val)=0;tclvalue(Plotselection2.val)=0;tclvalue(Plotselection3.val)=0;
      Disabled(c("SNPselection.Label", "SNPselection.port1", "SNPselection.port2"))
      Disabled(c("fixSNPsize.Label", "fixSNPsize.box", "fixSNPsize.ps"))
      Disabled(c("fixSNPprop.Label", "fixSNPprop.box", "fixSNPprop.ps"))
      Disabled(c("chromotest.port1", "chromotest.port2", "chromotest.LeftLabel", "chromotest.RightLabel", "chromotest.Label"))
      Disabled(c("DataVisual.port1", "DataVisual.port2", "DataVisual.LeftLabel", "DataVisual.Label", "DataVisual.RightLabel"))
      Disabled(c("Plotselection1.lab", "Plotselection2.lab", "Plotselection3.lab", "Plotselection1.port", "Plotselection2.port", "Plotselection3.port"))
      Disabled(c("QLOHintensity.Label", "QLOHintensity.box", "QLOHintensity.ps"))
      Disabled(c("Alphascaling.Label", "Alphascaling.box", "Alphascaling.ps"))
      Disabled(c("SNP.thinning.space", "SNP.thinning.label", "SNP.thinning.port1", "SNPthinning.LeftLabel", "SNPthinningsize.Label",
        "SNPthinningsize.box", "SNPthinningsize.Label2", "SNP.thinning.port2", "SNPthinning.RightLabel"))
      Disabled(c("associationtest.Label", "associationtest.port1", "associationtest.LeftLabel", "associationtest.RightLabel",
        "associationtest.Wilcoxon", "associationtest.Ftest", "associationtest.port1.1", "associationtest.port1.2",
        "associationtest.port2", "input.cov", "input.cov", "input.cov_end", "box.input3", "entry.Path3", "associationtest.port1.3",
        "associationtest.GEE"))
  })
  
  LOHintensityValue <- tclVar("1")
  LOHintensityValue.unclick <- tclVar("2")
  tkconfigure(LOHintensity.port1,variable=LOHintensityValue,value="1")
  tkconfigure(LOHintensity.port2,variable=LOHintensityValue,value="2")
  LOHintensity.LeftLabel <- tklabel(frame.para2,background=color,font=fonttitle,text="Yes")
  LOHintensity.RightLabel <- tklabel(frame.para2,background=color,font=fonttitle,text="No")
  LOHintensity.Label <- tklabel(frame.para2,background=color,font=fonttitle,text="     Estimate of LOH intensity:")
  tkgrid(LOHintensity.Label,LOHintensity.port1,LOHintensity.LeftLabel,LOHintensity.port2,LOHintensity.RightLabel)
  
  # SNP thinning #
  SNP.thinning.space<-tklabel(frame.para2.1,background=color,text="       ")
  SNP.thinning.label<-tklabel(frame.para2.1,background=color,font=fonttitle,text="    - SNP thinning:")
  SNP.thinning.port1<-tkradiobutton(frame.para2.1,background=color)
  SNP.thinning.port2<-tkradiobutton(frame.para2.1,background=color)
  
  SNPthinningvalue<- tclVar("2")
  tkconfigure(SNP.thinning.port1,variable=SNPthinningvalue,value="1")
  tkconfigure(SNP.thinning.port2,variable=SNPthinningvalue,value="2")
  
  SNPthinningsize.val <- tclVar("2")
  
  SNPthinning.LeftLabel <- tklabel(frame.para2.1,background=color,font=fonttitle,text="Yes")
  SNPthinning.RightLabel <- tklabel(frame.para2.1,background=color,font=fonttitle,text="No")
  SNPthinningsize.Label <- tklabel(frame.para2.1,background=color,font=fonttitle,text="(Select one every k SNPs: ")
  SNPthinningsize.box <- tkentry(frame.para2.1,width=7,textvariable=SNPthinningsize.val)
  SNPthinningsize.Label2 <- tklabel(frame.para2.1,background=color,font=fonttitle,text=", between 2 and 10 )")
  
  tkgrid(SNP.thinning.space,SNP.thinning.label,SNP.thinning.port1,SNPthinning.LeftLabel,SNPthinningsize.Label,SNPthinningsize.box,SNPthinningsize.Label2,SNP.thinning.port2,SNPthinning.RightLabel)
  
  # SNP selection in the window #
  
  SNPselection.Label <- tklabel(frame.para3,background=color,font=fonttitle,text="             - Window size: ")
  
  SNPselection.port1 <- tkradiobutton(frame.para3,background=color)
  SNPselection.port2 <- tkradiobutton(frame.para3,background=color)
  #SNPselection.port3 <- tkradiobutton(frame.para3,background=color)
  SNPselectionValue <- tclVar("1")
  tkconfigure(SNPselection.port1,variable=SNPselectionValue,value="1")
  tkconfigure(SNPselection.port2,variable=SNPselectionValue,value="2")
  
  fixSNPsize.val <- tclVar("50")
  fixSNPprop.val <- tclVar("0.05")
  
  fixSNPsize.Label <- tklabel(frame.para3,background=color,font=fonttitle,text="Fixed SNP number:     ")
  fixSNPsize.box <- tkentry(frame.para3,width=7,textvariable=fixSNPsize.val)
  fixSNPsize.ps <- tklabel(frame.para3,background=color,font=fonttitle,text=" , between 1 and the total number of SNPs                                ")
  tkgrid(SNPselection.Label,SNPselection.port1,fixSNPsize.Label,fixSNPsize.box,fixSNPsize.ps)
  
  fixSNPprop.Label <- tklabel(frame.para3,background=color,font=fonttitle,text="Fixed SNP proportion:")
  fixSNPprop.box <- tkentry(frame.para3,width=7,textvariable=fixSNPprop.val)
  fixSNPprop.ps <- tklabel(frame.para3,background=color,font=fonttitle,text=", between 0 and 1                                                                          ")
  
  fixsnp.space1<-tklabel(frame.para3,background=color,text="")
  tkgrid(fixsnp.space1,SNPselection.port2,fixSNPprop.Label,fixSNPprop.box,fixSNPprop.ps)      
  
  # Chromosome ranking test #
  
  chromotest.port1 <- tkradiobutton(frame.para4.1,background=color,command=function(){tclvalue(LOHintensityValue)=1;})
  chromotest.port2 <- tkradiobutton(frame.para4.1,background=color)
  chromotestValue <- tclVar("1")
  tkconfigure(chromotest.port1,variable=chromotestValue,value="1")
  tkconfigure(chromotest.port2,variable=chromotestValue,value="2")
  chromotest.LeftLabel <- tklabel(frame.para4.1,background=color,font=fonttitle,text="Yes")
  chromotest.RightLabel <- tklabel(frame.para4.1,background=color,font=fonttitle,text="No")
  chromotest.Label <- tklabel(frame.para4.1,background=color,font=fonttitle,text="     Chromosome ranking test:")
  tkgrid(chromotest.Label,chromotest.port1,chromotest.LeftLabel,chromotest.port2,chromotest.RightLabel)
  
  # Association test #
  associationcheck <- tclVar("1") ### check "Yes" or "No" for association test 20120106 20120130
  associationValue <- tclVar("1") ### check "Wilcoxon", "1": Wilcoxon
  ftestValue <- tclVar("0")
  geeValue <- tclVar("0")
  dir.cov <- tclVar("")
  #dir.cov1 <- tclVar("")
  associationtest.Label<- tklabel(frame.para4.1.1,background=color,font=fonttitle,text="     Association test:                   ")
  associationtest.LeftLabel <- tklabel(frame.para4.1.1,background=color,font=fonttitle,text="Yes")
  associationtest.RightLabel <- tklabel(frame.para4.1.1,background=color,font=fonttitle,text="No")
  associationtest.Label1 <- tklabel(frame.para4.1.2,background=color,font=fonttitle,text="                                                           -  ") 
  associationtest.Wilcoxon <- tklabel(frame.para4.1.2, background=color, font=fonttitle, text="Wilcoxon")
  associationtest.Ftest <- tklabel(frame.para4.1.2, background=color, font=fonttitle, text="Regression")
  associationtest.GEE <- tklabel(frame.para4.1.2, background=color, font=fonttitle, text="GEE")
  #associationtest.Label2 <- tklabel(frame.para4.1.3,background=color,font=fonttitle,text="                                                           -                             ") 
  #associationtest.GEE <- tklabel(frame.para4.1.3, background=color, font=fonttitle, text="GEE")
  associationtest.port1.1 <- tkcheckbutton(frame.para4.1.2, background=color)
  associationtest.port1.2 <- tkcheckbutton(frame.para4.1.2, background=color, command = function(){
    if(tclvalue(ftestValue) == "0" & tclvalue(geeValue) == "0"){
      Disabled(c("input.cov", "input.cov_end", "box.input3", "entry.Path3"))
    }else{
      Normal(c("input.cov", "input.cov_end", "box.input3", "entry.Path3"))
      tclvalue(geeValue) <- "0"
      #Disabled(c("input.cov1", "input.cov_end1", "box.input4", "entry.Path4"))
    }
  })
  associationtest.port1.3 <- tkcheckbutton(frame.para4.1.2, background = color, command = function(){
    if (tclvalue(geeValue) == "0" & tclvalue(ftestValue) == "0"){
      #Disabled(c("input.cov1", "input.cov_end1", "box.input4", "entry.Path4"))
      Disabled(c("input.cov", "input.cov_end", "box.input3", "entry.Path3"))
    }else{
      Normal(c("input.cov", "input.cov_end", "box.input3", "entry.Path3"))
      #Normal(c("input.cov1", "input.cov_end1", "box.input4", "entry.Path4"))
      tclvalue(ftestValue) <- "0"
      #Disabled(c("input.cov", "input.cov_end", "box.input3", "entry.Path3"))
    }
  })
  input.cov <- tklabel(frame.para4.1.2, background=color, font=fonttitle, text="(Covariates data directory:", state = "disabled")
  input.cov_end <- tklabel(frame.para4.1.2, background=color, font=fonttitle, text=")", state = "disabled")
  box.input3 <- tkbutton(frame.para4.1.2, text="...", state = "disabled", command=function(){
    tclvalue(dir.cov) <- tkgetOpenFile()
  })
  
  ### to check if users had set the covariate selection
  cova.select <- FALSE
  
  entry.Path3 <- tkentry(frame.para4.1.2, width="20", textvariable = dir.cov, state = "disabled", xscrollcommand = function(...){
    cova.select <<- FALSE
    ### if the directory is not a folder neither a excute file, normal the "Info" button
    if (file.exists(tclvalue(dir.cov))){
      if (!all(is.na(file.info(tclvalue(dir.cov)))) & !file.info(tclvalue(dir.cov))[, "isdir"]){
        #Normal("input.cov.info")
  ### interface for variable selection
        covpath.name <<- readLines(paste(tclvalue(dir.cov)), n = 1)
        covpath.name <<- strsplit(covpath.name, "\t", fixed = TRUE)[[1]]
        if (length(covpath.name) <= 6){
          msg <<- paste("Dear LOHAS users, please provide a valid phenotype-covariate file (please refer to Chapter 6.4 in user guide for variable nomenclature).")
          print(msg)
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
        }else{
          cova.select <<- FALSE
          cova_select <<- tktoplevel()
          tkwm.title(cova_select, "Covariate Select")
          choose_frame <<- tkframe(cova_select)
          var_frame <<- tkframe(choose_frame)
          var_scr <<- tkscrollbar(var_frame, repeatinterval = 5, command = function(...) tkyview(var_tl,...))
          var_tl <<- tklistbox(var_frame, height = 15, selectmode = "extended", yscrollcommand = function(...) tkset(var_scr,...), background = "white")
          tkgrid(tklabel(var_frame, text = "Variable"))
          tkgrid(var_tl, var_scr)
          tkgrid.configure(var_scr, rowspan = 4, sticky = "nsw")

          covpath.name <<- covpath.name[-(1:6)]
          model_cov <<- covpath.name
          tr.id <<- which(substr(covpath.name, 1, 1) %in% c("Q", "B"))
          cov.id <<- which(substr(covpath.name, 1, 1) %in% c("C", "D"))
          model_var_idx <<- cov.id
          model_cov_idx <<- tr.id
          model_var_idx <<- c()
          model_cov_idx <<- sort(c(cov.id, tr.id))

          #if (length(cov.id) > 0){
          #  for (x in cov.id){
          #  #for(x in ancova_var_idx){
          #    tkinsert(var_tl, "end", model_cov[x])
          #    #tkinsert(var_tl, "end", covpath.name[x])
          #  }
          #}else{}

          cova_frame <<- tkframe(choose_frame)
          cova_scr <<- tkscrollbar(cova_frame, repeatinterval = 5, command = function(...) tkyview(cova_tl,...))
          cova_tl <<- tklistbox(cova_frame, height = 15, selectmode = "extended", yscrollcommand = function(...) tkset(cova_scr, ...), background = "white")
          tkgrid(tklabel(cova_frame, text = "Covariate"))
          tkgrid(cova_tl, cova_scr)
          tkgrid.configure(cova_scr, rowspan = 4, sticky = "nsw")
          in_out <<- tkframe(choose_frame)

          if (length(c(cov.id, tr.id)) > 0){ # Default setting: all variables are put into model
          #if (length(tr.id) > 0){
            for (x in sort(c(cov.id, tr.id))){
          #  for (x in tr.id){
              tkinsert(cova_tl, "end", model_cov[x])
            }
          }else{}

          tkpack(tkbutton(in_out, text = "  >>  ", command = function(...){
            if(length(as.integer(tkcurselection(var_tl))) != 0){
              varIndex <<- as.integer(tkcurselection(var_tl))
              for(x in 1:length(varIndex)){
                tkdelete(var_tl, varIndex[x]-x+1)
              }
              model_cov_idx <<- sort(c(model_cov_idx, model_var_idx[varIndex+1]))
              for(x in varIndex){
                tkinsert(cova_tl, which(model_cov_idx==model_var_idx[x+1])-1, model_cov[model_var_idx[x+1]])
              }
              model_var_idx <<- model_var_idx[-(varIndex+1)]
            }
          }))

          tkpack(tklabel(in_out, text = "     "))
          tkpack(tkbutton(in_out, text = "  <<  ", command = function(...){
            if(length(as.integer(tkcurselection(cova_tl)))!=0){
              covIndex <<- as.integer(tkcurselection(cova_tl))
              for(x in 1:length(covIndex)){
                tkdelete(cova_tl, covIndex[x]-x+1)
              }
              model_var_idx <<- sort(c(model_var_idx, model_cov_idx[covIndex+1]))
              for(x in covIndex){
                tkinsert(var_tl, which(model_var_idx==model_cov_idx[x+1])-1, model_cov[model_cov_idx[x+1]])
              }
              model_cov_idx <<- model_cov_idx[-(covIndex+1)]
            }
          }))
          
          tkgrid(var_frame, in_out, cova_frame, padx = 10)
          tkgrid(choose_frame)
          # tkpack(in_out, side = "left")
          # tkpack(cova_frame, side = "left")
          fr_next <- tkframe(cova_select)
          
          tkpack(tkbutton(fr_next, text = "  Next  ", command = function(...){
            cov.select.ok <- 0
            if (!any(model_cov_idx %in% tr.id)){
              msg <<- "Dear LOHAS users, please choose at least one phenotype variable (please refer to Chapter 6.4 in user guide for variable nomenclature)."
              print(msg)
              #write(msg, paste(outpath, "/Log.txt", sep = ""))
              if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
              tkfocus(cova_select)
              cov.select.ok <- cov.select.ok+1
            }else{}
          
            if (!any(model_cov_idx %in% cov.id)){
              msg <<- "Dear LOHAS users, please choose at least one adjusted variable (please refer to Chapter 6.4 in user guide for variable nomenclature)."
              print(msg)
              #write(msg, paste(outpath, "/Log.txt", sep = ""))
              if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
              tkfocus(cova_select)
              cov.select.ok <- cov.select.ok+1
            }else{}
            
            if (cov.select.ok == 0){
              cova.select <<- TRUE
              ### add back the first 6 pedigree format columns
              cov.id <<- model_cov_idx[model_cov_idx %in% cov.id] + 6
              tr.id <<- model_cov_idx[model_cov_idx %in% tr.id] + 6
              tkdestroy(cova_select)
            }else{}
          }))
          
          tkgrid(fr_next)
          
          fr_info <- tkframe(cova_select)
          tkgrid(tklabel(fr_info, text = "Trait: \"Q\"uantitative ; \"B\"inary."), sticky = "w")
          tkgrid(tklabel(fr_info, text = "Adjusted variable: \"C\"ontinuous ; \"D\"iscrete."), sticky = "w")
          tkgrid(fr_info)
  ### END ### interface for variable selection
        }
      }else{
        #Disabled("input.cov.info")
      }
    }else{
      #Disabled("input.cov.info")
    }
  })
  
  #input.cov1 <- tklabel(frame.para4.1.3, background=color, font=fonttitle, text="(Covariates data directory:", state = "disabled")
  #input.cov_end1 <- tklabel(frame.para4.1.3, background=color, font=fonttitle, text=")", state = "disabled")
  #box.input4 <- tkbutton(frame.para4.1.3, text="...", command=function() tclvalue(dir.cov1) <- tkgetOpenFile(), state = "disabled")
  #entry.Path4 <- tkentry(frame.para4.1.3, width="20", textvariable=dir.cov1, state = "disabled")
  associationtest.port1 <- tkradiobutton(frame.para4.1.1,background=color, command=function(){
    if (tclvalue(groupValue) == "2"){ # two groups allow Wilcoxon test
      tclvalue(associationValue)="1";
      Normal(c("associationtest.port1.1", "associationtest.Wilcoxon", "associationtest.port1.2", "associationtest.Ftest", "associationtest.port1.3", "associationtest.GEE"))
    }else{
      Normal(c("associationtest.port1.2", "associationtest.Ftest", "associationtest.port1.3", "associationtest.GEE"))
    }
  })
  associationtest.port2 <- tkradiobutton(frame.para4.1.1,background=color, command=function(){
  tclvalue(associationValue)="2";
  tclvalue(ftestValue) <- "0"
  tclvalue(geeValue) <- "0"
  Disabled(c("associationtest.port1.1", "associationtest.Wilcoxon", "associationtest.port1.2", "associationtest.Ftest",
  "input.cov", "input.cov_end", "box.input3", "entry.Path3",  "associationtest.port1.3", "associationtest.GEE"))})
  #Disabled(c("associationtest.port1.1", "associationtest.Wilcoxon", "associationtest.port1.2", "associationtest.Ftest",
  #"input.cov", "input.cov_end", "box.input3", "entry.Path3",  "associationtest.port1.3", "associationtest.GEE", "input.cov1", "entry.Path4", "box.input4", "input.cov_end1"))})
  tkconfigure(associationtest.port1.1, variable=associationValue)
  tkconfigure(associationtest.port1.2, variable=ftestValue)
  tkconfigure(associationtest.port1.3, variable = geeValue)
  tkconfigure(associationtest.port1,variable=associationcheck,value="1")
  tkconfigure(associationtest.port2,variable=associationcheck,value="2")
  #associationtestwindow.Label<- tklabel(frame.para4.1.1,background=color,font=fonttitle,text="(")
  #associationtestwindow.LeftLabel <- tklabel(frame.para4.1.1,background=color,font=fonttitle,text="Single-window")
  #associationtestwindow.RightLabel <- tklabel(frame.para4.1.1,background=color,font=fonttitle,text="(Number of windows: ")
  
  #associationtestwindow.port1 <- tkradiobutton(frame.para4.1.1,background=color)
  #associationtestwindow.port2 <- tkradiobutton(frame.para4.1.1,background=color)
  #associationwindowValue <- tclVar("1")
  #tkconfigure(associationtestwindow.port1,variable=associationwindowValue,value="1")
  #tkconfigure(associationtestwindow.port2,variable=associationwindowValue,value="2")
  
  #associationtestwindowsize.val <- tclVar("1")
  #multiwindowsize.box <- tkentry(frame.para4.1.1,width=16,textvariable=associationtestwindowsize.val)
  #associationtestmultiwindow.Label <- tklabel(frame.para4.1.1,background=color,font=fonttitle,text=", between 1 and 10 )")
  
  #tkgrid(associationtest.Label,associationtest.port1,associationtest.LeftLabel,associationtestwindow.RightLabel,multiwindowsize.box,associationtestmultiwindow.Label,associationtest.port2,associationtest.RightLabel )
  tkgrid(associationtest.Label, associationtest.port1, associationtest.LeftLabel, associationtest.port2, associationtest.RightLabel )
  tkgrid(associationtest.Label1, associationtest.port1.1, associationtest.Wilcoxon, associationtest.port1.2, 
    associationtest.Ftest, associationtest.port1.3, associationtest.GEE, input.cov, entry.Path3, box.input3, 
    input.cov_end)
  #tkgrid(associationtest.Label1, associationtest.port1.1, associationtest.Wilcoxon, associationtest.port1.2, associationtest.Ftest, input.cov, entry.Path3, box.input3, input.cov_end)
  #tkgrid(associationtest.Label2, associationtest.port1.3, associationtest.GEE, input.cov1, entry.Path4, box.input4, input.cov_end1)
  #tkgrid(associationtest.Label, associationtest.port1, associationtest.LeftLabel , associationtest.port2, associationtest.RightLabel )
  
  # Regression test # load in covariates of F-test in regresion analysis 20120102
  
  #dir.cov <- tclVar("NA")
  #ftest.Label <- tklabel(frame.para4.3, background=color, font=fonttitle, text="    F-test in regression analysis:    ")
  #tkgrid(ftest.Label)
  #input.cov <- tklabel(frame.para4.3, background=color, font=fonttitle, text="          -Covariates data directory:    ")
  #box.input3 <- tkbutton(frame.para4.3, text="...", command=function() tclvalue(dir.cov) <- tkgetOpenFile())
  #entry.Path3 <- tkentry(frame.para4.3, width="50", textvariable=dir.cov)
  #tkgrid(input.cov, entry.Path3, box.input3)
  
  # Data visualzation #
  
  DataVisual.port1 <- tkradiobutton(frame.para4.2,background=color,command=function(){tclvalue(Plotselection1.val)=1;tclvalue(Plotselection2.val)=1;tclvalue(Plotselection3.val)=1;
      Normal(c("Plotselection1.lab", "Plotselection2.lab", "Plotselection3.lab", "Plotselection1.port", "Plotselection2.port", "Plotselection3.port"))
      Normal(c("QLOHintensity.Label", "QLOHintensity.box", "QLOHintensity.ps"))
      Normal(c("Alphascaling.Label", "Alphascaling.box", "Alphascaling.ps"))
   })
  DataVisual.port2 <- tkradiobutton(frame.para4.2,background=color,command=function(){tclvalue(Plotselection1.val)=0;tclvalue(Plotselection2.val)=0;tclvalue(Plotselection3.val)=0;
      Disabled(c("Plotselection1.lab", "Plotselection2.lab", "Plotselection3.lab", "Plotselection1.port", "Plotselection2.port", "Plotselection3.port"))
      Disabled(c("QLOHintensity.Label", "QLOHintensity.box", "QLOHintensity.ps"))
      Disabled(c("Alphascaling.Label", "Alphascaling.box", "Alphascaling.ps"))
  })
  DataVisualValue <- tclVar("1")
  DataVisualValue.unclick <- tclVar("2")
  
  tkconfigure(DataVisual.port1,variable=DataVisualValue,value="1")
  tkconfigure(DataVisual.port2,variable=DataVisualValue,value="2")
  DataVisual.LeftLabel <- tklabel(frame.para4.2,background=color,font=fonttitle,text="Yes")
  DataVisual.RightLabel <- tklabel(frame.para4.2,background=color,font=fonttitle,text="No")
  DataVisual.Label <- tklabel(frame.para4.2,background=color,font=fonttitle,text="     Data visualization:                ")
  tkgrid(DataVisual.Label,DataVisual.port1,DataVisual.LeftLabel,DataVisual.port2,DataVisual.RightLabel)
  
  Plotselection1.lab <- tklabel(frame.para5,background=color,font=fonttitle,text="Indivdual LOH plot")
  Plotselection2.lab <- tklabel(frame.para6,background=color,font=fonttitle,text="Combined LOH plot")
  Plotselection3.lab <- tklabel(frame.para7,background=color,font=fonttitle,text="LOH biplot")
  
  Plotselection1.port <- tkcheckbutton(frame.para5,background=color,command=function(){tclvalue(DataVisualValue)=1;tclvalue(LOHintensityValue)=1;})
  Plotselection2.port <- tkcheckbutton(frame.para6,background=color,command=function(){tclvalue(DataVisualValue)=1;tclvalue(LOHintensityValue)=1;})
  Plotselection3.port <- tkcheckbutton(frame.para7,background=color,command=function(){tclvalue(DataVisualValue)=1;tclvalue(LOHintensityValue)=1;})
  
  Plotselection1.val <- tclVar("1")
  Plotselection2.val <- tclVar("1")
  Plotselection3.val <- tclVar("1")
  
  tkconfigure(Plotselection1.port,variable=Plotselection1.val)
  tkconfigure(Plotselection2.port,variable=Plotselection2.val)
  tkconfigure(Plotselection3.port,variable=Plotselection3.val)
  
  QLOHintensity.val<- tclVar("0.9")
  Alphascaling.val<- tclVar("0")
  
  QLOHintensity.Label <- tklabel(frame.para6,background=color,font=fonttitle,text=" (Quantile of reference LOH intensity:")
  QLOHintensity.box <- tkentry(frame.para6,width=6,textvariable=QLOHintensity.val)
  QLOHintensity.ps <- tklabel(frame.para6,background=color,font=fonttitle,text=" , between 0 and 1 )")
  
  Alphascaling.Label <- tklabel(frame.para7,background=color,font=fonttitle,text=" (Alpha scaling: ")
  Alphascaling.box <- tkentry(frame.para7,width=6,textvariable=Alphascaling.val)
  Alphascaling.ps <- tklabel(frame.para7,background=color,font=fonttitle,text=" , between 0 and 1 )")
  
  Plotselection.space1<- tklabel(frame.para5,background=color,text="                                                      ")
  Plotselection.space2<- tklabel(frame.para6,background=color,text="                                                      ")
  Plotselection.space3<- tklabel(frame.para7,background=color,text="                                                      ")
  
  tkgrid(Plotselection.space1,Plotselection1.port,Plotselection1.lab)
  tkgrid(Plotselection.space2,Plotselection2.port,Plotselection2.lab,QLOHintensity.Label,QLOHintensity.box,QLOHintensity.ps)
  tkgrid(Plotselection.space3,Plotselection3.port,Plotselection3.lab,Alphascaling.Label,Alphascaling.box,Alphascaling.ps)
  
  ##### check selections are all right and run LOH GUI #####
  #figure.path<-paste(strsplit(LOHgui,"\\LOH_gui.R",fixed = TRUE),sep="")
  #figure.run.step1<-paste(figure.path,"\\run1.gif",sep="")
  #figure.run.step2<-paste(figure.path,"\\run2.gif",sep="")
  #figure.run.step3<-paste(figure.path,"\\run3.gif",sep="")
  #image1 <- tclVar()
  #tcl("image","create","photo",image1,file=figure.run.step1)
  
  ### LOH folder directory 20111227
  #LOH.file <- paste(strsplit(LOHgui,"\\LOH_gui_new110307.r",fixed = TRUE))
  LOH.file <- gsub(paste("/LOHAS_Interface.r", sep = ""), "", gsub("\\\\", "/", LOHAS.gui))
  
  #debug.path <- paste(strsplit(LOHgui,"\\LOH_gui_new110307.r",fixed = TRUE),"\\debug110307.r",sep="")
  # debug.path <- paste(LOH.file, "/debug", tempdate, ".r", sep = "")
  #Run.button <- tkbutton(LOH,text="  Run  ",image=image1,
  #			  command=function(){cat("Please wait a moment...\n")
  #						   source(debug.path)})
  
  #Run.button <- tkbutton(frame,text="  Run  ",
  #			  command=function(){ Match <<- FALSE; cat("Please wait a moment...\n")
  #			      if (file.exists(paste(gsub("\\\\", "/", paste(tclvalue(dir.input), "/", sep = "")), "/Inpath.txt", sep = "")) & file.exists(paste(gsub("\\\\", "/", paste(tclvalue(dir.input), "/", sep = "")), "/Datasetting.txt", sep = ""))){
  #                                source(debug.path)
  #                              }else{
  #                                tkmessageBox(title = "Warning", message = "Dear LOHAS users, please provide valid parameter setting files in the input directory.")
  #                                tkfocus(tn)
  #                              }
  #						   })
  
  Run.button <- tkbutton(frame, text = "  Run  ", command = function(){
    cat("Please wait a moment...\n")
    Match <<- FALSE
    Cluster <<- FALSE
    PlotAllinOne <<- TRUE
    PlotOneinOne <<- TRUE
    RV <<- FALSE
    sequencing <<- FALSE
    MemoryEnough <<- TRUE
    Deg <<- 2
    SN <<- 100000 # maximum number of markers each time in Regression test or GEE
    ###############################
    ### version: 1     ############
    ###############################
  
    log.file <<- paste("Beginning time for LOHAS: ", date(), sep = "")
    print(log.file)
  
    if(!Match){
      ###########################################################
      #####         alpha and chromosome number(1 ~ 23)     #####
      ###########################################################
  
      #for(x in c(1:23, "all")){ # chr23 are not analyzed 20130109
      for(x in c(1:22, "all")){
        assign(paste("chr.check", x, sep = ""), as.numeric(tclvalue(get(paste("chr", x,".val", sep = "")))), envir = .GlobalEnv)
      }
  
      if(chr.checkall == 1){
        #Chr.num <- 1:23
        Chr.num <<- 1:22
      }
      if(chr.checkall != 1){
        #Chr.num <- which(sapply(1:23, function(x) get(paste("chr.check", 1:23, sep="")[x]))==1)
        Chr.num <<- which(sapply(1:22, function(x) get(paste("chr.check", 1:22, sep = "")[x])) == 1)
      }
      ##############################################################################
      inpath <<- tclvalue(tkget(entry.Path))
      outpath <<- tclvalue(tkget(entry.Path2))
      ##############################################################################
      num.group <<- as.numeric(tclvalue(groupValue))
      ethnicity <<- patient <<- NULL
      if(num.group == 1) ethnicity <<- as.numeric(tclvalue(ethnicity.val))
      if(num.group == 2) patient <<- as.numeric(tclvalue(patient1.val))
      ##############################################################################
      chiptype <<- as.numeric(tclvalue(chipValue))
      chip1 <<- as.numeric(tclvalue(chip1.val))
      chip2 <<- as.numeric(tclvalue(chip2.val))
      chipMer <<- as.numeric(tclvalue(chip1_2.val))
      ##############################################################################
      description <<- (tclvalue(statsValue) == "1")
      doLOH <<- (tclvalue(LOHintensityValue) == "1")
      SNPthin <<- (tclvalue(SNPthinningvalue) == "1")
      SNPthin.v <<- as.numeric(tclvalue(SNPthinningsize.val))
      ##############################################################################
      SNPselectype <<- as.numeric(tclvalue(SNPselectionValue))
      propsize <<- ifelse(SNPselectype == "1", as.numeric(tclvalue(fixSNPsize.val)), as.numeric(tclvalue(fixSNPprop.val)))
      doRanktest <<- (tclvalue(chromotestValue) == "1")
      doAssotest <<- (tclvalue(associationValue) == "1" )
      doFtest <<- (tclvalue(ftestValue) == "1")
      doGEE <<- (tclvalue(geeValue) == "1")
      covpath <<- tclvalue(tkget(entry.Path3))
      ##############################################################################
      visual.data <<- (as.numeric(tclvalue(DataVisualValue)) == 1)
      plotselect1 <<- tclvalue(Plotselection1.val) == "1"
      plotselect2 <<- tclvalue(Plotselection2.val) == "1"
      plotselect3 <<- tclvalue(Plotselection3.val) == "1"
      QLOHselect <<- as.numeric(tclvalue(QLOHintensity.val))
      Alphaselect <<- as.numeric(tclvalue(Alphascaling.val))
      ##############################################################################
    } # !Match
    ###############################################################################
    ###############################################################################
    ### version 2 #################################################################
    ###############################################################################
    ###############################################################################
    if(Match){
      ###########################################################
      #####         alpha and chromosome number(1 ~ 23)     #####
      ###########################################################
  
      for(x in c(1:22, "all")){
        assign(paste("chr.check", x, sep=""), as.numeric(tclvalue(get(paste("chr", x,".val2", sep="")))), envir = .GlobalEnv)
      }
      if(chr.checkall==1){
        Chr.num <- 1:22
      }
      if(chr.checkall!=1){
        Chr.num <- which(sapply(1:22, function(x) get(paste("chr.check", 1:22, sep="")[x]))==1)
      }
  
      ##############################################################################
      inpath <- tclvalue(tkget(entry2.Path))
      outpath <- tclvalue(tkget(entry2.Path2))
      ##############################################################################
      num.group <- 1; ethnicity <- 2; patient <- NULL
      number.popu <- as.numeric(tclvalue(group2.val))
  
      Data.column <- read.table(paste(inpath,"/Datasetting.txt",sep=""), as.is=TRUE, header=FALSE, sep="\t")
  
      if(number.popu != ncol(Data.column)-1){
        msg <- "Please check the input file: the number of phase is not the same in datasetting file and the number of study groups you selected."
        if (length(args) <= 0) tkmessageBox(title=" Error ",message=msg,icon="error", type="ok")
        ok <- ok+1
      }
  
      ##############################################################################
      chiptype <- as.numeric(tclvalue(chipValue2))
      chip1 <- as.numeric(tclvalue(chip1.val2))
      chip2 <- as.numeric(tclvalue(chip2.val2))
      chipMer <- as.numeric(tclvalue(chip1_2.val2))
  
      ##############################################################################
      description <- (tclvalue(statsValue2)=="1")
      doLOH <- (tclvalue(LOHintensityValue2)=="1")
      SNPthin <- (tclvalue(SNPthinningvalue2)=="1")
      SNPthin.v <- as.numeric(tclvalue(SNPthinningsize.val2))
      ##############################################################################
      SNPselectype <- as.numeric(tclvalue(SNPselectionValue2))
      propsize <- ifelse(SNPselectype=="1", as.numeric(tclvalue(fixSNPsize.val2)), as.numeric(tclvalue(fixSNPprop.val2)))
      doRanktest <- (tclvalue(chromotestValue2)=="1")
      doAssotest <- (tclvalue(associationValue2)=="1" )
      doFtest <- (tclvalue(ftestValue2) == "1")
      ##############################################################################
      visual.data <- as.numeric(tclvalue(DataVisualValue2))
      plotselect1 <- tclvalue(Plotselection1.val2) == "1"
      plotselect2 <- tclvalue(Plotselection2.val2) == "1"
      plotselect3 <- tclvalue(Plotselection3.val2) == "1"
      QLOHselect <- as.numeric(tclvalue(QLOHintensity.val2))
      Alphaselect <- as.numeric(tclvalue(Alphascaling.val2))
      ##############################################################################
    } # Match
  
    ###############################################
    #####     Example 1-4 setting 20130109    #####
    ###############################################
    # Example 1
    if (toupper(inpath) == "EXAMPLE1"){
      inpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Input", sep = "")
      dir.create(inpath, showWarnings = FALSE)
      tempinpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Examples/Example1/", sep = "")
      write.table(tempinpath, paste(inpath, "/Inpath.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
      tempdatasetting <<- matrix(c("Example1", "1", "2", "3", "4", "", "", "", "1", "7"), ncol = 1)
      write.table(tempdatasetting, paste(inpath, "/Datasetting.txt", sep = ""),
        row.names = c("Population_name", "SNP_ID", "Chromosome", "Physical_position", "ABcall", "RS_number", "MAF_threshold", "MAF_weight", "Start_field", "End_field"),
        col.names = FALSE, quote = FALSE, sep = "\t")
      #rm(tempinpath, tempdatasetting)
    }else{}
    # Example 2
    if (toupper(inpath) == "EXAMPLE2"){
      inpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Input", sep = "")
      dir.create(inpath, showWarnings = FALSE)
      tempinpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Examples/Example2/", c("CEU/", "CHB/", "YRI/"), sep = "")
      write.table(tempinpath, paste(inpath, "/Inpath.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
      tempdatasetting <<- matrix(c(c("CEU", "CHB", "YRI"), rep(c("1", "2", "3", "4", "", "", "", "1", "7"), each = 3)), ncol = 3, byrow = TRUE)
      write.table(tempdatasetting, paste(inpath, "/Datasetting.txt", sep = ""),
        row.names = c("Population_name", "SNP_ID", "Chromosome", "Physical_position", "ABcall", "RS_number", "MAF_threshold", "MAF_weight", "Start_field", "End_field"),
        col.names = FALSE, quote = FALSE, sep = "\t")
      #rm(tempinpath, tempdatasetting)
    }else{}
    # Example 3
    if (toupper(inpath) == "EXAMPLE3"){
      inpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Input", sep = "")
      dir.create(inpath, showWarnings = FALSE)
      tempinpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Examples/Example3/", c("Case/", "Control/"), sep = "")
      write.table(tempinpath, paste(inpath, "/Inpath.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
        tempdatasetting <<- matrix(c(c("Case", "Control"), rep(c("1", "2", "3", "4", "5", "", "", "1", "9"), each = 2)), ncol = 2, byrow = TRUE)
      write.table(tempdatasetting, paste(inpath, "/Datasetting.txt", sep = ""),
        row.names = c("Population_name", "SNP_ID", "Chromosome", "Physical_position", "ABcall", "RS_number", "MAF_threshold", "MAF_weight", "Start_field", "End_field"),
        col.names = FALSE, quote = FALSE, sep = "\t")
      #rm(tempinpath, tempdatasetting)
    }else{}
    # Example 4
    if (toupper(inpath) == "EXAMPLE4"){
      inpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Input", sep = "")
      dir.create(inpath, showWarnings = FALSE)
      tempinpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Examples/Example4/", c("Case/", "Control/"), sep = "")
      write.table(tempinpath, paste(inpath, "/Inpath.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
      tempdatasetting <<- matrix(c(c("Case", "Control"), rep(c("1", "2", "3", "4", "", "0.05", "N", "1", "7"), each = 2)), ncol = 2, byrow = TRUE)
      write.table(tempdatasetting, paste(inpath, "/Datasetting.txt", sep = ""),
        row.names = c("Population_name", "SNP_ID", "Chromosome", "Physical_position", "ABcall", "RS_number", "MAF_threshold", "MAF_weight", "Start_field", "End_field"),
        col.names = FALSE, quote = FALSE, sep = "\t")
      #rm(tempinpath, tempdatasetting)
    }else{}
    # Example 5 PlotExample
    if (toupper(inpath) == "EXAMPLE5"){
      inpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Input", sep = "")
      dir.create(inpath, showWarnings = FALSE)
      tempinpath <<- paste(gsub("/Program$", "", gsub("\\\\", "/", LOH.file)), "/Examples/Example4/", c("Case/", "Control/"), sep = "")
      write.table(tempinpath, paste(inpath, "/Inpath.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
      tempdatasetting <<- matrix(c(c("Case", "Control"), rep(c("1", "2", "3", "4", "", "0.05", "Y", "1", "7"), each = 2)), ncol = 2, byrow = TRUE)
      write.table(tempdatasetting, paste(inpath, "/Datasetting.txt", sep = ""),
        row.names = c("Population_name", "SNP_ID", "Chromosome", "Physical_position", "ABcall", "RS_number", "MAF_threshold", "MAF_weight", "Start_field", "End_field"),
        col.names = FALSE, quote = FALSE, sep = "\t")
    }else{}
  
    dir.create(paste(outpath, sep = ""), showWarnings = FALSE, recursive = TRUE)
    #######################################################
    #     check for input and output directory            #
    #######################################################
    ok <<- 0
  
    if(inpath == ""){
      #msg <- "Please check the input directory."
      msg <<- "Dear LOHAS users, please provide a valid input directory name."
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon="error", type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }
    #if (ok == 0 & (!file.exists(paste(gsub("\\\\", "/", paste(tclvalue(dir.input), "/", sep = "")), "/Inpath.txt", sep = "")) | !file.exists(paste(gsub("\\\\", "/", paste(tclvalue(dir.input), "/", sep = "")), "/Datasetting.txt", sep = "")))){
    #  tkmessageBox(title = "Warning", message = "Dear LOHAS users, please provide valid parameter setting files in the input directory.")
    #  tkfocus(tn)
    #  ok <<- ok + 1
    #}else{}
    if (ok == 0 & (!file.exists(paste(gsub("\\\\", "/", paste(inpath, "/", sep = "")), "/Inpath.txt", sep = "")) | !file.exists(paste(gsub("\\\\", "/", paste(inpath, "/", sep = "")), "/Datasetting.txt", sep = ""))) & !plot.only){
      msg <<- "Dear LOHAS users, please provide valid parameter setting files in the input directory."
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ", message=msg, icon="error", type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }else{}
  
    if(ok == 0 & outpath == ""){
      #msg <- "Please check the output directory."
      msg <<- "Dear LOHAS users, please provide a valid output directory name."
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ", message=msg, icon="error", type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    if (ok == 0 & !(num.group %in% c(1, 2))){
      msg <<- paste("Dear LOHAS users, please provide a valid number of group(s), for example, ", ifelse(length(args) > 0, "--num.group=1", "num.group <- 1"), " indicates a one-group-study; and ", ifelse(length(args) > 0, "--num.group=2", "num.group <- 2"), " indicates a two-group-study.", sep = "")
      print(msg)
      write(msg, paste(outpath, "/Log.txt", sep = ""))
      if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    if (ok == 0 & num.group == 1){
      if (is.null(ethnicity)){
        msg <<- paste("Dear LOHAS users, please provide a valid setting of population(s) in single group study, for example, ", ifelse(length(args) > 0, "--ethnicity=1", "ethnicity <- 1"), " for single population; and ", ifelse(length(args) > 0, "--ethnicity=2", "ethnicity <- 2"), " for multiple populations.", sep = "")
        print(msg)
        write(msg, paste(outpath, "/Log.txt", sep = ""))
        if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
        tkfocus(tn)
        ok <<- ok+1
      }else{
        if (!(ethnicity %in% c(1, 2))){
          msg <<- paste("Dear LOHAS users, please provide a valid setting of population(s) in single group study, for example, ", ifelse(length(args) > 0, "--ethnicity=1", "ethnicity <- 1"), " for single population; and ", ifelse(length(args) > 0, "--ethnicity=2", "ethnicity <- 2"), " for multiple populations.", sep = "")
          print(msg)
          write(msg, paste(outpath, "/Log.txt", sep = ""))
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok+1
        }else{}
      }
    }
  
    if (ok == 0 & num.group == 2){
      if (is.null(patient)){
        msg <<- paste("Dear LOHAS users, please provide a valid indicator of case group in two groups study, for example, ", ifelse(length(args) > 0, "--patient=1", "patient <- 1"), " indicates the 1st group is the case group; and --patient=2 indicates the 2nd group is the case group.", sep = "")
        print(msg)
        write(msg, paste(outpath, "/Log.txt", sep = ""))
        if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
        tkfocus(tn)
        ok <<- ok+1
      }else{
        if (!(patient %in% c(1, 2))){
          msg <<- paste("Dear LOHAS users, please provide a valid indicator of case group in two groups study, for example, ", ifelse(length(args) > 0, "--patient=1", "patient <- 1"), " indicates the 1st group is the case group; and --patient=2 indicates the 2nd group is the case group.", sep = "")
          print(msg)
          write(msg, paste(outpath, "/Log.txt", sep = ""))
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok+1
        }else{}
      }
    }
    ##########################################
    #     check for selection of chips
    ##########################################
    if (ok == 0 & !(chiptype %in% c(1, 2, 3, 4))){
      msg <<- paste("Dear LOHAS users, please provide a valid type of gene chip, for example, ", ifelse(length(args) > 0, "--chiptype=1", "chiptype <- 1"), " indicates the Affy 100K; ", ifelse(length(args) > 0, "--chiptype=2", "chiptype <- 2"), " indicates the Affy 500K; ", ifelse(length(args) > 0, "--chiptype=3", "chiptype <- 3"), " indicates the Affy 6.0; and ", ifelse(length(args) > 0, "--chiptype=4", "chiptype <- 4"), " indicates the customized chip.", sep = "")
      print(msg)
      write(msg, paste(outpath, "/Log.txt", sep = ""))
      if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    if(ok == 0 & (chiptype %in% 1:2) & chip1 == 0 & chip2 == 0 & chipMer == 0 ){
      #msg <- "Please choose the chips."
      msg <<- "Dear LOHAS users, please specify which chips (chip1, chip2 and/or merged-chip) are analyzed."
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ",message=msg,icon="error",type="ok")
      tkfocus(tn)
      ok <<- ok+1
     }
  
    if(chiptype != 4){
      if(doLOH & ok == 0 & length(Chr.num) == 0){
        #msg <- "Please choose chromosomes."
        msg <<- "Dear LOHAS users, please specify which chromosomes are analyzed."
        print(msg)
        write(msg, paste(outpath,"/Log.txt",sep=""))
        if (length(args) <= 0) tkmessageBox(title=" Error ", message=msg, icon="error", type="ok")
        tkfocus(tn)
        ok <<- ok+1
      }
    }
    
    if (ok == 0 & plot.only){
      if (length(Chr.num) == 0){
        msg <<- "Dear LOHAS users, please specify which plots of chromosomes are drawn."
        print(msg)
        write(msg, paste(outpath,"/Log.txt",sep=""))
        if (length(args) <= 0) tkmessageBox(title=" Error ", message=msg, icon="error", type="ok")
        tkfocus(tn)
        ok <<- ok+1
      }else{}
    }else{}
  
    if (!doLOH & !plot.only & any(doRanktest, doAssotest, doFtest, doGEE, plotselect1, plotselect2, plotselect3)){ # must estimate LOH intensity before Association test and plots
      #msg <- "Please choose estimating LOH intensity to do association test or data visualization."
      msg <<- "Dear LOHAS users, please activate the estimation of LOH/LCSH intensity if you would like to run a chromosome ranking test, association test or data visualization."
      print(msg)
      write(msg, paste(outpath, "/Log.txt", sep = ""))
      if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
      tkfocus(tn)
      ok <<- ok + 1
    }else{
      if (num.group == 1){
        patient <<- NULL
        if (doAssotest){ # Wilcoxon only provides for case-control study
          #msg <- "Wilcoxon test only provides for case-control study."
          msg <<- "Dear LOHAS users, please choose a two-group analysis if you would like to run a Wilcoxon test."
          print(msg)
          write(msg, paste(outpath, "/Log.txt", sep = ""))
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok + 1
        }
        if (!is.null(ethnicity)){
          if (ethnicity == 2 & (doFtest | doGEE)){
            #msg <- "Association tests do not provide in one group multiple populations."
            msg <<- paste("Dear LOHAS users, please choose a one-group-single-population or two-group analysis if you would like to run a ", ifelse(doFtest, "regression", "GEE"), " analysis", sep = "")
            print(msg)
            write(msg, paste(outpath, "/Log.txt", sep = ""))
            if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
            tkfocus(tn)
            ok <<- ok + 1
          }
        }
      }else{
        ethnicity <<- NULL
      }
      if (all(doFtest, doGEE)){ # Regression and GEE models are mutually exclusive
        #msg <- "Regression model based on independent sample and GEE model is for family based sample."
        msg <<- "Dear LOHAS users, please only choose one analysis from regression or GEE according to your data attribute."
        print(msg)
        write(msg, paste(outpath, "/Log.txt", sep = ""))
        if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
        tkfocus(tn)
        ok <<- ok + 1
      }
    }
  
    ##############################################################################
    if (ok == 0 & !(SNPselectype %in% c(1, 2))){
      msg <<- paste("Dear LOHAS users, please provide a valid code of the method for determining a window size, for example, ", ifelse(length(args) > 0, "--SNPselectype=1", "SNPselectype <- 1"), " indicates that window size is given by a fixed SNP number; ", ifelse(length(args) > 0, "--SNPselectype=2", "SNPselectype <- 2"), " indicates that window size is defined by a fixed proportion of SNPs on a chromosome.", sep = "")
      print(msg)
      write(msg, paste(outpath, "/Log.txt", sep = ""))
      if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    if(ok == 0 & SNPselectype == 1 & propsize <= 1){
      #msg <- "Wrong number of SNPs."
      msg <<- "Dear LOHAS users, please input a numeral greater than 1 for the number of SNPs in a window."
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ",message=msg, icon="error",type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    if(ok == 0 & SNPselectype == 2 & (propsize <= 0 | propsize > 1)){
      #msg <- "Wrong proportion of SNPs."
      msg <<- paste("Dear LOHAS users, please input a numeral ", ifelse(propsize <= 0, "greater than 0", "not greater than 1"), " for the proportion of SNPs on a chromosome.", sep = "")
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ",message=msg, icon="error",type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }
    ################################################################################
    if(ok == 0 & plotselect2 == 1 & (QLOHselect <= 0 | QLOHselect >= 1)){
      #msg <- "Wrong Quantitle of reference LOH intensity."
      msg <<- paste("Dear LOHAS users, please input a numeral ", ifelse(QLOHselect <= 0, "greater than 0", "less than 1"), " for the quantile of reference LOH/LCSH intensity.", sep = "")
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ",message=msg,icon="error",type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    if(ok == 0 & plotselect3 == 1 & (Alphaselect < 0 | Alphaselect > 1)){
      #msg <- "Wrong Alpha scaling value."
      msg <<- paste("Dear LOHAS users, please input a numeral ", ifelse(Alphaselect < 0, "not less than 0", "not greater than 1"), " for the alpha scaling value.", sep = "")
      print(msg)
      write(msg, paste(outpath,"/Log.txt",sep=""))
      if (length(args) <= 0) tkmessageBox(title=" Error ",message=msg,icon="error",type="ok")
      tkfocus(tn)
      ok <<- ok+1
    }
  
    ################################################################################
    # check for variable names in covariate data( header name "Q", "B", "C", and "D" are required) 20130122
    if (ok == 0 & (doFtest | doGEE)){
      #if (!file.exists(covpath)){ # show error when covariate data was not found
      if (all(is.na(file.info(covpath))) | file.info(covpath)[, "isdir"]){
        #msg <- "Invalid path of covariates data."
        msg <<- "Dear LOHAS users, please provide a valid path of the phenotype-covariate file."
        print(msg)
        write(msg, paste(outpath, "/Log.txt", sep = ""))
        if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
        tkfocus(tn)
        ok <<- ok + 1
      }else{
        covpath.name <<- readLines(covpath, n = 1)
        covpath.name <<- strsplit(covpath.name, "\t", fixed = TRUE)[[1]]
        if (length(covpath.name) <= 6){
          msg <<- paste("Dear LOHAS users, please provide a valid phenotype-covariate file (please refer to Chapter 6.4 in user guide for variable nomenclature).")
          print(msg)
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok + 1
        }
        cov.data <<- read.table(covpath, header = TRUE, sep = "\t", as.is = TRUE)
        covpath.name <<- covpath.name[-(1:6)]
        if (doGEE & length(unique(cov.data[, 1])) <= 1){ # number of families should greater than 1 when doGEE
          #msg <- "The number of families must greater than 1 when running GEE analysis."
          msg <<- "Dear LOHAS users, cannot run a GEE analysis because your data conain only one family."
          print(msg)
          write(msg, paste(outpath, "/Log.txt", sep = ""))
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok + 1
        }
        if (!all(substr(covpath.name, 1, 1) %in% c("Q", "B", "C", "D"))){ # all covariate headers must be on of "Q", "B", "C", or "D"
          #msg <- "Invalid variable names in covariate data."
          msg <<- "Dear LOHAS users, please provide valid names of covariate variables (please refer to Chapter 6.4 in user guide for variable nomenclature)."
          print(msg)
          write(msg, paste(outpath,"/Log.txt",sep=""))
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok + 1
        }
        if (!any(substr(covpath.name, 1, 1) %in% c("Q", "B"))){ # there must be trait(s) in covariate data
          #msg <- "Please specify the trait variable(s)."
          msg <<- "Dear LOHAS users, please provide valid names of phenotype variables (please refer to Chapter 6.5 in user guide for variable nomenclature)."
          print(msg)
          write(msg, paste(outpath, "/Log.txt", sep = ""))
          if (length(args) <= 0) tkmessageBox(title = " Error ", message = msg, icon = "error", type = "ok")
          tkfocus(tn)
          ok <<- ok + 1
        }
        rm(covpath.name)
      }
    }
  
  
    ################################################################################
    if(ok == 0){
      #prog.path <- paste(strsplit(debug.path,"\\debug110307.r",fixed = TRUE),"\\LOHAS110307.r",sep="");
      write(log.file, paste(outpath,"/Log.txt",sep=""))
      prog.path <<- paste(LOH.file, "/LOHAS_MainFunction.r", sep = "")
      source(prog.path);
    }
    ################################################################################
  
  })
  
  
  
  ##### permutation the GUI #####
  tkgrid(frame.title)
  tkgrid.configure(frame.title,sticky="w")
  tkgrid(frame.up)
  tkgrid.configure(frame.up,sticky="w")
  tkgrid(frame.path)
  tkgrid.configure(frame.path,sticky="w")
  tkgrid(frame.group)
  tkgrid.configure(frame.group,sticky="w")
  tkgrid(frame.group.1)
  tkgrid.configure(frame.group.1,sticky="w")
  tkgrid(frame.group.1.1)
  tkgrid.configure(frame.group.1.1,sticky="w")
  #tkgrid(frame.group.1.2)
  #tkgrid.configure(frame.group.1.2,sticky="w")
  
  tkgrid(frame.chip)
  tkgrid.configure(frame.chip,sticky="w")
  tkgrid(frame.chip.1)
  tkgrid.configure(frame.chip.1,sticky="w")
  
  tkgrid(frame.chip.2)
  tkgrid.configure(frame.chip.2,sticky="w")
  tkgrid(frame.para)
  tkgrid.configure(frame.para,sticky="w")
  tkgrid(frame.para2)
  tkgrid.configure(frame.para2,sticky="w")
  tkgrid(frame.para2.1)
  tkgrid.configure(frame.para2.1,sticky="w")
  
  tkgrid(frame.para3)
  tkgrid.configure(frame.para3,sticky="w")
  tkgrid(frame.para4)
  tkgrid.configure(frame.para4,sticky="w")
  tkgrid(frame.para4.1)
  tkgrid.configure(frame.para4.1,sticky="w")
  tkgrid(frame.para4.1.1)
  tkgrid.configure(frame.para4.1.1,sticky="w")
  tkgrid(frame.para4.1.2)
  tkgrid.configure(frame.para4.1.2,sticky="w")
  #tkgrid(frame.para4.1.3)
  #tkgrid.configure(frame.para4.1.3,sticky="w")
  
  #tkgrid(frame.para4.3)
  #tkgrid.configure(frame.para4.3,sticky="w")
  
  tkgrid(frame.para4.2)
  tkgrid.configure(frame.para4.2,sticky="w")
  
  tkgrid(frame.para8)
  tkgrid.configure(frame.para8,sticky="w")
  tkgrid(frame.para5)
  tkgrid.configure(frame.para5,sticky="w")
  tkgrid(frame.para6)
  tkgrid.configure(frame.para6,sticky="w")
  tkgrid(frame.para7)
  tkgrid.configure(frame.para7,sticky="w")
  tkgrid(tklabel(frame,background=color,text=""))
  tkgrid(Run.button)
  tkgrid(frame)
  
  #tk2notetab.select(nb, "Main function")# default "Main functions"
  #tk2notetab.text(nb) # Text of the currently selected tab
  # no tab 20121030
  tkgrid(tn)
  tkgrid(tklabel(LOH,background=color,text=""))
  
  # tkwait.window(LOH)









