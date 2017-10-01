

#Chr.num <- NULL
#if(Cluster){

#inpath <- paste(inpath, "/", sep="")
#outpath <- paste(outpath, "/", sep="")

##### Index setting ######

####################################
####################################
#####   package "locfit"       #####
####################################

#source(paste(strsplit(debug.path,"\\debug110307.r",fixed = TRUE), "\\LOHfunction110307.r", sep=""))
#source(paste(LOH.file, "/LOHfunction", tempdate, ".r", sep = ""))
source(paste(LOH.file, "/LOHAS_Subfunction.r", sep = "")) ### Change into released name

### cancel the package path 20111227

#inpath <- paste(inpath, "/", sep="")
#outpath <- paste(outpath, "/", sep="")
######################################
#####  Beginning time of LOHAS   #####
######################################
# add the log file 20130121
write("Please wait a moment...", paste(outpath,"/Log.txt",sep=""), append = TRUE)
log.file <- paste("Beginning time for LOHAS: ", date(), sep = "")
print( log.file )
write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
print( paste("LOHAS is running. Please wait a moment for the results.", sep = "") )  # Processing reminder
write(paste("LOHAS is running. Please wait a moment for the results.", sep = ""), paste(outpath,"/Log.txt",sep=""), append = TRUE)

########################################
#####    Reading datasetting.txt   #####
########################################
########################################
#if(Cluster == FALSE){
#  setwd(inpath)
#}
if (!plot.only){ ### if (!plot.only) [1]
# the format of Datasetting.txt are fixed 20130109
# there are "9" rows about parameters setting
# 1. Population_name, 2. SNP_ID, 3. Chromosome, 4. Physical_position, 5. Genotype,
# 6. RS_number, 7. MAF_threshold, 8. Start_field, 9. End_field
Data.column <- read.table(paste(inpath, "/Datasetting.txt", sep = ""), as.is = TRUE, header = FALSE, sep = "\t", fill = TRUE, blank.lines.skip = FALSE, allowEscapes = TRUE)
number.popu <- ncol(Data.column)-1
if(Match) num.group <- number.popu
# stopifnot(number.popu ==  as.numeric(tclvalue(group2.val)))

population <- read.table(paste(inpath, "/Inpath.txt", sep = ""), stringsAsFactors = FALSE, allowEscapes = TRUE, header = FALSE, fill = TRUE)
# change the back slash into forward slash 20130109
for (i in 1:nrow(population)) population[i, 1] <- gsub("\\\\", "/", population[i, ])
inpath <- gsub("\\\\", "/", inpath)
outpath <- gsub("\\\\", "/", outpath)
# compulsorily insert "/" at the end of the path of genotype data to prevent users from forgeting typing "/" 20130109
for (i in 1:ncol(population)) {if (substr(population[i, 1], nchar(population[i, 1]), nchar(population[i, 1])) != "/") population[i, 1] <- paste(population[i, 1], "/", sep = "")}
Path <- lapply(1:number.popu, function(i) dir(population[i, ]))

if(chiptype == 4){
  # customized chip allows SNP information in arbitrary columns if user specified in "Datasetting.txt"
  # if the order of columns are Probe_set, Chromosome, Physical_position, Genotype, respectively, users can leave those column as blanks
  name.popu <- Data.column[1, -1]
  SnpCol.column <- as.numeric(Data.column[2, -1]);SnpCol.column[is.na(SnpCol.column)] <- 1
  ChrCol.column <- as.numeric(Data.column[3, -1]);ChrCol.column[is.na(ChrCol.column)] <- 2
  PosiCol.column <- as.numeric(Data.column[4, -1]);PosiCol.column[is.na(PosiCol.column)] <- 3
  CallCol.column <- as.numeric(Data.column[5, -1]);CallCol.column[is.na(CallCol.column)] <- 4
  RsCol.column <- as.numeric(Data.column[6, -1])
  Start_field <- as.numeric(Data.column[9, -1])
  End_field <- as.numeric(Data.column[10, -1])
  sequencing <- TRUE # for Coding()
  # examine whether there is MAF threshold, insert 0 if NOT
  if (!all(is.na(as.numeric(as.character(Data.column[7, -1]))))){
    # MAF threshold for sequencing data
    Rare_threshold <- as.numeric(Data.column[7, -1])
	MAF_weight <- toupper(as.character(Data.column[8, -1]))
    Rare_threshold[is.na(Rare_threshold)] <- 0
	MAF_weight[is.na(MAF_weight)] <- "N"
    #if (Rare_threshold[1] > 0 & Rare_threshold[1] <= 1){  # If threshold == 0, use original coding method
    if (any((Rare_threshold > 1) || (Rare_threshold < 0))){
      #MemoryEnough <- FALSE
      msg <- "Dear LOHAS users, the threshold of minor allele frequency should be between 0 and 1. Please reset the threshold in Datasetting.txt"
      if (length(args) <= 0){ 
        tkmessageBox(title=" Informaion ", message = msg, icon = "info", type = "ok")
        tkfocus(tn)
      }
      write(msg, paste(outpath,"/Log.txt",sep=""), append = TRUE)
    }else{}
  }else{
    #msg <- paste("The parameter setting has changed in LOHAS v2.0.")
    #msg <- paste(msg, paste("Please add the MAF threshold for rare-variants."))
    #msg <- paste(msg, paste("LOHAS v2.0 will set the threshold to 0 if you didn't provide it."))
    msg <- "Dear LOHAS users, a threshold of minor allele frequency for defining rare variants is not found in Datasetting.txt. A threshold of 0 is assumed for your analysis."
    print(msg)
    if (length(args) <= 0){
      tkmessageBox(title=" Informaion ", message = msg, icon = "info", type = "ok")
      tkfocus(tn)
    }
    write(msg, paste(outpath,"/Log.txt",sep=""), append = TRUE)
    Rare_threshold <- rep(0, (ncol(Data.column)-1))
	MAF_weight <- rep("N", (ncol(Data.column)-1))
  }
  ### read one sample from each group
  for(i in 1:number.popu){
    assign(paste("data_", i, sep=""), Read.User(i, 1, 1:23)) # chr23 are not analyzed 20130109
    #assign(paste("data_", i, sep=""), Read.User(i, 1, 1:22))
  }
  IDnumber <- sapply(1:number.popu, function(i) length(Path[[i]]))
  Label <- lapply(1:number.popu, function(i)
    # the strsplit() is used to avoid length different in sample ID
    # example: NA0001.txt and NA00001.txt
    # we should specify the End_field = 7, then the first sample ID would be "NA0001."
    # however, we don't want "." in sample ID
    # #####NOTE#####: hence, the
    sapply(1:IDnumber[i] , function(j) strsplit(substr(Path[[i]][j], Start_field[i], End_field[i]), ".", fixed = TRUE)[[1]][1])
  )
  Chr.num <- sort(unique(data_1[, 3]))
}

if(chiptype == 3){
  name.popu <- Data.column[1, -1]
  SnpCol.column <- as.numeric(Data.column[2, -1])
  ChrCol.column <- as.numeric(Data.column[3, -1]) # chromosome number is needed after combine into one file for each individual 20130109
  PosiCol.column <- as.numeric(Data.column[4, -1])
  CallCol.column <- as.numeric(Data.column[5, -1])
  RsCol.column <- as.numeric(Data.column[6, -1])
  Rare_threshold <- as.numeric(Data.column[7, -1])
  MAF_weight <- toupper(as.character(Data.column[8, -1]))
  Start_field <- as.numeric(Data.column[9, -1])
  End_field <- as.numeric(Data.column[10, -1])

  ### read one sample from each group
  #chr.order <- c(2,13,17,18,19,20,21,22,23,3,4,5,6,7,8,9,10,11,12,14,15,16,25,26)
  ### only  chromosome 1-22, X required 20111229
  #chr.order <- c(1, 12, 16:22, 2:11, 13:15, 23)
  ### chromosomes are combined to a single file for each individual 20130109
  #AFFY6.SNP <- read.table(paste(LOH.file, "/AFFY6_information.txt",sep = ""), as.is = T, header = T)
  IDnumber <- sapply(1:number.popu, function(i) length(Path[[i]]))
  Label <- lapply(1:number.popu, function(i)
    unique(sapply(1:IDnumber[i] , function(a) strsplit(substr(Path[[i]][a], Start_field[i], End_field[i]), ".", fixed = TRUE)[[1]][1]))
  )
  IDnumber <- sapply(1:number.popu, function(i) length(Label[[i]]))
  # chromosomes are combined to a single file for each individual 20130109
  #path <- lapply(1:number.popu, function(i){
  #  lapply(1:IDnumber[[i]], function(j)
  #    Path[[i]][which(substr(Path[[i]], Start_field[i], End_field[i])==Label[[i]][j])][chr.order])
  #})
  for(i in 1:number.popu){
    assign(paste("data_", i, sep=""), Read.Affy(i, 1, Chr.num))
  }
}

if(chiptype %in% 1:2){
  name.popu <- Data.column[1, -1]
  SnpCol.column <- as.numeric(Data.column[2, -1])
  ChrCol.column <- as.numeric(Data.column[3, -1])
  PosiCol.column <- as.numeric(Data.column[4, -1])
  CallCol.column <- as.numeric(Data.column[5, -1])
  RsCol.column <- as.numeric(Data.column[6, -1])
  Rare_threshold <- as.numeric(Data.column[7, -1])
  MAF_weight <- toupper(as.character(Data.column[8, -1]))
  Start_field <- as.numeric(Data.column[9, -1])
  End_field <- as.numeric(Data.column[10, -1])
  if(chipMer) Chip <- c(1, 2);
  if(!chipMer) Chip <- na.omit(c(ifelse(chip1, 1, NA), ifelse(chip2, 2, NA)))
  ### read one sample from each group  
  IDnumber <- sapply(1:number.popu, function(i) length(Path[[i]]))
  Label <- lapply(1:number.popu, function(i)
    unique(sapply(1:IDnumber[i] , function(a) strsplit(substr(Path[[i]][a], Start_field[i], End_field[i]), ".", fixed = TRUE)[[1]][1]))
  )
  IDnumber <- sapply(1:number.popu, function(i) length(Label[[i]]))
  path <- lapply(1:number.popu, function(i){
    lapply(1:IDnumber[[i]], function(j)
      Path[[i]][which(substr(Path[[i]], Start_field[i], End_field[i])==Label[[i]][j])])
  })
  for(i in 1:number.popu){
    assign(paste("data_", i, sep = ""), Read.500K(i, 1, Chr.num, TRUE, Chip))
  }
}

#########
interSNP <- NULL
for(x in Chr.num){
  SNP_1 <- as.character(data_1[which(data_1[, 3]==x), 1]) # if no as.character(), interSNP would be error when intersect. modify 20130128
  if(number.popu > 1){
    for(j in 2:number.popu){
      SNP_i <- as.character(get(paste("data_", j, sep=""))[which(get(paste("data_", j, sep=""))[, 3]==x), 1])
      SNP_1 <- intersect(SNP_1, SNP_i)
    }
  }
  interSNP <- c(interSNP, SNP_1)
}
}else{
  load(paste(inpath, "/Parameter.RData", sep = ""))
} ### END of if (!plot.only) [1]

####################################
#####    Data description     ######
####################################
####################################
### 20130117 modify indent in data description
### Data description folder is longer used, only provides "Data description.txt" under output folder
PlotAllinOne <- ifelse(length(Chr.num) == 1, FALSE, PlotAllinOne)
##PopGroup <- ifelse(num.group==1, ifelse(ethnicity==1, "","population "), "group ")
# popgroup <- ifelse(Match, "phase ", ifelse(num.group==1, "population ", "group "))
# PopGroup <- ifelse(Match, "Phase ", ifelse(num.group==1, "Population ", "Group "))
popgroup <- ifelse(num.group==1, "population ", "group ")
PopGroup <- ifelse(num.group==1, "Population ", "Group ")

if (!plot.only){ ### if (!plot.only) [2]

if(description){

  #dir.create(paste(outpath,"/Data description",sep=""), showWarnings = F)
  #DataDesc.dir <- paste(outpath,"/Data description",sep="")
  
  line1 <- paste("======================================")
  line2 <- paste("===        Data description        ===")
  line3 <- paste("======================================")
  line3.1 <- paste("Input/output path ---")
  line.datades <- paste("-    Data description: Yes")
  input.text <- paste("-    Input directory name: ", inpath, sep="")
  output.text <- paste("-    Output directory name: ", outpath, sep="")
  line.grouptitle <- paste("Study group ---", sep="")
  line4 <- paste("-    The number of study groups: ", num.group,
    ifelse(num.group==2, paste("; Patient group: ", ifelse(patient==1, "1st", "2nd"), " group", sep=""), ""), sep="")
  
  line5 <- paste("Data format ---",sep="")
  line5.text <- switch(chiptype,
    "1" = "Affy 100K",
    "2" = "Affy 500K",
    "3" = "Affy 6.0",
    "4" = "Customized SNP panel")
  line5.text2 <- NULL ## for Affy Chip one or two
  if(chiptype %in% 1:2) {
    PasteChip <- na.omit(c(ifelse(chip1, "Chip 1", NA), ifelse(chip2, "Chip 2", NA), ifelse(chipMer, "Merge 2 chips", NA)))
    line5.text2 <- paste("---", paste(PasteChip, collapse=" and "), sep="")
  }
  line.genechip <- paste("-    Gene chip: ", line5.text, line5.text2, sep="")
  line6 <- line7 <- line8 <- line9 <- line10 <- line11 <- line12 <- NULL 
  line.stattitle <- line.SNPthinning <- line.SNPthinning2 <- line.window <- line.chrranking <- line.assotest<- NULL
  line.space <- paste("         ",sep="")
  ###
  SNPnumber <- sapply(1:number.popu, function(i) nrow(get(paste("data_", i, sep=""))))
  
  ###
  interSNPnumber <- length(interSNP)
  ### sample number for each population
  if(!Match){
  for(i in 1:number.popu){
    if(number.popu > 1){
      line6 <- rbind(line6, paste("    o    The number of total SNPs in ", popgroup, i,": ", SNPnumber[i], sep=""))
      line8 <- rbind(line8, paste("-    The number of individuals in ", popgroup, i, ": ", IDnumber[i], sep=""))
    }
    if(number.popu == 1){
      line6 <- rbind(line6, paste("    o    The number of total SNPs: ", SNPnumber[i], sep=""))
      line8 <- rbind(line8, paste("-    The number of individuals: ", IDnumber[i], sep=""))
    }                                     
  }
  if(!(num.group == 1 & number.popu == 1))
    line7 <- paste("    o    The number of identical SNPs in ",
        ifelse(num.group == 1, "all populations: ", "both groups: "), interSNPnumber, sep="")
  }
  if(Match){
     line8 <- rbind(line8, paste("-    The number of individuals: ", IDnumber[i], sep=""))
     line6 <- rbind(line6, paste("    o    The number of total SNPs: ", SNPnumber[i], sep=""))
  }
  
  
  ###
  if(doLOH){
    line.loh <- paste("-    Estimate of LOH intensity: Yes")
    SNPnumber <- table(data_1[match(interSNP, data_1[, 1]), 3])[paste(Chr.num)]
    line11 <- paste("    o    The number of total SNPs after removing SNPs without chromosome information: ", sum(SNPnumber, na.rm = T), sep="")
    #line12 <- paste("The chromosome(s) studied in LOH analysis:",  ifelse(all(1:23 %in% Chr.num), "All Chromosomes", paste(Chr.num, collapse=", ")), sep=" ") # chr23 are not analyzed 20130109
    line12 <- paste("-    The chromosome(s) studied in LOH analysis:",  ifelse(all(1:22 %in% Chr.num), "All Chromosomes (1-22)", paste(Chr.num, collapse=", ")), sep=" ")
    ###
    if(SNPthin){
      ### maybe wrong!!  should computer this for each columns  
      SNPnumber.thin <- sum(round(SNPnumber/as.numeric(SNPthin.v), 0))
      line.SNPthinning <- paste("    o    The number of total SNPs used for LOH analysis with SNP thinning: ", SNPnumber.thin, sep="")
      line.SNPthinning2 <- paste("    o    SNP thinning: Yes---Select one every ", SNPthin.v," SNPs",sep="")
    }  
    if(1-SNPthin) {
      line.SNPthinning <- paste("    o    The number of total SNPs used for LOH analysis without SNP thinning: ", sum(SNPnumber), sep="")
      line.SNPthinning2 <- paste("    o    SNP thinning: No")
    }
    ###
    line.window <- paste("    o    Window size: Fixed SNP ", ifelse(SNPselectype == 1, "number", "proportion"),": ", propsize, sep="")
    line.stattitle <- paste("Statistical analysis ---", ifelse(RV, "Rare Variants", ""))
  }else{
    line.loh <- paste("-    Estimate of LOH intensity: No")
  }
  
  #if(Match | num.group==2){
    line.chrranking <- paste("-    Chromosome ranking test: ", ifelse(doRanktest, "Yes", "No"), sep="")
    #line.assotest <- paste("Association test: ", ifelse(doAssotest, "Yes", "No"), sep="")
  #}
  line.assotest <- paste("-    Association test: ", ifelse((doAssotest | doFtest | doGEE), "Yes", "No"), sep="")
  if (doAssotest | doFtest | doGEE){
    line.assotest <- c(line.assotest, paste("    o    Wilcoxon test: ", ifelse(doAssotest, "Yes", "No"), sep = ""))
    line.assotest <- c(line.assotest, paste("    o    Regression test: ", ifelse(doFtest, "Yes", "No"), sep = ""))
    line.assotest <- c(line.assotest, paste("    o    GEE test: ", ifelse(doGEE, "Yes", "No"), sep = ""))
    if (any(doFtest, doGEE)){
      line.assotest <- c(line.assotest, paste("    o    Covariate file path: ", covpath, sep = ""))
      # Show the variable names of Traits and Adjust covariates
      # Change into TAB delimited 20130109
      covpath.name <- read.table(covpath, header = TRUE, nrows = 1, sep = "\t");covpath.name <- colnames(covpath.name)[-(1:6)]
      if (!cova.select){
        line.assotest <- c(line.assotest, paste("          +    Trait", ifelse(length(which(substr(covpath.name, 1, 1) %in% c("Q", "B"))) == 1, "", "s"), ": ", paste(covpath.name[which(substr(covpath.name, 1, 1) %in% c("Q", "B"))], collapse = ", "), sep = ""))
        line.assotest <- c(line.assotest, paste("          +    Adjust covariate", ifelse(length(which(substr(covpath.name, 1, 1) %in% c("C", "D"))) == 1, "", "s"), ": ", paste(covpath.name[which(substr(covpath.name, 1, 1) %in% c("C", "D"))], collapse = ", "), sep = ""))
      }else{
        line.assotest <- c(line.assotest, paste("          +    Trait", ifelse(length(tr.id) == 1, "", "s"), ": ", paste(covpath.name[tr.id-6], collapse = ", "), sep = ""))
        line.assotest <- c(line.assotest, paste("          +    Adjust covariate", ifelse(length(cov.id) == 1, "", "s"), ": ", paste(covpath.name[cov.id-6], collapse = ", "), sep = ""))
      }
      if ("covpath.name" %in% ls()) rm(covpath.name)
    }else{}
  }
  line.plot1 <- line.plot2 <- line.plot3 <- NULL
  if (plotselect1 | plotselect2 | plotselect3){
    line.datavis <- paste("-    Data visualization: Yes")
    line.plot1 <- paste("    o    Individual LOH plot: ", ifelse(plotselect1, "Yes", "No"), sep = "")
    line.plot2 <- paste("    o    Combined LOH plot: ", ifelse(plotselect2, paste("Yes, quantile of reference LOH intensity: ", QLOHselect, sep = ""), "No"), sep = "")
    line.plot3 <- paste("    o    LOH biplot: ", ifelse(plotselect3, paste("Yes, alpha scaling: ", Alphaselect, sep = ""), "No"), sep = "")
  }else{
    line.datavis <- paste("-    Data visualization: No")
  }
  
  name <- c(line1, line2, line3, line3.1, input.text, output.text, line.space, line.grouptitle,
    line4, line8, line9, line.space, line5, line.genechip, line12, line6, line7, line10, line11,
    line.SNPthinning, line.space, line.stattitle, line.datades, line.loh, line.SNPthinning2, line.window, line.chrranking, line.assotest,
    line.space, line.datavis, line.plot1, line.plot2, line.plot3)
  description.name <- name
  #write.table(name, paste(DataDesc.dir,"/output.txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(name, paste(outpath,"/Data description.txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
  rm(name, line1, line2, line3, line3.1, input.text, output.text, line.space, line.grouptitle,
    line4, line8, line9, line5, line.genechip, line6, line7, line10, line11,
    line.SNPthinning, line12, line.stattitle, line.datades, line.loh, line.SNPthinning2, line.window, line.chrranking, line.assotest,
    line.datavis, line.plot1, line.plot2, line.plot3)
  
  }  ##description
  
  rm(list = paste("data_", 1:number.popu, sep = ""))
  gc(reset = TRUE)

}else{
  line1 <- paste("======================================")
  line2 <- paste("===   Original parameter setting   ===")
  line3 <- paste("======================================")
  line3.1 <- paste("Input/output path ---")
  input.text <- paste("-    Input directory name: ", inpath, sep="")
  output.text <- paste("-    Output directory name: ", outpath, sep="")
  line.space <- paste("         ",sep="")
  line.stattitle <- paste("Statistical analysis ---", ifelse(RV, "Rare Variants", ""))
  line.assotest <- paste("-    Association test: ", ifelse((doAssotest | doFtest | doGEE), "Yes", "No"), sep="")
  if (doAssotest | doFtest | doGEE){
    line.assotest <- c(line.assotest, paste("    o    Wilcoxon test: ", ifelse(doAssotest, "Yes", "No"), sep = ""))
    line.assotest <- c(line.assotest, paste("    o    Regression test: ", ifelse(doFtest, "Yes", "No"), sep = ""))
    line.assotest <- c(line.assotest, paste("    o    GEE test: ", ifelse(doGEE, "Yes", "No"), sep = ""))
    if (any(doFtest, doGEE)){
      line.assotest <- c(line.assotest, paste("    o    Covariate file path: ", covpath, sep = ""))
      # Show the variable names of Traits and Adjust covariates
      # Change into TAB delimited 20130109
      if (!cova.select){
        #covpath.name <- read.table(covpath, header = TRUE, nrows = 1, sep = "\t");covpath.name <- colnames(covpath.name)[-(1:6)]
        covpath.name <- strsplit(readLines(covpath, n = 1), "\t", fixed = TRUE)[[1]][-(1:6)]
        line.assotest <- c(line.assotest, paste("          +    Trait", ifelse(length(which(substr(covpath.name, 1, 1) %in% c("Q", "B"))) == 1, "", "s"), ": ", paste(covpath.name[which(substr(covpath.name, 1, 1) %in% c("Q", "B"))], collapse = ", "), sep = ""))
        line.assotest <- c(line.assotest, paste("          +    Adjust covariate", ifelse(length(which(substr(covpath.name, 1, 1) %in% c("C", "D"))) == 1, "", "s"), ": ", paste(covpath.name[which(substr(covpath.name, 1, 1) %in% c("C", "D"))], collapse = ", "), sep = ""))
      }else{
        covpath.name <- strsplit(readLines(covpath, n = 1), "\t", fixed = TRUE)[[1]]
        line.assotest <- c(line.assotest, paste("          +    Trait", ifelse(length(tr.id) == 1, "", "s"), ": ", paste(covpath.name[tr.id], collapse = ", "), sep = ""))
        line.assotest <- c(line.assotest, paste("          +    Adjust covariate", ifelse(length(cov.id) == 1, "", "s"), ": ", paste(covpath.name[cov.id], collapse = ", "), sep = ""))
      }
    }else{}
  }
  line.plot1 <- line.plot2 <- line.plot3 <- NULL
  if (plotselect1 | plotselect2 | plotselect3){
    line.datavis <- paste("-    Data visualization: Yes")
    line.plot1 <- paste("    o    Individual LOH plot: ", ifelse(plotselect1, "Yes", "No"), sep = "")
    line.plot2 <- paste("    o    Combined LOH plot: ", ifelse(plotselect2, paste("Yes, quantile of reference LOH intensity: ", QLOHselect, sep = ""), "No"), sep = "")
    line.plot3 <- paste("    o    LOH biplot: ", ifelse(plotselect3, paste("Yes, alpha scaling: ", Alphaselect, sep = ""), "No"), sep = "")
  }else{
    line.datavis <- paste("-    Data visualization: No")
  }
  temp.description.name <- c(description.name[1:3], line3.1, input.text, output.text, line.space,
    line.stattitle, line.assotest, line.space, line.datavis, line.plot1, line.plot2, line.plot3, 
    line.space, line.space, line.space, line1, line2, line3)
  description.name <- c(temp.description.name, description.name[-(1:3)])
  write.table(description.name, paste(outpath,"/Data description.txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")
  rm(temp.description.name, line1, line2, line3, line3.1, input.text, output.text, 
    line.space, line.stattitle, line.assotest, line.datavis, line.plot1, line.plot2, line.plot3)
} ### END of if (!plot.only) [2]

###########################################
#####    Main LOH analysis program    #####
###########################################
####################################
#####  Estimate LOH intensity  #####
####################################

dir.create(paste(outpath, "/Numeric results/", sep=""), showWarnings = F)
dir.create(pathnumLOH <- paste(outpath, "/Numeric results/", ifelse(num.group == 1, "LCSH", "LOH"), " intensity/", sep = ""), showWarnings = F)

LCSHLOH <- ifelse(RV, "RV", ifelse(num.group == 1, "LCSH", "LOH"))

if (!plot.only){ ### if (!plot.only) [3]
### femaleID
heter.rate <- function(Genotype){
  number.het <- sum(Genotype=="AB")
  number.homo <- sum((Genotype=="AA"|Genotype=="BB"))
  number.het/(number.het+number.homo)
}

if(any(Chr.num == 23)){
  femaleID <- list(NULL)
  for(i in 1:number.popu){
    Which <- which(sapply(1:IDnumber[i], function(j){
      Data <- ReadData(paste(chiptype), i, j, 23)
      # Data <- Data[match(interSNP, Data[,1]),]
      Data <- Data[which(Data[, "Probe"] %in% interSNP), ]
      Data <- Data[order(Data[, "pos"]), "call"]
      (heter.rate(Data) >= 0.1)
    }))
    femaleID[[i]] <- Label[[i]][Which]
    #femaleID[[i]] <- 1:IDnumber[i]
    rm(Which)
  }
}




#########################################################################
##### LOH and RV intensity estimation (coding by sapply)  2010-1230 #####
#########################################################################
### 2010-12-30
### for 100K/500K, Affy6.0, User

print("Intensity estimation ...")
write("Intensity estimation ...", paste(outpath,"/Log.txt",sep=""), append = TRUE)
dir.create(paste(outpath, "/Temp", sep = ""), showWarnings = FALSE) # create temp folder for partition genotype files 130611
if(chiptype %in% 1:2){
  if(chipMer)
  Chip12 <- list(1:2, ifelse(chip1, 1, NA), ifelse(chip2, 2, NA))
  if(!chipMer)
  Chip12 <- list(1:2, ifelse(chip1, 1, NA), ifelse(chip2, 2, NA))
  Chip12 <- lapply(which(!is.na(Chip12)), function(x) Chip12[[x]])
}


if(chiptype %in% 3:4){
  Chip12 <- 1
}

for(aa in 1:length(Chip12)){
  if(chiptype %in% 3:4) FILENAME <- NULL
  if(chiptype %in% 1:2) FILENAME <- paste("_Chip", paste(Chip12[[aa]], collapse=""), sep="")
##
for(chr.num in Chr.num){
  chr <- sprintf("%02.f", chr.num)
  eval(parse(text = paste("mafinfo.chip", aa, ".chr", chr, "<- list()", sep = ""))) # save the minor allel frequency of each population 20130418
  for(i in 1:number.popu){
    ptm <- proc.time()
    if (chiptype %in% 4){ # Customized chip: partition into 100000 SNPs a file for each population and save in outpath/Temp 130611
      for (j in 1:IDnumber[i]){
        DATA <- ReadData(paste(chiptype), i, j, chr.num)
        DATA <- DATA[which(DATA[, "Probe"] %in% interSNP), ]
        DATA <- DATA[, c("Probe", "RS_number", "pos", "call")]
        seppart <- unique(c(seq(0, nrow(DATA), by = SN), nrow(DATA)))
        if (j == 1) write.table(DATA[, 1:3], paste(outpath, "/Temp/popu", i, ".chr", chr, ".snpinfo", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
        for (sepi in 1:(length(seppart)-1)){
          write.table(matrix(DATA[(seppart[sepi]+1):seppart[sepi+1], "call"], nrow = 1), paste(outpath, "/Temp/popu", i, ".chr", chr, ".part", sepi, sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = ifelse(j == 1, FALSE, TRUE))
        }
      }
    rm(DATA);gc()
    }
    #if(chiptype %in% 3:4) Data <- as.matrix(CombRawData(i, chr.num))
    if(chiptype %in% 3) Data <- as.matrix(CombRawData(i, chr.num))
    if(chiptype %in% 1:2) Data <- as.matrix(CombRawData(i, chr.num, chip=Chip12[[aa]]))
    if (chiptype %in% 1:3){
      #Data[, -(1:3)] <- t(sapply(1:nrow(Data), function(k) Coding(Data[k, -(1:3)], Type=ifelse(RV, "RV", "LOH"), Sequencing = sequencing, Rare.threshold = ifelse(sequencing, Rare_threshold[i], NULL))))
      for (snpi in 1:nrow(Data)){
        tempdata <- Coding(Data[snpi, -(1:3)], Type=ifelse(RV, "RV", "LOH"), Sequencing = sequencing, Rare.threshold = Rare_threshold[i], MAF.weight = MAF_weight)
        Data[snpi, -(1:3)] <- tempdata$coding
        if (snpi == 1){
          eval(parse(text = paste("mafinfo.chip", aa, ".chr", chr, "[[", i, "]] <- ", tempdata$maf, sep = "")))
        }else{
          eval(parse(text = paste("mafinfo.chip", aa, ".chr", chr, "[[", i, "]] <- c(mafinfo.chip", aa, ".chr", chr, "[[", i, "]], ", tempdata$maf, ")", sep = "")))
        }
        rm(tempdata)
      }
      save(Data, file = paste(outpath, "/Temp/popu", i, ".chip", aa, ".chr", chr, ".coding.RData", sep = ""))
    }
    if (chiptype %in% 4){
      for (sepi in 1:(length(seppart)-1)){
        Data <- read.table(paste(outpath, "/Temp/popu", i, ".chr", chr, ".part", sepi, sep = ""), header = FALSE, sep = "\t", as.is = TRUE)
        Data <- t(Data)
		if(MAF_weight[i] == "Y") Weight <- matrix(1, nrow = nrow(Data), ncol = ncol(Data))
        for (snpi in 1:nrow(Data)){
           tempdata <- Coding(Data[snpi, ], Type=ifelse(RV, "RV", "LOH"), Sequencing = sequencing, Rare.threshold = ifelse(sequencing, Rare_threshold[i], NULL), MAF.weight = MAF_weight[i])
           Data[snpi, ] <- tempdata$coding
		   if(MAF_weight[i] == "Y") Weight[snpi,] <- tempdata$weight
           if (snpi == 1 & sepi == 1){
             eval(parse(text = paste("mafinfo.chip", aa, ".chr", chr, "[[", i, "]] <- ", tempdata$maf, sep = "")))
           }else{
             eval(parse(text = paste("mafinfo.chip", aa, ".chr", chr, "[[", i, "]] <- c(mafinfo.chip", aa, ".chr", chr, "[[", i, "]], ", tempdata$maf, ")", sep = "")))
           }
        }
        for (indi in 1:ncol(Data)){
			if(MAF_weight[i] == "Y"){
				write.table(cbind(Data[, indi], Weight[,indi]), paste(outpath, "/Temp/", paste(i, Label[[i]][indi], sep = "_"), ".chr", chr, ".txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, append = ifelse(sepi == 1, FALSE, TRUE))
			}else{
				write.table(cbind(Data[, indi]), paste(outpath, "/Temp/", paste(i, Label[[i]][indi], sep = "_"), ".chr", chr, ".txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE, append = ifelse(sepi == 1, FALSE, TRUE))
			}
		}
      }
    }
    
    if (chiptype %in% 1:3) load(paste(outpath, "/Temp/popu", i, ".chip", aa, ".chr", chr, ".coding.RData", sep = ""))
    ###
    if((chiptype %in% 1:3) & SNPthin){
      Data <- Data[seq(1, nrow(Data), by=SNPthin.v), ]
    }
    ###
    if(chr.num == Chr.num[1]) dir.create(paste(pathnumLOH, PopGroup, i, sep = ""), showWarnings = F)
    # save .RData of individual LOH intensity in temp, the individual LOH intensity .csv is no longer used
    # only combined LOH intensities are saved in "LOH intensity" folder
    # "LOH plot" folder not used 20130121
    if (!file.exists(paste(pathnumLOH, "temp", sep = ""))) dir.create(paste(pathnumLOH, "temp", sep = ""), showWarnings = F)
    for(j in 1:IDnumber[i]){
      #if(chr.num == Chr.num[1]) dir.create(paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], sep = ""), showWarnings = FALSE)
      if (chiptype %in% 1:3){
        ests.all <- Localfit(Data[, 1], Data[, 2], Data[,3+j], Data[, 3], Size2Prop(SNPselectype, propsize, nrow(Data)))
      }
      if (chiptype %in% 4){
        if (j == 1){
           snpinfo <- read.table(paste(outpath, "/Temp/popu", i, ".chr", chr, ".snpinfo", sep = ""), sep = "\t", as.is = TRUE)
           if (SNPthin) snpinfo <- snpinfo[seq(1, nrow(snpinfo), by=SNPthin.v), ]
        }
        Data <- read.table(paste(outpath, "/Temp/", paste(i, Label[[i]][j], sep = "_"), ".chr", chr, ".txt", sep = ""), header = FALSE, as.is = TRUE)
        if (SNPthin) Data <- Data[seq(1, nrow(Data), by=SNPthin.v), ]
		if(MAF_weight[i] == "Y"){
			ests.all <- Localfit(snpinfo[, 1], snpinfo[, 2], Data[, 1], snpinfo[, 3], Size2Prop(SNPselectype, propsize, nrow(Data)), weights = Data[,2])
		}else{
			ests.all <- Localfit(snpinfo[, 1], snpinfo[, 2], Data[, 1], snpinfo[, 3], Size2Prop(SNPselectype, propsize, nrow(Data)))
		}
      }
      ests.all[ests.all==""] <- "NA"   ###why???
      colnames(ests.all) <- c("Probe_Set", "RS_number", "Phy_Posi_(Mb)", paste(LCSHLOH, "_intensity", sep=""))
      if (!all(is.na(RsCol.column))){
      #if(chiptype %in% 1:3)
        #write.csv(ests.all, file = paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], "/Num_Popu",i,"_", Label[[i]][j], FILENAME, "_Chr", chr,".csv", sep=""), row.names=F, na="")
        partdata <- ests.all
        save(partdata, file = paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(Chip12[[aa]], collapse=""), sep = "")), "_Chr", chr, ".RData", sep = ""))
      }else{
        #write.csv(cbind(1:nrow(ests.all), ests.all[, -2]), file = paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], "/Num_Popu",i,"_", Label[[i]][j], FILENAME, "_Chr", chr,".csv", sep=""), row.names=F, na="")
        partdata <- cbind(1:nrow(ests.all), ests.all[, -2])
        save(partdata, file = paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(paste(Chip12[[aa]], collapse=""), collapse=""), sep = "")), "_Chr", chr, ".RData", sep = ""))
      }
      rm(ests.all, partdata); gc(reset = T)
    }
    ptm <- round(proc.time() - ptm, 1)
    log.file <- paste("Chromosome", chr, PopGroup, i, "finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
    print(log.file)
    write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
  }
  eval(parse(text = paste("save(mafinfo.chip", aa, ".chr", chr, ", file = \'", pathnumLOH, "temp/mafinfo.chip", aa, ".chr", chr, ".RData\')", sep = "")))
} ### Chromosome
##
} ### Chip12

}else{
  dir.create(paste(pathnumLOH, "/temp", sep = ""), recursive = TRUE, showWarnings = FALSE)
    if(chiptype %in% 1:2){
    if(chipMer)
    Chip12 <- list(1:2, ifelse(chip1, 1, NA), ifelse(chip2, 2, NA))
    if(!chipMer)
    Chip12 <- list(1:2, ifelse(chip1, 1, NA), ifelse(chip2, 2, NA))
    Chip12 <- lapply(which(!is.na(Chip12)), function(x) Chip12[[x]])
  }


  if(chiptype %in% 3:4){
    Chip12 <- 1
  }
} ### END of if (!plot.only) [3]

##########################################
##### Chromomsome ranking test #####  ####
##########################################
### only for two groups ### or match sample (more than two groups)
### 2010-08-27  ##

if (!plot.only){ ### if (!plot.only) [4]

if((!Match)& num.group == 2 & doRanktest){
  print("Chromomsome ranking test ...")
  write("Chromomsome ranking test ...", paste(outpath,"/Log.txt",sep=""), append = TRUE)
  ptm <- proc.time()
  dir.create(pathchrtest <- paste(outpath, "/Numeric results/", "Chromosome ranking test/",sep=""), showWarnings = F)
  for(aa in 1:length(Chip12)){
    if(chiptype %in% 3:4) FILENAME <- NULL
    if(chiptype %in% 1:2) FILENAME <- paste("_Chip", paste(Chip12[[aa]], collapse=""), sep="")
    Pvalue.chrTEST <- sapply(Chr.num, function(chr.num) round(RankTest(chr.num, interSNP, chip12=Chip12[[aa]]), 4))
    test.chr <- cbind(sprintf("%-15s", Chr.num), sprintf("%.4f", Pvalue.chrTEST))
    colnames(test.chr) <- c("Chromosome", "p-value")
    write.csv(test.chr, paste(pathchrtest, "Chromosome ranking test", FILENAME, ".csv",sep=""), row.names=F)
    rm(test.chr)
  }
  ptm <- round(proc.time() - ptm, 1)
  log.file <- paste("Chromosome ranking test finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
  print(log.file)
  write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
}

if(Match & doRanktest & num.group==2){
  ### the same as before
}

}else{} ### if (!plot.only) [4]

##############################
#####   Graphic result   #####
##############################

if(any(c(plotselect1, plotselect2, plotselect3))){
  dir.create(paste(outpath,"/Graphic results/",sep=""), showWarnings = F)
  #dir.create(pathnumIP <- paste(outpath,"/Numeric results/", LCSHLOH," plot/",sep=""), showWarnings = F)
}
ptm0 <- proc.time()
#################################################
#####   Individual and Combined LOH plots  #####
#################################################
#### 2010-08-30 ####
#####################
CPcolor <- c("cornflowerblue", "tomato3", "red2", "red4", "gray0")
color.biplot <- colors()[c(26,552,84,502,639,498,79,31,47,52,68,75,81,90,107,108,116,128,142,367,372,382,399,404,410,424,429,435,436,450,455,460,467,475,484,490,493,42,503,508,514,519,524,530,535,536,541,547,557,562,567,568,573,574,584,589,594,599,610,615,620,630,640)]
if(Match){
  color.biplot <- CPcolor
}
Matplot <- function(x, y, Legend=TRUE, Combined=TRUE, Match=FALSE, xlab="Physical position (Mb)", ylab=paste(LCSHLOH, "intensity"), ...){
### 101023
  if(Combined & Legend) par(fig=c(0.15, 0.9, 0, 1))
  matplot(x, y, ylim=c(0, 1), col=color.biplot, xlab=ifelse(Legend, xlab, ""), ylab=ifelse(Legend, ylab, ""),
    type="p", main=paste("Chromosome ",chr.num,sep=""), pch=20, ...)
  if(Combined & Legend){
    if(!Match & num.group!=2){
      par(mar=c(0,0,0,0), fig=c(0,1,0,1), new=T);
      plot(1:10, type="n");
      legend(1, 2.5, c(paste(LCSHLOH, "intensity"), paste("Threshold=", QLOHselect, sep="")), 
        col=c("blue", "red"), pch=16, cex=1.2, title=paste(ifelse(Combined, "Combined", "Individual"), LCSHLOH, "plot"))
    }
    if(!Match & num.group==2){
      par(mar=c(0,0,0,0), fig=c(0,1,0,1), new=T);
      plot(1:10, type="n");
      legend(1, 2.5, c(paste(LCSHLOH, "intensity"), paste(QLOHselect*100, "% QI", sep="")), 
        col=c("blue", "red"), pch=16, cex=1.2, title=paste(ifelse(Combined, "Combined", "Individual"), LCSHLOH, "plot"))   
    }
    if(Match){
      par(mar=c(0,0,0,0), fig=c(0,1,0,1), new=T);
      plot(1:10, type="n");
      legend(1, 2.5, paste("phase", 1:number.popu, sep=""),
        col=color.biplot[1:number.popu], pch=16, cex=1.2, title=paste(ifelse(Combined, "Combined", "Individual"), LCSHLOH, "plot"))         
    }
  }
  if(!Legend){
    mtext(xlab, 1, line=2, cex=0.6)
    mtext(ylab, 2, line=2, cex=0.6)
  }  
}

##############################################################

#if (!plot.only){ ### if (!plot.only) [5]

Type <- c(ifelse(plotselect1, FALSE, NA), ifelse(plotselect2, TRUE, NA))

if (!plot.only){
if(num.group==1) control <- 1:number.popu
if(num.group==2) control <- which(1:number.popu != patient)
if(Match) control <- 1
patient <- (1:number.popu) %w/o% control
}else{}

if(plotselect1)
  dir.create(pathplotI <- paste(outpath,"/Graphic results/Individual plot/",sep=""), showWarnings = F)
if(plotselect2)
  dir.create(pathplotC <- paste(outpath,"/Graphic results/Combined plot/",sep=""), showWarnings = F)
if(!Match)
  for(i in 1:number.popu) dir.create(paste(pathnumLOH, PopGroup, i, sep=""), showWarnings = F)

##############################################################

if(!Match){ ### not Match
if(any(doLOH, !all(is.na(Type)), doAssotest, doFtest, doGEE)){ #### The combined LOH intensity is also used in association test
#dir.create(paste(pathnumLOH, "temp", sep = ""), showWarnings = F) # save temp .RData
for(aa in 1:length(Chip12)){
  if(chiptype %in% 3:4) FILENAME <- NULL
  if(chiptype %in% 1:2) FILENAME <- paste("_Chip", paste(Chip12[[aa]], collapse=""), sep="")
  if (!plot.only){
    n.part <- c() # save the number of partitions of each chromosme
  }else{}

  QLOHlist <- NULL
  for(chr.num in Chr.num){
    proc.time()-> ptm
    chr <- sprintf("%02.f", chr.num)
    #parttemp <- read.csv(paste(pathnumLOH, PopGroup, 1, "/", Label[[1]][1], "/Num_Popu", 1,"_", Label[[1]][1], FILENAME, "_Chr", chr,".csv", sep=""), header = T)
    if (!plot.only){
      load(paste(pathnumLOH, "temp/Num_Popu", 1, "_", Label[[1]][1], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(Chip12[[aa]], collapse=""), sep = "")), "_Chr", chr, ".RData", sep = ""))
    }else{
      load(paste(inpath, "/Num_Popu", 1, "_", Label[[1]][1], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(Chip12[[aa]], collapse=""), sep = "")), "_Chr", chr, ".RData", sep = ""))
    }
    parttemp <- partdata;rm(partdata)
    if (!plot.only){
      n.part <- c(n.part, ceiling(nrow(parttemp)/SN));names(n.part)[length(n.part)] <- chr # number of partition files in chr
    }else{}
    part.skip <- unique(c(seq(0, nrow(parttemp), by = SN), nrow(parttemp)))
    rm(parttemp)
    for(i in c(control, patient)){
      Data <- CombData(i, chr.num, Chip12[[aa]])
      temp.combloh <- Data[, -(1:3)]
      colnames(temp.combloh) <- Label[[i]]
      
      if (!plot.only){
        for (k in 1:(length(part.skip)-1)){
          if (which(c(control, patient) == i) == 1){
            combloh <- temp.combloh[(part.skip[k]+1):part.skip[k+1], ]
            partdatainfo <- Data[(part.skip[k]+1):part.skip[k+1], 1:3]
            colnames(partdatainfo) <- c("Probe_Set", "RS_number", "Phy_Posi_(Mb)")
            save(partdatainfo, combloh, file = paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".RData", sep = ""))
            rm(partdatainfo, combloh)
          }else{
            load(paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".RData", sep = ""))
            combloh <- cbind(combloh, temp.combloh[(part.skip[k]+1):part.skip[k+1], ])
            save(partdatainfo, combloh, file = paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".RData", sep = ""))
            rm(partdatainfo, combloh);gc(reset = T)
          }
        }
      }else{}

      rm(temp.combloh)
      if(plotselect2){
        if(i %in% control){
          QLOH <- QLOHcontrol(num.group, chr.num, i, Data, QLOHselect)
        }
        Data <- cbind(Data, round(QLOH, 4))
        QLOHlist[[chr]] <- cbind(QLOHlist[[chr]], QLOH)
      }
      write.csv(Data, file = paste(pathnumLOH, PopGroup, i, "/Num_Popu",i, FILENAME, "_Chr", chr, ".csv",sep=""), row.names=F,quote=F)
      rm(Data); gc()
    }
    ptm <- round(proc.time() - ptm, 1)
    log.file <- paste("Chromosome", chr, "finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
    print(log.file)
    write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
  }

if(any(Chr.num == "23")){
  QLOHlist[["23"]] <- QLOHlist[["23"]][, c(seq(1, length=number.popu, by=2), seq(2, length=number.popu, by=2))]
}

####



for(Combined in Type){
  for(i in 1:number.popu){ 
    dir.create(paste(ifelse(Combined, pathplotC, pathplotI), PopGroup, i, sep=""), showWarnings = F)
  }
}
#    for(j in 1:IDnumber[i])
#      dir.create(paste(ifelse(Combined, pathplotC, pathplotI), PopGroup, i,"/", Label[[i]][j] ,sep=""))
if(PlotOneinOne){

for(chr.num in Chr.num){
  chr <- sprintf("%02.f", chr.num)
  proc.time() -> ptm
  for(i in c(control, patient)){
    if(chiptype %in% 3:4)
      Data <- read.csv(paste(pathnumLOH, PopGroup, i, "/Num_Popu",i,"_Chr", chr, ".csv",sep=""), header=TRUE, stringsAsFactors=FALSE)
    if(chiptype %in% 1:2)
      Data <- read.csv(paste(pathnumLOH, PopGroup, i, "/Num_Popu",i,"_Chip", paste(Chip12[[aa]], collapse=""), "_Chr", chr, ".csv",sep=""), header=TRUE, stringsAsFactors=FALSE)
    for(j in 1:IDnumber[[i]]){
      ### individual plot for each chromosome
      if(plotselect1){
        if(chr.num==Chr.num[1]) dir.create(paste(pathplotI, PopGroup, i,"/", Label[[i]][j] ,sep=""), showWarnings = F)
        plotname <- paste(PopGroup, i, "/", Label[[i]][j], "/IP_Popu", i, "_", Label[[i]][j], FILENAME, "_Chr", chr, ".jpeg", sep="") 
        jpeg(paste(pathplotI, plotname, sep=""), width = 1024, height = 768)    
        Matplot(Data[, 3], Data[, 3+j], TRUE, FALSE)
        dev.off()
      }
      ### combind plot for each chromosome
      if(plotselect2){
        if(chr.num==Chr.num[1]) dir.create(paste(pathplotC, PopGroup, i,"/", Label[[i]][j] ,sep=""), showWarnings = F)
          plotname <- paste(PopGroup, i, "/", Label[[i]][j], "/CP_Popu", i, "_", Label[[i]][j], FILENAME, "_Chr", chr, ".jpeg", sep="")
        q.column <- ifelse(chr.num==23, ifelse(!(Label[[i]][j] %in% femaleID[[i]]), number.popu+i, i), i)
        jpeg(paste(pathplotC, plotname, sep=""), width = 1024, height = 768)
        Matplot(Data[, 3], cbind(Data[, 3+j], QLOHlist[[chr]][, q.column]), TRUE)             
        dev.off()
      }        
    }
  }
  ptm <- round(proc.time() - ptm, 1)
  log.file <- paste("Chromosome", chr, "finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
  print(log.file)
  write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
}

} 


### All chromosome in one plot ###
if(PlotAllinOne){

for(i in 1:number.popu){
  ptm <- proc.time()
  for(j in 1:IDnumber[i]){
    for(Combined in na.omit(Type)){
      dir.create(paste(ifelse(Combined, pathplotC, pathplotI), PopGroup, i, "/AllChr/", sep=""), showWarnings = F)
      jpeg(paste(ifelse(Combined, pathplotC, pathplotI), PopGroup, i, "/AllChr/IP_Popu", i, "_", Label[[i]][j], FILENAME, "_AllChr.jpeg", sep=""), width = 1024, height = 768);
      ####
      par(mfcol = c(5,5), mar=c(4, 3, 2, 1))
      for(chr.num in Chr.num){
        chr <- sprintf("%02.f", chr.num)     
#          filename <- paste(PopGroup, i, "/", Label[[i]][j], "/Num_Popu", i, "_", Label[[i]][j], "_Chr", chr, ".txt", sep="")
#          Data <- read.table(paste(pathnumLOH, filename, sep=""), header=T)
        #filename <- paste(PopGroup, i, "/", Label[[i]][j], "/Num_Popu", i, "_", Label[[i]][j], FILENAME, "_Chr", chr, ".csv", sep="")
        #Data <- read.csv(paste(pathnumLOH, filename, sep=""))
        filename <- if (!plot.only){
                      paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(Chip12[[aa]], collapse=""), sep = "")), "_Chr", chr, ".RData", sep = "")
                    }else{
                      paste(inpath, "/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(Chip12[[aa]], collapse=""), sep = "")), "_Chr", chr, ".RData", sep = "")
                    }
        load(filename)
        Data <- partdata;rm(partdata)
        q.column <- ifelse(chr.num==23, ifelse(!(Label[[i]][j] %in% femaleID[[i]]), number.popu+i, i), i)
        if(!Combined)
          Matplot(Data[, 3], Data[, 4], FALSE, FALSE)
        if(Combined)
          Matplot(Data[, 3], cbind(Data[, 4], QLOHlist[[chr]][, q.column]), FALSE, TRUE)
        gc(); rm(Data); gc(reset=T)
      }
      if(length(Chr.num)<23){
        for(k in 1:(23-length(Chr.num))){
          plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
        }
      }
      if(!Combined){
        plot(1:10, main="", cex.main=3, type="n", lty=1, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
        legend(0.6, 10 , c("",""), title=paste(ifelse(Combined, "Combined", "Individual"), LCSHLOH, "plot", sep=" "), cex=2, bty="n")
      }
      if(Combined){
        plot(1:10, main=paste("Combined", LCSHLOH, "plot"), cex.main=2, lty=1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
        legend(1.5, 10 , c("LOH intensity", paste(QLOHselect*100,"% QI",sep="")), col=c("blue", "red"),pch=16, cex=2,bty="n")
      }
      ####
      dev.off()
    }
  } 
  ptm <- round(proc.time() - ptm, 1)
  log.file <- paste("PopGroup", i, "finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
  print(log.file)
  write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
}

} ### All chromosome in one plot ###

} ### Chip12

} #### Type
} ### not Match End
#}else{} ### END of if (!plot.only) [5]
#########################################################
#########################################################


if(Match){ ### Match
  Matrix <- matrix(4:(3 + sum(IDnumber)), ncol = number.popu, byrow=TRUE)
  for(chr.num in Chr.num){    
    proc.time() -> ptm
    chr <- sprintf("%02.f", chr.num)
    DATA <- CombData(1, chr.num)
    Data <- as.data.frame(matrix(NA, nrow=nrow(DATA), ncol=3+sum(IDnumber)))
    Data[, 1:ncol(DATA)] <- DATA
    colnames(Data)[1:ncol(DATA)] <- colnames(DATA)
    for(i in 2:number.popu){
      DATA <- CombData(i, chr.num)[, -(1:3)]
      colnames(Data)[sum(!is.na(Data[1, ]))+(1:ncol(DATA))] <- colnames(DATA)
      Data[, sum(!is.na(Data[1, ]))+(1:ncol(DATA))] <- DATA
    }
    Data <- Data[, c(1:3, as.vector(Matrix))]
    write.csv(Data, file = paste(pathnumLOH, "Num_Chr", chr, ".csv",sep=""), row.names=F,quote=F)
    ptm <- round(proc.time() - ptm, 1)
    log.file <- paste("Chromosome", chr, "finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
    print(log.file)
    write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
    ### One chromosome in one plot ###
    Data <- Data[, c(1:3, as.vector(Matrix))]
    if(PlotOneinOne){ ### One Chromosome in one plot
    ### individual plot for each chromosome
    if(plotselect1){
      for(i in 1:number.popu) dir.create(paste(pathplotI, PopGroup, i, sep=""), showWarnings = F)     
      a <- 0
      for(i in 1:number.popu){
        for(j in 1:IDnumber[[i]]){
          if(chr.num==Chr.num[1]) dir.create(paste(pathplotI, PopGroup, i,"/", Label[[i]][j] ,sep=""), showWarnings = F)
          a <- a+1
          plotname <- paste(PopGroup, i, "/", Label[[i]][j], "/IP_Popu", i, "_", Label[[i]][j], "_Chr", chr, ".jpeg", sep="") 
          jpeg(paste(pathplotI, plotname, sep=""), width = 1024, height = 768)    
          Matplot(Data[, 3], Data[, 3+a], TRUE, FALSE)
          dev.off()
        }      
      }
    }
    if(plotselect2){
      ID <- lapply(1:IDnumber[1], function(j) sapply(1:number.popu, function(i) paste(Label[[i]][j], sep="")))
      ID <- sapply(ID, paste, collapse="_")
      for(i in 1:length(ID)) dir.create(paste(pathplotC, ID[i], sep=""), showWarnings = F)
      for(i in 1:IDnumber[1]){
        plotname <- paste(ID[i], "/CP_", ID[i], "_Chr", chr, ".jpeg", sep="") 
        jpeg(paste(pathplotC, plotname, sep=""), width = 1024, height = 768)    
        Matplot(Data[, 3], Data[, Matrix[, i]], TRUE, TRUE, TRUE)
        dev.off()
      }
    }
  } ### One Chromosome in one plot
  rm(Data); gc()
  }

### All chromosome in one plot ###
if(PlotAllinOne){
if(plotselect1){
for(i in 1:number.popu){
  ptm <- proc.time()
  for(j in 1:IDnumber[i]){
    dir.create(paste(pathplotI, PopGroup, i, "/AllChr/", sep=""), showWarnings = F)
    jpeg(paste(pathplotI, PopGroup, i, "/AllChr/IP_Popu", i, "_", Label[[i]][j], FILENAME, "_AllChr.jpeg", sep=""), width = 1024, height = 768);
    ####
    par(mfcol = c(5,5), mar=c(4, 3, 2, 1))
    for(chr.num in Chr.num){
      chr <- sprintf("%02.f", chr.num)     
      #filename <- paste(PopGroup, i, "/", Label[[i]][j], "/Num_Popu", i, "_", Label[[i]][j], FILENAME, "_Chr", chr, ".csv", sep="")
      #Data <- read.csv(paste(pathnumLOH, filename, sep=""))
      load(paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], ifelse(chiptype %in% 3:4, "", paste("_Chip", paste(Chip12[[aa]], collapse=""), sep = "")), "_Chr", chr, ".RData", sep = ""))
      Data <- partdata;rm(partdata)
      Matplot(Data[, 3], Data[, 4], FALSE, FALSE)
      gc(); rm(Data); gc(reset=T)
    }
    if(length(Chr.num)<23){
      for(k in 1:(23-length(Chr.num))){
      plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
      }
    }    
    plot(1:10, main="", cex.main=3, type="n", lty=1, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    legend(0.6, 10 , c("",""), title=paste("Individual", LCSHLOH, "plot", sep=" "), cex=2, bty="n")
    ####
    dev.off()
  } 
  ptm <- round(proc.time() - ptm, 1)
  log.file <- paste("PopGroup", i, "finished: user,", ptm[1], "seconds; elapsed,", ptm[3], "seconds")
  print(log.file)
  write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
}  
} #### Individual plot End

if(plotselect2){ #### Combinded plot
  dir.create(paste(pathplotC, "AllChr/", sep=""), showWarnings = F)
  ID <- lapply(1:IDnumber[1], function(j) sapply(1:number.popu, function(i) paste(Label[[i]][j], sep="")))      
  ID <- sapply(ID, paste, collapse="_")  
  for(i in 1:IDnumber[1]){
    jpeg(paste(pathplotC, "AllChr/", ID[i], FILENAME, "_AllChr.jpeg", sep=""), width = 1024, height = 768)
    ####
    par(mfcol = c(5,5), mar=c(4, 3, 2, 1))
    for(chr.num in Chr.num){
      chr <- sprintf("%02.f", chr.num)
      Data <- read.csv(paste(pathnumLOH, "Num_Chr", chr, ".csv",sep=""))
      Data <- Data[, c(1:3, as.vector(Matrix))]
      plotname <- paste(ID[i], "/CP_", ID[i], "_Chr", chr, ".jpeg", sep="") 
      Matplot(Data[, 3], Data[, Matrix[, i]], FALSE, TRUE, TRUE)
    }
    if(length(Chr.num) < 23){
      for(k in 1:(23-length(Chr.num))){
        plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
      }
    }
    plot(1:10, main="", cex.main=3, type="n", lty=1, xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    legend(0.6, 10 , paste("phase", 1:number.popu, sep=""), col=color.biplot[1:number.popu], title=paste("Combinded", LCSHLOH, "plot", sep=" "), cex=2, bty="n", pch=20)
    ####
    dev.off()       
  }
} #### Combinded plot
} #### All Chromosome in One plot

} ### Match


#####################################################
##### Association test and -log10 pvalue plot   #####
##################################################### 
###########################
#####    Two groups    ####
###########################
##### With CB plot #####

#if (!plot.only){ ### if (!plot.only) [6]
#if(!Match & (doAssotest || doFtest) & num.group == 2){
if(!Match & (doAssotest || doFtest || doGEE)){  # Add GEE test and association test can be used when num.group == 1
print("Association test...")
write("Association test...", paste(outpath,"/Log.txt",sep=""), append = TRUE)
dir.create(pathnumAT <- paste(outpath,"/Numeric results/Association test/",sep=""), showWarnings = F)

for(aa in 1:length(Chip12)){
  if(chiptype %in% 3:4) FILENAME <- NULL
  if(chiptype %in% 1:2) FILENAME <- paste("_Chip", paste(Chip12[[aa]], collapse=""), sep="")
  
  if (0){
  ##### Partition LOH intensity (./Numerical results/temp)#####
  if (doFtest || doGEE){
    n.part <- c()
    if (!file.exists(paste(pathnumLOH, "temp", sep = ""))) dir.create(paste(pathnumLOH, "temp", sep = ""), showWarnings = F)
    for(chr.num in (Chr.num %w/o% 23)){
      chr <- sprintf("%02.f", chr.num)
      parttemp <- read.csv(paste(pathnumLOH, PopGroup, 1, "/", Label[[1]][1], "/Num_Popu",1,"_", Label[[1]][1], FILENAME, "_Chr", chr,".csv", sep=""), header = T)
      n.part <- c(n.part, ceiling(nrow(parttemp)/SN));names(n.part)[length(n.part)] <- chr # number of partition files in chr
      part.skip <- unique(c(seq(0, nrow(parttemp), by = SN), nrow(parttemp)))
      rm(parttemp)
      for (k in 1:(length(part.skip)-1)){
        combloh <- c()
        partdatainfo <- read.csv(paste(pathnumLOH, PopGroup, 1, "/", Label[[1]][1], "/Num_Popu",1,"_", Label[[1]][1], FILENAME, "_Chr", chr,".csv", sep=""), header = F, skip = (1+part.skip[k]), nrows = (part.skip[k+1]-part.skip[k]));partdatainfo <- partdatainfo[, -4]
        for(i in 1:number.popu){
          if (k == 1){
            for(j in 1:IDnumber[i]){
              #partdata <- read.csv(paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], "/Num_Popu",i,"_", Label[[i]][j], FILENAME, "_Chr", chr,".csv", sep=""), header = F, skip = (1+part.skip[k]), nrows = (part.skip[k+1]-part.skip[k]))
              partdata <- read.csv(paste(pathnumLOH, PopGroup, i, "/", Label[[i]][j], "/Num_Popu",i,"_", Label[[i]][j], FILENAME, "_Chr", chr,".csv", sep=""), header = T)
              combloh <- cbind(combloh, partdata[(1+part.skip[k]):part.skip[k+1], 4])
              save(partdata, file = paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], FILENAME, "_Chr", chr, ".RData", sep = ""))
              rm(partdata);gc()
            }
          }else{
            for(j in 1:IDnumber[i]){
              load(paste(pathnumLOH, "temp/Num_Popu", i, "_", Label[[i]][j], FILENAME, "_Chr", chr, ".RData", sep = ""))
              combloh <- cbind(combloh, partdata[(1+part.skip[k]):part.skip[k+1], 4])
              rm(partdata);gc()
            }
          }
        }
        #write.table(partdatainfo, paste(pathnumLOH, "temp/chr", chr.num, ".info.part", k, ".txt", sep = ""), row.names = F, col.names = c("Probe_Set", "RS_number", "Phy_Posi_(Mb)"), quote = F, sep = "\t")
        #write.table(combloh, paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".txt", sep = ""), row.names = F, col.names = unlist(Label), quote = F, sep = "\t")
        colnames(partdatainfo) <- c("Probe_Set", "RS_number", "Phy_Posi_(Mb)")
        colnames(combloh) <- unlist(Label)
        save(partdatainfo, combloh, file = paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".RData", sep = ""))
        rm(partdatainfo, combloh);gc();gc(reset = T)
      }
    }
  }
  }else{} # END of if (0) move to

for(chr.num in (Chr.num %w/o% 23)){
  chr <- sprintf("%02.f", chr.num)
  if(!any(!all(is.na(Type)), doAssotest, doFtest, doGEE)){
    fitcontrol <- CombData(control, chr.num, chip12=ifelse(chiptype %in% 3:4, 3, paste(Chip12[[aa]], collapse="")))
    fitcase <- CombData(patient, chr.num, chip12=ifelse(chiptype %in% 3:4, 3, paste(Chip12[[aa]], collapse="")))
    QLOH <- QLOHcontrol(num.group, chr.num, control, fitcontrol, QLOHselect)
    fitcontrol <- cbind(fitcontrol, QLOH)
    fitcase <- cbind(fitcase, QLOH)
    colnames(fitcase)[4:(dim(fitcase)[2]-1)] <- sapply(4:(dim(fitcase)[2]-1), function(k) substr(colnames(fitcase)[k], 3, nchar(colnames(fitcase)[k]) ) )  # set the individual id equal to that in the covariate data 20120220
    colnames(fitcontrol)[4:(dim(fitcontrol)[2]-1)] <- sapply(4:(dim(fitcontrol)[2]-1), function(k) substr(colnames(fitcontrol)[k], 3, nchar(colnames(fitcontrol)[k]) ) )
  }
  if(any(!all(is.na(Type)), doAssotest, doFtest, doGEE)){
    fitcase <- read.csv(paste(pathnumLOH, PopGroup, patient, "/Num_Popu", patient, FILENAME, "_Chr", chr, ".csv",sep=""), header=TRUE, stringsAsFactors=FALSE)
    fitcontrol <- read.csv(paste(pathnumLOH, PopGroup, control, "/Num_Popu", control, FILENAME, "_Chr", chr, ".csv",sep=""), header=TRUE, stringsAsFactors=FALSE)
    colnames(fitcase)[4:(dim(fitcase)[2]-1)] <- sapply(4:(dim(fitcase)[2]-1), function(k) substr(colnames(fitcase)[k], 4, nchar(colnames(fitcase)[k]) ) )
    colnames(fitcontrol)[4:(dim(fitcontrol)[2]-1)] <- sapply(4:(dim(fitcontrol)[2]-1), function(k) substr(colnames(fitcontrol)[k], 4, nchar(colnames(fitcontrol)[k]) ) )
  }
  snp <- fitcase[, 3]

  ########## Wilcoxon Test ##########################################################
  if (doAssotest){
    m1 <- sapply(1:nrow(fitcase), function(j)
    median(as.numeric(c(fitcase[j, 4:(ncol(fitcase)-1)], fitcontrol[j, 4:(ncol(fitcontrol)-1)]))))
    xx <- rep(0,(ncol(fitcase) -  4));
    yy <- rep(0,(ncol(fitcontrol) - 4));
    r <- nonzero(snp, 0)
    ##### proportion of cases and controls larger than quantile of controls in CB plot #####
    prop.case <- sapply(1:nrow(fitcase), function(k) sum(fitcase[k, 4:(ncol(fitcase)-1)]> fitcase[k,ncol(fitcase)] )/(ncol(fitcase)-4) )
    prop.control <- sapply(1:nrow(fitcontrol), function(k) sum(fitcontrol[k, 4:(ncol(fitcontrol)-1)]> fitcontrol[k, ncol(fitcontrol)] )/(ncol(fitcontrol)-4) )
    ################################################################################
    p_value <- c()
    p_value_columnnames <- c()
    ##### Window size 1 #####
    windowtest.size <- 1
    p_w1 <- c()
    for(k in 1:length(snp)){
      for(j in 1:(ncol(fitcase)-4)){
        xx[j] <- KL(fitcase[k, j+3], m1[k], r[k], "sum")
      }
      for(j in 1:(ncol(fitcontrol)-4)){
        yy[j] <- KL(fitcontrol[k, j+3],m1[k],r[k], "sum")
      }
      p_w.temp <- wilcox.test(yy, xx, alternative="less")
      p_w1 <- c(p_w1, abs(-log10(p_w.temp$p.value)))
    }
    p_value <- cbind(p_value, sprintf("%-9s",sprintf("%.4f",p_w1)))
    p_value_columnnames <- c(p_value_columnnames, "-log10P(#W=1)")
  }


  ############ F test ################################################################
  if (doFtest){
    ### case: 1, control: 0
    covpath <- gsub("\\\\", "/", covpath)
    covariateData <- read.table(covpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Define the format to be similar to linkage format
    # Family_ID, Individual_ID, Father_ID, Mother_ID, Sex, Disease_status
    # the above six column is required either Regression or GEE test
    # And the additional columns are used in model
    # ( Use the first alphabet of the column name to determine the type of variable )
    # "Q" for quantitative trait, eg. QSBP
    # "B" for binary trait, eg. Bdisease
    # "C" for quantitative covariate, eg. Cage
    # "D" for discrete covariate, eg. Dgender
    covariate.name <- as.character(covariateData[, 2]) # 2nd column is Individual_ID
    covariate.column <- colnames(covariateData)
    
    if (!cova.select){
      statepos <- which(substr(covariate.column, 1, 1) %in% c("Q", "B")) # find the covariate be tested in our hypothesis (H0: beta_state = 0 vs H1: beta_state =\= 0)
      modelpos <- 7:ncol(covariateData) # independent variables in the model
      if (length(statepos) == 0){
        stop(paste("There is no quantitative/binary trait."))
      }else{
        modelpos <- modelpos[!(modelpos %in% statepos)]  # the covariates be tested in our hypothesis are removed and put it to the last in aov()
      }
    }else{
      statepos <- tr.id
      modelpos <- cov.id
    }
    
    # the discrete covariates should be changed into dummy variables
    discretepos <- which(substr(covariate.column, 1, 1) %in% c("B", "D")); discretepos <- discretepos[discretepos != 6]

    if (length(discretepos) > 0){
      if (FALSE){
        modelpos <- modelpos[!(modelpos %in% discretepos)]
        for (discreteposi in discretepos){
          unique.discretepos <- sort(unique(covariateData[, discreteposi]))
          # define the smallest value to be the reference state
          for (uniquei in unique.discretepos[-1]){
            dummyvar <- rep(0, nrow(covariateData))
            dummyvar[which(covariateData[, discreteposi] == uniquei)] <- 1
            covariateData <- cbind(covariateData, dummyvar);colnames(covariateData)[ncol(covariateData)] <- paste(colnames(covariateData)[discreteposi], "_", uniquei, sep = "")
            rm(dummyvar)
          }
        }
        # add the dummy variable positions in covariateData to modelpos
        modelpos <- c(modelpos, ((length(covariate.column)+1):ncol(covariateData)))
        covariate.column <- colnames(covariateData) # because the number of column increases
      }else{}
      for (discreteposi in discretepos) covariateData[, discreteposi] <- as.factor(covariateData[, discreteposi])
    }else{}

    if (0){  ##### UNUSED
    # "ID" and "state" are specific colnames for individual id and the covariate be test in our hypothesis (H0: beta_state = 0 vs H1: beta_state =\= 0)
    # the default columns for "ID" and "state" are 1st and 2nd respectively
    idpos <- ifelse(any(colnames(covariateData) == "ID"), which(colnames(covariateData) == "ID")[1], 1)
    covariate.name <- as.character(covariateData[, idpos])
    statepos <- ifelse(any(colnames(covariateData) == "state"), which(colnames(covariateData) == "state")[1], 2)
    #covariateData[, 1] <- sapply(1:dim(covariateData)[1], function(k) ifelse(covariateData[k, 2] == 1, substr(covariateData[k, 1], Start_field[-control], End_field[-control]), substr(covariateData[k, 1], Start_field[control], End_field[control]) ) )
    }else{}  ##### UNUSED END

    #for (chr.num in Chr.num){
      #CombLI <- cbind(as.matrix(fitcase[, 4:(dim(fitcase)[2]-1)]), as.matrix(fitcontrol[, 4:(dim(fitcontrol)[2]-1)])) #  column: sample ID, row: marker
      #CombLI <- t(CombLI)  ### combined LOH intensity column: marker, row: sample ID
      #for (iid in 1:dim(CombLI)[1]){  ### Align the Individual ID between covariateData and CombLI
      #  iid.1 <- which(covariateData[, 1] == rownames(CombLI)[iid])
      #  if (length(iid.1) != 1){
      #    cat("Errors in Individual ID\n")
      #  }else{}
      #  iid.2 <- covariateData[iid, ]
      #  covariateData[iid, ] <- covariateData[iid.1, ]
      #  covariateData[iid.1, ] <- iid.2
      #}
      #rm(list = c("iid.1", "iid.2", "iid"))

      #chr <- sprintf("%02.f", chr.num)
      total.f.pvalue <- c()  # combind partitions in a chromosome
      for (k in 1:n.part[chr]){
        #CombLI <- read.table(paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".txt", sep = ""), header = T, sep = "\t")
        #CombLI.info <- read.table(paste(pathnumLOH, "temp/chr", chr.num, ".info.part", k, ".txt", sep = ""), header = T, sep = "\t")
        if (!plot.only){
          load(paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".RData", sep = ""))
        }else{
          load(paste(inpath, "/chr", chr.num, ".part", k, ".RData", sep = ""))
        }
        CombLI.info <- partdatainfo
        CombLI <- combloh
        rm(partdatainfo, combloh)
        if (k == 1){
          CombLI.name <- colnames(CombLI)
          # align the ID in covariate data and LOH intensity data because their order may be different
          # use pmatch() function to modify covariate data
          # since the length of ID can be different, and an error occurs when aligning "NA001" and "NA0012" for example
          # we should align longer ID first
          covariate.name.length <- nchar(covariate.name)
          alignpos <- rep(NA, length(covariate.name))
          for (align.i in sort(unique(covariate.name.length), decreasing = T)){
            alignpos[which(covariate.name.length == align.i)] <- pmatch(covariate.name[which(covariate.name.length == align.i)], CombLI.name)
          }
          covariateData <- covariateData[if(any(is.na(alignpos))){!is.na(alignpos)}else{T}, ]
          covariate.name <- covariateData[, 2]
        }else{}
        CombLI <- CombLI[, CombLI.name %in% covariate.name]
        CombLI <- CombLI[, alignpos]
        CombLI.name <- colnames(CombLI)
        CombLI <- t(matrix(as.numeric(as.matrix(CombLI)), ncol = length(CombLI.name)))
        tf.pvalue <- c() # p-value matrix
        for (stateposi in statepos){
          #f.aov <- aov(as.formula(paste("CombLI ~ ", paste("as.numeric(covariateData[, ", paste(c((1:ncol(covariateData))[-c(idpos, statepos)], statepos), "])", sep = ""), collapse = "+"), sep = "") ) )  ### fit an analysis of variance model
          f.aov <- aov(as.formula(paste("CombLI ~ ", paste("as.numeric(covariateData[, ", paste(c(modelpos, stateposi), "])", sep = ""), collapse = "+"), sep = "") ) )  ### fit an analysis of variance model
          f.summ <- summary(f.aov)
          f.pvalue <- sapply(1:length(f.summ), function(k) f.summ[[k]][paste("as.numeric(covariateData[, ", stateposi, "])", sep = ""), "Pr(>F)"])
          tf.pvalue <- cbind(tf.pvalue, sprintf("%-9s", sprintf("%.4f", -log10(f.pvalue) )))
          colnames(tf.pvalue)[ncol(tf.pvalue)] <- paste("Regression_", covariate.column[stateposi], "_-log10(pv)", sep = "")
          rm(f.aov, f.summ, f.pvalue);gc(reset = T)
        }
        #write.csv(cbind(CombLI.info, f.pvalue), paste(pathnumAT,"Regression_test", FILENAME, "_Chr",chr,".csv",sep=""), row.names = F, col.names = ifelse(k == 1,  c(colnames(CombLI.info), "-Log10(Pvalue)"), F),  quote  = F, append = ifelse(k == 1, F, T))
        #write.csv(cbind(CombLI.info, tf.pvalue), paste(pathnumAT,"Regression_test", FILENAME, "_Chr",chr,".csv",sep=""), row.names = F, col.names = ifelse(k == 1,  c(colnames(CombLI.info), paste("-log10p.", covariate.column[statepos], sep = "")), F),  quote  = F, append = ifelse(k == 1, F, T))
        #rm(tf.pvalue);gc(reset = T)
      #}  ### END of for (chr.num in Chr.num)
      total.f.pvalue <- rbind(total.f.pvalue, tf.pvalue)
      rm(tf.pvalue);gc(reset = T)
    }
  }

############ GEE ################################################################
  if (doGEE){
    ### case: 1, control: 0
    covpath <- gsub("\\\\", "/", covpath)
    covariateData <- read.table(covpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # make sure that Family_ID are change to numerical form
    covariateData[, 1] <- as.factor(covariateData[, 1])
    covariateData[, 1] <- as.numeric(covariateData[, 1])
    # Define the format to be similar to linkage format
    # Family_ID, Individual_ID, Father_ID, Mother_ID, Sex, Disease_status
    # the above six column is required either Regression or GEE test
    # And the additional columns are used in model
    # ( Use the first alphabet of the column name to determine the type of variable )
    # "Q" for quantitative trait, eg. QSBP
    # "B" for binary trait, eg. Bdisease
    # "C" for quantitative covariate, eg. Cage
    # "D" for discrete covariate, eg. Dgender
    covariate.name <- as.character(covariateData[, 2]) # 2nd column is Individual_ID
    covariate.column <- colnames(covariateData)

    if (!cova.select){
      statepos <- which(substr(covariate.column, 1, 1) %in% c("Q", "B")) # find the covariate be tested in our hypothesis (H0: beta_state = 0 vs H1: beta_state =\= 0)
      modelpos <- 7:ncol(covariateData) # independent variables in the model
      if (length(statepos) == 0){
        stop(paste("There is no quantitative/binary trait."))
      }else{
        modelpos <- modelpos[!(modelpos %in% statepos)]  # the covariates be tested in our hypothesis are removed and put it to the last in aov()
      }
    }else{
      statepos <- tr.id
      modelpos <- cov.id
    }
    
    # the discrete covariates should be changed into dummy variables
    discretepos <- which(substr(covariate.column, 1, 1) %in% c("B", "D")); discretepos <- discretepos[discretepos != 6]

    if (length(discretepos) > 0){
      if (FALSE){
        modelpos <- modelpos[!(modelpos %in% discretepos)]
        for (discreteposi in discretepos){
          unique.discretepos <- sort(unique(covariateData[, discreteposi]))
          # define the smallest value to be the reference state
          for (uniquei in unique.discretepos[-1]){
            dummyvar <- rep(0, nrow(covariateData))
            dummyvar[which(covariateData[, discreteposi] == uniquei)] <- 1
            covariateData <- cbind(covariateData, dummyvar);colnames(covariateData)[ncol(covariateData)] <- paste(colnames(covariateData)[discreteposi], "_", uniquei, sep = "")
            rm(dummyvar)
          }
        }
        # add the dummy variable positions in covariateData to modelpos
        modelpos <- c(modelpos, (length(covariate.column)+1):ncol(covariateData))
        covariate.column <- colnames(covariateData) # because the number of column increases
      }else{}
      ### change discrete variable as factor
      for (discreteposi in discretepos) covariateData[, discreteposi] <- as.factor(covariateData[, discreteposi])
    }else{}

    if (0){ ##### UNUSED
    # "ID" and "state" are specific colnames for individual id and the covariate be test in our hypothesis (H0: beta_state = 0 vs H1: beta_state =\= 0)
    # the default columns for "ID" and "state" are 1st and 2nd respectively
    idpos <- ifelse(any(colnames(covariateData) == "ID"), which(colnames(covariateData) == "ID")[1], 1)
    covariate.name <- as.character(covariateData[, idpos])
    statepos <- ifelse(any(colnames(covariateData) == "state"), which(colnames(covariateData) == "state")[1], 2)
    }else{} ##### UNUSED END

    #for (chr.num in Chr.num){
      #chr <- sprintf("%02.f", chr.num)
      total.g.pvalue <- c()  # combind partitions in a chromosome
      for (k in 1:n.part[chr]){
        #CombLI <- read.table(paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".txt", sep = ""), header = T, sep = "\t")
        #CombLI.info <- read.table(paste(pathnumLOH, "temp/chr", chr.num, ".info.part", k, ".txt", sep = ""), header = T, sep = "\t")
        if (!plot.only){
          load(paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".RData", sep = ""))
        }else{
          load(paste(inpath, "/chr", chr.num, ".part", k, ".RData", sep = ""))
        }
        CombLI.info <- partdatainfo
        CombLI <- combloh
        rm(partdatainfo, combloh)
        if (k == 1){
          CombLI.name <- colnames(CombLI)
          # align the ID in covariate data and LOH intensity data because their order may be different
          # use pmatch() function to modify LOH intensity data because the family number in covariate data should be ordered
          # since the length of ID can be different, and an error occurs when aligning "NA001" and "NA0012" for example
          # we should align longer ID first
          covariate.name.length <- nchar(covariate.name)
          alignpos <- rep(NA, length(covariate.name))
          for (align.i in sort(unique(covariate.name.length), decreasing = T)){
            alignpos[which(covariate.name.length == align.i)] <- pmatch(covariate.name[which(covariate.name.length == align.i)], CombLI.name)
          }
          covariateData <- covariateData[if(any(is.na(alignpos))){!is.na(alignpos)}else{T}, ]
          family.order <- order(covariateData[, 1])
          covariateData <- covariateData[family.order, ]  # the Family_ID should be ordered
          covariate.name <- covariateData[, 2]
          alignpos <- na.omit(alignpos)[family.order]
        }
        CombLI <- CombLI[, CombLI.name %in% covariate.name]
        #CombLI <- CombLI[, family.order]  # the Individual_ID in CombLI aligns to the Individual_ID in covariate data
        CombLI <- CombLI[, alignpos]  # the Individual_ID in CombLI aligns to the Individual_ID in covariate data
        CombLI.name <- colnames(CombLI)
        CombLI <- matrix(as.numeric(as.matrix(CombLI)), ncol = length(CombLI.name))

        tg.pvalue <- c() # p-value matrix
        for (stateposi in statepos){
          tempg.pvalue <- rep(NA, nrow(CombLI))
          for (snpi in 1:nrow(CombLI)){
            geedataframe <- data.frame(as.numeric(CombLI[snpi, ]), as.data.frame(covariateData[, modelpos]), as.numeric(covariateData[, stateposi]), covariateData[, 1]) # dependent, independent, trait variable, and family number
            names(geedataframe) <- c("LOHintensity", covariate.column[c(modelpos, stateposi)], "Family_ID")
            gee.test <- geem(as.formula(paste("LOHintensity ~ ", paste(covariate.column[c(modelpos, stateposi)], collapse = "+"), sep = "")), id = Family_ID, data = geedataframe, family = gaussian, corstr = "exchangeable")
            gee.sum <- summary(gee.test)
            tempg.pvalue[snpi] <- -log10(tail(gee.sum$p, 1))
            rm(geedataframe, gee.test, gee.sum);gc(reset = T)
          }
          tg.pvalue <- cbind(tg.pvalue, sprintf("%.4f", tempg.pvalue))
          colnames(tg.pvalue)[ncol(tg.pvalue)] <- paste("GEE_", covariate.column[statepos[stateposi]], "_-log10(pv)", sep = "")
          rm(tempg.pvalue);
        }
        total.g.pvalue <- rbind(total.g.pvalue, tg.pvalue)
        rm(tg.pvalue);gc(reset = T)
        #write.csv(cbind(CombLI.info, tg.pvalue), paste(pathnumAT,"GEE", FILENAME, "_Chr",chr,".csv",sep=""), row.names = F, col.names = ifelse(k == 1,  c(colnames(CombLI.info), "-Log10(Pvalue)"), F),  quote  = F, append = ifelse(k == 1, F, T))
      }
    #}  ### END of for (chr.num in Chr.num)
  }
  #################################################################################
  #max.len.rs <-max(nchar(as.character(fitcase[,2])))+1
  #fitcase[,2]<-sprintf(paste("%-",max.len.rs,"s",sep=""),fitcase[,2] )
  if (doAssotest){
    all.matrix <- cbind(fitcase[,1:3],
      sapply(1:nrow(fitcase), function(k) mean( unlist( fitcase[k, 4:(ncol(fitcase)-1)]) ) ),
      sapply(1:nrow(fitcontrol), function(k) mean( unlist( fitcontrol[k, 4:(ncol(fitcontrol)-1)]))),
      sapply(1:nrow(fitcase), function(k) median( unlist( fitcase[k, 4:(ncol(fitcase)-1)]) ) ),
      sapply(1:nrow(fitcontrol), function(k) median( unlist( fitcontrol[k, 4:(ncol(fitcontrol)-1)]) ) ),
      sapply(1:nrow(fitcase), function(k) sd( unlist( fitcase[k, 4:(ncol(fitcase)-1)]) ) ),
      sapply(1:nrow(fitcontrol), function(k) sd( unlist( fitcontrol[k,4:(ncol(fitcontrol)-1)]) ) ))
    colnames(all.matrix) <- c("Probe_Set", "RS_number", "Phy_Posi_(Mb)","Mean_case","Mean_cont","Median_case","Median_cont","Std_case","Std_cont")

    all.matrix <- cbind(all.matrix, prop.case, prop.control, p_value)
    colnames(all.matrix)[(dim(all.matrix)[2]-2):dim(all.matrix)[2]] <- c("Prop_case","Prop_cont", p_value_columnnames)
  }else{  #  Association test matrix when single group, i.e., no Wilcoxon test
    all.matrix <- c()
    for (k in 1:n.part[chr]){
      if (!plot.only){
        load(paste(pathnumLOH, "temp/chr", chr.num, ".part", k, ".RData", sep = ""))
      }else{
        load(paste(inpath, "/chr", chr.num, ".part", k, ".RData", sep = ""))
      }
      #all.matrix <- rbind(all.matrix, read.table(paste(pathnumLOH, "temp/chr", chr.num, ".info.part", k, ".txt", sep = ""), header = T, sep = "\t"))
      all.matrix <- rbind(all.matrix, partdatainfo)
      rm(partdatainfo, combloh);gc()
    }
  }
  if (doFtest){
    all.matrix <- cbind(all.matrix, total.f.pvalue)
    #rm(total.f.pvalue);
    gc(reset = T)
  }else{}
  if (doGEE){
    all.matrix <- cbind(all.matrix, total.g.pvalue)
    #rm(total.g.pvalue);
    gc(reset = T)
  }else{}
  write.csv(all.matrix, file=paste(pathnumAT,"Association_test", FILENAME, "_Chr",chr,".csv",sep=""), row.names=F)


  #################################################################################
  ##### P-value Plot #####
  if(PlotOneinOne){
    xlim.temp <- c(min(snp), max(snp))
    color.pvalueplot <- colors()[c(26,552,84,502,639,498,79,31,47,52)]
    pathnumPP <- paste(outpath, "/Graphic results/Association test/", sep = "")
    dir.create(pathnumPP, showWarnings = FALSE)
    if (doAssotest){
      jpeg(paste(pathnumPP, "Wilcoxon_test", FILENAME, "_Chr", chr,".jpeg",sep=""), width = 1024, height = 768)
      #par(fig=c(0.2,0.95,0,0.95), mar=c(5, 4, 5, 2))
      par(fig=c(0.05,0.95,0,0.95), mar=c(5, 4, 5, 2)) # chang setting because the legend is removed
      plot(snp[1:length(snp)],as.numeric(p_value[1:(length(snp)),1]) , col=color.pvalueplot[1],ylim=c(0,10),xlim=xlim.temp, xlab="Physical position (Mb)",ylab=expression(-log[10] (p) ),type="h",main=paste("Chromosome ",chr.num,sep=""))
      par(oma=c(0,0,0,0),mar=c(0,0,0,0),fig=c(0,1,0,1), new=T)
      plot(1:10,type="n")
      # 2012-12-27 Remove the legend in p-value plot
      #ch <- paste("-LP(#W=1)",sep="")
      #legend(1, 3.5, ch, col=color.pvalueplot[1], pch=16, cex=1.2,title=expression(-Log[10](P)))
      dev.off()
    }else{}
    if (doFtest){
      for (stateposi in 1:ncol(total.f.pvalue)){
        jpeg(paste(pathnumPP, "Regression_test", FILENAME, "_Chr", chr, "_", substr(colnames(total.f.pvalue)[stateposi], 12, (nchar(colnames(total.f.pvalue)[stateposi])-11)), ".jpeg",sep=""), width = 1024, height = 768)
        #par(fig=c(0.2,0.95,0,0.95), mar=c(5, 4, 5, 2))
        par(fig=c(0.05,0.95,0,0.95), mar=c(5, 4, 5, 2)) # chang setting because the legend is removed
        plot(snp[1:length(snp)],as.numeric(total.f.pvalue[1:(length(snp)), stateposi]) , col=color.pvalueplot[1],ylim=c(0,10),xlim=xlim.temp, xlab="Physical position (Mb)",ylab=expression(-log[10] (p) ),type="h",main=paste("Chromosome ",chr.num,sep=""))
        par(oma=c(0,0,0,0),mar=c(0,0,0,0),fig=c(0,1,0,1), new=T)
        plot(1:10,type="n")
        #ch <- expression(paste("-LP", (beta[(substr(colnames(total.f.pvalue)[stateposi], 12, (nchar(colnames(total.f.pvalue)[stateposi])-11)))] != 0), sep = ""))
        # 2012-12-27 Remove the legend in p-value plot
        #ch <- bquote(-LP(beta[.(substr(colnames(total.f.pvalue)[stateposi], 12, (nchar(colnames(total.f.pvalue)[stateposi])-11)))] != 0))
        #legend(1, 3.5, ch, col=color.pvalueplot[1], pch=16, cex=1.2,title=expression(-Log[10](P)))
        dev.off()
      }
    }else{}
    if (doGEE){
      for (stateposi in 1:ncol(total.g.pvalue)){
        jpeg(paste(pathnumPP, "GEE", FILENAME, "_Chr", chr, "_", substr(colnames(total.g.pvalue)[stateposi], 5, (nchar(colnames(total.g.pvalue)[stateposi])-11)), ".jpeg",sep=""), width = 1024, height = 768)
        #par(fig=c(0.2,0.95,0,0.95), mar=c(5, 4, 5, 2))
        par(fig=c(0.05,0.95,0,0.95), mar=c(5, 4, 5, 2)) # chang setting because the legend is removed
        plot(snp[1:length(snp)],as.numeric(total.g.pvalue[1:(length(snp))]) , col=color.pvalueplot[1],ylim=c(0,10),xlim=xlim.temp, xlab="Physical position (Mb)",ylab=expression(-log[10] (p) ),type="h",main=paste("Chromosome ",chr.num,sep=""))
        par(oma=c(0,0,0,0),mar=c(0,0,0,0),fig=c(0,1,0,1), new=T)
        plot(1:10,type="n")
        #ch <- expression(paste("-LP", (beta[substr(colnames(total.g.pvalue)[stateposi], 5, (nchar(colnames(total.g.pvalue)[stateposi])-11))] != 0), sep = ""))
        # 2012-12-27 Remove the legend in p-value plot
        #ch <- bquote(-LP(beta[.(substr(colnames(total.g.pvalue)[stateposi], 5, (nchar(colnames(total.g.pvalue)[stateposi])-11)))] != 0))
        #legend(1, 3.5, ch, col=color.pvalueplot[1], pch=16, cex=1.2,title=expression(-Log[10](P)))
        dev.off()
      }
    }else{}
  }
  ##########################################
  gc()
  #rm(fitcase, fitcontrol, snp, mu.norm, m11, m1, r, prop.control, prop.case, p_value, all.matrix)
  gc()
  ##########################################
}  ### END of for(chr.num in (Chr.num %w/o% 23))
if(PlotAllinOne){
  all.matrix.name <- colnames(all.matrix)
  for (iid in 1:3){
    #if (c(doAssotest, doFtest, doGEE)[iid]){
    for (all.matrix.namei in which(!is.na(pmatch(substr(all.matrix.name, 1, 3), c("-lo", "Reg", "GEE")[iid])))){
      jpeg(paste(pathnumPP, c("Wilcoxon_test", "", "")[iid], ifelse(iid == 1, "", all.matrix.name[all.matrix.namei]), FILENAME, "_AllChr.jpeg",sep=""), width = 1024, height = 768)
      par(mfcol =c(5,5), mar=c(4, 3, 2, 1))
      for(chr.num in (Chr.num %w/o% 23)){
        chr <- sprintf("%02.f", chr.num)
        data.temp <- read.csv(paste(pathnumAT, "Association_test", FILENAME, "_Chr",chr,".csv",sep=""))
        xlim.temp <- c(min(data.temp[, 3]), max(data.temp[, 3]))
        ydata.temp <- as.vector(data.temp[1:nrow(data.temp), all.matrix.namei])
        xdata.temp <- as.vector(data.temp[1:nrow(data.temp), 3])
        plot(xdata.temp, ydata.temp , col=color.pvalueplot[1],xlim=xlim.temp, ylim=c(0,10), xlab="",ylab="", type="h", main=paste("Chromosome ",chr.num,sep=""))
        mtext("Physical position (Mb)", 1, line=2, cex=0.6)
        mtext(expression(-log[10] (p) ), 2, line=2, cex=0.6)
        ##############
        gc()
      #  rm(data.temp)
        gc(reset=TRUE)
        ##############
      }
      if(length(Chr.num)<22){
        for(i in 1:(22-length(Chr.num %w/o% 23))){
          plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
        }
      }
      plot(1:10,main=expression(-Log[10](P)),cex.main=2,lty=1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
      plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
      par(oma=c(0,0,0,0),mar=c(0,0,0,0),fig=c(0,1,0,1), new=T);
      plot(1:10,type="n");
      #ifelse(iid == 1, legend(8.8, 6, paste("-LP(#W=1)",sep=""), col=color.pvalueplot[1],bty="n", pch=16, cex=2), legend(8.8, 6, expression(paste("-LP", (beta[state] != 0), sep = "")), col=color.pvalueplot[1],bty="n", pch=16, cex=2));
      if (iid == 1){
        legend(8.8, 6, paste("-LP(#W=1)",sep=""), col=color.pvalueplot[1],bty="n", pch=16, cex=2)
      }else{
        legend(8.8, 6, expression(paste("-LP", (beta != 0), sep = "")), col=color.pvalueplot[1],bty="n", pch=16, cex=2)
      }
      dev.off()
    }
   #}else{}
 }
}
} ### END of for(aa in 1:length(Chip12))
  #delete.file <- dir(paste(pathnumLOH, "temp/", sep = ""), pattern = ".txt", full.names = T)
  #for (deletei in 1:length(delete.file)) file.remove(delete.file[deletei])
  ptm <- round(proc.time() - ptm, 1)
  log.file <- paste("Association test finished: user, ", ptm[1], " seconds; elapsed, ", ptm[3], " seconds", sep = "")
  print(log.file)
  write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
  # delete the files in temp directory
  #unlink(paste(pathnumLOH, "temp/", sep = ""), recursive = T)
} ### END of if(!Match & (doAssotest || doFtest || doGEE))

#}else{} ### END of if (!plot.only) [6]
#############################################
############  biplot   ######################
#############################################
if(plotselect3){

##############################################
##### Functions of new speturm biplot    #####
#################################################################
##### Physical position information of p arm and qarm from  #####
##### the annotation file from Affymetrix 2008/11/25        #####
#################################################################
##### Affy 100K, 500K #####
# 100K: Chr 13, 14, 15, 22 only has q-arm #
# 500K: Chr 13, 14, 15, 21, 22 only has q-arm #


options(digits=10)
##### Chip12 #####
chip12_pq_100K<- c(120.928505,89.771345,90.193711,48.613152,46.287844,58.843839,57.846410,43.462675,43.559955,38.584925,51.383437,34.290302,18.425192,19.387587,19.852603,35.003380,22.048100,15.073063,24.194887,26.155905,10.039984,15.276762,57.543114)
chip12_pq_500K<- c(120.939802,89.750715,90.330595,48.667700,46.349901,58.835467,57.882795,43.776595,46.875500,39.064552,51.320589,34.341801,18.155568,19.362325,18.846941,34.996986,22.069633,14.961178,24.160374,26.200755,14.131220,15.258423,58.075602)
##### Chip1 (Hind, Nsp) #####
chip1_pq_100K<- c(120.928505,89.563119,90.193711,48.315779,45.749379,58.843839,57.616519,43.462675,43.549615,39.114808,51.207947,34.290302,18.425192,19.387587,19.852603,35.003380,21.950207,15.018265,24.194887,25.597055,10.019671,15.276762,57.543114)
chip1_pq_500K<- c(120.227505,89.750715,90.330595,48.336433,46.349901,58.835467,57.882795,43.776595,38.737064,39.064552,51.258029,33.928129,18.110262,19.362325,18.846941,34.996986,22.009077,14.951870,24.146152,26.170699,14.131220,14.490036,58.075602)
##### Chip2 (Xbar, Sty) #####
chip2_pq_100K<- c(120.902013,89.771345,90.045737,48.613152,46.287844,58.724642,57.846410,43.312864,43.559955,38.584925,51.383437,34.007171,18.321079,19.285288,19.208413,34.493676,22.048100,15.073063,24.164050,26.155905,10.039984,14.919559,57.458364)
chip2_pq_500K<- c(120.939802,88.911124,90.191897,48.667700,45.770840,58.761209,57.834342,43.227364,46.875500,38.961358,51.320589,34.341801,18.155568,19.272965,18.809644,34.823053,22.069633,14.961178,24.160374,26.200755,13.636494,15.258423,57.283339)

##### Affy 6.0 #####
# Affy 6.0: Chr 13, 14, 15, 22 only has q-arm #

chip_pq_Affy6<- c(120.992603,91.649657,90.366641,49.278866,46.419092,58.885344 ,58.023925,43.898071,46.765644,39.090896,51.420212,34.727104,17.943628,19.272965,18.331687,35.063218,22.159777,15.205871,24.313463,26.237925,10.189119,14.435070,58.513210,11.241114)

######################


dir.create(pathnumB <- paste(outpath, "/Numeric results/", LCSHLOH, " biplot/", sep=""), showWarnings = F)
dir.create(pathplotB <- paste(outpath , "/Graphic results/Biplot/", sep=""), showWarnings = F)
###############################################################
Mean <- function(x){sum(x)/length(x)}
Format <- function(x, width, ...){format(x, width = width, justify="right", ...)}
Sprintf <- function(x, width, ...) {sprintf(paste("% ", width, ".4f", sep=""), x, ...)}

################################################################
#group.name <- sprintf("%-10s", rep(as.matrix(name.popu), IDnumber))
#Popu.name <- sprintf(paste("%-", max(End_field - Start_field+1), "s",sep="" ), unlist(Label))
group.name <- sprintf("%-10s", unlist(lapply(c(control, patient), function(k){rep(as.matrix(name.popu)[1, k], IDnumber[[k]])})))
Popu.name <- sprintf(paste("%-", max(End_field - Start_field+1), "s",sep="" ), unlist(lapply(c(control, patient), function(k){Label[[k]]})))
################################################################

for(aa in 1:length(Chip12)){
  if(chiptype %in% 3:4) FILENAME <- NULL
  if(chiptype %in% 1:2) FILENAME <- paste("_Chip", paste(Chip12[[aa]], collapse=""), sep="")

if (PlotAllinOne){
jpeg(paste(pathplotB, "BP", FILENAME, "_AllChr.jpeg",sep=""), width = 1024, height = 768)
par(oma=c(0,1,0,1), fig=c(0,0.9,0,1), mar=c(3.3, 3.2, 4.5, 3.3), mfcol =c(5,5))
}else{}
###
#ind_name <- as.character(sequence(IDnumber))
#TextColor <- rep(color.biplot[1:length(IDnumber)], IDnumber)
ind_name <- as.character(unlist(lapply(c(control, patient), function(k){sequence(IDnumber[[k]])})))
TextColor <- unlist(lapply(c(control, patient), function(k){rep(color.biplot[k], IDnumber[[k]])}))
###
for(chr.num in Chr.num){
  chr <- sprintf("%02.f", chr.num)
  #if (sequencing){
  if (TRUE){
  #if (chiptype == 4){ # allow all customize chip, not just for sequencing data
  # Yu-Tin's SVD version (start)
  totalmAFchr <- 0
  for (k in 1:n.part[chr]){
    if (!plot.only){
      load(paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".RData", sep = ""))
    }else{
      load(paste(inpath, "/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".RData", sep = ""))
    }
    mAFchr <- as.matrix(combloh)
    MeanAF <- apply(mAFchr, 1, function(x) mean(as.numeric(x), na.rm = T))
    mAFchr <- matrix(as.numeric(mAFchr), nrow = length(MeanAF)) - MeanAF
    assign(paste("mAFchr", k, sep = ""), mAFchr)
    #save(paste("mAFchr", k, sep = ""), file = paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".mAFchr.RData", sep = ""))
    eval(parse(text = paste("save(mAFchr", k, ", file = \"", pathnumLOH, "/temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".mAFchr.RData\")", sep = "")))
    eval(parse(text = paste("rm(mAFchr", k, ")", sep = "")))
    totalmAFchr <- totalmAFchr + t(mAFchr)%*%mAFchr
    rm(mAFchr, MeanAF, combloh);gc(reset = T)
  }
  pcK <- 2
  evA_value <- eigen(totalmAFchr)$values
  evA_value <- round(evA_value*1000000)/1000000
  V <- eigen(totalmAFchr/nrow(totalmAFchr))$vectors; #V=round(V*1000)/1000
  u <- V[, 1:2]
  transposeV <- t(V)
  d <- sqrt(evA_value)
  ##### numberic output 1
  Length1 <- c(8, 14, 7, 18)
  prop.variation.name <- c("Variable", "Singular_value", "PEV_(%)", "Cumulative_PEV_(%)")
  PEV <- d^2/sum(d^2)
  cumulative.PEV <- sapply(1:length(PEV), function(x) sum(PEV[1:x]))
  prop.variation.info <- cbind(Format(paste("SV", 1:length(d), sep=""), Length1[1]),
    Format(Sprintf(d^2, ""), Length1[2]), Sprintf(PEV, Length1[3]), Sprintf(cumulative.PEV, Length1[4]))
  ##### numberic output 2
  diagEV <- diag(d)
  inverse_diagEV <- diag(1/d)
  inverse_diagEV[is.infinite(inverse_diagEV)] <- 0
  v <- NULL
  for(k in 1:n.part[chr]){
    load(paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".mAFchr.RData", sep = ""))
    v <- rbind(v, get(paste("mAFchr", k, sep = ""))%*%u)
    rm(list = paste("mAFchr", k, sep = ""))
  }
  v[, 1] <- v[, 1]*inverse_diagEV[1, 1]
  v[, 2] <- v[, 2]*inverse_diagEV[2, 2]
  lambda <- diag(d) #individual x individual
  gc(reset = T)
  lambda <- lambda[1:pcK, 1:pcK]
  # C: dimension = Individual_number x 2
  # G: dimension = SNP_number x 2
  # if(alpha == 1){
  C <- u
  G <- v
  if(Alphaselect == 1){
    G <- G%*%lambda
    #C <- C
  }else if(Alphaselect == 0){
    G <- G*sqrt(nrow(C)) # C*{# individuals}
    C <- t(lambda%*%t(C))/sqrt(nrow(C))
  }
  Length2 <- c(10, max(c(nchar(Popu.name), 4)), 7, 7)
  ind.component <- cbind(group.name,
    Format(Popu.name, Length2[2]),
    Format(sprintf("%.4f", C[, 1]), Length2[3]),
    Format(sprintf("%.4f", C[, 2]), Length2[4]))
  ind.component.name <- c("Population", Format("Name", Length2[2]), Format("Ind_C1", Length2[3]), Format("Ind_C2", Length2[4]))
  ###
  Length3 <- c(13, 13, 7, 7)
  if (!plot.only){
    load(paste(pathnumLOH, "temp/Num_Popu", 1, "_", Label[[1]][1], FILENAME, "_Chr", chr, ".RData", sep = ""))
  }else{
    load(paste(inpath, "/Num_Popu", 1, "_", Label[[1]][1], FILENAME, "_Chr", chr, ".RData", sep = ""))
  }
  marker.component <- cbind(as.character(partdata[, 1]), as.character(rep(NA, nrow(partdata))), Format(partdata[, 3], Length3[2]), Sprintf(G[, 1], Length3[3]), Sprintf(G[, 2], Length3[4]))
  marker.component.name <- c("Probe_Set", "RS_number", "Phy_Posi_(Mb)" , "Mar_C1", "Mar_C2")
  CG <- list(C = G, G = C)
  dd <- cbind(partdata[, 1], rep(NA, nrow(partdata)), partdata[, 3])
  rm(partdata, C, G)
  ###
  # Yu-Tin's SVD version (end)
  }else{}

  #if (!sequencing){
  if (FALSE){
  #if (chiptype != 4){
    #dd <- CombData(1, chr.num, chip12=ifelse(chiptype %in% 3:4, 3, paste(Chip12[[aa]], collapse="")))
    ### temporary management 20121227
    #if(number.popu > 1){
    #  for(i in 2:number.popu){
    #     dd <- cbind(dd, CombData(i, chr.num, chip12=ifelse(chiptype %in% 3:4, 3, paste(Chip12[[aa]], collapse="")))[, -(1:3)])
    #  }
    #}else{}
    # load .RData in temp folder then combine to "dd"
    dd <- c()
    for (k in 1:n.part[chr]){
      load(paste(pathnumLOH, "temp/chr", chr.num, ifelse(chiptype %in% 3:4, "", FILENAME), ".part", k, ".RData", sep = ""))
      #dd <- rbind(dd, cbind(as.matrix(partdatainfo), as.matrix(combloh)))
      dd <- rbind(dd, cbind(partdatainfo, combloh))
      rm(partdatainfo, combloh);gc(reset = T)
    }
    P  <- dd[, -(1:3)] - apply(dd[, -(1:3)], 1, Mean)
    S1  <-  svd.GB(as.matrix(P));
    ##### numberic output 1
    Length1 <- c(8, 14, 7, 18)
    prop.variation.name <- c("Variable", "Singular_value", "PEV_(%)", "Cumulative_PEV_(%)")
    PEV <- (S1$d)^2/sum((S1$d)^2)
    cumulative.PEV <- sapply(1:length(PEV), function(x) sum(PEV[1:x]))
    prop.variation.info <- cbind(Format(paste("SV", 1:length(S1$d), sep=""), Length1[1]),
      Format(Sprintf((S1$d)^2, ""), Length1[2]), Sprintf(PEV, Length1[3]), Sprintf(cumulative.PEV, Length1[4]))
    ##### numberic output 2
    CG <- svd.CG(S1, Alphaselect)
    Length2 <- c(10, max(c(nchar(Popu.name), 4)), 7, 7)
    ind.component <- cbind(group.name,
      Format(Popu.name, Length2[2]),
      Format(sprintf("%.4f",CG$G[,1]), Length2[3]),
      Format(sprintf("%.4f",CG$G[,2]), Length2[4]))
    ind.component.name <- c("Population", Format("Name", Length2[2]), Format("Ind_C1", Length2[3]), Format("Ind_C2", Length2[4]))
    ###
    Length3 <- c(13, 13, 7, 7)
    marker.component <- cbind(as.character(dd[,1]), as.character(dd[,2]), Format(dd[,3], Length3[2]), Sprintf(CG$C[,1], Length3[3]), Sprintf(CG$C[,2], Length3[4]))
    marker.component.name <- c("Probe_Set", "RS_number", "Phy_Posi_(Mb)" , "Mar_C1", "Mar_C2")
    ###
  }else{}  # END of if (!sequencing)

  save.info.dir <-  paste(pathnumB, "Num", FILENAME, "_Chr", chr, ".csv",sep="")
  if(chiptype == 4){
    outtable <- rbind(prop.variation.name, prop.variation.info, NULL, ind.component.name, ind.component, NULL, marker.component.name[-2], marker.component[, -2])
  }else{
    outtable <- rbind(c(prop.variation.name, NA), cbind(prop.variation.info, NA), NULL, c(ind.component.name, NA), cbind(ind.component, NA), NULL, marker.component.name, marker.component)
  }
  write.csv(outtable, file=save.info.dir, na = "", row.names=FALSE)
  ##############################################################################
  jpeg(paste(pathplotB, "BP", FILENAME, "_Chr",chr, ".jpeg", sep=""), width = 1024, height = 768);
  par(fig=c(0.125,0.875,0,0.95), mar=c(5, 4, 5, 2));
  dd3 <- as.numeric(dd[, 3])
  BIPLOT(CG$C, CG$G, row.col=marker.color(nrow(CG$C)), col.col=TextColor,
    ind_name_list = ind_name,
    row_name_list = round(dd3*10)/10,
    title = list(paste("Chromosome ",chr.num,sep=""),cex=1.5, col="green"))
  par(oma=c(0,0,0,0), mar=c(0,0,0,0), fig=c(0,1,0,1), new=T);
  plot(1:10,type="n");
  if(length(name.popu)<=32){
    legend(1, 2.5+0.2*{number.popu-2}, name.popu, col=color.biplot[c(1:number.popu)], pch=16, cex=1.2,title=paste(LCSHLOH, "biplot"));
  }
  if(length(name.popu) > 32 &length(name.popu) <= 60){
    legend(0.8, 2.5+0.18*{{number.popu/2}-2}, name.popu, col=color.biplot[c(1:72)],ncol=2, pch=16, cex=1,title=paste(LCSHLOH, "biplot"));
  }
  par(oma=c(0,0,0,0),fig=c(0.9,1,0,1), new=T);
  ##############################################################################
  ## p-arm and q-arm
 	if(chiptype %in% 1:3){
    ChipType <- switch(chiptype,
      "1" = "100K",
      "2" = "500K",
      "3" = "Affy6")
    ChipType <- paste("chip", ifelse(chiptype %in% 1:2, Chip12[[aa]], ""), "_pq_", sep="", ChipType)
    chip12_pq <- get(ChipType)
    if(chiptype %in% c(1, 3)) Q.Chr.num <- c(13, 14, 15, 22)
    if(chiptype %in% 2) Q.Chr.num <- c(13, 14, 15, 21, 22)
    if(chr.num %in% Q.Chr.num){
      qarm <- dd3[which.max(diff(dd3))] # start of qarm
      if(abs(chip12_pq[chr.num] - qarm)> 1){
        qarm <- chip12_pq[chr.num]
      }
    }
    else{
      Q <- which.max(diff(dd3))
      parm <- dd3[Q] # end of parm
      qarm <- dd3[Q+1] # start of qarm
      if(abs(chip12_pq[chr.num]-parm)>1){
        parm<- dd3[sum(dd3 < chip12_pq[chr.num])]
        qarm<- dd3[sum(dd3 < chip12_pq[chr.num])+1]
      }
    } 	 	
    if(chr.num %in% Q.Chr.num){demo.newcolor.q(qarm, dd3)}
   	else{demo.newcolor.pq(parm, qarm, dd3)}
 	}
  ##
  if(chiptype %in% 4) demo.color()
  dev.off()
  ################################
  BIPLOT(CG$C, CG$G, row.col=marker.color(nrow(CG$C)), col.col=TextColor,
    ind_name_list = ind_name,
    row_name_list = round(dd3*10)/10,
    title = list(paste("Chromosome ",chr.num,sep=""),cex=1.1, col="green"), All=TRUE)
  gc()
  rm(dd, dd3, Length1, Length2, Length3, CG, ind.component, marker.component)
  gc(reset=TRUE)	
}

if(length(Chr.num) < 23){
  for(i in 1:(23 - length(Chr.num)) ){
    plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  }
}

if (PlotAllinOne){
plot(1:10,main=paste(LCSHLOH, "biplot"), cex.main=2,lty=1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
plot(1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
par(oma=c(0,0,0,0), mar=c(0,0,0,0), fig=c(0,1,0,1), new=T);
plot(1:10,type="n");
if(number.popu <= 12){
  legend(8.8, 4, name.popu, col=color.biplot[c(1:number.popu)],bty="n", pch=16, cex=2);
}
if(number.popu > 12 & number.popu <= 32){
  legend(8.6, 4.1, name.popu, col=color.biplot[c(1:number.popu)],ncol=2,bty="n", pch=16, cex=1.5);
}
if(number.popu > 32 & number.popu <= 48){
  legend(8.3, 4.1, name.popu, col=color.biplot[c(1:number.popu)],ncol=3,bty="n", pch=16, cex=1.5);
}
if(number.popu > 48 & number.popu <= 60){
  legend(8.4, 4.1, name.popu, col=color.biplot[c(1:number.popu)],ncol=3,bty="n", pch=16, cex=1.3);
}

}else{}
dev.off()
} ###Chip_12

}

# delete the files in temp directory
if (file.exists(paste(outpath, "/Temp", sep = ""))){
  remove.files <- dir(paste(outpath, "Temp/", sep = ""), full.names = T)
  rm(remove.files)
}else{}
#if (length(remove.files) > 0) for (i in 1:length(remove.files)) file.remove(remove.files[i])
#unlink(paste(pathnumLOH, "temp", sep = ""), recursive = TRUE)
#if (file.exists(paste(pathnumLOH, "temp/", sep = "")) & version$os == "mingw32") shell(paste("rmdir ", gsub("/", "\\\\", paste(pathnumLOH, "temp/", sep = "")), sep = ""))
#rm(remove.files)

###################################
#####  Ending time of LOHAS   #####
###################################
#tcl("image","create","photo",image1,file=figure.run.step3)
log.file <- paste("Ending time for LOHAS: ", date(), sep = "")
print( log.file )
write(log.file, paste(outpath,"/Log.txt",sep=""), append = TRUE)
print(paste("Computation of LOHAS is finished.", sep = ""))
write(paste("Computation of LOHAS is finished.", sep = ""), paste(outpath,"/Log.txt",sep=""), append = TRUE)

if (any(names(warnings()) != "")){
  write("", paste(outpath,"/Log.txt",sep=""), append = TRUE)
  write("Warning messages:", paste(outpath,"/Log.txt",sep=""), append = TRUE)
  for (warningsi in 1:length(warnings()))
    write(names(warnings())[warningsi], paste(outpath,"/Log.txt",sep=""), append = TRUE)
}

### save parameter setting .RData
if (!plot.only){
  save.para <- na.omit(c("name.popu", "SnpCol.column", "ChrCol.column", "PosiCol.column", "CallCol.column", "RsCol.column",
                         "Rare_threshold", "Start_field", "End_field", "chiptype", if(chiptype %in% 1:2){"Chip"}else{NA}, 
                         "IDnumber", "Label", "n.part", "number.popu", "num.group", "control", "patient", "description.name"))
  #save(name.popu, SnpCol.column, ChrCol.column, PosiCol.column, CallCol.column, RsCol.column, Rare_threshold, Start_field, End_field, if(chiptype %in% 1:2){Chip}else{i},
  #IDnumber, Label, Chip12, n.part, number.popu, num.group, control, patient, file = paste(pathnumLOH, "temp/Parameter.RData", sep = ""))
  save(list = save.para, file = paste(pathnumLOH, "temp/Parameter.RData", sep = ""))
}else{}
# for release gui
if (length(args) <= 0) tkmessageBox(title = " Information ", message = "Computation of LOHAS is finished!", icon = "info", type = "ok")

#rm(list = ls());gc()