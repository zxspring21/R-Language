ptm <- proc.time()
library(pspline)
subfunction = paste(pathway_program, "/ALOHA_subfunction.r", sep = ""); source(subfunction)
pathway_input = data_input
pathway_output = data_output
MainName = unlist(strsplit(pathway_input, "/"))
MainName = MainName[length(MainName)]
pathway_output = paste(pathway_output, "/", MainName, sep = "")
dir.create(pathway_output, showWarnings = F) 
output_GR = paste(pathway_output, "/Graphical result", sep = "")
dir.create(output_GR, showWarnings = F)
output_NR = paste(pathway_output, "/Numerical result", sep = "")
dir.create(output_NR, showWarnings = F)
if(studygroup_choice == 2){
	output_AILOH = paste(output_NR, "/AF_AI_LOH", sep = "")
	dir.create(output_AILOH, showWarnings = F)
}
Lf.txt = paste(pathway_output, "/Log.txt", sep = "")
cat("###-----------------------------------------------### \n", file = Lf.txt, append = F, sep = "")
cat("###                ALOHA Log file                 ### \n", file = Lf.txt, append = T, sep = "")
cat("###-----------------------------------------------### \n", file = Lf.txt, append = T, sep = "")
pathway_sample_list = paste(pathway_output, "/Sample list", sep = "")
dir.create(pathway_sample_list, showWarnings = F)
cat("--------------------- Input data checking --------------------- \n")
cat("The input data checking is starting. \n")
expr_data = data_format_checking(pathway_input, studygroup_choice, Lf.txt, pathway_sample_list)
cat(paste(expr_data, "\n", sep = ""))
cat("The input data checking is finished. \n")
sample_list_file = list.files(pathway_sample_list, full.names = T)
sample_list_name = list.files(pathway_sample_list, full.names = F)  
sample_list_name = strsplit(sample_list_name, ".txt")
sample_list_name = unlist(sample_list_name)
nSL = length(sample_list_file) 
for(k in 1:nSL){ 
	A = read.delim(sample_list_file[k], stringsAsFactors = F)
	n = nrow(A) 
	gender = A[, 3]
	n_female = sum(gender == "F")
	n_male = n - n_female
	gender_statment = paste(n_male, " males and ", n_female, " females", sep = "")
	cat("There are ", gender_statment, " in ", sample_list_name[k], "\n", sep = "")
	if(studygroup_choice == 2 & sample_list_name[k] == "Control" & (n_male == 0 | n_female == 0)){
		  cat("\nALOHA warning message: There is only one gender in control study group.\n")
    }
}
if(studygroup_choice == 2) 
{
	if(AFref_choice[1] == 1) 
	{
		pathway_DB = paste(pathway, "DATABASE", sep = "")
		chip_file = list.files(pathway_DB, full.names = T)[chipsize] 
		popu_file = list.files(chip_file, full.names = T) 
		popu_file = popu_file[c(1, 5, 2, 4, 3)] 
		popu_file = popu_file[ipopu]
		pathway_AF = paste(popu_file, "/AF_ref", sep = "")
	} else { 
		calculating_AFref = paste(output_NR, "/User_provided_AFref" ,sep = "")
		dir.create(calculating_AFref, showWarnings = F)
		NULL_SNP = AFref_checking(AFref_input, calculating_AFref) 
		if(!is.null(NULL_SNP)){
			pathway_empty_SNP = paste(output_NR, "/Null_SNP_list.txt", sep = "")
			write.table(NULL_SNP, pathway_empty_SNP, row.names = F , quote = F, sep = "\t")
		}
		pathway_AF = calculating_AFref
		cat("ALOHA message: The AF reference of user provided is suitable for ALOHA calculation.\n")
	}
	AFref_file = list.files(pathway_AF, full.names = T)
	cat("---------------- AI/LOH calculation ---------------- \n")
	cat("AI and LOH detection ... \n")
	cat("    \n", file = Lf.txt, append = T, sep = "")
	cat("   ~AI and LOH detection are preparing. \n", file = Lf.txt, append = T, sep = "")
	SL_file = list.files(pathway_sample_list, full.names = T)
	SL_name = list.files(pathway_sample_list, full.names = F)
	popufile = list.files(pathway_input, full.names = T) 
	popuname = list.files(pathway_input, full.names = F)
	for(k in 1:length(popufile)){ 
		NR_CC = paste(output_AILOH, "/", popuname[k], sep = "") 
		dir.create(NR_CC, showWarnings = F)
		ind_file = list.files(popufile[k], full.names = T)
		ind_name = list.files(popufile[k], full.names = F)
		sample_table = read.delim(SL_file[k], stringsAsFactors = F)		
		for(i in 1:length(ind_name)){ 
			log_message = paste("       ~AI and LOH detection of ", ind_name[i], " are preparing. \n", sep = "")
			cat(log_message, file = Lf.txt, append = T, sep = "")
			cat(ind_name[i], "\n")
			output_ind = paste(NR_CC, "/", ind_name[i], sep = "")
			dir.create(output_ind, showWarnings = F)
            gender = sample_table[sample_table[, 2] == ind_name[i], 3]
			chr_file = list.files(ind_file[i], full.names = T)
			chr_name = list.files(ind_file[i], full.names = F)
			if(gender == "F") AFref_chr = AFref_file[c(1:23)]
			if(gender == "M") AFref_chr = AFref_file[c(1:22, 24)]
			for(j in 1:length(chr_file)){ 
				AF = read.delim(chr_file[j], stringsAsFactors = F)
				ref = read.delim(AFref_chr[j], stringsAsFactors = F) 
				AF[[1]] = as.numeric(sapply(strsplit(as.character(AF[[1]]), "SNP_A-"), function(x) x[2]));
				colnames(AF)[1] = "Probe_set" 
				colnames(ref)[1] = "Probe_set"
				AF = AF[, -c(2, 5)] 
				AFref = ref[, -c(2:3)] 
				AF.merge = merge(AF, AFref, by = "Probe_set")
				AF.merge = Double_pp_checking(AF.merge)
				nSNP = nrow(AF.merge)
				AF.ind = AF.merge[, 1:4]
				AF.ref = AF.merge[, -c(3, 4)] 
				AIchr = AI(AF.ind, AF.ref, a, nSNP) 
				LOHchr = LOH(AF.ind, AF.ref, a, nSNP)
				stopifnot(all(AIchr$Probe_set == LOHchr$Probe_set))
				AI_LOH = cbind(AIchr, LOHchr[5])
				AI_LOH = AI_LOH[order(AI_LOH$Phy_position),]
				pathway_chr_result = paste(output_ind, "/", chr_name[j], sep = "")
				write.table(AI_LOH, pathway_chr_result, col.names = T, row.names = F, quote = F, sep = "\t")
			} 
			log_message = paste("       ~AI and LOH detection of ", ind_name[i], " are finished. \n", sep = "")
			cat(log_message, file = Lf.txt, append = T, sep = "")
		} 
	} 
	cat("   ~AI and LOH detection are finished. \n", file = Lf.txt, append = T, sep = "")
	cat("AI and LOH detection are finished. \n")
	cat("The calculation of AI and LOH proportion are preparing. \n")
	log_message = paste("  \n", sep = ""); cat(log_message, file = Lf.txt, append = T, sep = "")
	log_message = paste("   ~The calculation of AI and LOH proportion in are preparing. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	NR_file = list.files(output_AILOH, full.names = T)
	nNRfile = length(NR_file) 
	for(k in 1:nNRfile) 
	{
		ind_file = list.files(NR_file[k], full.names = T)
		ind_name = list.files(NR_file[k], full.names = F)
		nind = length(ind_name)
		for(i in 1:nind) 
		{
			chr_file = list.files(ind_file[i], full.names = T)
			chr_name = list.files(ind_file[i], full.names = F)
			nchr = length(chr_name)
			for(j in 1:nchr) 
			{
				AI_LOH = read.delim(chr_file[j], stringsAsFactors = F)
				index.AI = AI_LOH[, 5]
				AIprop = WindowSlidingProp(index.AI, w)
				index.LOH = AI_LOH[, 6]
				LOHprop = WindowSlidingProp(index.LOH, w)
				AIprop = as.numeric(AIprop)
				AIprop = round(AIprop*100000)/100000
				AIprop = format(AIprop, digits = 10)
				LOHprop = as.numeric(LOHprop)
				LOHprop = round(LOHprop*100000)/100000
				LOHprop = format(LOHprop, digits = 10)
				AI_LOH = cbind(AI_LOH, AIprop, LOHprop)
				write.table(AI_LOH, chr_file[j], row.names = F, quote = F, sep = "\t")
			}
		}
	}
	log_message = paste("   ~ The calculation of AI and LOH proportion are finished. \n\n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	cat("The calculation of AI and LOH proportion are finished. \n")
	cat("Replace AI/LOH proportion of control group by quantile... \n")
	log_message = paste("   ~Replacement of AI and LOH proportion in control group are preparing. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	pathway_Control = paste(output_NR, "/AF_AI_LOH/Control", sep = "")
	ind_file = list.files(pathway_Control, full.names = T)
	ind_name = list.files(pathway_Control, full.names = F)
	pathway_Ref = paste(output_NR, "/Control_AILOH-Ref", sep = "")
	dir.create(pathway_Ref, showWarnings = F)
	control_sample_list = paste(pathway_sample_list, "/Control.txt", sep = "")
	control_table = read.delim(control_sample_list, colClasses = c(NULL, NULL, "character", NULL))
	gender = as.character(control_table[, 3])
	f_ind_group = ind_file[gender == "F"]
	nFemale = nrow(f_ind_group)
	m_ind_group = ind_file[gender == "M"]
	nMale = nrow(m_ind_group)
	auto_chr = NULL
	chr_name = c(list.files(ind_file[1], full.names = F), "Chr_24.txt") 
	AI_LOH_prop_control = lapply(1:24, function(i){
		gc(verbose = F, TRUE)
		if(i < 23){
			temp_ind_file = ind_file
			ichr = i
		} else if(i == 23){ 
			temp_ind_file = f_ind_group 
			ichr = i
		} else if(i == 24){ temp_ind_file = m_ind_group; ichr = 23 } 
		if(length(temp_ind_file) == 0){
			chr = NULL
			chr = matrix(rep(NA, 4*10), 10, 4)
			chr = data.frame(chr)
			colnames(chr)[1] = "Probe_set"
			colnames(chr)[2] = "Phy_position"
			colnames(chr)[3] = "qAI"
			colnames(chr)[4] = "qLOH"
			return(chr)
		} else {
			chr_file = sapply(temp_ind_file, function(m) list.files(m, full.names = T)[ichr]) 
			ind_chr = lapply(chr_file, function(j) read.delim(j, colClasses = c(rep("NULL", 6), rep("numeric", 2))))
			names(ind_chr) = NULL
			chr_AI = sapply(ind_chr, function(x) x$AIprop) 
			chr_LOH = sapply(ind_chr, function(x) x$LOHprop) 
			chr_info = read.delim(chr_file[1], colClasses = c(rep("integer", 2), rep("NULL", 6)))
			qAI = apply(chr_AI, 1, quantile, q_AILOH)
			qLOH = apply(chr_LOH, 1, quantile, q_AILOH)
			chr = cbind(chr_info, qAI, qLOH)
		}
	})
	auto_chr = NULL
	for(autosome in 1:22) auto_chr = rbind(auto_chr, AI_LOH_prop_control[[autosome]][, c("qAI", "qLOH")])
	female_chr = AI_LOH_prop_control[[23]][, c("qAI", "qLOH")]
	male_chr = AI_LOH_prop_control[[24]][, c("qAI", "qLOH")]
	upper.AILOH.prop = sapply(list(autosome = auto_chr, female = female_chr, male = male_chr), function(x) {
		apply(x, 2, function(y) {
			m = mean(y)
			s = sd(y)
			m + 6*s 
		})
	})
	for(i in 1:24){ 
		chr_AILOH_prop = AI_LOH_prop_control[[i]]
		AI_prop = chr_AILOH_prop$qAI
		LOH_prop = chr_AILOH_prop$qLOH
		if(i < 23){
			AI_upper = upper.AILOH.prop["qAI", "autosome"]
			LOH_upper = upper.AILOH.prop["qLOH", "autosome"]
		} else if(i == 23){
			AI_upper = upper.AILOH.prop["qAI", "female"]
			LOH_upper = upper.AILOH.prop["qLOH", "female"]
		} else {
			AI_upper = upper.AILOH.prop["qAI", "male"]
			LOH_upper = upper.AILOH.prop["qLOH", "male"]
		}
		AI_prop[AI_prop > AI_upper] = AI_upper
		LOH_prop[LOH_prop > LOH_upper] = LOH_upper
		chr_AILOH_prop$qAI = AI_prop
		chr_AILOH_prop$qLOH = if(i != 24){ LOH_prop } else NA
		pathway_output_chr = paste(pathway_Ref, "/", chr_name[i], sep = "")
		write.table(chr_AILOH_prop, pathway_output_chr, col.names = T, row.names = F, sep = "\t", quote = F)
	}
	log_message = paste("   ~Replacement of AI and LOH proportion in control group are finished. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	log_message = paste("  \n", sep = ""); cat(log_message, file = Lf.txt, append = T, sep = "")
	log_message = paste("   ~AI and LOH graphical plotting are preparing. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	cat("Plotting graphical result ... \n")
	ref_file = list.files(pathway_Ref, full.name = T)
	ref_name = list.files(pathway_Ref, full.name = F)
	pathway_case = paste(output_NR, "/AF_AI_LOH/Case", sep = "")
	ind_file = list.files(pathway_case, full.names = T)
	ind_name = list.files(pathway_case, full.names = F)
	nind = length(ind_name)
	Fref = ref_file[c(1:23)]; Mref = ref_file[c(1:22, 24)]
	output_AFAILOH_figure = paste(output_GR, "/AILOH figure", sep = "")
	dir.create(output_AFAILOH_figure, showWarnings = F)
	output_AFAILOH_figure_indbase = paste(output_AFAILOH_figure, "/Each individual", sep = "") 
	dir.create(output_AFAILOH_figure_indbase, showWarnings = F) 
	case_sample_list = paste(pathway_sample_list, "/Case.txt", sep = "")
	case_table = read.delim(case_sample_list, stringsAsFactors = F)
	case_gender = case_table[, 3]
	for(i in 1:nind){ 
		gender = case_gender[i]
		if(gender == "F") ref_chr = Fref
		if(gender == "M") ref_chr = Mref
		chr_file = list.files(ind_file[i], full.names = T) 
		chr_name = list.files(ind_file[i], full.names = F)
		nchr = length(chr_name)
		whole_chr = NULL;
		for(j in 1:nchr){ 
			ind_chr = read.delim(chr_file[j], stringsAsFactors = F)
			ref = read.delim(ref_chr[j], stringsAsFactors = F, colClasses = c("integer", "NULL", "numeric", "numeric"))
			if(j == 23 & gender == "M"){
				ind_chr[, 8] = NA
				write.table(ind_chr, chr_file[j], row.names = F, quote = F, sep = "\t")
			}
			stopifnot(all(ind_chr$Probe_set == ref$Probe_set))
			chr = merge(ind_chr, ref, by = "Probe_set")
			if(!nrow(chr)){ 
				chr = cbind(ind_chr, qAI = NA, qLOH = NA)                                                   
			}
			chr = cbind(Chr = j, chr)
			chr = chr[order(chr$Phy_position),]
			whole_chr = rbind(whole_chr, chr)
		}
		whole_chr$Probe_set = NULL  
		whole_chr$Phy_position = whole_chr$Phy_position/1000000
		pp.temp = PPcum(whole_chr[, c("Chr", "Phy_position")])
		pp = pp.temp[[1]]
		pp_bound = pp.temp[[2]]
		ind_info = whole_chr
		ind_info[, 2] = pp[, 2]
		Ref = ind_info[, c("Chr", "qAI", "qLOH")] 
		plot_list = c(1:3)
		nplot = length(plot_list)
		pathway_output_ind = paste(output_AFAILOH_figure_indbase, "/", ind_name[i], ".png", sep = "")
		png(filename = pathway_output_ind, width = 1680, height = 1050)
		layout(matrix(1:nplot, nplot, 1))
		par(mai = c(0.25, 0.5, 0.25, 0.2))
		def_set = rep(10, 23)  
		for(iplot in plot_list)
		{
			if(iplot == 1)
			{
				AF_table = ind_info[, c(2, 3, 4, 5, 1)] 
				yregion = AF_figure(AF_table)
				label_setting = c("Allele frequency", "AF")
				figure_public_setting(pp_bound, yregion, label_setting)
			} else if(iplot == 2){
				AI_table = cbind(ind_info[, c(2, 7, 5)], Ref[, c(2, 1)]) 
				yregion = AILOH_figure(AI_table, def_set)
				label_setting = c("Allelic imbalance", "AI proportion")
				figure_public_setting(pp_bound, yregion, label_setting)
			} else {
				LOH_table = cbind(ind_info[, c(2, 8, 6)], Ref[,c(3, 1)]) 
				yregion = AILOH_figure(LOH_table, def_set)
				label_setting = c("Loss of heterozygosity", "LOH proportion")
				figure_public_setting(pp_bound, yregion, label_setting)
			}
		}
		dev.off()
	}
	log_message = paste("   ~AI and LOH graphical plotting are finished. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	cat("Combined all individuals' AI/LOH figure plotting is preparing..\n")      
	prop_AILOH=paste(output_NR,"/Prop_AILOH",sep="")
	dir.create(prop_AILOH,showWarnings=F)
	prop_Ref=paste(prop_AILOH,"/Ref_prop",sep="")
	dir.create(prop_Ref,showWarnings=F)
	prop_Case=paste(prop_AILOH,"/Case_prop",sep="")
	dir.create(prop_Case,showWarnings=F)                
	ref_file=list.files(pathway_Ref,full.name=T)
	ref_name=list.files(pathway_Ref,full.name=F)      
	nchr=length(ref_name)
	for(i in 1:nchr){
		ref=read.table(ref_file[i],head=T,sep="\t")
		if(sum(as.numeric(is.na(ref)))!=40){
			temp_pp=ref[,2];
			temp_ref_AI=ref[,3]
			temp_ref_LOH=ref[,4]                      
			if(sum(is.na(temp_ref_AI))==0){ref_AI<-sm.spline(temp_pp, temp_ref_AI, df=def_set);ref_AI=unlist(ref_AI[3])}
			if(!sum(is.na(temp_ref_AI))==0){ref_AI=NA}            
			if(sum(is.na(temp_ref_LOH))==0){ref_LOH<-sm.spline(temp_pp, temp_ref_LOH, df=def_set);ref_LOH=unlist(ref_LOH[3])}
			if(!sum(is.na(temp_ref_LOH))==0){ref_LOH=NA}            
			ref=cbind(ref[,c(1,2)],ref_AI,ref_LOH)            
		} 
		colnames(ref)[3]="Ref_AI_sm"         
		colnames(ref)[4]="Ref_LOH_sm"
		sm_result=paste(prop_Ref,"/",ref_name[i],sep="")
		write.table(ref,sm_result,col.names=T,row.names=F,quote=F,sep="\t")                   
	}
	ref_file=list.files(prop_Ref,full.names=T)
	Male_ref=ref_file[c(1:22,24)]
	Female_ref=ref_file[c(1:23)]        
	ind_file=list.files(pathway_case,full.names=T)
	ind_name=list.files(pathway_case,full.names=F)
	nind=length(ind_name)
	case_table=paste(pathway_sample_list,"/Case.txt",sep="")
	case_table=read.table(case_table,head=T,sep="\t",colClasses=c(NULL,NULL,"character",NULL))
	case_table=case_table[,c(2,3)]
	ind_list=cbind(ind_name,ind_file)
	colnames(ind_list)[1]="Ind_ID"
	ind_list=merge(ind_list,case_table,by="Ind_ID")                  
	whole_AI=NULL;
	whole_LOH=NULL;        
	for(i in 1:nind){
		chr_file=list.files(as.character(ind_list[i,2]),full.names=T)
		gender=as.character(ind_list[i,3])
		if(gender=="F"){Ref_file=Female_ref}
		if(gender=="M"){Ref_file=Male_ref}          
		nchr=length(chr_file)
		whole_pp=NULL;
		for(j in 1:nchr){
			ref=read.table(Ref_file[j],head=T,sep="\t")
			if(sum(as.numeric(is.na(ref)))==40){
				chr=read.table(chr_file[j],head=T,sep="\t")
				pp=cbind(j,chr[,2],2,2) 
				colnames(pp)[1]="j"
				colnames(pp)[2]="Phy_position"              
				colnames(pp)[3]="AI"
				colnames(pp)[4]="LOH"              
			}
			if(sum(as.numeric(is.na(ref)))!=40){
				ref=ref[,-2]
				chr=read.table(chr_file[j],head=T,sep="\t")
				chr=chr[c(1,2,7,8)]
				temp_pp=chr[,2]
				temp_AI=chr[,3]
				temp_LOH=chr[,4]
				if(sum(is.na(temp_AI))==0){ind_AI<-sm.spline(temp_pp, temp_AI, df=def_set);ind_AI=unlist(ind_AI[3])}
				if(!sum(is.na(temp_AI))==0){ind_AI=NA}                
				if(sum(is.na(temp_LOH))==0){ind_LOH<-sm.spline(temp_pp, temp_LOH, df=def_set);ind_LOH=unlist(ind_LOH[3])}
				if(!sum(is.na(temp_LOH))==0){ind_LOH=NA}                                
				ind=cbind(chr[,c(1,2)],ind_AI,ind_LOH) 
				colnames(ind)[3]="ind_AI_sm"         
				colnames(ind)[4]="ind_LOH_sm"
				chr=merge(ind,ref,by="Probe_set")
				chr=sort.append(chr,2)
				chr=cbind(j,chr)
				pp=chr[,c(1,3)];              
				AI=chr[,4]-chr[,6]
				AI[AI>0]=1;AI[AI<=0]=0
				AI[is.na(AI)==T]=2
				LOH=chr[,5]-chr[,7]  
				LOH[LOH>0]=1;LOH[LOH<=0]=0                      
				LOH[is.na(LOH)==T]=2
				pp=cbind(pp,AI,LOH) 
			}                      
			whole_pp=rbind(whole_pp,pp)
		}
		whole_AI=cbind(whole_AI,whole_pp[,3])
		whole_LOH=cbind(whole_LOH,whole_pp[,4])          
		if(i==nind){
			pp=PPcum(whole_pp[,c(1,2)])[[1]]
			pp_bound=PPcum(whole_pp[,c(1,2)])[[2]]
		}
	}
	gender=c("chr","pp",ind_list[,3])
	AI=cbind(pp,whole_AI);AI=rbind(gender,AI)      
	LOH=cbind(pp,whole_LOH);LOH=rbind(gender,LOH)
	GR_ALL=paste(output_AFAILOH_figure,"/All samples",sep="")              
	dir.create(GR_ALL,showWarnings=F)
	AI_GR=paste(GR_ALL,"/Abnorm_AI.png",sep="")                      
	png(filename=AI_GR,width=1680,height=1200)
	All_sample_AILOH_figure(AI,pp_bound,"Allelic imbalance")
	dev.off()
	LOH_GR=paste(GR_ALL,"/Abnorm_LOH.png",sep="")
	png(filename=LOH_GR,width=1680,height=1200)        
	All_sample_AILOH_figure(LOH,pp_bound,"Loss of heterozygosity")        
	dev.off()
	cat("Combined all individuals' AI/LOH figure plotting is finished..\n")  
} 
log_message = paste("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
log_message = paste(" \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
log_message = paste("   ~Biplot calculation is preparing. \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
cat("--------------------- Biplot calcaulation --------------------- \n")
inpath <- pathway_input
outpath <- pathway_output
alpha <- a_biplot
outpath = paste(outpath, "/Biplot_Temp", sep = "")
dir.create(outpath, showWarnings = F)
dir.create(paste(outpath, "/temp", sep = ""), showWarnings = F)
if(studygroup_choice == 2){
	popu_file = list.files(pathway_input, full.names = T)
	popu_name = list.files(pathway_input, full.names = F)
	popu_name = popu_name[c(1, 2)]
}
if(studygroup_choice == 1) popu_name = list.files(inpath, full.names = F)
TarDir = popu_name
LEGEND = popu_name
AllCol <- NULL
AllIndivIndex <- NULL
AllIndivNum <- NULL
AllIndivName <- NULL
COL <- c("red", "blue", "orange", "cyan", "green", "purple", "chocolate", "cadetblue", "magenta", "olivedrab")
SNP_table = NULL
for(chr in 1:23){ 
	chrname <- ifelse(chr < 10, paste("0", chr, sep = ""), chr)
	log_message = paste("       ~The biplot calculation of chromosome ", chrname, " is preparing. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
	AFchr = NULL
	chr_snp = NULL
	for(i in 1:length(TarDir)){ 
		nowpath <- paste(inpath, TarDir[i], sep = "/")
		IDlist <- dir(nowpath)
		ID.check <- file.info(paste(nowpath, "/", IDlist, sep = ""))$isdir
		IDlist <- IDlist[ID.check]
		if(chr == 1){
			AllCol <- c(AllCol, COL[i])
			AllIndivNum <-c(AllIndivNum, length(IDlist))
			AllIndivName <-c(AllIndivName, IDlist)
			AllIndivIndex <- c(AllIndivIndex, 1:length(IDlist))
		}
		OutAF = NULL
		for(ind in 1:length(IDlist)){ 
			nowID <- IDlist[ind]
			A = read.delim(paste(nowpath, "/", nowID, "/Chr_", chrname, ".txt", sep = ""))
			A = A[, c(1, 6)] 
			AF = A[[2]]
			OutAF = cbind(OutAF, AF)
		} 
		probe_set = as.character(A[[1]])
		nSNP = nrow(A)
		chr_snp = c(chr_snp, nSNP)
		MeanAF <- apply(OutAF, 1, function(x) mean(x, na.rm = T))
		tmp <- sapply(1:nrow(OutAF), function(ii){
			nowAF <- OutAF[ii, ]
			OutAF[ii, which(is.na(nowAF))] <<- MeanAF[ii]
		})
		rm("tmp")
		OutAF = cbind(probe_set, OutAF) 
		head_name = c("probe_set", IDlist)
		colnames(OutAF) = head_name
		if(i == 1) AFchr = OutAF
		if(i != 1) AFchr = merge(AFchr, OutAF, by = "probe_set")
	} 
	AFchr = AFchr[, -1] 
	AFchr = as.matrix(AFchr)
	AFchr = matrix(as.numeric(AFchr), nrow(AFchr), ncol(AFchr))
	rm("OutAF")
	NAnum = apply(AFchr, 1, function(x) length(which(is.na(x))))
	posi = which(NAnum == ncol(AFchr))
	if(length(as.matrix(posi)) == 0) mAFchr = AFchr
	if(length(as.matrix(posi)) != 0) mAFchr = AFchr[-posi, ]
	rm("NAnum", "AFchr")
	MeanAF = apply(mAFchr, 1, function(x) mean(x, na.rm = T))
	mAFchr = mAFchr - MeanAF
	rm("MeanAF")
	K <- 2
	N <- sum(AllIndivNum)
	S1 <- svd(t(mAFchr)) 
	u <- S1$u 
	v <- S1$v 
	lambda <- diag(S1$d) 
	u <- u[, 1:K]
	v <- v[, 1:K]
	lambda <- lambda[1:K, 1:K]
	if(alpha == 1){
		C <- u%*%lambda
		G <- v
	} else {
		C <- sqrt(N)*u
		G <- t(lambda%*%(t(v))/sqrt(N))
	}
	write.table(S1$d, paste(outpath, "/temp/lambdaforBP_chr", chrname, ".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
	write.table(C, paste(outpath, "/temp/CforBP_chr", chrname, ".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
	write.table(G, paste(outpath, "/temp/GforBP_chr", chrname, ".txt", sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
	rm("S1", "u", "v", "lambda")
	BPMG.chr(AllIndivIndex, AllIndivNum, AllIndivName, AllCol, chr, LEGEND, outpath)
	cat("Calculation of AF biplot for chromosome ", chr, " is finished.\n", sep = "")
	SNP_table = rbind(SNP_table, chr_snp) 
	log_message = paste("       ~The biplot calculation of chromosome ", chrname, " is finished. \n", sep = "")
	cat(log_message, file = Lf.txt, append = T, sep = "")
} 
BPMG.WG(AllIndivIndex, AllIndivNum, AllIndivName, AllCol, LEGEND, outpath)
bp_file = list.files(outpath, pattern = ".tiff", full.names = T)
if(length(bp_file)){
	bp_name = list.files(outpath, pattern = ".tiff", full.names = F)
	GR_biplot = paste(output_GR, "/Biplot", sep = "")
	dir.create(GR_biplot, showWarnings = F)
	GR_bp = paste(GR_biplot, "/", bp_name, sep = "")
	for(j in 1:24) file.rename(bp_file[j], GR_bp[j])
}
bp_temp_NR = paste(outpath, "/temp", sep = "")
output_BP_NR = paste(output_NR, "/Biplot", sep = "")
dir.create(output_BP_NR, showWarnings = F)
NR_bp = list.files(bp_temp_NR, full.names = T)
popu_file = list.files(inpath, full.names = T)
npopu = length(popu_name)
ind_name = NULL
for(j in 1:npopu){
	popu_ind_name = list.files(popu_file[j])
	popu_ind_name = paste(popu_name[j], "_", popu_ind_name, sep = "")
	ind_name = c(ind_name, popu_ind_name)
}
ind_file = list.files(popu_file[1], full.names = T)[1]
chr_file = list.files(ind_file[1], full.names = T)
for(k in 1:23){
	chip_info = read.delim(chr_file[k])[, c(1, 3)] 
	C_table = read.delim(NR_bp[k], header = F)
	G_table = read.delim(NR_bp[k + 23], header = F)
	L_table = read.delim(NR_bp[k + 46], header = F)
	L_table = unlist(L_table)
	A = list(G_table, C_table, L_table)
	if(k < 10){
		ichr = paste("0", k, sep = "")
	} else ichr = k
	chr_txt = paste(output_BP_NR, "/Chr_", ichr, ".txt", sep = "") 
	Biplot_NR_manage(chr_txt, A, ind_name, chip_info)
}
unlink(outpath, recursive = T)
log_message = paste("   ~Biplot calculation is finished. \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
log_message = paste("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
log_message = paste(" \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
log_message = paste("   ~Data description is preparing. \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
cat("Data description!!! \n")
data_descrption = paste(output_NR, "/Data description.txt", sep = "")
single_line = c("-----------------------------------------------------")
double_line = c("=========================================================")
title_part = paste("===============      Data description     ===============", sep = "")
title_part = rbind(double_line, title_part, double_line)
write.table(title_part, file = data_descrption, row.names = F, quote = F, col.names = F, append = F, sep = "\t")
write.table(c(" "), file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table("1. Input/output path -", file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
DD_pathway_input = paste("   Input directory name: ", pathway_input, sep = "")
DD_pathway_output = paste("   Output directory name: ", pathway_output, sep = "")
write.table(DD_pathway_input, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table(DD_pathway_output, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
DD_num_group = paste("   The number of study groups: ", length(list.files(pathway_input)), sep = "")
write.table(DD_num_group, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
for(i in 1:length(TarDir)){
	DD_numSNP_popu = paste("   The number of SNPs in ", TarDir[i], ": ", sum(SNP_table[, i], na.rm = T), sep = "")
	write.table(DD_numSNP_popu, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
}
efficient_SNPs = sum(SNP_table[, 1], na.rm = T)
DD_num_effSNP_popu = paste("   The number of identical SNPs in all populations: ", efficient_SNPs, sep = "")
write.table(DD_num_effSNP_popu, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
DD_popu_file = list.files(pathway_input, full.names = T)
for(i in 1:length(DD_popu_file)){ 
	indname = list.files(DD_popu_file[i], full.names = F)
	nind = length(indname)
	DD_num_ind = paste("   The number of individuals in ", TarDir[i], ": ", nind, sep = "")
	write.table(DD_num_ind, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
}
write.table(c(" "),file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table("2. AF reference -", file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
if(sum(AFref_choice) == 0) AFref_resource_option = "NA"
if(AFref_choice[1] == 1) AFref_resource_option = "Database"
if(AFref_choice[2] == 1) AFref_resource_option = "User provided"
temp_AFref_resource = paste("   Reference: ", AFref_resource_option, sep = "");
write.table(temp_AFref_resource, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
chip_list = c("Affymetrix 100K", "Affymetrix 500K")
popu_list = c("Asia (CHB+JPT)", "YRI", "CEU", "TWN", "Combined") 
if(AFref_choice[1] == 1){
	temp_population = paste("   ALOHA database: ", popu_list[ipopu], sep = "");
	write.table(temp_population, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
	temp_chiptype = paste("   Genome-wide SNP chip: ", chip_list[chipsize], sep = "")
	write.table(temp_chiptype, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
}
if(AFref_choice[2] == 1){
	AFref_exp = paste("   Input directory name :  ", AFref_input, sep = "");
	write.table(AFref_exp, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
}
write.table(c(" "),file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table("3. AI/LOH analysis -", file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
temp_Alpha = paste("   Confidence level: ", ifelse(studygroup_choice == 1, NA, 1-a), sep = "")
write.table(temp_Alpha, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
temp_window = paste("   Window size: ", ifelse(studygroup_choice == 1, NA, w), sep = "")
write.table(temp_window, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
temp_qAILOH = paste("   Upper bound (Quantile) of reference: ", ifelse(studygroup_choice == 1, NA, q_AILOH), sep = "")
write.table(temp_qAILOH, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table(c(" "),file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table("4. AF biplot -", file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
temp_Alpha_biplot = paste("   Alpha scaling: ", a_biplot, sep = "")
write.table(temp_Alpha_biplot, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
write.table(c(" "), file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
log_message = paste("   ~Data description is finished. \n", sep = "")
cat(log_message, file = Lf.txt, append = T, sep = "")
u = proc.time() - ptm
sec = u[3]
hour = sec/3600
minute = (hour - floor(hour))*60
second = (minute - floor(minute))*60
write.table(" ", file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
time_stance = paste("Elapsed time: ", floor(hour), "-H,  ", floor(minute), "-M,  ", floor(second), "-S ", sep = "")
write.table(time_stance, file = data_descrption, row.names = F, quote = F, col.names = F, append = T, sep = "\t")
cat("Computation of ALOHA is finished. \n")

