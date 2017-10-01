studygroup_choice = as.numeric(tclvalue(inputdata_format.statsValue));  
data_input = tclvalue(pathway_input); 
if(data_input == "") stop("ALOHA error message: The directory of data input is incorrect.")
data_output = tclvalue(pathway_output)
if(data_input %in% c("Test1", "Test2")){
	if(!file.exists(data_output)){
		data_output = paste(pathway, "OUTPUT/Test_Example_Output", sep = "")
		dir.create(data_output, showWarnings = F)
	}
	cat("The output path is ", data_output, "\n")
	data_input = paste(pathway, "EXAMPLE/", data_input, sep = "")
} else if(!file.exists(data_output)){
	stop("ALOHA error message: The directory of result output is incorrect.")
}
group_choice = as.numeric(tclvalue(groupValue))
chipsize_choice = as.numeric(tclvalue(chipValue))
AFref_choice = as.numeric(c(tclvalue(database.val), tclvalue(userprovided.val)))
AFref_input = tclvalue(pathway_AFrefinput) 
if(studygroup_choice == 2 & AFref_choice[2] == 1 & AFref_input == "") 
	stop("ALOHA error message: The directory of user-provided AF reference is incorrect.")
chipsize = chipsize_choice      
ipopu = group_choice
if(studygroup_choice == 2){
	a = as.numeric(tclvalue(Alphascaling.val))
	if(a < 0.9 | a > 1) {
	  stop("ALOHA error message: The inputted confidence level is incorrect (the value should range from 0.9 to 1).")
	}
	a = 1-a 
	w = as.numeric(tclvalue(windowsize.val))
	if(w < 10) {
	  stop("ALOHA error message: The inputted window size is incorrect (the value should be an integer >= 10).")
	}
	q_AILOH = as.numeric(tclvalue(qRef.val))
	if(q_AILOH < 0.9 | q_AILOH > 1) {
	  stop("ALOHA error message: The inputted upper bound of reference is incorrect (the value should range from 0.9 to 1).")
	}
}
a_biplot = as.numeric(tclvalue(biplot.Alphascaling.val))
if(!a_biplot %in% c(0, 1)) stop("ALOHA error message: The inputted alpha scaling is incorrect (the value should be either 0 or 1).")
main_function = paste(pathway_program, "/Main_function.r", sep="")
source(main_function)
