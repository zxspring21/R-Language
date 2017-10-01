
require(tcltk)
#addTclPath("C:/Tcl/lib") 
tclRequire("BWidget")
t1<-tktoplevel()

tktitle(t1)<-"Real data analysis"


frametitle<-tkframe(t1)
frameintro<-tkframe(frametitle,relief="ridge",borderwidth=5)
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontbody <- tkfont.create(size=11)
fontbody2 <-tkfont.create(size=3)

tkgrid(tklabel(frameintro,text="Welcome to use KBAT",font=fontHeading),sticky="w")
tkgrid(tklabel(frameintro,text="",font=fontbody2))
tkgrid(tklabel(frameintro,text="    KBAT (Kernel-based association test) is a convenient analysis tool for disease gene association mapping. Several",font=fontbody),sticky="w")
tkgrid(tklabel(frameintro,text="powerful association tests in KBAT are developed based on the concept of p-value combination and sliding window.",font=fontbody),sticky="w")
tkgrid(tklabel(frameintro,text="The methods provide systematic genome-wide or candidate-region searches for disease susceptibility genes.",font=fontbody),sticky="w")
tkgrid(tklabel(frameintro,text="Numerical/graphic results are outputted together to provide insight into the disease-marker association in study regions.",font=fontbody),sticky="w")
tkgrid(tklabel(frameintro,text="",font=fontbody2))

tkgrid(tklabel(frameintro,text="Reference: ",font=fontbody),sticky="w")
tkgrid(tklabel(frameintro,text="Hsin-Chou Yang, Hsin-Yi Hsieh & Cathy SJ Fann. (2008) KBAT: Kernel-based association test. Genetics in press.",font=fontbody),sticky="w")



tkgrid(frameintro)
tkgrid(frametitle)
tkgrid(tklabel(t1,text=""))


framebody<-tkframe(t1)
framepara<-tkframe(framebody)



Path1<-tclVar("C:\\\\KBAT\\\\Real\\\\Input")
entry.Path<-tkentry(framepara,width="23",textvariable=Path1,xscrollcommand=function(...)
												   {
												    path <- paste(tclvalue(Path1),"\\map.txt",sep="")
												    if(tclvalue(Path1)!="" & file.exists(path))
												    {
												     x <- tclVar(nrow(read.table(path)))
												     tkconfigure(entry.Num,textvariable=x)
												     tkconfigure(entry.End,textvariable=x)
												    }
												   })
pathlab<-tklabel(framepara,text="Directory of data input:")
pathbut<-tkbutton(framepara,text="...",
			 command=function() { tclvalue(Path1) <- tkchooseDirectory(initialdir=tclvalue(Path1))
						    path <- paste(tclvalue(Path1),"\\map.txt",sep="")
						    if(tclvalue(Path1)!="" & file.exists(path))
						    {
						     x <- tclVar(nrow(read.table(path)))
						     tkconfigure(entry.Num,textvariable=x)
						     tkconfigure(entry.End,textvariable=x)
						    }
						  })
tkgrid(pathlab,entry.Path,pathbut)
tkgrid.configure(pathlab,sticky="w")
tkgrid.configure(entry.Path,sticky="w")
tkgrid.configure(pathbut,sticky="w")


Pathout<-tclVar("C:\\\\KBAT\\\\Real\\\\Output")
entry.Pathout<-tkentry(framepara,width="23",textvariable=Pathout)
pathoutlab<-tklabel(framepara,text="Directory of results output:")
pathoutbut<-tkbutton(framepara,text="...",
			     command=function() tclvalue(Pathout) <- tkchooseDirectory())
tkgrid(pathoutlab,entry.Pathout,pathoutbut)
tkgrid.configure(pathoutlab,sticky="w")
tkgrid.configure(entry.Pathout,sticky="w")
tkgrid.configure(pathoutbut,sticky="w")


Num1<-tclVar("123")
entry.Num<-tkentry(framepara,width="23",textvariable=Num1)
numlab<-tklabel(framepara,text="Total number of SNPs:")
tkgrid(numlab,entry.Num)
tkgrid.configure(numlab,sticky="w")
tkgrid.configure(entry.Num,sticky="w")


St1<-tclVar("1")
entry.St<-tkentry(framepara,width="23",textvariable=St1)
stlab<-tklabel(framepara,text="The first marker of study region:")
tkgrid(stlab,entry.St)
tkgrid.configure(stlab,sticky="w")
tkgrid.configure(entry.St,sticky="w")


End1<-tclVar("123")
entry.End<-tkentry(framepara,width="23",textvariable=End1)
endlab<-tklabel(framepara,text="The last marker of study region:")
tkgrid(endlab,entry.End)
tkgrid.configure(endlab,sticky="w")
tkgrid.configure(entry.End,sticky="w")


welab<-tklabel(framepara,text="Weighting procedure:")
weightlab <- c("Distance","LD","LD and/or distance")
comboBox1 <- tkwidget(framepara,"ComboBox",editable=FALSE,values=weightlab)
tkgrid(welab,comboBox1)
tkgrid.configure(welab,sticky="w")


ldlab<-tklabel(framepara,text="Data format of LD information:")
ldtype <- c("Not available","LD measure","Genotype data")
comboBox2 <- tkwidget(framepara,"ComboBox",editable=FALSE,values=ldtype)
tkgrid(ldlab,comboBox2)
tkgrid.configure(ldlab,sticky="w")


adawslab<-tklabel(framepara,text="Determination of bandwidth/window size:")
adawstype<- c("Bandwidth","Window")
comboBox3 <- tkwidget(framepara,"ComboBox",editable=FALSE,values=adawstype)
tkgrid(adawslab,comboBox3)
tkgrid.configure(adawslab,sticky="w")


adblab<-tklabel(framepara,text="Bandwidth or m (window size=2m+1):   ")
band1<-tclVar("")
entry.band<-tkentry(framepara,width="23",textvariable=band1)
fixeg<-tklabel(t1,text="(e.g., 1, 3, 5)")
tkgrid(adblab,entry.band,fixeg)
tkgrid.configure(adblab,sticky="w")
tkgrid.configure(entry.band,sticky="w")
tkgrid.configure(fixeg,sticky="w")

##Truncation Threshold##
tuo1<-tclVar("1")
entry.tuo<-tkentry(framepara,width="23",textvariable=tuo1)
tuolab<-tklabel(framepara,text="Truncation threshold (Theta):")
tuoeg<-tklabel(framepara,text="(e.g., 0.05, 0.1, 1)")
tkgrid(tuolab,entry.tuo,tuoeg)
tkgrid.configure(tuolab,sticky="w")
tkgrid.configure(entry.tuo,sticky="w")
tkgrid.configure(tuoeg,sticky="w")


scr <- tkscrollbar(t1, repeatinterval=6, command=function(...)tkyview(t2,...))
t2<-tklistbox(framepara,height=4,width="23",selectmode="multiple",yscrollcommand=function(...)tkset(scr,...),background="white")
statlab<-tklabel(framepara,text="Statistic:")
tkgrid(statlab,t2,scr)
tkgrid.configure(statlab,sticky="w")
tkgrid.configure(scr,rowspan=4,sticky="W")
fruits <-c("SLM","MPM","PPM","WPPM-PD","WPPM-LD","WPPM-PDLD","KBAT-PD","KBAT-PDLD")
for (i in (1:8))
{
    tkinsert(t2,"end",fruits[i])
}
tkselection.set(t2,0)


B1<-tclVar("1000")
entry.B<-tkentry(framepara,width="23",textvariable=B1)
blab<-tklabel(framepara,text="Number of Monte Carlo replications:")
Beg<-tklabel(framepara,text="(Between 10 and 10,000)")
tkgrid(blab,entry.B,Beg)
tkgrid.configure(Beg,sticky="w")
tkgrid.configure(blab,sticky="w")
tkgrid.configure(entry.B,sticky="w")


xlabe1<-tclVar("Position")
entry.xlabe<-tkentry(framepara,width="23",textvariable=xlabe1)
xllab<-tklabel(framepara,text="Label of the horizontal axis:")
tkgrid(xllab,entry.xlabe)
tkgrid.configure(xllab,sticky="w")
tkgrid.configure(entry.xlabe,sticky="w")


tkgrid(framepara)
tkgrid(framebody)
tkpack(framepara,side="right",anchor="w")
tkgrid(tklabel(t1,text=""))


KBAT.path<-gsub("\\KBAT_R.R|/KBAT_R.R|\\KBAT_R.R|/KBAT_R.R","",KBAT)
sae<-paste(KBAT.path,"\\Program\\Save.R",sep="")
OK.but <- tkbutton(t1,text="  Run  ",command=function(){
 fruitChoice<-fruits[as.numeric(tkcurselection(t2))+1]
cat("Please wait a while, KBAT is checking the input information...\n")
 source(sae)}
)
oklab<-tklabel(t1,text="")
tkgrid(OK.but)


tkfocus(t1)



