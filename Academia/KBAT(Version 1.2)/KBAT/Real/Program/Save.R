

a<-matrix(0,14,1)
PATH1<-tkget(entry.Path)
PATH1<-tclvalue(PATH1)
a[1,]<-paste("PATH1<-c('",PATH1,"')",sep="")

PATHout<-tkget(entry.Pathout)
PATHout<-tclvalue(PATHout)
a[2,]<-paste("PATHout<-c('",PATHout,"')",sep="")

Num<-tkget(entry.Num)
Num<-tclvalue(Num)
Num<-as.numeric(Num)
a[3,]<-paste("Num<-c(Num)",sep="")

St<-tkget(entry.St)
St<-tclvalue(St)
St<-as.numeric(St)
a[4,]<-paste("St<-c(St)",sep="")

End<-tkget(entry.End)
End<-tclvalue(End)
End<-as.numeric(End)
a[5,]<-paste("End<-c(End)",sep="")

We1<-as.character(tclvalue(tkget(comboBox1)))
{if(We1=="Distance"){
Weight<-1
}
else{
{if(We1=="LD"){
Weight<-2
}else{Weight<-3}}
}}
a[6,]<-paste("Weight<-c(Weight)",sep="")

LD1<-as.character(tclvalue(tkget(comboBox2)))
if(LD1=="LD measure") LDdata<-1 else if(LD1=="Genotype data") LDdata<-2 else LDdata<-0
a[7,]<-paste("LDdata<-c(LDdata)",sep="")

borw<-as.character(tclvalue(tkget(comboBox3)))
borw2<-ifelse(borw=="Window",1,2)
a[8,]<-paste("borw2<-c(",borw2,")",sep="")

ws<-tkget(entry.band)
op<-tclvalue(ws)
a[9,]<-paste("op<-c(",op,")",sep="")

tuo<-tkget(entry.tuo)
Tuo.size<-tclvalue(tuo)
a[10,]<-paste("Tuo.size<-c(",Tuo.size,")",sep="")

B<-tkget(entry.B)
B<-tclvalue(B)
B<-as.numeric(B)
a[11,]<-paste("B<-c(B)",sep="")

Xlabe<-tkget(entry.xlabe)
Xlabe<-tclvalue(Xlabe)
a[12,]<-paste("Xlabe<-c('",Xlabe,"')",sep="")

fruitChoice<-fruits[as.numeric(tkcurselection(t2))+1]
Stat<-as.matrix(fruitChoice)
a[13,]<-"Stat<-Stat"

a[14,]<-"bonf<-0.05/Num"

write.table(a,paste(KBAT.path,"\\Program\\","Inf.R",sep=""),quote=F,row.names=F,col.names=F)
infpath<-paste(KBAT.path,"\\Program\\","Inf.R",sep="")
source(infpath)



cerr<-0
if(St>Num|St<=0){
 msg <- paste("The first marker of data should be between 1 and ", Num,".",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

if(End>Num|End<=0){
 msg <- paste("The end marker of data should be between 1 and ", Num,".",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

{if(borw2==1){
O1<-(op)>((Num-1)/2)
}else{O1<-c("FALSE")}}
O2<-(op)<=0

if(any(O1=="TRUE")|any(O2=="TRUE")){
 msg <- paste("The window size should be between 1 and ",((Num-1)/2),".",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

T1<-(Tuo.size)>1
T2<-(Tuo.size)<0
if(any(T1=="TRUE")|any(T2=="TRUE")){
 msg <- paste("Truncation threshold should be between 0 and 1.",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

if(Weight==2){
list.stat<-c("WPPM-PD","WPPM-PDLD","KBAT-PD","KBAT-PDLD")
check.stat<-is.element(Stat,list.stat)
if(any(check.stat)=="TRUE"){
 msg <- paste("Calculation of ",Stat," requares the map information.",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}}

if(Weight==1){
list.stat<-c("WPPM-LD","WPPM-PDLD","KBAT-PDLD")
check.stat<-is.element(Stat,list.stat)
if(any(check.stat)=="TRUE"){
 msg <- paste("Calculation of ",Stat," requares the LD information.",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}}

if(LDdata!=0 & Weight==1){
 msg <- paste("LD information will not be used although 'LD measure' is selected in 'Data format of LD information'.",sep="")
 tkmessageBox(title="KBAT warning message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

if(LDdata==0 & Weight!=1){
 msg <- paste("Data format of LD information (LD measure or Genotype data)should be selected.",sep="")
 tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

if(nrow(Stat)==0){
 msg <- paste("Please select at least one test statistic.",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

if(B>10000){
 msg <- paste("The number of Monte Carlo replications should be between 10 and 10,000.",sep="")
   tkmessageBox(title="KBAT error message",message=msg,icon="warning",type="ok")
cerr<-cerr+1
}

{if(cerr==0){
read<-paste(KBAT.path,"\\Program\\ReadFile.R",sep="") 
source(read)}}