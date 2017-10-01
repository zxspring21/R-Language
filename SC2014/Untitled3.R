library(xlsx)
data <- read.xlsx("/Users/portran/Desktop/Academia/PATIENT1000_0703.xlsx",1) 
#classes<-sapply(data,class)
data.age<-data[,c(3)]
data.diagnoise <- data[ ,c(6)]
x<- c(young=0, adult=0,old=0,l=0,ll=0)
y<- c(a=0,b=0,c=0,d=0,e=0)
for(j in 1:1000){
  if(data.diagnoise[i]==703){
    y[1]=y[1]+1
  }else if(data.diagnoise[i]==7051){
    y[2]=y[2]+1
  }else if((data.diagnoise[i]==5714)){
    y[3]=y[3]+1
  }else if((data.diagnoise[i]==702)){
    y[4]=y[4]+1
  }else{ #7041
    y[5]=y[5]+1
  }
    
}
for(i in 1:1000){
  if(data.age[i]<19){
    x[1]=x[1]+1
  }else if (data.age[i]<=64){
    x[2]=x[2]+1
  }else{
    x[3]=x[3]+1
  }
}
hist(x,main="people and diagonise", pro=T)
#plot(x,y,main = "people and diagonise")
#classes
#head(data,1)
#head(data)
#d<-scan()

#summary(mydata)
#sapply(mydata,sd)
#xtabs(BIRTH,data=mydata)

