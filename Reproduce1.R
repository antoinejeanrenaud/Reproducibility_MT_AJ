########
# PACKAGES to use
library(flexwinn)
library(tictoc)
library(sarima)
library(dplyr)
library(splines)
library(ggplot2)
library(ggpubr)

##########
#data generation

# Data generation with trend
#param: number of change point
datagen<-function(change=5){
  data<-numeric(0)
  datat<-numeric(0)
  datan<-numeric(0)
  v<-numeric(0)
  for (i in 1:(change+1)){
    mean<-runif(1,min=-2,max=2)
    max<-runif(1,1,2)
    trend<-sample(c(1,-1),1)
    sdfac<-runif(1,min=0.5,max=1)
    n<-sample(50:150,1)
    m<-sample(c(1,2,3,4),1)
    sdbeta<-sqrt((2*5)/((2+5)^2*(2+5+1)))
    if (m==1){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rnorm(n)#rt(n,10)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+trend*seq(0,max,length.out=n))
      datan<-c(datan,new+mean+trend*seq(0,max,length.out=n))

    }
    if (m==2){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+trend*seq(0,sqrt(max),length.out=n)^2)
      datan<-c(datan,new+mean+trend*seq(0,sqrt(max),length.out=n)^2)
    }
    if (m==3){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+trend*sqrt(seq(0,max,length.out=n)))
      datan<-c(datan,new+mean+trend*sqrt(seq(0,max,length.out=n)))

    }
    if (m==4){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean)
      datan<-c(datan,new+mean)

    }
    if (m==5){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+max*sin(seq(0,1,length.out=n)))
      datan<-c(datan,new+mean+max*sin(seq(0,10,length.out=n)))

    }
    v<-c(v,n)

  }
  v<-cumsum(v)
  v<-v[-length(v)]
  return(list(data,datat,v,datan))
}

######
#data  generation with no trend

datagen.no.trend<-function(change=5){
  data<-numeric(0)
  datat<-numeric(0)
  datan<-numeric(0)
  v<-numeric(0)
  for (i in 1:(change+1)){
    mean<-runif(1,min=-4,max=4)
    sdfac<-runif(1,min=0.6,max=1)
    n<-sample(50:150,1)

      new<-(rbeta(n,2,5)-2/7)/0.16#rt(n,10)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean)
      datan<-c(datan,new+mean)



      v<-c(v,n)
    }



  v<-cumsum(v)
  v<-v[-length(v)]
  return(list(data,datat,v,datan))
}


########
#Cost of segment calculator with fixed knots
#utility function
fksplinecost<-function(data,knots,index1=1,index2=length(data)){
  size<-length(data)
  if (size==1){
    return(0)
  }
  if (size<5){
    sd<-sd(data)
    mu<-mean(data)
    neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi)))
    if (sd==0){
      return(0)
    }else{
    return(neglog)
    }
  }
  cov<-index1:index2
  newknots<-knots[knots<index2& knots>index1]
  spline<-lm(data~ns(cov,knots = newknots,intercept=TRUE))
  # col<-ncol(ns(cov,knots = newknots,intercept=TRUE))
  # if (col>4) print(col)
  mu<-spline$fitted.values
  sd<-sd(data-mu)
  if (sd==0){
    return(0)
  }
  neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi)))
  if (size<15){
  return(neglog) }else{
    return(neglog)
    }

}

########
#Estimating changepoints in distribution with fixed knots

fkPELT<-function(data,knots){
n<-length(data)
f<-numeric(length = n+1)
f[1]<- -3*log(300)
cp<-rep(list(numeric()),n+1)
R<-rep(list(numeric()),n)
R[[1]]<-0


for (t in 1:n){
  m<-length(R[[t]])
  neglog<-numeric(m)
  newR<-R[[t]]
  for (r in 1:m){
    neglog[r]<-fksplinecost(data[(newR[r]+1):t],knots=knots,index1=newR[r]+1,index2=t)
  }
  newneglog<-neglog
  
  newm<-length(newR)
  stat<-numeric(newm)
  for (r in 1:newm){
    tot<-f[newR[r]+1]+3*log(300)+newneglog[r]
    stat[r]<-tot

  }
  f[t+1]<-min(stat)
  t1<-newR[which.min(stat)]
  cp[[t+1]]<-c(cp[[t1+1]],t1)
  for (r in 1:m){

    tot<-f[R[[t]][r]+1]+log(300)+neglog[r]
    if (tot<=f[t+1]& t<n){

      R[[t+1]]<-c(R[[t+1]],R[[t]][r])
    }
  }
  if (t<n&t>39){
    R[[t+1]]<-c(R[[t+1]],t-19)}
}

cp[[n+1]]<-cp[[n+1]][-1]
return(cp[[n+1]])
}




########
# Test if find correctly changepoints
test<-function(nchange=6){
x<-1:200
y<-1:150
#data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(146,mean=4,sd=1),rnorm(500,mean = sin(0.2*seq(1,125,by=0.25)),sd=0.5),rnorm(300,mean=seq(0.01,3,by=0.01)))
#data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(200,mean = 0.0003*x^2,sd=1),rnorm(300,mean=0.5,sd=1.1),rnorm(300,mean = -0.5,sd=1.5),rnorm(100,mean = 0.5,sd=0.67))
#data<-c(rnorm(450,mean = 0,sd=1),rnorm(150,mean=0.5,sd=1.1),rnorm(300,mean = 0,sd=1.5),rnorm(100,mean = 0.5,sd=1))
#data<-rnorm(1000)
tic()
data<-datagen(nchange)
par(mfrow=c(2,1))
untransformed<-ts(data[[2]])
distorted<-ts(data[[1]])
plot(untransformed,main="data untransformed")
plot(distorted,main="data transformed with estimated changepoints")
knots<-seq(1,length(data[[1]]),length.out=floor(length(data[[1]])/60)+2)
knots<-knots[-c(1,length(knots))]
knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
changepoints<-fkPELT(data[[1]],knots = knots)
for (b in data[[3]]){

  abline(v=b,col="blue")

}
correct<-0
for (a in changepoints){
  if (sum(a==data[[3]])>0){
    correct<-1+correct
  abline(v=a,col="red")}else{
    abline(v=a,col="green")
  }
}

print(changepoints)
print(correct/length(data[[3]]))
par(mfrow=c(1,1))
toc()
}
set.seed(167)
test(10)
####RED: TRUE POSITIVES
####BLUE: FALSE NEGATIVES
####GREEN: FALSE POSITIVES

##############
# TEST 2
test2<-function(){
  x<-1:200
  y<-1:150
  data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(146,mean=4,sd=1),rnorm(500,mean = sin(0.2*seq(1,125,by=0.25)),sd=0.5),rnorm(300,mean=seq(0.01,3,by=0.01)))
  #data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(200,mean = 0.0003*x^2,sd=1),rnorm(300,mean=0.5,sd=1.1),rnorm(300,mean = -0.5,sd=1.5),rnorm(100,mean = 0.5,sd=0.67))
  #data<-c(rnorm(450,mean = 0,sd=1),rnorm(150,mean=0.5,sd=1.1),rnorm(300,mean = 0,sd=1.5),rnorm(100,mean = 0.5,sd=1))
  #data<-rnorm(1000)
  tic()
  plot(ts(data),main="data")
  plot(ts(data),main="data with estimated changepoints")
  knots<-seq(1,length(data),length.out=floor(length(data)/60)+2)
  knots<-knots[-c(1,length(knots))]
  knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
  changepoints<-fkPELT(data,knots = knots)

  for (a in changepoints){

      abline(v=a,col="red")
  }

  print(changepoints)

  par(mfrow=c(1,1))
  toc()
}
test2()
############

############
# Test limitation for perfect accuracy
limitsd<-function(nrep=10){
  ratio<-seq(1,3,by=0.2)
  precision<-numeric(length = length(ratio))
  for (j in 1:length(ratio)) {


  preci<-vector(length = nrep)
  for (i in 1:nrep){
  preci[i]<-FALSE
  print(paste(i,j,sep="_"))
  dat<-c(rnorm(500),rnorm(500)+ratio[j])
  knots<-seq(1,length(dat),length.out=floor(length(dat)/60)+2)
  knots<-knots[-c(1,length(knots))]
  knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
  ch<-fkPELT(dat,knots = knots)
  if (sum(ch==500)>0){
    preci[i]<-TRUE
  }
  }
  precision[j]<-mean(preci)
  }

  plot(ratio,precision)
}

limitsd(nrep=100)

# Test limitation for approximate accuracy
limitsd2<-function(nrep=10){
  ratio<-seq(1,3,by=0.2)
  precision<-numeric(length = length(ratio))
  for (j in 1:length(ratio)) {


    preci<-vector(length = nrep)
    for (i in 1:nrep){
      print(paste(i,j,sep="_"))
      preci[i]<-FALSE
      dat<-c(rnorm(500),rnorm(500)+ratio[j])
      knots<-seq(1,length(dat),length.out=floor(length(dat)/60)+2)
      knots<-knots[-c(1,length(knots))]
      knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
      ch<-fkPELT(dat,knots = knots)
      if (sum(497<=ch & ch<=503)>0){
        preci[i]<-TRUE
      }
    }
    precision[j]<-mean(preci)
  }

  plot(ratio,precision)
}

limitsd2(nrep=100)
#######
#indvidual test for limit depending on add: change in mean
add<-2
dat<-c(rnorm(500),rnorm(500)+add)
knots<-seq(1,length(dat),length.out=floor(length(dat)/60)+2)
knots<-knots[-c(1,length(knots))]
knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
ch<-fkPELT(dat,knots = knots)
data<-ts(dat)
plot(data)
for (i in ch){
abline(v=ch,col="red")
}










