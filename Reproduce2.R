library(flexwinn)
library(ggplot2)
library(ggpubr)
library(tictoc)
## data generation
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
    n<-sample(200:250,1)
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

####compare perfect win and flexwinn
corch<-numeric(5)
corwinn<-numeric(5)
for (i in 1:5){
data<-datagen(10)
met1<-data[[1]]
changepoints<-data[[3]]
result1<-winn(as.data.frame(met1),changepoints,graph=FALSE)
result2<-flexwinn(as.data.frame(met1),graph=FALSE)
corwinn[i]<-cor(result1[[1]],data[[2]])
corch[i]<-cor(result2[[1]],data[[2]])

}

mean(corch)
mean(corwinn)
##########
####
#comparing different accuracy in plates winn vs flexwinn
accuracyplates<-function(ndata=10,nrep=10){
  accuracy<-seq(0,1,by=0.1)
  sizeacc<-length(accuracy)
  msen<-numeric(nrep)
  msewinn<-numeric(nrep)
  corn<-numeric(nrep)
  corwinn<-numeric(nrep)
  mmsewinn<-numeric(sizeacc)
  sdmsewinn<-numeric(sizeacc)
  mmsen<-numeric(sizeacc)
  sdmsen<-numeric(sizeacc)
  mcorwinn<-numeric(sizeacc)
  sdcorwinn<-numeric(sizeacc)
  mcorn<-numeric(sizeacc)
  sdcorn<-numeric(sizeacc)
  passwn<-numeric(nrep)
  passwn.ch<-numeric(nrep)
  m.passwn<-numeric(sizeacc)
  m.passwn.ch<-numeric(sizeacc)
  for (acc in 1:length(accuracy)){
    set.seed(423)
    #set.seed(333)
    #set.seed(123)

    for (h in 1:nrep){
      data<-datagen(ndata)
      true<-data[[3]]
      worsech<-true
      for (pr in 1:length(worsech)){
        good<-rbinom(1,1,p=accuracy[acc])
        if (accuracy[acc]==1|accuracy[acc]==0){
          rbinom(1,1,p=0.5)
        }
        if (good==0){
          worsech[pr]<-worsech[pr]+sample(c(-10:-1,1:10),1,replace = TRUE)
        }
        if (good==1){
          sample(c(-10:-1,1:10),1,replace = TRUE)
        }

      }
      met1<-data[[1]]
      met1<-as.data.frame(met1)
      resultwinn<-winn(met1,worsech)
      corwinn[h]<-cor(resultwinn[[1]],data[[2]],method = "spearman")
      msewinn[h]<-mean((resultwinn[[1]]-data[[2]])^2)



      print(paste(h,acc,sep="_"))

      resultwinn.ch<-flexwinn(met1)

      normalized<-resultwinn.ch[[1]]
      #passwn.ch[h]<-resultwinn.ch$summary.transf["pass.wn.2"][[1]]
      #passwn.ch[h]<-normalized[[2]]


      corn[h]<-cor(normalized,data[[2]],method = "spearman")

      msen[h]<-mean((normalized-data[[2]])^2)


    }
    mmsewinn[acc]<-mean(msewinn)
    sdmsewinn[acc]<-sd(msewinn)
    mmsen[acc]<-mean(msen)
    sdmsen[acc]<-sd(msen)
    mcorwinn[acc]<-mean(corwinn)
    sdcorwinn[acc]<-sd(corwinn)
    mcorn[acc]<-mean(corn)
    sdcorn[acc]<-sd(corn)
  }


  meanmse<-c(mmsewinn,mmsen)
  sdmse<-c(sdmsewinn,sdmsen)
  table1<-cbind.data.frame(mse=meanmse,sdmse,accuracy=rep(accuracy,2),type=c(rep("winn",sizeacc),rep("flexwinn",sizeacc)))
  meancor<-c(mcorwinn,mcorn)
  sdcor<-c(sdcorwinn,sdcorn)
  table2<-cbind.data.frame(correlation=meancor,sdcor,accuracy=rep(accuracy,2),type=c(rep("winn",sizeacc),rep("flexwinn",sizeacc)))


  plotmse<-ggplot(data=table1,aes(x=accuracy,y=mse,colour=type)) +
    geom_point()+ggtitle(paste("mse with untransformed data /unaccurate ", "plates",sep=""))+
    geom_errorbar(aes(ymin=mse-sdmse/sqrt(nrep)*qt(0.975,nrep), ymax=mse+sdmse/sqrt(nrep)*qt(0.975,nrep)))
  plotcor<-ggplot(data=table2,aes(x=accuracy,y=correlation,colour=type)) +
    geom_point()+ggtitle(paste("spearman correlation with untransformed data /unaccurate ","plates",sep=""))+
    geom_errorbar(aes(ymin=correlation-sdcor/sqrt(nrep)*qt(0.975,nrep), ymax=correlation+sdcor/sqrt(nrep)*qt(0.975,nrep)))

  ggarrange(plotmse,plotcor,nrow=2)

}

set.seed(423)
accuracyplates(n=10,nrep=100)
#########
#Compare methods
#compare winn and flexwinn in different ways
#method 1: inaccuracy
#method 2: too many 
#method 3: missing
#method 0: perfect
comparemethods<-function(n1=5,n2=10,nrep=10,method=1){
  name<-c("perfect","unaccurate","too many","change in")
  n1<-n1
  n2<-n2
  msen<-numeric(nrep)
  msewinn<-numeric(nrep)
  corn<-numeric(nrep)
  corwinn<-numeric(nrep)
  mmsewinn<-numeric(n2-n1+1)
  sdmsewinn<-numeric(n2-n1+1)
  mmsen<-numeric(n2-n1+1)
  sdmsen<-numeric(n2-n1+1)
  mcorwinn<-numeric(n2-n1+1)
  sdcorwinn<-numeric(n2-n1+1)
  mcorn<-numeric(n2-n1+1)
  sdcorn<-numeric(n2-n1+1)
  passwn<-numeric(nrep)
  passwn.ch<-numeric(nrep)
  m.passwn<-numeric(n2-n1+1)
  m.passwn.ch<-numeric(n2-n1+1)
  for (ndata in n1:n2){

    for (h in 1:nrep){
      data<-datagen(ndata)
      true<-data[[3]]
      worsech<-true
      if (method==1){
        for (pr in 1:length(worsech)){
          good<-rbinom(1,1,p=0.5)
          if (good==0){
            worsech[pr]<-worsech[pr]+sample(c(-10:-1,1:10),1,replace = TRUE)

          }

        }
      }
      if (method==2){
        numberbad2<-sample(1:floor(length(true)/2),1)
        bad2<-1:(length(true)-1)#sample(1:(length(true)-1),numberbad2)
        for (j in 1:length(bad2)){
          add<-floor((worsech[bad2[j]]+worsech[bad2[j]+1])/2)
          worsech<-c(worsech,add)
        }
      }
      if (method==3){
        numberbad3<-sum(rbinom(length(true),1,0.2))
        bad3<-sample(1:length(true),numberbad3)
        worsech<-worsech[-bad3]

      }
      worsech<-sort(worsech)
      met1<-data[[1]]
      met1<-as.data.frame(met1)
      resultwinn<-winn(met1,worsech,graph = TRUE)
      corrected.df<-resultwinn[[1]]
      corwinn[h]<-cor(corrected.df,data[[2]],method = "spearman")
      msewinn[h]<-mean((corrected.df-data[[2]])^2)




      print(paste(ndata,h,sep="_"))
      resultwinn.ch<-flexwinn(met1,graph = TRUE)


      normalized<-resultwinn.ch[[1]]



      corn[h]<-cor(normalized,data[[2]],method = "spearman")

      msen[h]<-mean((normalized-data[[2]])^2)


    }

    mmsewinn[ndata-(n1-1)]<-mean(msewinn)
    sdmsewinn[ndata-(n1-1)]<-sd(msewinn)
    mmsen[ndata-(n1-1)]<-mean(msen)
    sdmsen[ndata-(n1-1)]<-sd(msen)
    mcorwinn[ndata-(n1-1)]<-mean(corwinn)
    sdcorwinn[ndata-(n1-1)]<-sd(corwinn)
    mcorn[ndata-(n1-1)]<-mean(corn)
    sdcorn[ndata-(n1-1)]<-sd(corn)
  }


  number<-(n1:n2)
  meanmse<-c(mmsewinn,mmsen)
  sdmse<-c(sdmsewinn,sdmsen)
  table1<-cbind.data.frame(mse=meanmse,sdmse,number_chpoints=rep(number,2),type=c(rep("winn",n2-n1+1),rep("flexwinn",n2-n1+1)))
  meancor<-c(mcorwinn,mcorn)
  sdcor<-c(sdcorwinn,sdcorn)
  table2<-cbind.data.frame(correlation=meancor,sdcor,number_chpoints=rep(number,2),type=c(rep("winn",n2-n1+1),rep("flexwinn",n2-n1+1)))


  plotmse<-ggplot(data=table1,aes(x=number_chpoints,y=mse,colour=type)) +
    geom_point()+ggtitle(paste("mse with untransformed data /",name[method+1]," plates",sep=""))+
    geom_errorbar(aes(ymin=mse-sdmse/sqrt(nrep)*qt(0.975,nrep-1), ymax=mse+sdmse/sqrt(nrep)*qt(0.975,nrep-1)))
  plotcor<-ggplot(data=table2,aes(x=number_chpoints,y=correlation,colour=type)) +
    geom_point()+ggtitle(paste("spearman correlation with untransformed data /",name[method+1]," plates",sep=""))+
    geom_errorbar(aes(ymin=correlation-sdcor/sqrt(nrep)*qt(0.975,nrep-1), ymax=correlation+sdcor/sqrt(nrep)*qt(0.975,nrep-1)))

  ggarrange(plotmse,plotcor,nrow=2)

}

#set.seed(913)
set.seed(135) #too many nrep 100
#set.seed(455)#unaccurate
#set.seed(124)#nrep100
#set.seed(366)#missing
comparemethods(n1=10,n2=10,nrep=1,method=2)

######
#Missing plates
#compare different probabbilities in missing plates
missingplates<-function(ndata=10,nrep=10){
  missing<-seq(0,0.5,by=0.1)
  sizemiss<-length(missing)
  msen<-numeric(nrep)
  msewinn<-numeric(nrep)
  corn<-numeric(nrep)
  corwinn<-numeric(nrep)
  mmsewinn<-numeric(sizemiss)
  sdmsewinn<-numeric(sizemiss)
  mmsen<-numeric(sizemiss)
  sdmsen<-numeric(sizemiss)
  mcorwinn<-numeric(sizemiss)
  sdcorwinn<-numeric(sizemiss)
  mcorn<-numeric(sizemiss)
  sdcorn<-numeric(sizemiss)
  passwn<-numeric(nrep)
  passwn.ch<-numeric(nrep)
  m.passwn<-numeric(sizemiss)
  m.passwn.ch<-numeric(sizemiss)
  for (mis in 1:length(missing)){
    #set.seed(243)
    #set.seed(333)
    set.seed(123)

    for (h in 1:nrep){
      data<-datagen(ndata)
      true<-data[[3]]
      worsech<-true
      miss.sample<-as.logical(rbinom(length(worsech),1,p=missing[mis]))
      if (missing[mis]==0){rbinom(length(worsech),1,p=0.5)}
      worsech<-worsech[!miss.sample]
      met1<-data[[1]]
      met1<-as.data.frame(met1)
      resultwinn<-winn(met1,worsech)
      corrected.df<-resultwinn[[1]]
      corwinn[h]<-cor(corrected.df,data[[2]])
      msewinn[h]<-mean((corrected.df-data[[2]])^2)




      print(paste(mis,h,sep="_"))
      resultwinn.ch<-flexwinn(met1)


      normalized<-resultwinn.ch[[1]]



      corn[h]<-cor(normalized,data[[2]])

      msen[h]<-mean((normalized-data[[2]])^2)


    }

    mmsewinn[mis]<-mean(msewinn)
    sdmsewinn[mis]<-sd(msewinn)
    mmsen[mis]<-mean(msen)
    sdmsen[mis]<-sd(msen)
    mcorwinn[mis]<-mean(corwinn)
    sdcorwinn[mis]<-sd(corwinn)
    mcorn[mis]<-mean(corn)
    sdcorn[mis]<-sd(corn)
  }


  meanmse<-c(mmsewinn,mmsen)
  sdmse<-c(sdmsewinn,sdmsen)
  table1<-cbind.data.frame(mse=meanmse,sdmse,missing=rep(missing,2),type=c(rep("winn",sizemiss),rep("flexwinn",sizemiss)))
  meancor<-c(mcorwinn,mcorn)
  sdcor<-c(sdcorwinn,sdcorn)
  table2<-cbind.data.frame(correlation=meancor,sdcor,missing=rep(missing,2),type=c(rep("winn",sizemiss),rep("flexwinn",sizemiss)))


  plotmse<-ggplot(data=table1,aes(x=missing,y=mse,colour=type)) +
    geom_point()+ggtitle(paste("mse with untransformed data /missing ", "plates",sep=""))+
    geom_errorbar(aes(ymin=mse-sdmse/sqrt(nrep)*qt(0.975,nrep), ymax=mse+sdmse/sqrt(nrep)*qt(0.975,nrep)))
  plotcor<-ggplot(data=table2,aes(x=missing,y=correlation,colour=type)) +
    geom_point()+ggtitle(paste("spearman correlation with untransformed data /missing ","plates",sep=""))+
    geom_errorbar(aes(ymin=correlation-sdcor/sqrt(nrep)*qt(0.975,nrep), ymax=correlation+sdcor/sqrt(nrep)*qt(0.975,nrep)))

  ggarrange(plotmse,plotcor,nrow=2)

}
set.seed(123)
missingplates(ndata = 10,nrep = 100)


