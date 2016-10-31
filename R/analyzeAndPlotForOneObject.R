#analysis 
#Function to check for whole numbers
is.whole <- function(x) 
  is.numeric(x) && floor(x)==x 

analyzeAndPlotForOneObject<-function(normalizeddata,dataname,projectname)
{
{
  datafile<-normalizeddata
  if(class(datafile) == "character")
  {
    getrawdata<-as.matrix((read.table(datafile,header=F,sep="\t",stringsAsFactors=F)))    
  }else if(class(datafile)=="data.frame")
  {
    getrawdata<-as.matrix(datafile)
  }else if(class(datafile)=="matrix")
  {
    getrawdata<-datafile  
  }
  currentjob<-projectname
  
  if(is.na(currentjob[2])){currentjob[2]=currentjob[1]}
  
  #getrawdata<-normalizeddata
  methodnames<-dataname
  getEDdata<-getrawdata[1,]
    
  filterED<-as.numeric(getEDdata[-which(getEDdata<1)])
  filterrawdata<-getrawdata[,-(1:(length(getEDdata)-length(filterED)))]
  colnames(filterrawdata)<-getrawdata[2,-(1:(length(getEDdata)-length(filterED)))]
  filterrawdata<-(as.matrix((filterrawdata[-(1:2),])))
  class(filterrawdata)<-"numeric"
  methodlist<-filterrawdata
  
  #getrawdata<-normalizeddata[[3]]
  #filterrawdata<-normalizeddata[[4]]
  
  #HKflag<-normalizeddata[[6]]
  
  pooledvarmem<-vector()
  Medianofsamples<-vector()
  datamadsamplesmem<-vector()
  datalabels<-vector()
  
  avgcvmem<-matrix(nrow=length(levels(as.factor(unlist(filterED)))),ncol=length(methodlist),byrow=T)
  avgmadmem<-matrix(nrow=length(levels(as.factor(unlist(filterED)))),ncol=length(methodlist),byrow=T)
  avgvarmem<-matrix(nrow=length(levels(as.factor(unlist(filterED)))),ncol=length(methodlist),byrow=T)
  #datamadpeptmem<-matrix(nrow=nrow(filterrawdata),ncol=length(methodlist),byrow=T)
  
  
  anpvalue<-anfdr<-kwpvalue<-kwfdr<-vector()
  
    datastore<-(methodlist)
    
    templabel<-methodnames
    
    datamadsamples<-NA
    datamadpeptides<-NA
    pooledvar<-NA
    dataskew<-NA
    datacv<-NA
    tempcv<-NA
    tempskew<-NA
    datakurt<-NA
    pvobj1<-0
    pvobj2<-0
    
    
    datalabels<-c(datalabels,templabel)
    #MAD rows
    x<-1
    z<-1
    y<-1
    flag<-1
    flag1<-1
    count<-0
    meti<-1
    madmem<-matrix(nrow=nrow(datastore),ncol=length(levels(as.factor(unlist(filterED)))),byrow=T)
    tempcv<-vector()
    varmem<-vector()
    tempvar<-vector()
    nonmissingmat<-vector()
    
    
    for(i in 1:length(filterED))
    {
      if(x!=filterED[i] || i==length(filterED))
      {
        y<-i-1
        if(i==length(filterED))
        {
          y<-i
        }
        if(flag==1)
        {
          count<-count+1
          madmem[,count]<-apply(datastore[,z:y],1,function(x) {mad(x,na.rm=T)})
          nonmissingmat<-(apply(datastore[,z:y],1,function(x) {((sum(!is.na(x))))}))-1
          tempvar<-nonmissingmat*apply(datastore[,z:y],1,function(x) {var(x,na.rm=TRUE)})
        }
        if(flag==2)
        {  
          count<-count+1
          madmem[,count]<-apply(datastore[,z:y],1,function(x) {mad(x,na.rm=T)})  
          nonmissingmat<-(apply(datastore[,z:y],1,function(x) {((sum(!is.na(x))))}))-1
          tempvar<-nonmissingmat*apply(datastore[,z:y],1,function(x) {var(x,na.rm=TRUE)})
        }
        varmem<-c(varmem,((sum(tempvar,na.rm=T))/(sum(nonmissingmat,na.rm=T))))
        z<-i
        x<-filterED[i]
        flag=2;
      }
    }
    avgvarmem<-varmem
    temmadmatsum<-apply(madmem,2,mean,na.rm=T)
    avgmadmem<-temmadmatsum
    
    tempcvmat<-matrix(nrow=nrow(datastore),ncol=length(levels(as.factor(unlist(filterED)))),byrow=T)
    
    for(i in 1:nrow(datastore))
    {
      
      tempcv<-numSummary(datastore[i,],statistics=c("cv"),groups=unlist(filterED))
      tempcvmat[i,]<-tempcv$table
    }
    
    temcvmatsum<-apply(tempcvmat,2,mean,na.rm=T)
    avgcvmem<-((temcvmatsum*100))
    
    
    
    #ANOVA
    an<-rowSums(is.na(datastore))
    datastoretmp<-datastore[an<(ncol(datastore)/2),]
    anpvalue<-cbind(anpvalue,apply(datastoretmp,1,function(x) summary(aov(unlist(x)~filterED))[[1]][[5]][1]))
    anfdr<-cbind(anfdr,p.adjust(anpvalue[,meti],method="BH"))
    
    #Kruskal Wallis
    kwpvalue<-cbind(kwpvalue,apply(datastoretmp,1,function(x) kruskal.test(unlist(x)~filterED,na.action="na.exclude")[[3]][1]))
    kwfdr<-cbind(kwfdr,p.adjust(kwpvalue[,meti],method="BH")) 
  
  
  
}
print("Finished analysis, preparing plots and report....")

#plotting
{
  palette(c("red","green","blue","orange","darkgray","blueviolet","darkslateblue","darkviolet","gray","bisque4","brown","cadetblue4","darkgreen","darkcyan","darkmagenta","darkgoldenrod4","coral1"))
  jobdir<-paste(getwd(),sep="")
  pdf(file=paste(jobdir,"/Norm_report-",currentjob[1],".pdf",sep=""),paper="a4r",width=0,height=0)    
  
  def.par<-par(no.readonly=T)
  par(mfrow=c(4,1))
  textplot(paste("Project Name: ",currentjob[2],sep=""),halign="center",valign="center",cex=2.5)
  plot(0:10,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
  textplot("Normalization report for one datafile by Normalyzer",halign="center",valign="center",cex=1.5)
  textplot(paste("Report created on: ",Sys.Date(),sep=""),halign="center",valign="center",cex=1.5)
  pageno=2
  
  
  
  #CV
{
  tout<-rbind(c(1,2,3),c(4))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(2,2,3,2),xpd=NA)
  boxplot(avgcvmem,main="PCV - Intragroup",names=c(methodnames),border="red",density=20,cex=0.3,cex.axis=0.9,las=2,frame.plot=F)
  stripchart(as.data.frame(avgcvmem),vertical=T,cex=0.4,las=2,pch=20,add=T,col="darkgray")
  boxplot(avgmadmem,main="PMAD - Intragroup",names=c(methodnames),border="red",density=20,cex=0.3,cex.axis=0.9,las=2,frame.plot=F)
  stripchart(as.data.frame(avgmadmem),vertical=T,cex=0.4,las=2,pch=20,add=T,col="darkgray")
  boxplot(avgvarmem,main="PEV - Intragroup",names=c(methodnames),border="red",density=20,cex=0.3,cex.axis=0.9,las=2,frame.plot=F)
  stripchart(as.data.frame(avgvarmem),vertical=T,cex=0.4,las=2,pch=20,add=T,col="darkgray")
  mtext("Data Summary",side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
  pageno=pageno+1
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0)    
}
  
 
  #MA plots
{
  tout<-rbind(c(1,2,3,4,5),c(6,7,8,9,10),c(11,12,13,14,15))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(3,2,3,2),xpd=NA)
  pc=0
  
    datastore<-as.matrix(methodlist)
    tempcolname<-colnames(datastore)
    #datastore<-datastore[,1:sum(filterED==1)]
    datastore1<-datastore[!is.na(datastore[,1]),]
  
    for(i in 1:ncol(datastore1))
    {
      if(is.whole(i/15))
      {
        mtext("MA plots",side=3,adj=0,outer=TRUE,col="gray",line=1)
        mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
        mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0) 
        mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
        pageno=pageno+1
      }
      limma::plotMA(datastore1,array=i,xlab="",ylab="",main=paste(tempcolname[i],methodnames))
    }
  mtext("MA plots",side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0) 
  #mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
   
  
}
  #Scatterplots
{
  tout<-rbind(c(1,2,3,4),c(5,6,7,8),c(9,10,11,12))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(3,2,3,2),xpd=NA)
  pc=0
 
  pairs(datastore,cex=0.01)
  
  mtext("Scatterplots",side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0)  
  mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
  pageno=pageno+1  
}
  
  #qqplot
{
  tout<-rbind(c(1,2,3,4,5),c(6,7,8,9,10),c(11,12,13,14,15))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(3,2,3,2),xpd=NA)
  pc=0
   
    datastore<-(methodlist)
    tempcolname<-colnames(datastore)
  for(i in 1:ncol(datastore1))
  {
    if(is.whole(i/15))
    {
      mtext("Q-Q plots",side=3,adj=0,outer=TRUE,col="gray",line=1)
      mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
      mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0) 
      mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
      pageno=pageno+1
    }
    qqnorm(datastore[,i],main=paste(tempcolname[i],methodnames),xlab="",ylab="")
    #pc=pc+1
  }
  
  mtext("Q-Q plots",side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0) 
  #mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
}
  
  #boxplot
{
  tout<-rbind(c(1,2),c(3,4),c(5,6))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(3,2,3,2),xpd=NA)
  mindata<-1000
  maxdata<-0
  
    tempmin<-min(methodlist,na.rm=T)
    tempmax<-max(methodlist,na.rm=T)
    if(tempmin<mindata)
    {
      mindata<-tempmin
    }
    if(tempmax>maxdata)
    {
      maxdata<-tempmax
    }
    
    par(mar=c(5,1,1,1))
    boxplot(methodlist,cex=0.1,cex.axis=0.7,las=2,names=colnames(methodlist),main="boxplot",col=(filterED),outcol="lightgray",ylim=c((mindata-1),(maxdata+1)))
  
  mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
  pageno=pageno+1    
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0)  
}
   #RLEplots
   
    deviations = methodlist - rowMedians(methodlist,na.rm=T)
    boxplot(deviations,outcol="lightgray",cex=0.1,cex.axis=0.7,las=2,names=colnames(methodlist),main="RLE Plot",col=(filterED))
 
  #density plots
{
   
    datastore<-(methodlist)
    tempd<-density(datastore[,1],na.rm=T)
    plot(density(datastore[,1],na.rm=T),xlab="",ylab="",ylim=c(min(tempd$y),max(tempd$y)*1.5),main="density plot",lty=2,lwd=1,col="darkgray")
    for(j in 2:ncol(datastore))
    {
      lines(density(datastore[,j],na.rm=T),,lty=2,lwd=1,col="darkgray")
    }
}

  #mds PLOT
     
    datastore<-(methodlist)
    d<-dist(scale(t(na.omit(datastore)),center=TRUE,scale=TRUE))
    fit<-cmdscale(d,eig=TRUE,k=2)
    x<-fit$points[,1]
    y<-fit$points[,2]
    plot(x,y,type="n", main="mds Plot",xlab="",ylab="")
    text((fit$points[,1]),(fit$points[,2]),col=filterED,labels=filterED)
    
  #meanSDplot
   
    datastore<-(methodlist)
    meanSdPlot(datastore,xlab="",ylab="",main="meanSD Plot")
    
  #check correlation
{
  par(mfrow=c(1,1))
  tempcorsum<-vector()
  tempcorsumspear<-vector()
  corsumspear<-vector()
  avgcorsum<-matrix(nrow=((length(levels(as.factor(filterED)))*2))-1,ncol=length(methodlist),byrow=T)
  avgcorsumspear<-matrix(nrow=((length(levels(as.factor(filterED)))*2))-1,ncol=length(methodlist),byrow=T)
  corsum<-vector()
  
    flag1<-1
    datastore<-as.matrix(methodlist)
    for(j in 1:length(filterED))
    {
      if(j!= length(filterED) && flag1==filterED[j] && flag1==filterED[j+1])
      {
        z<-j+1
        test<-cor(as.matrix(datastore[,j:z]),use="complete.obs",method="pearson")
        tempcorsum<-c(tempcorsum,test[2])
        test<-cor(as.matrix(datastore[,j:z]),use="complete.obs",method="spearman")
        tempcorsumspear<-c(tempcorsumspear,test[2])
        
      }
      else
      {
        corsum<-c(corsum,mean(tempcorsum))
        tempcorsum<-vector()
        corsumspear<-c(corsumspear,mean(tempcorsumspear))
        tempcorsumspear<-vector()
        if(j!=length(filterED))
        {
          z<-j+1
          flag1<-filterED[j]
          test<-cor(as.matrix(datastore[,j:z]),use="complete.obs",method="pearson")
          tempcorsum<-c(tempcorsum,test[2])
          test<-cor(as.matrix(datastore[,j:z]),use="complete.obs",method="spearman")
          tempcorsumspear<-c(tempcorsumspear,test[2])
        }
      }
    }
    avgcorsum<-corsum
    corsum<-vector()
    avgcorsumspear<-corsumspear
    corsumspear<-vector()
  
  tout<-rbind(c(1,2),c(3))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(2,2,3,2),xpd=NA)
  abc<-boxplot(avgcorsum,main="Pearson correlation - Intragroup",names=c(methodnames),border="red",density=20,cex=0.3,cex.axis=0.9,las=2)
  stripchart(as.data.frame(avgcorsum),vertical=T,cex=0.4,las=2,pch=20,add=T,col="darkgreen")
  #grid(nx=NULL,ny=10,col="lightgray")
  abc<-boxplot(avgcorsumspear,main="Spearman correlation - Intragroup",names=c(methodnames),border="red",density=20,cex=0.3,cex.axis=0.9,las=2)
  stripchart(as.data.frame(avgcorsumspear),vertical=T,cex=0.4,las=2,pch=20,add=T,col="darkgreen")
  #grid(nx=NULL,ny=10,col="lightgray")
  mtext("Correlation plots",side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
  pageno=pageno+1    
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0)    
}
  
  #dendrograms
{
  tout<-rbind(c(1,2,3,4),c(5,6,7,8),c(9,10,11,12))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(2,2,3,2),xpd=NA)
  colt<-(c("red","green","blue","orange","darkgray","blueviolet","darkslateblue","darkviolet","gray","bisque4","brown","cadetblue4","darkgreen","darkcyan","darkmagenta","darkgoldenrod4","coral1"))
  
  
    temp<-scale(t(na.omit(methodlist)),center=TRUE,scale=TRUE)
    hc<-hclust(dist(temp),"ave")
    plot(as.phylo(hc),main=methodnames,cex=0.5,tip.color=colt[filterED])
    axisPhylo(side=1)
  
  mtext(paste("Dendrograms - Built from",ncol(temp),"variables containing non-missing data",sep=" "),side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
  pageno=pageno+1    
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0)    
}
  
  #DE plots
{
  tout<-rbind(c(1,2,3),c(4))
  layout(tout)
  par(mar=c(2,2,2,1),oma=c(2,2,3,2),xpd=NA)
  barplot(colSums(anfdr<0.05),main="ANOVA",names=c(methodnames),border="red",density=20,cex=0.5,cex.axis=0.9,las=2,ylab="No. of Variables with FDR<0.05")
  barplot(colSums(kwfdr<0.05),main="Kruskal Wallis",names=c(methodnames),border="red",density=20,cex=0.5,cex.axis=0.9,las=2,ylab="No. of Variables with FDR<0.05")
  mtext("Differential Expression",side=3,adj=0,outer=TRUE,col="gray",line=1)
  mtext(paste("Page ",pageno,sep=""),side=1,adj=1,outer=T,col="gray",line=0)
  pageno=pageno+1
  mtext(paste("Project: ",currentjob[2],sep=""),side=1,outer=T,col="gray",line=0)
  mtext("Normalyzer Report",side=1,adj=0,outer=T,col="gray",line=0)    
}
  
  dev.off()
  
  }
}