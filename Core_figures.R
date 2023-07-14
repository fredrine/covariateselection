#### Load libraries and specify folders ----
library(splines)
library(pbs)
library(MGLM)
library(fields)
library(lmtest)
library(MASS)
library(mvtnorm)
library(ggbiplot)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(sparseMVN)
library(rgl)
library(fields)
library(tidyverse)
library(RColorBrewer)
library(dendextend)
library(rmatio)
library(R.matlab)
library(EnvStats)
library(rhdf5)
library(GenBinomApps)
library(ggvenn)

figure_folder = "/Users/fredrine/Documents/CalciumData/Figures/" # Folder in which the figures are saved
functions_folder = "/Users/fredrine/Documents/CalciumData/Scripts/" # Folder in which the functions_model_selection.R script lies
processed_data_folder = "/Users/fredrine/Documents/CalciumData/Processed_data/" # Folder in which the results from the initial processing lie (binarized cell activity, covariate matrix with splined versions of the covariates)
result_folder = "/Users/fredrine/Documents/CalciumData/Result_data/" # Folder in which the results from running the forward selection methods lie

#### Load functions ----
source(paste(functions_folder,"functions_model_selection.R",sep=""))
#### Figure plotting functions ----
# Figure 1A
Nonsense_corr = function(cex=2.5){
  parlwd = par("lwd")
  par(lwd=cex)
  correlation_coefficients_iid = readRDS(paste(figure_folder,"correlation_coefficients_iid.RDS",sep=""))
  correlation_coefficients_ac_a = readRDS(paste(figure_folder,"correlation_coefficients_ac_a.RDS",sep=""))
  correlation_coefficients_ac_b = readRDS(paste(figure_folder,"correlation_coefficients_ac_b.RDS",sep=""))
  sample_size = 200
  opaq = 1
  #par(mfrow=c(1,2),mar=c(5,5,5,5),lwd=3)
  
  #breaks = seq(-1,1,by=0.025)
  breaks = seq(-1,1,by=0.04)
  rs = seq(-1,1,length.out = 1000)
  ds = (1-rs^2)^((sample_size-4)/2)/beta(1/2,1/2*(sample_size-2))
  yl = c(0,max(c(max(ds,max(hist(correlation_coefficients_iid,breaks=breaks,plot=F)$density)))))
  h=hist(correlation_coefficients_ac_b,xlab="",ylim=yl,lwd=cex,main="", col=rgb(230/255,159/255,0/255,opaq),lty=1,breaks=breaks,freq=F,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,xaxt="n")
  title(xlab="Pearson correlation coefficient",cex=cex,cex.lab=cex,line=3.5)
  axis(side=1, at=seq(-1,1,by=0.5),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(-1,1,by=0.5),tick=F, labels=c("-1.0","-0.5","0.0","0.5","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2*5/3)
  #h=hist(correlation_coefficients_iid,xlab="Correlation coefficient",ylim=yl,main="", col=rgb(0/255,158/255,115/255,opaq),lty=1,breaks=breaks,freq=F,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
  hist(correlation_coefficients_ac_a, col=rgb(86/255,180/255,233/255,opaq),lty=1, add=T,breaks=breaks,freq=F)
  hist(correlation_coefficients_iid, col=rgb(0/255,158/255,115/255,opaq),lty=1, add=T,breaks=breaks,freq=F)
  #hist(correlation_coefficients_ac_b, col=rgb(240/255,228/255,66/255,opaq),lty=1, add=T,breaks=breaks,freq=F)
  lines(rs,ds,col="red",lty=1,lwd=cex*2)
  par(xpd=T)
  #legend(0.05,6.5,legend=c("Theoretical distribution","I.I.D.","Small autocorrelation","Large autocorrelation"),
  if (FALSE){
    legend(-0.03,1.159476*yl[2],legend=c("Theoretical distribution","I.I.D.","Small autocorrelation","Large autocorrelation"),
           col = c("red", NA, NA,NA),
           lty = c(1,NA,NA,NA),
           lwd = c(2*cex,cex,cex,cex),
           density=c(0,100,100,100),
           fill = c("red", "#009E73","#56B4E9","#E69F00"),
           border = c(NA,"black","black","black"),
           cex = cex,
           x.intersp = c(2.37,1,1,1),
           bty=c("n","n","n","n"),
           xjust=c(0.08,0,0,0)
    )
  }
  par(xpd=F,lwd=parlwd)
}

# Figure 1B
Nonsense_corr_pval = function(cex=2.5){
  parlwd = par("lwd")
  par(lwd=cex)
  opaq=1
  pvals_iid = readRDS(paste(figure_folder,"pvals_iid.RDS",sep=""))
  pvals_ac_a = readRDS(paste(figure_folder,"pvals_ac_a.RDS",sep=""))
  pvals_ac_b = readRDS(paste(figure_folder,"pvals_ac_b.RDS",sep=""))
  
  #breaks = seq(0,1,by=0.0125)
  breaks = seq(0,1,by=0.02)
  yl = c(0,3)
  hist(pvals_ac_a,xlab="",ylim=yl,main="", col=rgb(86/255,180/255,233/255,opaq),lwd=cex,lty=1,breaks=breaks,freq=F,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,xaxt="n")
  title(xlab="p-value (from correlation test)",cex=cex,cex.lab=cex,line=3.5)
  axis(side=1, at=seq(0,1,by=0.2),tick=T, labels=rep("",6),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,1,by=0.2),tick=F, labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2)
  hist(pvals_ac_b, col=rgb(230/255,159/255,0/255,opaq),lty=1, add=T,breaks=breaks,freq=F)
  hist(pvals_iid, col=rgb(0/255,158/255,115/255,opaq),lty=1, add=T,breaks=breaks,freq=F)
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  #legend("topright",legend=c("Theoretical distribution","I.I.D.","Small autocorrelation","Large autocorrelation"),col=c("black",rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),pch=c(NA,15,15,15),lwd=c(2,NA,NA,NA),lty=c(1,1,1,1),cex=cex)
  #legend("topright",legend=c("Theoretical distribution","","",""),col=c("black",NA,NA,NA),lwd=c(cex*2,NA,NA,NA),lty=c(2,1,1,1),cex=cex)
  #legend("topright",legend=c("                        ","I.I.D.","Small autocorrelation","Large autocorrelation"),bg="transparent",col=c(NA,NA,NA,NA),fill=c(NA,"#D55E00","#56B4E9","#F0E442"),lwd=c(NA,NA,NA,NA),lty=c(0,1,1,1),cex=cex)
  par(xpd=T)
  #legend(1.05/2*1,1.193975*yl[2],legend=c("Theoretical distribution","I.I.D.","Small autocorrelation","Large autocorrelation"),
  #legend(0.55/2*1*0.9,1.159476*yl[2]*0.9,legend=c("Theoretical distribution","I.I.D.","Small autocorrelation","Large autocorrelation"),
  legend(0.07,3.3,legend=c("Expected distribution (IID assumption)","IID data","Weakly autocorrelated data","Strongly autocorrelated data"),
         col = c("red", NA, NA,NA),
         lty = c(1,NA,NA,NA),
         lwd = c(2*cex,cex,cex,cex),
         density=c(0,100,100,100),
         fill = c("red", "#009E73","#56B4E9","#E69F00"),
         border = c(NA,"black","black","black"),
         cex = cex,
         x.intersp = c(0.72,-0.7,-0.7,-0.7),
         bty=c("n","n","n","n"),
         xjust=c(0.06,0,0,0)
  )
  par(xpd=F,lwd=parlwd)
}

# Figure 1C
Nonsense_dev = function(cex=2.5){
  parlwd = par("lwd")
  par(lwd=cex*1.2)
  
  lrs = readRDS(paste(figure_folder,"lrs.RDS",sep=""))
  
  opaq = 1
  lrs_ = lrs
  lrs_[which(lrs >= 30)] = NA
  #breaks = seq(0,30,by=0.375)
  breaks = seq(0,30,by=0.6)
  xs = seq(0,30,length.out = 1000)
  ds = dchisq(xs,1)
  yl = c(0,max(hist(lrs_,plot=F)$density))
  h=hist(lrs_[4,],xlab="",ylim=yl,main="", col=rgb(230/255,159/255,0/255,opaq),lty=1,lwd=cex,breaks=breaks,freq=F,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,xaxt="n")
  title(xlab="Deviance (Poisson GLM)",cex=cex,cex.lab=cex,line=3.5)
  axis(side=1, at=seq(0,30,by=5),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,30,by=5),tick=F, labels=seq(0,30,by=5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2*yl[2]/3)
  hist(lrs_[2,], col=rgb(86/255,180/255,233/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  hist(lrs_[1,], col=rgb(0/255,158/255,115/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  hist(lrs_[3,], col=rgb(240/255,228/255,66/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  lines(xs,ds,col="red",lty=1,lwd=cex*2)
  par(xpd=T)
  #legend(1.05/2*30,1.193975*yl[2],legend=c("Theoretical distribution","I.I.D.","Autocorrelated w. AR","Autocorrelated w.o. AR"),
  if (FALSE){
    legend(0.97/2*30,1.159476*yl[2],legend=c("Theoretical distribution","I.I.D., no missing variable","I.I.D, missing variable","Autocorr., no missing variable","Autocorr., missing variable"),
           col = c("red", NA, NA,NA,NA),
           lty = c(1,NA,NA,NA,NA),
           lwd = c(2*cex,cex,cex,cex,cex),
           density=c(0,100,100,100,100),
           fill = c("red", "#009E73","#56B4E9","#F0E442","#E69F00"),
           border = c(NA,"black","black","black","black"),
           cex = cex,
           #x.intersp = c(2.45,1,1,1),
           bty=c("n","n","n","n","n"),
           x.intersp = c(2.37,1,1,1,1),
           xjust=c(0.08,0,0,0,0)
    )
  }
  
  par(xpd=F,lwd=parlwd)
}

# Figure 1D
Nonsense_dev_pval = function(cex=2.5){
  parlwd = par("lwd")
  par(lwd=cex)
  LRpvals = readRDS(paste(figure_folder,"LRpvals.RDS",sep=""))
  opaq = 1
  #breaks = seq(0,1,by=0.0125)
  breaks = seq(0,1,by=0.02)
  yl = c(0,3)
  h=hist(LRpvals[4,],xlab="",ylim=yl,main="", col=rgb(230/255,159/255,0/255,opaq),lty=1,lwd=cex,breaks=breaks,freq=F,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,xaxt="n")
  title(xlab="p-value (from likelihood-ratio test)",cex=cex,cex.lab=cex,line=3.5)
  axis(side=1, at=seq(0,1,by=0.2),tick=T, labels=rep("",6),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,1,by=0.2),tick=F, labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2)
  
  hist(LRpvals[2,], col=rgb(86/255,180/255,233/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  hist(LRpvals[1,], col=rgb(0/255,158/255,115/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  hist(LRpvals[3,], col=rgb(240/255,228/255,66/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  #h=hist(LRpvals[1,],xlab="p-value",ylim=yl,main="", col=rgb(0/255,158/255,115/255,opaq),lty=1,breaks=breaks,freq=F,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
  #hist(LRpvals[2,], col=rgb(86/255,180/255,233/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  #hist(LRpvals[3,], col=rgb(240/255,228/255,66/255,opaq),lty=1,breaks=breaks, add=T,freq=F)
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  par(xpd=T)
  #legend(1.05/2*1,1.193975*yl[2],legend=c("Theoretical distribution","I.I.D.","Autocorrelated w. AR","Autocorrelated w.o. AR"),
  legend(0.07,3.3,legend=c("Expected distribution (GLM assumptions)","IID data, no missing variable","IID data, missing variable","Autocorrelated data, no missing variable","Autocorrelated data, missing variable"),
         col = c("red", NA, NA,NA,NA),
         lty = c(1,NA,NA,NA,NA),
         lwd = c(2*cex,cex,cex,cex,cex),
         density=c(0,100,100,100,100),
         fill = c("red", "#009E73","#56B4E9","#F0E442","#E69F00"),
         border = c(NA,"black","black","black","black"),
         cex = cex,
         #x.intersp = c(2.45,1,1,1),
         x.intersp = c(0.73,-0.7,-0.7,-0.7,-0.7),
         bty=c("n","n","n","n"),
         xjust=c(0.057,0,0,0,0)
  )
  par(xpd=F,lwd=parlwd)
}

# Figure 1E
CV_error_rates = function(cex=2.5){
  parlwd = par("lwd")
  par(lwd=cex)
  
  fakeresmat = readRDS(paste(result_folder,"CV_variants_no_effect.RDS",sep=""))
  
  N = 12000
  fn = c(20)
  cs = c(2,4,8,15,30,60,150,300,600)
  ns = c(TRUE,FALSE)
  RANDOMFOLDS = FALSE #Random folds when chunksize = 1
  RANDOMFOLDSONLAST = TRUE #The above is true for the last variant
  
  fold_nums = c()
  chunksizes = c()
  neighbor_skips = c()
  
  for (f in fn){
    for (c in cs){
      for (s in ns){
        fold_nums = c(fold_nums,f)
        chunksizes = c(chunksizes,N/f/c)
        neighbor_skips = c(neighbor_skips,s)
      }
    }
  }
  
  fold_nums = c(fold_nums,20)
  chunksizes = c(chunksizes,1)
  neighbor_skips = c(neighbor_skips,FALSE)
  
  num_variants = length(fold_nums)
  
  num_trials = dim(fakeresmat)[1]/num_variants
  
  plot_error_rates = function(fake_resmat,num_variants,fold_nums,chunksizes,neighbor_skips,num_trials,whichrandom=0,acfs=NULL){
    if (whichrandom == 0){
      whichrandom = num_variants
    }
    FakeResultDF = data.frame(fake_resmat)
    colnames(FakeResultDF) = c("candidate","pval","CV_score","num_folds","chunksize","skip_neighbor","variant")
    
    errorrates = array(0,num_variants)
    errors = array(0,num_variants)
    lowerconf = array(0,num_variants)
    upperconf = array(0,num_variants)
    
    for (i in 1:num_variants){
      vals = FakeResultDF$CV_score[which(FakeResultDF$variant == i)]
      errors[i] = length(which(vals > 0))
      errorrates[i] = length(which(vals > 0))/num_trials
      ci = clopper.pearson.ci(errors[i],num_trials,0.05,"two.sided")
      lowerconf[i] = ci$Lower.limit
      upperconf[i] = ci$Upper.limit
    }
    
    if (TRUE){
      cols = rep("#56B4E9",num_variants)
      cols[which(neighbor_skips)] = "#009E73"
      cols[whichrandom] = "gold"
      xaxis = sort(unique(log(chunksizes)/log(10)))
      xaxsec = sort(unique((chunksizes)))
      x = log(chunksizes)/log(10)
      x[which(neighbor_skips)] = x[which(neighbor_skips)]+0.01
      x[which(!neighbor_skips)] = x[which(!neighbor_skips)]-0.01
      x[whichrandom] = x[whichrandom]-0.08
    }
    nonrandom = setdiff(1:num_variants,whichrandom)
    
    xl = c(-0.15,2.5)
    plot(x[nonrandom],errorrates[nonrandom],bg=cols[nonrandom],pch=22,lwd=cex*1.2,cex=cex*1.2,xaxt="n",yaxt="n",cex.lab=cex,cex.axis=cex,cex.main=cex,cex.sub=cex,ylab="Error rate of CV scheme",xlab="",xaxt="n",main="",ylim=c(0,1),xlim=xl)
    if (!is.null(acfs)){
      lines(log(1:300)/log(10),acfs,col="darkgray")
    }
    title(xlab="Block size in bins (12000 bins in total)",cex=cex,cex.lab=cex,line=3.5)
 
    
    for (i in nonrandom){
      lines(c(x[i],x[i]),c(lowerconf[i],upperconf[i]),lwd=cex*1.5)
      lines(c(x[i]-0.03,x[i]+0.03),c(lowerconf[i],lowerconf[i]),lwd=cex*1.5)
      lines(c(x[i]-0.03,x[i]+0.03),c(upperconf[i],upperconf[i]),lwd=cex*1.5)
    }
    points(x[nonrandom],errorrates[nonrandom],bg=cols[nonrandom],pch=22,cex=cex*1.2,lwd=cex*1.2,cex.lab=cex,cex.axis=cex,cex.main=cex,cex.sub=cex)
    axis(side=1, at=xaxis,tick=T, labels=rep("",9),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    axis(side=1, at=xaxis,tick=F, labels=xaxsec,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2*1/3)
    
    axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels=c("0.0",0.2,0.4,0.6,0.8,"1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    
    lines(c(-0.13,-0.13),c(lowerconf[whichrandom],upperconf[whichrandom]),lwd=cex*1.5)
    lines(c(-0.13-0.03,-0.13+0.03),c(lowerconf[whichrandom],lowerconf[whichrandom]),lwd=cex*1.5)
    lines(c(-0.13-0.03,-0.13+0.03),c(upperconf[whichrandom],upperconf[whichrandom]),lwd=cex*1.5)
    points(-0.13,errorrates[whichrandom],bg="#E69F00",pch=22,cex=cex*1.2,lwd=cex*1.2,cex.lab=cex,cex.axis=cex,cex.main=cex,cex.sub=cex)

    legend(1.1,1.07,legend=c("No blocking/skipping","Blocking","Blocking + skipping"),
           col = c(NA, NA, NA),
           lty = c(NA,NA,NA),
           lwd = c(cex,cex,cex,cex),
           density=c(100,100,100),
           fill = c("#E69F00","#56B4E9","#009E73"),
           border = c("black","black","black"),
           cex = cex,
           x.intersp = c(-0.7,-0.7,-0.7),
           bty=c("n","n","n"),
           xjust=c(0,0,0)
    )
  }
  plot_error_rates(fakeresmat,num_variants,fold_nums,chunksizes,neighbor_skips,num_trials)
  par(lwd=parlwd)
}

# Figure 1F
SR_error_rates = function(cex=2.5){
  parlwd = par("lwd")
  par(lwd=cex)
  fakeresmat = readRDS(paste(result_folder,"CV_SR_variants_no_effect.RDS",sep=""))
  
  N = 12000
  fn = c(20)
  cs = c(2,4,8,15,30,60,150,300,600)
  ns = c(TRUE,FALSE)
  RANDOMFOLDS = FALSE #Random folds when chunksize = 1
  RANDOMFOLDSONLAST = TRUE #The above is true for the last variant
  fold_nums = c()
  chunksizes = c()
  neighbor_skips = c()
  
  for (f in fn){
    for (c in cs){
      for (s in ns){
        fold_nums = c(fold_nums,f)
        chunksizes = c(chunksizes,N/f/c)
        neighbor_skips = c(neighbor_skips,s)
      }
    }
  }
  
  fold_nums = c(fold_nums,20)
  chunksizes = c(chunksizes,1)
  neighbor_skips = c(neighbor_skips,FALSE)
  num_variants = length(fold_nums)
  num_trials = dim(fakeresmat)[1]/num_variants
  
  plot_error_rates = function(fake_resmat,num_variants,fold_nums,chunksizes,neighbor_skips,num_trials,whichrandom=0){
    if (whichrandom == 0){
      whichrandom = num_variants
    }
    FakeResultDF = data.frame(fake_resmat)
    colnames(FakeResultDF) = c("candidate","pval","CV_score","num_folds","chunksize","skip_neighbor","variant")
    
    errorrates = array(0,num_variants)
    errors = array(0,num_variants)
    lowerconf = array(0,num_variants)
    upperconf = array(0,num_variants)
    
    for (i in 1:num_variants){
      pvals = FakeResultDF$pval[which(FakeResultDF$variant == i)]
      errors[i] = length(which(pvals <= 0.05))
      errorrates[i] = length(which(pvals <= 0.05))/num_trials
      ci = clopper.pearson.ci(errors[i],num_trials,0.05,"two.sided")
      lowerconf[i] = ci$Lower.limit
      upperconf[i] = ci$Upper.limit
    }
    
    if (TRUE){
      cols = rep("#56B4E9",num_variants)
      cols[which(neighbor_skips)] = "#009E73"
      cols[whichrandom] = "gold"
      xaxis = sort(unique(log(chunksizes)/log(10)))
      xaxsec = sort(unique((chunksizes)))
      x = log(chunksizes)/log(10)
      x[which(neighbor_skips)] = x[which(neighbor_skips)]+0.01
      x[which(!neighbor_skips)] = x[which(!neighbor_skips)]-0.01
      x[whichrandom] = x[whichrandom]-0.08
    }
    nonrandom = setdiff(1:num_variants,whichrandom)
    yl = c(0,0.35)
    xl = c(-0.15,2.5)
    plot(x[nonrandom],errorrates[nonrandom],bg=cols[nonrandom],pch=22,lwd=cex*1.2,cex=cex*1.2,xaxt="n",yaxt="n",cex.lab=cex,cex.axis=cex,cex.main=cex,cex.sub=cex,ylab="Error rate of signed-rank test",xlab="",main="",ylim=yl,xlim=xl)
    title(xlab="Block size in bins (12000 bins in total)",cex=cex,cex.lab=cex,line=3.5)
    
    for (i in nonrandom){
      lines(c(x[i],x[i]),c(lowerconf[i],upperconf[i]),lwd=cex*1.5)
      lines(c(x[i]-0.03,x[i]+0.03),c(lowerconf[i],lowerconf[i]),lwd=cex*1.5)
      lines(c(x[i]-0.03,x[i]+0.03),c(upperconf[i],upperconf[i]),lwd=cex*1.5)
    }
    points(x[nonrandom],errorrates[nonrandom],bg=cols[nonrandom],pch=22,lwd=cex*1.2,xaxt="n",cex=cex*1.2,cex.lab=cex,cex.axis=cex,cex.main=cex,cex.sub=cex)
    #axis(side=1, at=xaxis, labels=xaxsec,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    axis(side=1, at=xaxis,tick=T, labels=rep("",9),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    axis(side=1, at=xaxis,tick=F, labels=xaxsec,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2*0.35/3)
    axis(side=2, at=c(0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35), labels=c("0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    abline(h=0.05,col="red",lty=1,lwd=cex*2)
    legend(1.1,1.07*yl[2],legend=c("Desired error rate","Blocking","Blocking + skipping"),
           col = c("red", NA, NA),
           lty = c(1,NA,NA),
           lwd = c(2*cex,cex,cex,cex),
           density=c(0,100,100,100),
           fill = c("red","#56B4E9","#009E73"),
           border = c(NA,"black","black"),
           cex = cex,
           x.intersp = c(0.72,-0.7,-0.7),
           bty=c("n","n","n"),
           xjust=c(0.097,0,0)
    )
    
  }
  
  plot_error_rates(fakeresmat,num_variants,fold_nums,chunksizes,neighbor_skips,num_trials)
  
  par(lwd=parlwd)
}

# Figure 1G
Multiple_comparisons = function(cex=2.5){
  parlwd=par("lwd")
  par(lwd=cex)
  x = 1:70
  r1 = 0.05
  y1 = 1-(1-r1)^x
  r2 = 0.02
  y2 = 1-(1-r2)^x
  r3 = 0.01
  y3 = 1-(1-r3)^x
  plot(x,y1,col="#009E73",cex=cex,xaxt="n",yaxt="n",lwd=cex*2,cex.main=cex,cex.lab=cex,cex.axis=cex,xlab="",ylab="Experimentwise error rate",ylim=c(0,1),type="l",xlim=c(1,70),xaxt="n")
  title(xlab="Number of comparisons",cex=cex,cex.lab=cex,line=3.5)
  axis(side=1, at=seq(0,70,by=10),tick=T, labels=rep("",8),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,70,by=10),tick=F, labels=seq(0,70,by=10),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.2*1/3)
  
  points(x,y1,col="#009E73",cex=cex*2/3,pch=19)
  lines(x,y2,col="#56B4E9",cex=cex,lwd=cex*2,cex.main=cex,cex.lab=cex,cex.axis=cex)
  points(x,y2,col="#56B4E9",cex=cex*2/3,pch=15)
  lines(x,y3,col="#E69F00",cex=cex,lwd=cex*2,cex.main=cex,cex.lab=cex,cex.axis=cex)
  points(x,y3,col="#E69F00",cex=cex*2/3,pch=17)
  #axis(side=1, at=c(0,10,20,30,40,50,60,70), labels=c("0","10","20","30","40","50","60","70"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=c(0.0,0.2,0.4,0.6,0.8,1.0), labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  
  legend(1.4,1.07,legend=c("Error rate of a single test:","0.05","0.03","0.01"),
         col = c("white","#009E73","#56B4E9","#E69F00"),
         pch = c(NA,19,15,17),
         lty = c(NA,1,1,1),
         lwd = c(NA,cex*2,cex*2,cex*2),
         cex = c(cex,cex,cex,cex),
         x.intersp = c(-3,1,1,1),
         bty=c("n","n","n","n"),
         xjust=c(0,0,0,0)
  )
  par(lwd=parlwd)
}

# Figure 2
Hypothetical_p_distributions = function(cex=3){
  pvals_single = readRDS(paste(figure_folder,"pvals_single.RDS",sep=""))
  pvals_max = readRDS(paste(figure_folder,"pvals_max.RDS",sep=""))
  pvals_max_corr = readRDS(paste(figure_folder,"pvals_max_corr.RDS",sep=""))
  
  parlwd=par("lwd")
  par(lwd=cex)
  breaks = seq(0,1,by=0.04)
  yl = c(0,3)
  
  hist(pvals_max,breaks=breaks,ylim=yl,main="Invalid test (overly permissive)",xlab="p-value",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,freq=F,xaxt="n")
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  axis(1,at=seq(0,1,by=0.2),tick=T,labels=rep("",6),pos=-0.05,cex.axis=cex)
  axis(1,at=seq(0,1,by=0.2),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),tick=F,pos=-0.1,cex.axis=cex,col="white")
  
  hist(pvals_max_corr,breaks=breaks,ylim=yl,xlab="p-value",main="Valid test",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,freq=F,xaxt='n')
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  axis(1,at=seq(0,1,by=0.2),tick=T,labels=rep("",6),pos=-0.05,cex.axis=cex)
  axis(1,at=seq(0,1,by=0.2),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),tick=F,pos=-0.1,cex.axis=cex,col="white")
  legend("topright",legend=c("Theoretical distribution"),lty=1,col="red",lwd=cex*2,cex=cex,bty="n")
  
  hist(1-pvals_max,breaks=breaks,ylim=yl,xlab="p-value",main="Conservative test",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,freq=F,xaxt='n')
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  axis(1,at=seq(0,1,by=0.2),tick=T,labels=rep("",6),pos=-0.05,cex.axis=cex)
  axis(1,at=seq(0,1,by=0.2),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),tick=F,pos=-0.1,cex.axis=cex,col="white")
  
  par(lwd=parlwd)
}

# Figure 3
Simulated_false_results = function(cex=3,line=-0.1,logscale=T){
  alpha = 0.05
  
  newlabs = c()
  newlabs[1] = expression(paste("CV"[""]))
  newlabs[2] = expression(paste("mSR"["MaxT"]))
  newlabs[3] = expression(paste("mSRR"["MaxT"]))
  newlabs[4] = expression(paste("CS"["Bonf"]))
  newlabs[5] = expression(paste("SR"[""]))
  newlabs[6] = expression(paste("SR"["Bonf"]))
  
  methods1 = c("CV","mSR","mSRR","CS")
  methods2 = c("CV","SR","SR_B")
  FakeMat_ = readRDS(paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))
  FakeMat2_ = readRDS(paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))
  
  covnames = c("B (1D)","A (2D)","C (1D)")
  
  plot_combined_results_better_bars(FakeMat_,FakeMat2_,methods1,methods2,newlabs,alpha=alpha,cex=cex,covnames=covnames,hline=T,logscale=logscale,line=line)
}

# Figure 4
Simulated_true_results = function(cex=3,line=-0.1,logscale=T){
  alpha = 0.05
  
  newlabs = c()
  newlabs[1] = expression(paste("CV"[""]))
  newlabs[2] = expression(paste("mSR"["MaxT"]))
  newlabs[3] = expression(paste("mSRR"["MaxT"]))
  newlabs[4] = expression(paste("CS"["Bonf"]))
  newlabs[5] = expression(paste("SR"[""]))
  newlabs[6] = expression(paste("SR"["Bonf"]))
  
  methods1 = c("CV","mSR","mSRR","CS")
  methods2 = c("CV","SR","SR_B")

  TrueMat_ = readRDS(paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
  TrueMat2_ = readRDS(paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
  
  covnames = c("B (1D)","A (2D)","C (1D)")
  
  plot_combined_results_better_bars(TrueMat_,TrueMat2_,methods1,methods2,newlabs,alpha=alpha,cex=cex,covnames=covnames,logscale = logscale,line=line)
}

# Figure 5
Calcium_false_results = function(cex=3,line=-0.1,logscale=T){
  alpha = 0.05
  
  newlabs = c()
  newlabs[1] = expression(paste("CV"[""]))
  newlabs[2] = expression(paste("mSR"["MaxT"]))
  newlabs[3] = expression(paste("mSRR"["MaxT"]))
  newlabs[4] = expression(paste("CS"["Bonf"]))
  newlabs[5] = expression(paste("SR"[""]))
  newlabs[6] = expression(paste("SR"["Bonf"]))
  
  methods1 = c("CV","mSR","mSRR","CS")
  methods2 = c("CV","SR","SR_B")
  Fake_ = readRDS(paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))
  Fake2_ = readRDS(paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))
  
  plot_combined_results_better_bars(Fake_,Fake2_,methods1,methods2,newlabs,alpha=alpha,cex=cex,covnames=c("HD","Pos","Spe"),hline = T,logscale = logscale,line=line)
  
}

# Figure 6
Calcium_true_results = function(cex=3,line=-0.1,logscale=F){
  alpha = 0.05
  
  newlabs = c()
  newlabs[1] = expression(paste("CV"[""]))
  newlabs[2] = expression(paste("mSR"["MaxT"]))
  newlabs[3] = expression(paste("mSRR"["MaxT"]))
  newlabs[4] = expression(paste("CS"["Bonf"]))
  newlabs[5] = expression(paste("SR"[""]))
  newlabs[6] = expression(paste("SR"["Bonf"]))
  
  methods1 = c("CV","mSR","mSRR","CS")
  methods2 = c("CV","SR","SR_B")
  True_ = readRDS(paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
  True2_ = readRDS(paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))

  plot_combined_results_better_bars(True_,True2_,methods1,methods2,newlabs,alpha=alpha,cex=cex,covnames=c("HD","Pos","Spe"),logscale = logscale,line=line)
}

# Figure 7
Ratemaps = function(cex=3){
  trackdata2 = readRDS(paste(processed_data_folder,"trackdata2.RDS",sep=""))
  Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
  gridandposCS = readRDS(paste(result_folder,"gridandposCS.RDS",sep=""))
  posnotgridCS = readRDS(paste(result_folder,"posnotgridCS.RDS",sep=""))
  gridnotposCS = readRDS(paste(result_folder,"gridnotposCS.RDS",sep=""))
  
  
  Xvec=trackdata2$headX[1:12000]
  Yvec=trackdata2$headY[1:12000]
  Xvec2=trackdata2$HD[1:12000]
  Yvec2=trackdata2$bodyspeed[1:12000]
  
  meanrates = apply(Y1,2,mean)/(1/7.25)
  colorbar = hcl.colors(100,palette="SunsetDark")
  
  plot_rms_from_list_normalised(X=Xvec,Y=Yvec,Y1,gridandposCS[1:9],main=paste("Grid cells classified as position (average rate of events: ",round(mean(meanrates[gridandposCS])*1000)/1000,")",sep=""),min_occ=0.1,binning=40,smoothing_sd=rep(1,4),noaxis=T,cex=cex,mpar=NULL,cbar=colorbar)
  plot_rms_from_list_normalised(X=Xvec,Y=Yvec,Y1,posnotgridCS[1:9],main=paste("Non-grid cells classified as position (average rate of events: ",round(mean(meanrates[posnotgridCS])*1000)/1000,")",sep=""),min_occ=0.1,binning=40,smoothing_sd=rep(1,4),noaxis=T,cex=cex,mpar=NULL,cbar=colorbar)
  plot_rms_from_list_normalised(X=Xvec,Y=Yvec,Y1,gridnotposCS[1:9],main=paste("Grid cells not classified as position (average rate of events: ",round(mean(meanrates[gridnotposCS])*1000)/1000,")",sep=""),min_occ=0.1,binning=40,smoothing_sd=rep(1,4),noaxis=T,cex=cex,mpar=NULL,cbar=colorbar)
}

# Figure 7B
Colorbar = function(cex=3){
  par(mfrow=c(3,1))
  plot.new()
  plot(rep(0.45,100),1:100,col=colorbar,pch=15,cex=cex*10,xlab="",ylab="",bty="n",xaxt="n",yaxt="n",xlim=c(0.45,1.1),ylim=c(1,100))
  text(0.62,y=97,pos=4,labels="95th",cex=cex)
  text(0.62,y=85,pos=4,labels="percentile",cex=cex)
  text(0.62,y=2,pos=4,labels="0 events/sec",cex=cex)
  plot.new()
}

# Figure 8
Venn_diagrams = function(cex=3){
  grid_cells = readRDS(paste(processed_data_folder,"gridcells.RDS",sep=""))
  
  hdtunedCS = readRDS(paste(result_folder,"hdtunedCS.RDS",sep=""))
  postunedCS = readRDS(paste(result_folder,"postunedCS.RDS",sep=""))
  spetunedCS = readRDS(paste(result_folder,"spetunedCS.RDS",sep=""))
  
  hdtunedmSRR = readRDS(paste(result_folder,"hdtunedmSRR.RDS",sep=""))
  postunedmSRR = readRDS(paste(result_folder,"postunedmSRR.RDS",sep=""))
  spetunedmSRR = readRDS(paste(result_folder,"spetunedmSRR.RDS",sep=""))
  
  hdtunedSR = readRDS(paste(result_folder,"hdtunedSR.RDS",sep=""))
  postunedSR = readRDS(paste(result_folder,"postunedSR.RDS",sep=""))
  spetunedSR = readRDS(paste(result_folder,"spetunedSR.RDS",sep=""))
  
  
  xCS = list("HD"=hdtunedCS,"Speed"=spetunedCS,"Position"=postunedCS,"Grid (Zong)"=grid_cells)
  xmSRR = list("HD"=hdtunedmSRR,"Speed"=spetunedmSRR,"Position"=postunedmSRR,"Grid (Zong)"=grid_cells)
  xSR = list("HD"=hdtunedSR,"Speed"=spetunedSR,"Position"=postunedSR,"Grid (Zong)"=grid_cells)
  
  bluecol = "#56B4E9"
  orangecol = "#D55E00"
  yellowcol = "#F0E442"
  graycol = "#999999"
  
  size = cex*2
  
  vennCS = ggvenn(xCS, stroke_size = cex/2,text_size=size*1.25,set_name_size=size*1.35,digits=0,fill_alpha = 0.7,
                  fill_color = c(orangecol,yellowcol,bluecol, graycol))  + ggtitle(expression('CS'['Bonf'])) + theme(plot.title = element_text(hjust = 0.5,vjust=6,size=size*6))
  
  
  vennSR = ggvenn(xmSRR, stroke_size = cex/2,text_size=size*1.25,set_name_size=size*1.35,digits=0,fill_alpha = 0.7,
                  fill_color = c(orangecol,yellowcol,bluecol, graycol)) + ggtitle(expression('mSRR'['MaxT'])) + theme(plot.title = element_text(hjust = 0.5,vjust=6,size=size*6))
  
  
  vennWT = ggvenn(xSR, stroke_size = cex/2,text_size=size*1.25,set_name_size=size*1.35,digits=0,fill_alpha = 0.7,
                  fill_color = c(orangecol,yellowcol,bluecol, graycol)) + ggtitle('SR') + theme(plot.title = element_text(hjust = 0.5,vjust=6,size=size*6))
  
  grid.arrange(vennCS,vennSR,vennWT,ncol=3)
}

# Figure 9
Blocking_illustration = function(cex=3){
  
  parlwd = par("lwd")
  par(lwd=cex)
  lwd = cex
  bluecol = "#56B4E9"
  greencol = "#009E73"
  #orangecol = "#D55E00"
  yellowcol = "#F0E442"
  
  nfolds = 8
  xl = c(nfolds+4,4*nfolds*10)
  yl = c(6,114)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  
  plot(0,0,col="white",xlim=xl,ylim=yl,xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
  
  xlefts = seq(1,3*nfolds*10,by=10)
  xrights = seq(9,3*nfolds*10,by=10)
  labelxs = seq(5,3*nfolds*10,by = 10)
  ybottoms = c(2,42,62,82,102)
  ytops = c(18,58,78,98,118)
  labelys = c(10,50,70,90,110)
  numbers = rep(1:nfolds,3)
  cols = list()
  cols[[1]] = rep(c(bluecol,yellowcol,rep(greencol,nfolds-3),yellowcol),3)
  cols[[2]] = rep(c(bluecol,yellowcol,rep(greencol,nfolds-3),yellowcol),3)
  cols[[3]] = rep(c(bluecol,yellowcol,rep(greencol,nfolds-3),yellowcol),3)
  cols[[4]] = rep(c(bluecol,yellowcol,rep(greencol,nfolds-3),yellowcol),3)
  cols[[5]] = rep(c(bluecol,yellowcol,rep(greencol,nfolds-3),yellowcol),3)
  cols[[1]][nfolds:(3*nfolds)] = cols[[5]][1:(2*nfolds+1)]
  cols[[1]][1:(nfolds-1)] = cols[[5]][(2*nfolds+2):(3*nfolds)]
  cols[[4]][2:(3*nfolds)] = cols[[5]][1:(3*nfolds-1)]
  cols[[4]][1] = cols[[5]][(3*nfolds)]
  cols[[3]][3:(3*nfolds)] = cols[[5]][1:(3*nfolds-2)]
  cols[[3]][1:2] = cols[[5]][(3*nfolds-1):(3*nfolds)]
  cols[[2]][4:(3*nfolds)] = cols[[5]][1:(3*nfolds-3)]
  cols[[2]][1:3] = cols[[5]][(3*nfolds-2):(3*nfolds)]
  
  for (i in 1:length(xlefts)){
    for (j in 1:length(ybottoms)){
      rect(xlefts[i],ybottoms[j],xrights[i],ytops[j],density=100,col=cols[[j]][i],border="black",lty=1,lwd=lwd)
      text(labelxs[i],labelys[j],labels=paste(numbers[i]),cex=cex)
    }
  }
  points(c(3*nfolds*5-nfolds,3*nfolds*5,3*nfolds*5+nfolds),rep((30),3),pch=19,col="black",cex=2)
  
  
  legend((3*nfolds*10),yl[2]*1.1,legend=c("Block from the \ntraining set","Skipped block","Block from the \ntest set","Number = the fold the \nblock is assigned to"),
         col = c(NA, NA, NA,NA),
         lty = c(NA,NA,NA,NA),
         lwd = c(cex*2,cex*2,cex*2,cex*2,NA),
         density=c(100,100,100,0),
         fill = c(greencol,yellowcol,bluecol,"white"),
         border = c("black","black","black","white"),
         cex = cex,
         x.intersp = c(-1,-1,-1,-2.8),
         y.intersp = c(2,2,2,2),
         bty=c("n","n","n","n"),
         xjust=c(0,0,0,0)
  )
  par(lwd=parlwd)
}

# Figure 10
Cyclical_illustration = function(cex=3){
  parlwd = par("lwd")
  par(lwd=cex)
  lwd=cex
  colorbar = hcl.colors(100,palette="SunsetDark")
  
  xl = c(0,100)
  yl = c(0,90)
  
  
  plot(0,0,col="white",xlim=xl,ylim=yl,xlab="",ylab="",main="",xaxt="n",bty="n",yaxt="n")
  
  xlefts = seq(0,99,by=1)
  xrights = 1+xlefts
  ybottoms = c(15,45,75)
  ytops = c(30,60,90)
  cols = list()
  cols[[1]] = colorbar
  cols[[2]] = colorbar
  cols[[3]] = colorbar
  cols[[1]][31:100] = cols[[3]][1:70]
  cols[[1]][1:30] = cols[[3]][71:100]
  cols[[2]][69:71] = "white"
  
  
  for (i in 1:length(xlefts)){
    for (j in 1:length(ybottoms)){
      rect(xlefts[i],ybottoms[j],xrights[i],ytops[j],density=100,col=cols[[j]][i])
    }
  }
  text(0,67,"Original data",pos=4,cex=cex)
  text(0,7,"Cyclically shifted data",pos=4,cex=cex)
  
  lwd = 4
  lines(c(71,71),c(43,41),col="black",lwd=lwd)
  lines(c(100,100),c(43,41),col="black",lwd=lwd)
  lines(c(71,100),c(41,41),col="black",lwd=lwd)
  lines(c(85.5,85.5),c(41,34),col="black",lwd=lwd)
  lines(c(0,85.5),c(34,34),col="black",lwd=lwd)
  lines(c(0,0),c(41,34),col="black",lwd=lwd)
  lines(c(0,-0.5),c(41,38),col="black",lwd=lwd)
  lines(c(0,0.5),c(41,38),col="black",lwd=lwd)
  
  par(lwd=parlwd)
}

# Figure 11
Multiple_corrections = function(cex=3){
  pvals_single = readRDS(paste(figure_folder,"pvals_single.RDS",sep=""))
  pvals_max = readRDS(paste(figure_folder,"pvals_max.RDS",sep=""))
  pvals_max_corr = readRDS(paste(figure_folder,"pvals_max_corr.RDS",sep=""))
  
  parlwd=par("lwd")
  par(lwd=cex)
  breaks = seq(0,1,by=0.04)
  yl = c(0,3)
  ncov = 5
  
  hist(pvals_max,breaks=breaks,ylim=yl,main="Uncorrected",xlab="p-value",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,freq=F,xaxt="n")
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  axis(1,at=seq(0,1,by=0.2),tick=T,labels=rep("",6),pos=-0.05,cex.axis=cex)
  axis(1,at=seq(0,1,by=0.2),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),tick=F,pos=-0.1,cex.axis=cex,col="white")
  
  pvals_bonfcorr = pvals_max*ncov
  pvals_bonfcorr[which(pvals_bonfcorr > 1)] = 1
  hist(pvals_bonfcorr,breaks=breaks,ylim=yl,xlab="p-value",main="Bonferroni corrected",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,freq=F,xaxt='n')
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  axis(1,at=seq(0,1,by=0.2),tick=T,labels=rep("",6),pos=-0.05,cex.axis=cex)
  axis(1,at=seq(0,1,by=0.2),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),tick=F,pos=-0.1,cex.axis=cex,col="white")
  
  hist(pvals_max_corr,breaks=breaks,ylim=yl,xlab="p-value",main="MaxT corrected",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,freq=F,xaxt='n')
  lines(c(0,1),c(1,1),col="red",lty=1,lwd=cex*2)
  axis(1,at=seq(0,1,by=0.2),tick=T,labels=rep("",6),pos=-0.05,cex.axis=cex)
  axis(1,at=seq(0,1,by=0.2),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),tick=F,pos=-0.1,cex.axis=cex,col="white")
  legend("topright",legend=c("Theoretical distribution"),lty=1,col="red",lwd=cex*2,cex=cex,bty="n")
  
  par(lwd=parlwd)
}

# Figure 12
Data_similarities = function(cex=3){
  parlwd = par("lwd")
  par(lwd=cex)
  
  set.seed(1)
  N = 12000
  sublength = 2500
  sd = 20
  bx = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  by = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  
  trackdata2 = readRDS(paste(processed_data_folder,"trackdata2.RDS",sep=""))
  headX = trackdata2$headX[1:N]
  headY = trackdata2$headY[1:N]
  
  yl = c(-44,44)
  plot(headX[1:sublength],type="l",xlab="",ylab="",cex.main=cex,cex.lab=cex,cex.axis=cex,ylim=yl,xaxt="n",yaxt="n")
  title(xlab="Time bins",cex=cex,cex.lab=cex,line=4)
  title(ylab="X",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,sublength,by=500),tick=T, labels=rep("",sublength/500+1),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,sublength,by=500),tick=F, labels=seq(0,sublength,by=500),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-44-0.2*88/3)
  axis(side=2, at=seq(-40,40,by=20),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(-40,40,by=20),tick=F, labels=seq(-40,40,by=20),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-44-0.2*88/3)
  
  yl = c(-0.1,1)
  par(cex.main=cex*1.2)
  a=acf(headX,lag.max=200,plot=F)
  
  par(lwd=cex/2)
  plot(a,main="Calcium data",cex.main=cex,cex.lab=cex,cex.axis=cex,xaxt="n",yaxt="n",xlab="",ylab="",ylim=yl)
  par(lwd=cex)
  title(xlab="Lag",cex=cex,cex.lab=cex,line=4)
  title(ylab="ACF",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,200,by=50),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,200,by=50),tick=F, labels=seq(0,200,by=50),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.12-0.2*1/3)
  axis(side=2, at=seq(0,1,by=0.2),tick=T, labels=rep("",6),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,by=0.2),tick=F, labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-0.33-0.2*0.66/3)
  
  par(cex.main=1.2)
  
  yl = c(-44,44)
  xl = c(-44,44)
  plot(headX[1:sublength],headY[1:sublength],type="l",xlab="",ylab="",xlim=xl,ylim=yl,cex.main=cex,cex.lab=cex,cex.axis=cex,xaxt="n",yaxt="n")
  title(xlab="X",cex=cex,cex.lab=cex,line=4)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(-40,40,by=20),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(-40,40,by=20),tick=F, labels=seq(-40,40,by=20),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-44-0.2*88/3)
  axis(side=2, at=seq(-40,40,by=20),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(-40,40,by=20),tick=F, labels=seq(-40,40,by=20),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-44-0.2*88/3)
  
  
  
  yl = c(-0.33,0.33)
  plot(bx[1:sublength],type="l",xlab="",ylab="",cex.main=cex,cex.lab=cex,cex.axis=cex,ylim=yl,xaxt="n")
  title(xlab="Time bins",cex=cex,cex.lab=cex,line=4)
  title(ylab="X",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,sublength,by=500),tick=T, labels=rep("",sublength/500+1),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,sublength,by=500),tick=F, labels=seq(0,sublength,by=500),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.33-0.2*0.66/3)
  axis(side=2, at=seq(-0.3,0.3,by=0.1),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(-0.3,0.3,by=0.1),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-0.33-0.2*0.66/3)
  
  yl = c(-0.1,1)
  par(cex.main=cex*1.2)
  a=acf(bx,lag.max=200,plot=F)
  par(lwd=cex/2)
  plot(a,main="Simulated data",cex.main=cex,cex.lab=cex,cex.axis=cex,xaxt="n",yaxt="n",xlab="",ylab="",ylim=yl)
  par(lwd=cex)
  title(xlab="Lag",cex=cex,cex.lab=cex,line=4)
  title(ylab="ACF",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,200,by=50),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,200,by=50),tick=F, labels=seq(0,200,by=50),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.12-0.2*1/3)
  axis(side=2, at=seq(0,1,by=0.2),tick=T, labels=rep("",6),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,by=0.2),tick=F, labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-0.33-0.2*0.66/3)
  
  par(cex.main=1.2)
  
  yl = c(-0.33,0.33)
  xl = c(-0.33,0.33)
  plot(bx[1:sublength],by[1:sublength],type="l",xlab="",ylab="",xlim=xl,ylim=yl,cex.main=cex,cex.lab=cex,cex.axis=cex,xaxt="n",yaxt="n")
  title(xlab="X",cex=cex,cex.lab=cex,line=4)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(-0.3,0.3,by=0.1),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(-0.3,0.3,by=0.1),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.33-0.2*0.66/3)
  axis(side=2, at=seq(-0.3,0.3,by=0.1),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(-0.3,0.3,by=0.1),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-0.33-0.2*0.66/3)
  
  
  par(lwd=parlwd)
}

# Figure 13
Simulation_illustration = function(cex=3){
  set.seed(1)
  parlwd=par("lwd")
  par(lwd=cex)
  b_proportion = 0.5
  N = 12000
  sd = 20
  bx = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  by = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  h = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  a = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  c = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  
  
  baserate = 0.03
  scale = 0.25
  rad = 0.06
  center = 0.15
  bpart = dnorm(bx,mean=center,sd=rad)*dnorm(by,mean=center,sd=rad)*(rad^2*2*pi) + dnorm(bx,mean=-center,sd=rad)*dnorm(by,mean=-center,sd=rad)*(rad^2*2*pi)
  maxb = max(bpart)
  if (maxb > 1){
    bpart = bpart/maxb
  }
  
  rad = 0.06
  center = 0.1
  hpart = dnorm(h,mean=center,sd=rad)*(rad*sqrt(2*pi))
  
  
  p_combined = baserate + scale*(bpart*b_proportion+hpart*(1-b_proportion))
  
  
  p = p_combined
  
  
  colorbar = hcl.colors(100,palette="SunsetDark")
  
  smoothingsd = c(0,0,1,1)
  y = rbinom(N,1,p)
  rm = ratemapXY(bpart,p,y,X=bx,Y=by,binning=40,smoothingsd = smoothingsd,binsize=1)
  image.plot(rm$rm1,main="p(x,y)",col=colorbar,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,axis.args=list(cex.axis=cex),xaxt='n',yaxt='n')
  title(xlab="X",cex=cex,cex.lab=cex,line=3.7)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,1,by=1/6),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,1,by=1/6),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.03)#,pos=-0.3-0.2*0.6/3)
  axis(side=2, at=seq(0,1,by=1/6),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,by=1/6),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  
  parmar = par("mar")
  par(mar=rep(parmar[1],4))
  plot(h,hpart,main="p(h)",ylab="",xlab="",cex.main=cex,cex.lab=cex,cex.axis=cex,xaxt="n",yaxt="n")
  title(xlab="h",cex=cex,cex.lab=cex,line=3.7)
  axis(side=1, at=seq(-0.3,0.3,by=0.1),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(-0.3,0.3,by=0.1),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.06)
  axis(side=2, at=seq(0,1,by=0.2),tick=T, labels=rep("",6),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,by=0.2),tick=F, labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  
  par(mar=parmar)
  
  image.plot(rm$rm2,main="p(x,y,h)",col=colorbar,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,axis.args=list(cex.axis=cex),xaxt='n',yaxt='n')
  title(xlab="X",cex=cex,cex.lab=cex,line=3.7)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,1,by=1/6),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,1,by=1/6),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.03)#,pos=-0.3-0.2*0.6/3)
  axis(side=2, at=seq(0,1,by=1/6),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,by=1/6),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  
  
  image.plot(rm$rm3,main="Smoothed sample y",col=colorbar,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,axis.args=list(cex.axis=cex),xaxt='n',yaxt='n')
  title(xlab="X",cex=cex,cex.lab=cex,line=3.7)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,1,by=1/6),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,1,by=1/6),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.03)#,pos=-0.3-0.2*0.6/3)
  axis(side=2, at=seq(0,1,by=1/6),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,by=1/6),tick=F, labels=c("-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  
  par(lwd=parlwd)
}

# Figure 14
Calcium_illustration = function(cex=3){
  parlwd=par("lwd")
  par(lwd=cex)
  trackdata2 = readRDS(paste(processed_data_folder,"trackdata2.RDS",sep=""))
  Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
  
  sublength=2000
  headX = trackdata2$headX
  headY = trackdata2$headY
  
  fakeys_ = as.matrix(Y1[1:12000,1:50])
  trueys_ = as.matrix(Y1[1:12000,51:100])
  
  truedat_ = trackdata2[1:12000,c(3,1,2,9)]
  fakedat_ = trackdata2[13051:25050,c(3,1,2,9)]
  image(trueys_[1:sublength,],main="Cell activity",cex=cex,cex.lab=cex,cex.axis=cex,cex.main=cex, xaxt="n", yaxt="n",xlab="",xaxt="n")
  title(xlab="Time bins",cex=cex,cex.lab=cex,line=4)
  axis(side=1, at=seq(0,sublength,by=500)/sublength,tick=T, labels=rep("",sublength/500+1),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,sublength,by=500)/sublength,tick=F, labels=seq(0,sublength,by=500),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.02)
  
  
  yl = c(-pi,pi)
  plot(truedat_[1:sublength,1],main="Head direction",type="l",xlab="",ylab="",cex.main=cex,cex.lab=cex,cex.axis=cex,ylim=yl,xaxt="n",yaxt="n")
  title(xlab="Time bins",cex=cex,cex.lab=cex,line=4)
  title(ylab="X",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,sublength,by=500),tick=T, labels=rep("",sublength/500+1),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,sublength,by=500),tick=F, labels=seq(0,sublength,by=500),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-pi-0.35)
  axis(side=2, at=seq(-3,3,by=1),tick=T, labels=rep("",7),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(-3,3,by=1),tick=F, labels=seq(-3,3,by=1),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-44-0.2*88/3)
  
  
  yl = c(-44,44)
  xl = c(-44,44)
  plot(headX[1:sublength],headY[1:sublength],type="l",xlab="",ylab="",xlim=xl,ylim=yl,cex.main=cex,cex.lab=cex,cex.axis=cex,xaxt="n",yaxt="n")
  title(xlab="X",cex=cex,cex.lab=cex,line=4)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(-40,40,by=20),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(-40,40,by=20),tick=F, labels=seq(-40,40,by=20),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-44-4.9)
  axis(side=2, at=seq(-40,40,by=20),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(-40,40,by=20),tick=F, labels=seq(-40,40,by=20),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)#,pos=-44-0.2*88/3)
  
  
  
  rm = ratemapXY(trueys_[,2],X=truedat_[,2],Y=truedat_[,3],binning=40,binsize=1/7.5,smoothingsd=1)
  quant = quantile(rm$rm1,0.95,na.rm=T)
  rm$rm1[which(rm$rm1 > quant)] = quant
  image(rm$rm1/quant,main="Smoothed ratemap",cex=cex,cex.lab=cex,cex.axis=cex,cex.main=cex,cex.sub=cex, xaxt="n", yaxt="n",col =colorbar)
  title(xlab="X",cex=cex,cex.lab=cex,line=3.8)
  title(ylab="Y",cex=cex,cex.lab=cex,line=3.2)
  axis(side=1, at=seq(0,1,length.out=5),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=1, at=seq(0,1,length.out=5),tick=F, labels=seq(-40,40,length.out=5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex,pos=-0.025)#,pos=-0.3-0.2*0.6/3)
  axis(side=2, at=seq(0,1,length.out=5),tick=T, labels=rep("",5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
  axis(side=2, at=seq(0,1,length.out=5),tick=F, labels=seq(-40,40,length.out=5),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)

  
  par(lwd=parlwd)
}

#### FIG1AB Nonsense correlations ----
set.seed(1)
sample_size = 200
num_samples = 10000
smoothing_parameter_a = 4 # TRY 2 or 3
smoothing_parameter_b = 12

# IID case
correlation_coefficients_iid = array(0,num_samples)
pvals_iid = array(0,num_samples)
for (i in 1:num_samples){
  sample1 = rnorm(sample_size)
  sample2 = rnorm(sample_size)
  correlation_coefficients_iid[i] = cor(sample1,sample2)
  pvals_iid[i] = cor.test(sample1,sample2,alternative="two.sided")$p.value
}

# Autocorrelated case
correlation_coefficients_ac_a = array(0,num_samples)
correlation_coefficients_ac_b = array(0,num_samples)
pvals_ac_a = array(0,num_samples)
pvals_ac_b = array(0,num_samples)
start = proc.time()[[3]]
for (i in 1:num_samples){
  samples_a = generate_X(sample_size,smoothing_kernel = "gaussian",smoothing_parameter = smoothing_parameter_a,ncov=2)
  samples_b = generate_X(sample_size,smoothing_kernel = "gaussian",smoothing_parameter = smoothing_parameter_b,ncov=2)
  correlation_coefficients_ac_a[i] = cor(samples_a[,1],samples_a[,2])
  correlation_coefficients_ac_b[i] = cor(samples_b[,1],samples_b[,2])
  pvals_ac_a[i] = cor.test(samples_a[,1],samples_a[,2],alternative="two.sided")$p.value
  pvals_ac_b[i] = cor.test(samples_b[,1],samples_b[,2],alternative="two.sided")$p.value
  print_progress(i,num_samples,start)
}

saveRDS(correlation_coefficients_iid,paste(figure_folder,"correlation_coefficients_iid.RDS",sep=""))
saveRDS(pvals_iid,paste(figure_folder,"pvals_iid.RDS",sep=""))

saveRDS(correlation_coefficients_ac_a,paste(figure_folder,"correlation_coefficients_ac_a.RDS",sep=""))
saveRDS(pvals_ac_a,paste(figure_folder,"pvals_ac_a.RDS",sep=""))

saveRDS(correlation_coefficients_ac_b,paste(figure_folder,"correlation_coefficients_ac_b.RDS",sep=""))
saveRDS(pvals_ac_b,paste(figure_folder,"pvals_ac_b.RDS",sep=""))


#### FIG1CD Nonsense LRT for GLM ----
nsim=10000
set.seed(1)
smoothing_parameter_a = 0
smoothing_parameter_b = 4
smoothing_parameter_c = 8
smoothing_kernel="gaussian"
n=200
coef1=0.5
coef2=0.25
base=5
LRpvals = matrix(0,nrow=4,ncol=nsim)
lrs = matrix(0,nrow=4,ncol=nsim)
for (i in 1:nsim){
  if (i %% round(nsim/10) == 0){
    print(paste(i,"/",nsim,sep=""))
  }
  X_a = generate_X(n,smoothing_kernel,0)
  X_b = generate_X(n,smoothing_kernel,0)
  X_c = generate_X(n,smoothing_kernel,smoothing_parameter_c)
  X_d = generate_X(n,smoothing_kernel,smoothing_parameter_c)
  Y_a = generate_Y(n,X_a,coef1=coef1,coef2=0,base=base)
  Y_b = generate_Y(n,X_b,coef1=coef1,coef2=coef2,base=base)
  Y_c = generate_Y(n,X_c,coef1=coef1,coef2=0,base=base)
  Y_d = generate_Y(n,X_d,coef1=coef1,coef2=coef2,base=base)
  
  fits_a = list()
  fits_a$fit1 = glm(Y_a~X_a[,1],family=poisson)
  fits_a$fit13 = glm(Y_a~X_a[,1]+X_a[,3],family=poisson)
  fits_b = list()
  fits_b$fit1 = glm(Y_b~X_b[,1],family=poisson)
  fits_b$fit13 = glm(Y_b~X_b[,1]+X_b[,3],family=poisson)
  fits_c = list()
  fits_c$fit1 = glm(Y_c~X_c[,1],family=poisson)
  fits_c$fit13 = glm(Y_c~X_c[,1]+X_c[,3],family=poisson)
  fits_d = list()
  fits_d$fit1 = glm(Y_d~X_d[,1],family=poisson)
  fits_d$fit13 = glm(Y_d~X_d[,1]+X_d[,3],family=poisson)
  
 
  lrt_a = lrtest(fits_a$fit13,fits_a$fit1)
  lrt_b = lrtest(fits_b$fit13,fits_b$fit1)
  lrt_c = lrtest(fits_c$fit13,fits_c$fit1)
  lrt_d = lrtest(fits_d$fit13,fits_d$fit1)
  
  lrs[1,i] = -2*diff(lrt_a$LogLik)
  lrs[2,i] = -2*diff(lrt_b$LogLik)
  lrs[3,i] = -2*diff(lrt_c$LogLik)
  lrs[4,i] = -2*diff(lrt_d$LogLik)
  
  LRpvals[1,i] = lrt_a$`Pr(>Chisq)`[2] #IID, no hidden
  LRpvals[2,i] = lrt_b$`Pr(>Chisq)`[2] #IID, but hidden variable
  LRpvals[3,i] = lrt_c$`Pr(>Chisq)`[2] #Autocorr, no hidden
  LRpvals[4,i] = lrt_d$`Pr(>Chisq)`[2] #Autocorr, and hidden
  
}

saveRDS(lrs,paste(figure_folder,"lrs.RDS",sep=""))
saveRDS(LRpvals,paste(figure_folder,"LRpvals.RDS",sep=""))



#### FIG1E Cross-validation scheme error rates ----
set.seed(1)

num_trials = 300
N = 12000
sd = 20
bprop = 0.5
alpha = 0.05
scale = 0.25 

hdknots_sim = 5
posknots_sim = 2
speknots_sim = 5
knots = c(hdknots_sim,posknots_sim,posknots_sim,speknots_sim)

ncov = 3
covinds_sim = list(`1`=1:(hdknots_sim+1),`2`=(1+hdknots_sim+1):(1+hdknots_sim+(1+posknots_sim)^2),`3`=(2+hdknots_sim+(1+posknots_sim)^2):(2+hdknots_sim+(1+posknots_sim)^2+speknots_sim))


if (TRUE){
  cv.config = config.default()
  cv.config$alpha = alpha
  cv.config$family = "binomial"
} # Set configs

fn = c(20)
cs = c(2,4,8,15,30,60,150,300,600)
ns = c(TRUE,FALSE)
RANDOMFOLDS = FALSE #Random folds when chunksize = 1
RANDOMFOLDSONLAST = TRUE #The above is true for the last variant

fold_nums = c()
chunksizes = c()
neighbor_skips = c()

for (f in fn){
  for (c in cs){
    for (s in ns){
      fold_nums = c(fold_nums,f)
      chunksizes = c(chunksizes,N/f/c)
      neighbor_skips = c(neighbor_skips,s)
    }
  }
}

fold_nums = c(fold_nums,20)
chunksizes = c(chunksizes,1)
neighbor_skips = c(neighbor_skips,FALSE)

num_variants = length(fold_nums)

fakeresmat = matrix(NA,nrow=num_trials*num_variants,ncol=7)
variant = rep(1:num_variants,num_trials)
fakeresmat[,7] = variant

acfs = matrix(NA,nrow=num_trials,ncol=300)
acfsy = matrix(NA,nrow=num_trials,ncol=300)

start = proc.time()[[3]]
for (i in 1:num_trials){
  fakedat = make_simulated_data(N,sd,ncov=ncov,scale=scale) 
  acfs[i,] = acf(fakedat$a,lag.max=300,plot=F)$acf[2:301]
  acfsy[i,] = acf(fakedat$y,lag.max=300,plot=F)$acf[2:301]
  
  fakeX = make_design_matrix(fakedat[,2:dim(fakedat)[2]],c(),c(),c(2),knots)
  
  fakeresmat[((i-1)*num_variants+1):(i*num_variants),1:6] = forward_selection_CV_variants(Y=fakedat$y,X=fakeX,
                                                                                                     cov_inds=covinds_sim,fold_nums,chunksizes,neighbor_skips,RANDOMFOLDS=RANDOMFOLDS,config=cv.config,RANDOMFOLDSONLAST=RANDOMFOLDSONLAST)[,1:6]
  
  print_progress(i,num_trials,start)
}


saveRDS(fakeresmat,paste(result_folder,"CV_variants_no_effect.RDS",sep=""))
#
#### FIG1F SR error rate for different blocksizes ----
set.seed(1)

num_trials = 300
N = 12000
sd = 20 
bprop = 0.5
alpha = 0.05
scale = 0.25 

hdknots_sim = 5
posknots_sim = 2
speknots_sim = 5
knots = c(hdknots_sim,posknots_sim,posknots_sim,speknots_sim)

ncov = 3
covinds_sim = list(`1`=1:(hdknots_sim+1),`2`=(1+hdknots_sim+1):(1+hdknots_sim+(1+posknots_sim)^2),`3`=(2+hdknots_sim+(1+posknots_sim)^2):(2+hdknots_sim+(1+posknots_sim)^2+speknots_sim))


if (TRUE){
  cv.config = config.default()
  cv.config$alpha = alpha
  cv.config$family = "binomial"
} # Set configs

fn = c(20)
cs = c(2,4,8,15,30,60,150,300,600)
ns = c(TRUE,FALSE)
RANDOMFOLDS = FALSE #Random folds when chunksize = 1
RANDOMFOLDSONLAST = TRUE #The above is true for the last variant

fold_nums = c()
chunksizes = c()
neighbor_skips = c()

for (f in fn){
  for (c in cs){
    for (s in ns){
      fold_nums = c(fold_nums,f)
      chunksizes = c(chunksizes,N/f/c)
      neighbor_skips = c(neighbor_skips,s)
    }
  }
}

fold_nums = c(fold_nums,20)
chunksizes = c(chunksizes,1)
neighbor_skips = c(neighbor_skips,FALSE)

num_variants = length(fold_nums)

fakeresmat = matrix(NA,nrow=num_trials*num_variants,ncol=7)
variant = rep(1:num_variants,num_trials)
fakeresmat[,7] = variant

acfs = matrix(NA,nrow=num_trials,ncol=300)
acfsy = matrix(NA,nrow=num_trials,ncol=300)

start = proc.time()[[3]]
for (i in 1:num_trials){
  fakedat = make_simulated_data(N,sd,ncov=ncov,scale=scale) 
  acfs[i,] = acf(fakedat$a,lag.max=300,plot=F)$acf[2:301]
  acfsy[i,] = acf(fakedat$y,lag.max=300,plot=F)$acf[2:301]
  
  fakeX = make_design_matrix(fakedat[,2:dim(fakedat)[2]],c(),c(),c(2),knots)
  
  fakeresmat[((i-1)*num_variants+1):(i*num_variants),1:6] = forward_selection_CV_variants(Y=fakedat$y,X=fakeX,
                                                                                                     cov_inds=covinds_sim,fold_nums,chunksizes,neighbor_skips,RANDOMFOLDS=RANDOMFOLDS,config=cv.config,RANDOMFOLDSONLAST=RANDOMFOLDSONLAST,BUTSIGNEDRANK=TRUE)[,1:6]
  
  print_progress(i,num_trials,start)
}


saveRDS(fakeresmat,paste(result_folder,"CV_SR_variants_no_effect.RDS",sep=""))
#



#### FIG2 / FIG 11 Multiple correction stuff ----
num_trials = 10000
ncov = 5
#m = 30
m = 100
sd = 0
nperm = 999

pvals_single = array(0,num_trials)
pvals_max = array(0,num_trials)
pvals_max_corr = array(0,num_trials)

start = proc.time()[[3]]
for (i in 1:num_trials){
  y = rnorm(m)
  x = matrix(rnorm(m*ncov),nrow=m,ncol=ncov)
  teststats = array(0,ncov)
  pvals = array(0,ncov)
  #e = rnorm(m,sd=sd)
  
  for (j in 1:ncov){
    #teststats[j] = mean(x[,j]-y)
    teststats[j] = signrank(x[,j]-y)
  }
  
  permstat = matrix(0,nrow=nperm,ncol=ncov)
  for (k in 1:nperm){
    u = round(runif(m))*2-1
    
    for (j in 1:ncov){
      #w = wilcox.test(x=y,y=x[,j]+e,paired=TRUE,alternative = "less")
      #teststats[j] = w$statistic
      #pvals[j] = w$p.value
      #permstat[k,j] = mean(u*(x[,j]-y),na.rm=T)
      permstat[k,j] = signrank(u*(x[,j]-y))
    }
  }
  maxineach = apply(permstat,1,max)
  
  for (j in 1:ncov){
    pvals[j] = (length(which(permstat[,j] >= teststats[j]))+1)/(nperm+1)
  }
  
  pvals_max_corr[i] = (length(which(maxineach >= teststats[which.max(teststats)]))+1)/(nperm+1)
  
  pvals_single[i] = pvals[1]
  pvals_max[i] = min(c(1,pvals[which.max(teststats)]))
  
  print_progress(i,num_trials,start)
}



saveRDS(pvals_single,paste(figure_folder,"pvals_single.RDS",sep=""))
saveRDS(pvals_max,paste(figure_folder,"pvals_max.RDS",sep=""))
saveRDS(pvals_max_corr,paste(figure_folder,"pvals_max_corr.RDS",sep=""))

#
#### Saving jpegs ----
basewidth = 12
baseheight = 7
mar = c(5,5,5,5)
#mar = c(4.5,5,1,0.5)


jpeg(paste(figure_folder,"jpegs/blank.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
plot.new()
dev.off()
#Singles

# Fig 1A
jpeg(paste(figure_folder,"jpegs/Nonsense_corr.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Nonsense_corr()
dev.off()

# Fig 1B
jpeg(paste(figure_folder,"jpegs/Nonsense_corr_pval.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Nonsense_corr_pval()
dev.off()

# Fig 1C
jpeg(paste(figure_folder,"jpegs/Nonsense_dev.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Nonsense_dev()
dev.off()

# Fig 1D
jpeg(paste(figure_folder,"jpegs/Nonsense_dev_pval.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Nonsense_dev_pval()
dev.off()

# Fig 1E
jpeg(paste(figure_folder,"jpegs/CV_error_rates.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
CV_error_rates()
dev.off()

# Fig 1F
jpeg(paste(figure_folder,"jpegs/SR_error_rates.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
SR_error_rates()
dev.off()

# Fig 1G
jpeg(paste(figure_folder,"jpegs/Multiple_comparisons.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Multiple_comparisons()
dev.off()

# Fig 2
jpeg(paste(figure_folder,"jpegs/Hypothetical_p_distributions.jpeq",sep=""),
     units="in", width=basewidth*2, height=baseheight, res=300)
par(mfrow=c(1,3),mar=mar)
Hypothetical_p_distributions()
dev.off()


# Figure 3, Simulated false results
jpeg(paste(figure_folder,"jpegs/Simulated_false_results_logscale.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9*2, res=300)
par(mfrow=c(2,1),mar=mar)
Simulated_false_results(logscale=T)
dev.off()
jpeg(paste(figure_folder,"jpegs/Simulated_false_results.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9*2, res=300)
par(mfrow=c(2,1),mar=mar)
Simulated_false_results(logscale=F)
dev.off()

# Figure 4, Simulated true results
jpeg(paste(figure_folder,"jpegs/Simulated_true_results_logscale.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9, res=300)
par(mfrow=c(1,1),mar=mar)
Simulated_true_results(logscale=T)
dev.off()
jpeg(paste(figure_folder,"jpegs/Simulated_true_results.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9, res=300)
par(mfrow=c(1,1),mar=mar)
Simulated_true_results(logscale=F)
dev.off()

# Figure 5, Calcium false results
jpeg(paste(figure_folder,"jpegs/Calcium_false_results_logscale.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9*2, res=300)
par(mfrow=c(2,1),mar=mar)
Calcium_false_results(logscale=T)
dev.off()
jpeg(paste(figure_folder,"jpegs/Calcium_false_results.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9*2, res=300)
par(mfrow=c(2,1),mar=mar)
Calcium_false_results(logscale=F)
dev.off()

# Figure 6, Calcium true results
jpeg(paste(figure_folder,"jpegs/Calcium_true_results_logscale.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9, res=300)
par(mfrow=c(1,1),mar=mar)
Calcium_true_results(logscale=T)
dev.off()
jpeg(paste(figure_folder,"jpegs/Calcium_true_results.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*0.9, res=300)
par(mfrow=c(1,1),mar=mar)
Calcium_true_results(logscale=F)
dev.off()

# Fig 7
jpeg(paste(figure_folder,"jpegs/Ratemaps.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*1.7, res=300)
par(mfrow=c(1,3),mar=c(5, 5, 5, 5))
Ratemaps()
dev.off()

jpeg(paste(figure_folder,"jpegs/Colorbar.jpeq",sep=""),
     units="in", width=basewidth/3, height=baseheight*1.8, res=300)
par(mar=c(4, 3, 4, 3))
Colorbar()
dev.off()

# Fig 8
jpeg(paste(figure_folder,"jpegs/Venn_diagrams.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*2*0.9, res=300)
par(mar=mar)
Venn_diagrams()
dev.off()

# Fig 9
jpeg(paste(figure_folder,"jpegs/Blocking_illustration.jpeq",sep=""),
     units="in", width=basewidth*2, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Blocking_illustration()
dev.off()

# Fig 10
jpeg(paste(figure_folder,"jpegs/Cyclical_illustration.jpeq",sep=""),
     units="in", width=basewidth*2, height=baseheight, res=300)
par(mfrow=c(1,1),mar=mar)
Cyclical_illustration()
dev.off()

# Fig 11
jpeg(paste(figure_folder,"jpegs/Multiple_corrections.jpeq",sep=""),
     units="in", width=basewidth*2, height=baseheight, res=300)
par(mfrow=c(1,3),mar=mar)
Multiple_corrections()
dev.off()

# Fig 12
jpeg(paste(figure_folder,"jpegs/Data_similarities.jpeq",sep=""),
     units="in", width=basewidth*2, height=baseheight*1.5, res=300)
par(mfrow=c(2,3),mar=c(6,6,5,5))
Data_similarities()
dev.off()

# Fig 13
jpeg(paste(figure_folder,"jpegs/Simulation_illustration.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*1.16, res=300)
par(mfrow=c(1,4),mar=c(5,6,5,11))
Simulation_illustration()
dev.off()

# Fig 14
jpeg(paste(figure_folder,"jpegs/Calcium_illustration.jpeq",sep=""),
     units="in", width=basewidth*3, height=baseheight*1.22, res=300)
par(mfrow=c(1,4),mar=c(5,6,5,6))
Calcium_illustration()
dev.off()

# Two doubles
jpeg(paste(figure_folder,"jpegs/Nonsense_corr_double2.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight*2, res=600)
par(mfrow=c(2,1),mar=mar)
Nonsense_corr()
Nonsense_corr_pval()
dev.off()

jpeg(paste(figure_folder,"jpegs/Nonsense_dev_double.jpeq",sep=""),
     units="in", width=basewidth, height=baseheight*2, res=300)
par(mfrow=c(2,1),mar=mar)
Nonsense_dev()
Nonsense_dev_pval()
dev.off()

