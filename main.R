#### Info ----

# Main script for loading data from Zong et. al. (2022), 
# processing it to the format we desire,
# and running different forward selection procedures to compare their properties

#### Libraries ----
library(splines)
library(pbs)
library(MGLM)
library(fields)
library(lmtest)
library(MASS)
library(mvtnorm)
library(ggbiplot)
library(fields)
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
library(umap)
library(glmnet)
library(EnvStats)
library(rhdf5)
library(GenBinomApps)
library(ggvenn)
library(vioplot)

#result_folder = "/Users/fredrine/Documents/CalciumData/Result_data/"
result_folder = "/Users/fredrine/Documents/CalciumData/Result_data_revision/"
#result_folder = "/Users/fredrine/Desktop/testing/"
functions_folder = "/Users/fredrine/Documents/CalciumData/Scripts/"
processed_data_folder = "/Users/fredrine/Documents/CalciumData/Processed_data/"
#processed_data_folder = "/Users/fredrine/Desktop/testing/"
data_folder = "/Users/fredrine/Documents/CalciumData/Prepped data/"
#### Load functions ----
source(paste(functions_folder,"functions_model_selection.R",sep=""))
#### Load data for rat 97045, date 20210317----
SAVE = FALSE

animal =  97045
date = "20210317" 
numbins = 26100
numtrackbins = 52200
filepath = paste(data_folder,animal,"/",date,"/","NAT.mat",sep="")
h5ls(filepath) # Find what to use as name below:
NAT = h5read(filepath,name="/#refs#/b")
trackdata = NAT[,c(2,3,4,7,8,9,6,5,10)]
filteredevents = seq(16,dim(NAT)[2],by=4)
neural = matrix(NA,nrow=dim(NAT)[1]/2,ncol=length(filteredevents))
for (i in 1:length(filteredevents)){
  neural[,i] = BinMean(NAT[,filteredevents[i]],2,na.rm=T)
}
NAT = NULL
neural= data.frame(neural)
colnames(neural) = paste("cell",1:(dim(neural)[2]),sep="")
neuralsum = apply(neural,2,sum)

otherfilepath = paste(data_folder,animal,"/",date,"/","NeuronInformation.mat",sep="")
neuroninfo = read.mat(otherfilepath) # Find what to use as name below:
repeats = unlist(neuroninfo$NeuronInformation$RepeatCell)

keepers = sort(setdiff(1:(dim(neural)[2]),union(repeats,which(is.na(neural[1,])))))
newneural = neural[,keepers]
for (i in 1:dim(newneural)[2]){
  newneural[which(newneural[,i] != 0),i] = 1
}

colnames(newneural) = paste("cell",keepers,sep="")

if (TRUE){
  gridcells = c(2,3	,4	,5	,6	,7,	9,	11,	14,	15,	17,	18,	21,	22,	25,	26,	28,	30,	33,	34,	40,	43,	45,	46,	
                47,	48,	49,	50,	51,	52,	53,	54,	55,	58,	60,	62,	63,	65,	67,	75,	76,	79,	81,	83,	84,
                85,	86,	89,	90,	91,	97,	98,	100,	101,	105,	109,	110,	112,	114,	116,	117,	123,	124,	127,	128,	130,	
                131,	132,	139,	140,	141,	142,	149,	150,	151,	154,	156,	160,	162,	165,	166,	173,	178,	181,	184,	188,	
                191,	192,	197,	198,	199,	200,	201,	203,	206,	207,	208,	210,	213,	214,	215,	216,	218,	221,	224,	225,	
                226,	227,	228,	232,	234,	235,	237,	239,	241,	242	,243	,244	,245	,246	,247,	248,	251,	252,	253,	254	,
                255,	256,	257,	258,	263,	265,	266,	267,	268,	269,	270,	272,	273,	274,	275,	276,	277,	278	,281,	282,	
                283,	285,	286,	288,	291,	292,	295,	296,	297,	298	,300	,305	,306,	307,	308,	309,	311,	312,	314,	315	,
                317,	322,	326,	327,	330,	331,	333,	334,	335,	338,	340,	346,	347,	348,	349,	350,	353,	354,	356,	358,	
                363,	365,	369,	370,	371,	374,	378,	382,	385,	388,	389,	399,	409,	410,	430,	432,	435,	437,	438,	441,	
                444,	445 )
}## Grid cells for 97045, 20210317

trackdata = data.frame(trackdata)
colnames(trackdata) = c("headX","headY","HD","bodyX","bodyY","BD","HBangle","headspeed","bodyspeed")
trackdata$HD = trackdata$HD/180*pi - pi
trackdata$BD = trackdata$BD/180*pi - pi

# Increasing the binsize for the tracked data for the two halves separately
trackdata2 = rbind(make_data(NULL,1/7.25,covariates=trackdata[1:26100,],cov_binsize=1/14.5,onlycount=FALSE,periodiccovs = c(3,6),max_T=NULL)[2:10],
                   make_data(NULL,1/7.25,covariates=trackdata[26101:52200,],cov_binsize=1/14.5,onlycount=FALSE,periodiccovs = c(3,6),max_T=NULL)[2:10])

# Saving the trackingdata
if (SAVE){
  saveRDS(trackdata2,paste(processed_data_folder,"trackdata2.RDS",sep=""))
}
#### Load data for rat 97045, date 20210314----
SAVE = FALSE

animal =  97045
date = "20210314" 
numbins = 26100
numtrackbins = 52200
filepath = paste(data_folder,animal,"/",date,"/","NAT.mat",sep="")
h5ls(filepath) # Find what to use as name below:
NAT = h5read(filepath,name="/#refs#/b")
trackdata = NAT[,c(2,3,4,7,8,9,6,5,10)]
filteredevents = seq(16,dim(NAT)[2],by=4)
neural = matrix(NA,nrow=dim(NAT)[1]/2,ncol=length(filteredevents))
for (i in 1:length(filteredevents)){
  neural[,i] = BinMean(NAT[,filteredevents[i]],2,na.rm=T)
}
NAT = NULL
neural= data.frame(neural)
colnames(neural) = paste("cell",1:(dim(neural)[2]),sep="")
neuralsum = apply(neural,2,sum)

otherfilepath = paste(data_folder,animal,"/",date,"/","NeuronInformation.mat",sep="")
neuroninfo = read.mat(otherfilepath) # Find what to use as name below:
repeats = unlist(neuroninfo$NeuronInformation$RepeatCell)

keepers = sort(setdiff(1:(dim(neural)[2]),union(repeats,which(is.na(neural[1,])))))
newneural = neural[,keepers]
for (i in 1:dim(newneural)[2]){
  newneural[which(newneural[,i] != 0),i] = 1
}

colnames(newneural) = paste("cell",keepers,sep="")

trackdata = data.frame(trackdata)
colnames(trackdata) = c("headX","headY","HD","bodyX","bodyY","BD","HBangle","headspeed","bodyspeed")
trackdata$HD = trackdata$HD/180*pi - pi
trackdata$BD = trackdata$BD/180*pi - pi

# Increasing the binsize for the tracked data for the two halves separately
trackdata2 = rbind(make_data(NULL,1/7.25,covariates=trackdata[1:26100,],cov_binsize=1/14.5,onlycount=FALSE,periodiccovs = c(3,6),max_T=NULL)[2:10],
                   make_data(NULL,1/7.25,covariates=trackdata[26101:52200,],cov_binsize=1/14.5,onlycount=FALSE,periodiccovs = c(3,6),max_T=NULL)[2:10])

# Saving the trackingdata
if (SAVE){
  saveRDS(trackdata2,paste(processed_data_folder,"trackdata2_other_date.RDS",sep=""))
}
#### Constructing the design matrix used in the GLM ----
SAVE = FALSE

hdknots = 6
posknots = 4
speknots = 5
datalength = 12000
halfwaypoint = 13050
spikethresh = 0.02 # Discarding cells with events in fewer than 2% of the bins (event rate < 0.15)

X1 = make_design_matrix(trackdata2[1:datalength,c(3,1,2,9)],c(4),c(1),c(2),c(hdknots,posknots,posknots,speknots))
X2 = make_design_matrix(trackdata2[(halfwaypoint+1):(halfwaypoint+datalength),c(3,1,2,9)],c(4),c(1),c(2),c(hdknots,posknots,posknots,speknots))
covinds = list(`1`=1:hdknots,`2`=(hdknots+1):(hdknots+(1+posknots)^2),`3`=(1+hdknots+(1+posknots)^2):dim(X1)[2])

Y1 = newneural[1:datalength,]
Y2 = newneural[(halfwaypoint+1):(halfwaypoint+datalength),]

spikemean1 = apply(Y1,2,mean)
spikemean2 = apply(Y2,2,mean)

keepers2 = setdiff(1:dim(Y1)[2],union(which(spikemean1<spikethresh),which(spikemean2<spikethresh)))

Y1 = Y1[,keepers2]
Y2 = Y2[,keepers2]
colnames(Y1) = colnames(newneural)[keepers2]
colnames(Y2) = colnames(newneural)[keepers2]

cellstokeep = keepers[keepers2]

gridcellskept = which(cellstokeep %in% gridcells)

# Saving the design matrices and the set of grid cells
if (SAVE){
  saveRDS(gridcellskept,paste(processed_data_folder,"gridcells.RDS",sep=""))
  #saveRDS(Y1,paste(processed_data_folder,"Y1_other_date.RDS",sep=""))
  saveRDS(Y1,paste(processed_data_folder,"Y1.RDS",sep=""))
  #saveRDS(Y2,paste(processed_data_folder,"Y2_other_date.RDS",sep=""))
  saveRDS(Y2,paste(processed_data_folder,"Y2.RDS",sep=""))
  saveRDS(X1,paste(processed_data_folder,"X1.RDS",sep=""))
  saveRDS(X2,paste(processed_data_folder,"X2.RDS",sep=""))
  saveRDS(covinds,paste(processed_data_folder,"covinds.RDS",sep=""))
  saveRDS(cellstokeep,paste(processed_data_folder,"cellstokeep.RDS",sep=""))
}
#### Make and save simulated data ----
SAVE = T

set.seed(1)
numcells = 300
N = 12000
Y1 = make_simulated_data_from_real_covariates(num_cells=numcells)
Y2 = make_simulated_data_from_real_covariates(num_cells=numcells, HIDDEN_ONLY=T)

hdknots_sim = 4
posknots_sim = 2
speknots_sim = 3
knots = c(hdknots_sim,posknots_sim,posknots_sim,speknots_sim)
covinds_sim = list(`1`=1:hdknots_sim,`2`=(hdknots_sim+1):(hdknots_sim+(1+posknots_sim)^2),`3`=(1+hdknots_sim+(1+posknots_sim)^2):(1+hdknots_sim+(1+posknots_sim)^2+speknots_sim))

trackdata2 = readRDS(paste(processed_data_folder,"trackdata2.RDS",sep=""))
X = make_design_matrix(trackdata2[1:N,c(3,1,2,9)],c(4),c(1),c(2),c(hdknots_sim,posknots_sim,posknots_sim,speknots_sim))



if (SAVE){
  saveRDS(Y1,paste(processed_data_folder,"Simulated_from_real_cov_Y1.RDS",sep=""))
  saveRDS(Y2,paste(processed_data_folder,"Simulated_from_real_cov_Y2.RDS",sep=""))
  saveRDS(X,paste(processed_data_folder,"Simulated_from_real_cov_X.RDS",sep=""))
  saveRDS(covinds_sim,paste(processed_data_folder,"covinds_sim_new.RDS",sep=""))
}

#### Make and save config-files ----

config.20.folds = config.default()
config.20.folds = prepare_folds(config.20.folds)

config.10.folds = config.default()
config.10.folds$nfolds=10
config.10.folds$SKIP_NEIGHBOR=F
config.10.folds = prepare_folds(config.10.folds)

saveRDS(config.20.folds,paste(processed_data_folder,"config_20_folds.RDS",sep=""))
saveRDS(config.10.folds,paste(processed_data_folder,"config_10_folds.RDS",sep=""))

#### Forward selection on simulated data, 20 folds  ----
#hdknots = 4
#posknots = 2
#speknots = 3
#fam = "binomial" 
#N = 12000
SAVE = F

X = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_X.RDS",sep=""))
Y1 = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_Y2.RDS",sep=""))
covinds_sim = readRDS(paste(processed_data_folder,"covinds_sim_new.RDS",sep=""))

set.seed(1)
numcells = dim(Y1)[2]
config.20.folds = readRDS(paste(processed_data_folder,"config_20_folds.RDS",sep=""))

methods = c("CV","signrankmax","signrankmax2","cyclicskip","lrtbonf")
bettermethodnames =c("CV","mSR","mSRR","CS","LRT")
TrueResultMat = array(0, dim=c(length(methods) * 2 + 1, config.20.folds$maxstep,numcells))
FakeResultMat = array(0, dim=c(length(methods) * 2 + 1, config.20.folds$maxstep,numcells))

checkpoints = c(50,100,150,200,250)
start = proc.time()[[3]]
for (i in 1:numcells){
  TrueResultMat[,,i] = forward_selection(Y=Y1[,i],X=X,cov_inds=covinds_sim,methods=methods,config=config.20.folds)
  FakeResultMat[,,i] = forward_selection(Y=Y2[,i],X=X,cov_inds=covinds_sim,methods=methods,config=config.20.folds,FAKE=T)
  
  print_progress(i,numcells,start)
  if (i %in% checkpoints){
    plot_results(TrueResultMat[,,1:i],bettermethodnames,0.05,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
    plot_results(FakeResultMat[,,1:i],bettermethodnames,0.05,hline=T,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
  }
}

if (SAVE){
  saveRDS(TrueResultMat,paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
  saveRDS(FakeResultMat,paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))
}

#TrueResultMat = readRDS(paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
#FakeResultMat = readRDS(paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))

plot_results(TrueResultMat,bettermethodnames,0.05,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
plot_results(FakeResultMat,bettermethodnames,0.05,hline=T,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))


#### Forward selection on simulated data, 10 folds  ----
#hdknots = 4
#posknots = 2
#speknots = 3
#fam = "binomial" 
#N = 12000
SAVE = T

X = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_X.RDS",sep=""))
Y1 = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_Y2.RDS",sep=""))
covinds_sim = readRDS(paste(processed_data_folder,"covinds_sim_new.RDS",sep=""))

set.seed(1)
numcells = dim(Y1)[2]
config.10.folds = readRDS(paste(processed_data_folder,"config_10_folds.RDS",sep=""))

methods = c("CV","signrank")
bettermethodnames =c("CV","SR")
TrueResultMat2 = array(0, dim=c(length(methods) * 2 + 1, config.10.folds$maxstep,numcells))
FakeResultMat2 = array(0, dim=c(length(methods) * 2 + 1, config.10.folds$maxstep,numcells))


checkpoints = c(50,100,150,200,250)
start = proc.time()[[3]]
for (i in 1:numcells){
  TrueResultMat2[,,i] = forward_selection(Y=Y1[,i],X=X,cov_inds=covinds_sim,methods=methods,config=config.10.folds)
  FakeResultMat2[,,i] = forward_selection(Y=Y2[,i],X=X,cov_inds=covinds_sim,methods=methods,config=config.10.folds,FAKE=T)
  
  print_progress(i,numcells,start)
  
  if (i %in% checkpoints){
    plot_results(TrueResultMat2[,,1:i],bettermethodnames,0.05,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
    plot_results(FakeResultMat2[,,1:i],bettermethodnames,0.05,hline=T,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
  }
}

if (SAVE){
  saveRDS(TrueResultMat2,paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
  saveRDS(FakeResultMat2,paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))
}

#TrueResultMat2 = readRDS(paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
#FakeResultMat2 = readRDS(paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))

plot_results(TrueResultMat2,bettermethodnames,0.05,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
plot_results(FakeResultMat2,bettermethodnames,0.05,hline=T,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))


#### Forward selection on calcium data, 20 folds ----
#hdknots = 6
#posknots = 4
#speknots = 5
#fam = "binomial" 
#datalength = 12000
#spikethresh = 0.02

SAVE = T

Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Y2.RDS",sep=""))
X1 = readRDS(paste(processed_data_folder,"X1.RDS",sep=""))
X2 = readRDS(paste(processed_data_folder,"X2.RDS",sep=""))
covinds = readRDS(paste(processed_data_folder,"covinds.RDS",sep=""))

numcells = dim(Y1)[2]
methods = c("CV","signrankmax","signrankmax2","cyclicskip","lrtbonf")
bettermethodnames = c("CV","mSR","mSRR","CS","LRT")

config.20.folds = readRDS(paste(processed_data_folder,"config_20_folds.RDS",sep=""))

TrueResults = array(0, dim=c(length(methods) * 2 + 1, config.20.folds$maxstep,numcells))
FakeResults = array(0, dim=c(length(methods) * 2 + 1, config.20.folds$maxstep,numcells))

set.seed(1)

checkpoints = c(50,100,150,200)
start = proc.time()[[3]]
for (i in 1:numcells){
  tryCatch({
    TrueResults[,,i] = forward_selection(Y=Y1[,i],X=X1,cov_inds=covinds,methods=methods,
                                         config=config.20.folds)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    FakeResults[,,i] = forward_selection(Y=Y2[,i],X=X1,cov_inds=covinds,methods=methods,
                                         config=config.20.folds,FAKE=T)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  print_progress(i,numcells,start)
  
  if (i %in% checkpoints){
    plot_results(TrueResults[,,1:i],bettermethodnames,0.05,cex=2,nstep=3)
    plot_results(FakeResults[,,1:i],bettermethodnames,0.05,hline = T,cex=2,nstep=3)
  }
}

if (SAVE){
  saveRDS(TrueResults,paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
  saveRDS(FakeResults,paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))
}


#### Plotting results from calcium data, 20 folds ----
TrueResults = readRDS(paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
FakeResults = readRDS(paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))

alpha = 0.05
bettermethodnames = c("CV","mSR","mSRR","CS","LRT")


plot_results(TrueResults,bettermethodnames,alpha,cex=2,nstep=3)
plot_results(FakeResults,bettermethodnames,alpha,hline = T,cex=2,nstep=3)

#
#### Forward selection on calcium data, 10 folds ----
#hdknots = 6
#posknots = 4
#speknots = 5
#fam = "binomial" 
#datalength = 12000
#spikethresh = 0.02

SAVE = T

Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Y2.RDS",sep=""))
X1 = readRDS(paste(processed_data_folder,"X1.RDS",sep=""))
X2 = readRDS(paste(processed_data_folder,"X2.RDS",sep=""))
covinds = readRDS(paste(processed_data_folder,"covinds.RDS",sep=""))

numcells = dim(Y1)[2]
methods = c("CV","signrank")
bettermethodnames = c("CV","SR")

config.10.folds = readRDS(paste(processed_data_folder,"config_10_folds.RDS",sep=""))

TrueResults2 = array(0, dim=c(length(methods) * 2 + 1, config.10.folds$maxstep,numcells))
FakeResults2 = array(0, dim=c(length(methods) * 2 + 1, config.10.folds$maxstep,numcells))

set.seed(1)

checkpoints = c(50,100,150,200)
start = proc.time()[[3]]
for (i in 1:numcells){
  tryCatch({
    TrueResults2[,,i] = forward_selection(Y=Y1[,i],X=X1,cov_inds=covinds,methods=methods,
                                                    config=config.10.folds)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    FakeResults2[,,i] = forward_selection(Y=Y2[,i],X=X1,cov_inds=covinds,methods=methods,
                                                    config=config.10.folds,FAKE=T)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  print_progress(i,numcells,start)
  
  if (i %in% checkpoints){
    plot_results(TrueResults2[,,1:i],bettermethodnames,0.05,cex=2,nstep=3)
    plot_results(FakeResults2[,,1:i],bettermethodnames,0.05,hline = T,cex=2,nstep=3)
  }
}

if (SAVE){
  saveRDS(TrueResults2,paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
  saveRDS(FakeResults2,paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))
}

#### Plotting results from calcium data, 10 folds ----
TrueResults2 = readRDS(paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
FakeResults2 = readRDS(paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))

alpha = 0.05
bettermethodnames = c("CV","SR")#,"SR_B")

plot_results(TrueResults2,bettermethodnames,alpha,cex=2,nstep=3)
plot_results(FakeResults2,bettermethodnames,alpha,hline = T,cex=2,nstep=3)

#
#### Make tuning tables ----
grid_cells = readRDS(paste(processed_data_folder,"gridcells.RDS",sep=""))

resmat = readRDS(paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
csind = 5
msrrind = 4
srind = 3
postunedCS = union(union(which(resmat[1,1,] == 2 & resmat[csind,1,] <= 0.05),which(resmat[1,2,] == 2 & resmat[csind,2,] <= 0.05)),which(resmat[1,3,] == 2 & resmat[csind,3,] <= 0.05))
hdtunedCS = union(union(which(resmat[1,1,] == 1 & resmat[csind,1,] <= 0.05),which(resmat[1,2,] == 1 & resmat[csind,2,] <= 0.05)),which(resmat[1,3,] == 1 & resmat[csind,3,] <= 0.05))
spetunedCS = union(union(which(resmat[1,1,] == 3 & resmat[csind,1,] <= 0.05),which(resmat[1,2,] == 3 & resmat[csind,2,] <= 0.05)),which(resmat[1,3,] == 3 & resmat[csind,3,] <= 0.05))
postunedmSRR = union(union(which(resmat[1,1,] == 2 & resmat[msrrind,1,] <= 0.05),which(resmat[1,2,] == 2 & resmat[msrrind,2,] <= 0.05)),which(resmat[1,3,] == 2 & resmat[msrrind,3,] <= 0.05))
hdtunedmSRR = union(union(which(resmat[1,1,] == 1 & resmat[msrrind,1,] <= 0.05),which(resmat[1,2,] == 1 & resmat[msrrind,2,] <= 0.05)),which(resmat[1,3,] == 1 & resmat[msrrind,3,] <= 0.05))
spetunedmSRR = union(union(which(resmat[1,1,] == 3 & resmat[msrrind,1,] <= 0.05),which(resmat[1,2,] == 3 & resmat[msrrind,2,] <= 0.05)),which(resmat[1,3,] == 3 & resmat[msrrind,3,] <= 0.05))

resmat = readRDS(paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
postunedSR = union(union(which(resmat[1,1,] == 2 & resmat[srind,1,] <= 0.05),which(resmat[1,2,] == 2 & resmat[srind,2,] <= 0.05)),which(resmat[1,3,] == 2 & resmat[srind,3,] <= 0.05))
hdtunedSR = union(union(which(resmat[1,1,] == 1 & resmat[srind,1,] <= 0.05),which(resmat[1,2,] == 1 & resmat[srind,2,] <= 0.05)),which(resmat[1,3,] == 1 & resmat[srind,3,] <= 0.05))
spetunedSR = union(union(which(resmat[1,1,] == 3 & resmat[srind,1,] <= 0.05),which(resmat[1,2,] == 3 & resmat[srind,2,] <= 0.05)),which(resmat[1,3,] == 3 & resmat[srind,3,] <= 0.05))

alltunedCS = intersect(hdtunedCS,intersect(postunedCS,spetunedCS))

tuningtableCS = matrix(0,nrow=4,ncol=4)
tuningtableCS[1,1] = length(hdtunedCS)
tuningtableCS[2,2] = length(postunedCS)
tuningtableCS[3,3] = length(spetunedCS)
tuningtableCS[4,4] = length(grid_cells)
tuningtableCS[1,2] = tuningtableCS[2,1] = length(intersect(hdtunedCS,postunedCS))
tuningtableCS[1,3] = tuningtableCS[3,1] = length(intersect(hdtunedCS,spetunedCS))
tuningtableCS[1,4] = tuningtableCS[4,1] = length(intersect(hdtunedCS,grid_cells))
tuningtableCS[2,3] = tuningtableCS[3,2] = length(intersect(postunedCS,spetunedCS))
tuningtableCS[2,4] = tuningtableCS[4,2] = length(intersect(postunedCS,grid_cells))
tuningtableCS[3,4] = tuningtableCS[4,3] = length(intersect(spetunedCS,grid_cells))

tuningtablemSRR = matrix(0,nrow=4,ncol=4)
tuningtablemSRR[1,1] = length(hdtunedmSRR)
tuningtablemSRR[2,2] = length(postunedmSRR)
tuningtablemSRR[3,3] = length(spetunedmSRR)
tuningtablemSRR[4,4] = length(grid_cells)
tuningtablemSRR[1,2] = tuningtablemSRR[2,1] = length(intersect(hdtunedmSRR,postunedmSRR))
tuningtablemSRR[1,3] = tuningtablemSRR[3,1] = length(intersect(hdtunedmSRR,spetunedmSRR))
tuningtablemSRR[1,4] = tuningtablemSRR[4,1] = length(intersect(hdtunedmSRR,grid_cells))
tuningtablemSRR[2,3] = tuningtablemSRR[3,2] = length(intersect(postunedmSRR,spetunedmSRR))
tuningtablemSRR[2,4] = tuningtablemSRR[4,2] = length(intersect(postunedmSRR,grid_cells))
tuningtablemSRR[3,4] = tuningtablemSRR[4,3] = length(intersect(spetunedmSRR,grid_cells))

tuningtableSR = matrix(0,nrow=4,ncol=4)
tuningtableSR[1,1] = length(hdtunedSR)
tuningtableSR[2,2] = length(postunedSR)
tuningtableSR[3,3] = length(spetunedSR)
tuningtableSR[4,4] = length(grid_cells)
tuningtableSR[1,2] = tuningtableSR[2,1] = length(intersect(hdtunedSR,postunedSR))
tuningtableSR[1,3] = tuningtableSR[3,1] = length(intersect(hdtunedSR,spetunedSR))
tuningtableSR[1,4] = tuningtableSR[4,1] = length(intersect(hdtunedSR,grid_cells))
tuningtableSR[2,3] = tuningtableSR[3,2] = length(intersect(postunedSR,spetunedSR))
tuningtableSR[2,4] = tuningtableSR[4,2] = length(intersect(postunedSR,grid_cells))
tuningtableSR[3,4] = tuningtableSR[4,3] = length(intersect(spetunedSR,grid_cells))

colnames(tuningtableCS) = colnames(tuningtableSR) = colnames(tuningtablemSRR) = c("HD","Position","Speed","Grid")
rownames(tuningtableCS) = rownames(tuningtableSR) = rownames(tuningtablemSRR) = c("HD","Position","Speed","Grid")

tuningdfCS = data.frame(tuningtableCS)
tuningdfmSRR = data.frame(tuningtablemSRR)
tuningdfSR = data.frame(tuningtableSR)


gridandposCS = intersect(grid_cells,postunedCS)
posnotgridCS = setdiff(postunedCS,grid_cells)
gridnotposCS = setdiff(grid_cells,postunedCS)
hdspeenotposCS = setdiff(intersect(hdtunedCS,spetunedCS),postunedCS)
hdnotposCS = setdiff(hdtunedCS,postunedCS)
hdspeeposCS = intersect(intersect(hdtunedCS,spetunedCS),postunedCS)
hdspeedCS = intersect(hdtunedCS,spetunedCS)
onlyHDCS = setdiff(setdiff(hdtunedCS,spetunedCS),postunedCS)
onlySpeedCS = setdiff(setdiff(spetunedCS,hdtunedCS),postunedCS)

saveRDS(gridandposCS,paste(result_folder,"gridandposCS.RDS",sep=""))
saveRDS(posnotgridCS,paste(result_folder,"posnotgridCS.RDS",sep=""))
saveRDS(gridnotposCS,paste(result_folder,"gridnotposCS.RDS",sep=""))

saveRDS(tuningdfCS,paste(result_folder,"tuningdfCS.RDS",sep=""))
saveRDS(tuningdfmSRR,paste(result_folder,"tuningdfmSRR.RDS",sep=""))
saveRDS(tuningdfSR,paste(result_folder,"tuningdfSR.RDS",sep=""))

saveRDS(hdtunedCS,paste(result_folder,"hdtunedCS.RDS",sep=""))
saveRDS(postunedCS,paste(result_folder,"postunedCS.RDS",sep=""))
saveRDS(spetunedCS,paste(result_folder,"spetunedCS.RDS",sep=""))

saveRDS(hdtunedmSRR,paste(result_folder,"hdtunedmSRR.RDS",sep=""))
saveRDS(postunedmSRR,paste(result_folder,"postunedmSRR.RDS",sep=""))
saveRDS(spetunedmSRR,paste(result_folder,"spetunedmSRR.RDS",sep=""))

saveRDS(hdtunedSR,paste(result_folder,"hdtunedSR.RDS",sep=""))
saveRDS(postunedSR,paste(result_folder,"postunedSR.RDS",sep=""))
saveRDS(spetunedSR,paste(result_folder,"spetunedSR.RDS",sep=""))


#### Error rates and power estimates and confidence intervals ----


# Simulated data
ntrials = 300
alpha = 0.05

methods = c("CV","LRT","SR","mSR","mSRR","CS")
error_rates_and_CI = matrix(NA,nrow=4,ncol=length(methods))
rownames(error_rates_and_CI) = c("Errors","Error_rate","Lower_bound","Upper_bound")
colnames(error_rates_and_CI) = methods
power_and_CI = matrix(NA,nrow=4,ncol=length(methods))
rownames(power_and_CI) = c("Successes","Power","Lower_bound","Upper_bound")
colnames(power_and_CI) = methods

TrueResultMat1 = readRDS(paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
FakeResultMat1 = readRDS(paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))
TrueResultMat2 = readRDS(paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
FakeResultMat2 = readRDS(paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))

# Dimension of result matrices is (number of methods + 1) X (number of steps) X (number of cells)

m_ind = c(2,6,3,4,5)
e_ind = c(1,2,4,5,6)
for (i in 1:5){
  m = m_ind[i]
  e = e_ind[i]
  error_rates_and_CI[1,e] = length(which(FakeResultMat1[m,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,e] = error_rates_and_CI[1,e]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,e],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,e] = ci$Lower.limit
  error_rates_and_CI[4,e] = ci$Upper.limit
  
  
  power_and_CI[1,e] = sum(c(
    length(which(TrueResultMat1[1,1,1:100] == 2 & TrueResultMat1[m,1,1:100] <= alpha)) + length(which(TrueResultMat1[1,2,1:100] == 2 & TrueResultMat1[m,2,1:100] <= alpha)) + length(which(TrueResultMat1[1,3,1:100] == 2 & TrueResultMat1[m,3,1:100] <= alpha)),
    length(which(TrueResultMat1[1,1,101:200] == 1 & TrueResultMat1[m,1,101:200] <= alpha)) + length(which(TrueResultMat1[1,2,101:200] == 1 & TrueResultMat1[m,2,101:200] <= alpha)) + length(which(TrueResultMat1[1,3,101:200] == 1 & TrueResultMat1[m,3,101:200] <= alpha)),
    length(which(TrueResultMat1[1,1,201:300] == 3 & TrueResultMat1[m,1,201:300] <= alpha)) + length(which(TrueResultMat1[1,2,201:300] == 3 & TrueResultMat1[m,2,201:300] <= alpha)) + length(which(TrueResultMat1[1,3,201:300] == 3 & TrueResultMat1[m,3,201:300] <= alpha))
  ))
  # Position (covariate number 2) included means a success
  power_and_CI[2,e] = power_and_CI[1,e]/ntrials
  ci = clopper.pearson.ci(power_and_CI[1,e],ntrials,alpha,CI="two.sided")
  power_and_CI[3,e] = ci$Lower.limit
  power_and_CI[4,e] = ci$Upper.limit
}

if (TRUE){
  error_rates_and_CI[1,3] = length(which(FakeResultMat2[3,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,3] = error_rates_and_CI[1,3]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,3],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,3] = ci$Lower.limit
  error_rates_and_CI[4,3] = ci$Upper.limit
  
  power_and_CI[1,3] = sum(c(
    length(which(TrueResultMat2[1,1,1:100] == 2 & TrueResultMat2[3,1,1:100] <= alpha)) + length(which(TrueResultMat2[1,2,1:100] == 2 & TrueResultMat2[3,2,1:100] <= alpha)) + length(which(TrueResultMat2[1,3,1:100] == 2 & TrueResultMat2[3,3,1:100] <= alpha)),
    length(which(TrueResultMat2[1,1,101:200] == 1 & TrueResultMat2[3,1,101:200] <= alpha)) + length(which(TrueResultMat2[1,2,101:200] == 1 & TrueResultMat2[3,2,101:200] <= alpha)) + length(which(TrueResultMat2[1,3,101:200] == 1 & TrueResultMat2[3,3,101:200] <= alpha)),
    length(which(TrueResultMat2[1,1,201:300] == 3 & TrueResultMat2[3,1,201:300] <= alpha)) + length(which(TrueResultMat2[1,2,201:300] == 3 & TrueResultMat2[3,2,201:300] <= alpha)) + length(which(TrueResultMat2[1,3,201:300] == 3 & TrueResultMat2[3,3,201:300] <= alpha))
  ))
  # Position (covariate number 2) included means a success
  power_and_CI[2,3] = power_and_CI[1,3]/ntrials
  ci = clopper.pearson.ci(power_and_CI[1,3],ntrials,alpha,CI="two.sided")
  power_and_CI[3,3] = ci$Lower.limit
  power_and_CI[4,3] = ci$Upper.limit
}

saveRDS(error_rates_and_CI,paste(result_folder,"Error_rates_with_CI_Simulated_data.RDS",sep=""))
saveRDS(power_and_CI,paste(result_folder,"Power_estimates_with_CI_Simulated_data.RDS",sep=""))

error_rates_and_CI = round(error_rates_and_CI*1000)/1000
power_and_CI = round(power_and_CI*1000)/1000

t(error_rates_and_CI)
t(power_and_CI)

# Calcium data
alpha = 0.05
grid_cells = readRDS(paste(processed_data_folder,"gridcells.RDS",sep=""))
ntrials = 249
ngridcells = length(grid_cells)

methods = c("CV","LRT","SR","mSR","mSRR","CS")
error_rates_and_CI = matrix(NA,nrow=4,ncol=length(methods))
rownames(error_rates_and_CI) = c("Errors","Error_rate","Lower_bound","Upper_bound")
colnames(error_rates_and_CI) = methods
power_and_CI = matrix(NA,nrow=4,ncol=length(methods))
rownames(power_and_CI) = c("Successes","Power","Lower_bound","Upper_bound")
colnames(power_and_CI) = methods

TrueResultMat1 = readRDS(paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
FakeResultMat1 = readRDS(paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))
TrueResultMat2 = readRDS(paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
FakeResultMat2 = readRDS(paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))

# Dimension of result matrices is (number of methods + 1) X (number of steps) X (number of cells)

m_ind = c(2,6,3,4,5)
e_ind = c(1,2,4,5,6)
for (i in 1:5){
  m = m_ind[i]
  e = e_ind[i]
  error_rates_and_CI[1,e] = length(which(FakeResultMat1[m,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,e] = error_rates_and_CI[1,e]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,e],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,e] = ci$Lower.limit
  error_rates_and_CI[4,e] = ci$Upper.limit
  
  
  power_and_CI[1,e] = length(which(TrueResultMat1[1,1,grid_cells] == 2 & TrueResultMat1[m,1,grid_cells] <= alpha)) + length(which(TrueResultMat1[1,2,grid_cells] == 2 & TrueResultMat1[m,2,grid_cells] <= alpha)) + length(which(TrueResultMat1[1,3,grid_cells] == 2 & TrueResultMat1[m,3,grid_cells] <= alpha))
  # Position included for a grid cell is considered a success
  power_and_CI[2,e] = power_and_CI[1,e]/ngridcells
  ci = clopper.pearson.ci(power_and_CI[1,e],ngridcells,alpha,CI="two.sided")
  power_and_CI[3,e] = ci$Lower.limit
  power_and_CI[4,e] = ci$Upper.limit
}

if (TRUE){
  error_rates_and_CI[1,3] = length(which(FakeResultMat2[3,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,3] = error_rates_and_CI[1,3]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,3],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,3] = ci$Lower.limit
  error_rates_and_CI[4,3] = ci$Upper.limit
  
  power_and_CI[1,3] = length(which(TrueResultMat2[1,1,grid_cells] == 2 & TrueResultMat2[3,1,grid_cells] <= alpha)) + length(which(TrueResultMat2[1,2,grid_cells] == 2 & TrueResultMat2[3,2,grid_cells] <= alpha)) + length(which(TrueResultMat2[1,3,grid_cells] == 2 & TrueResultMat2[3,3,grid_cells] <= alpha))
  # Position included for a grid cell is considered a success
  power_and_CI[2,3] = power_and_CI[1,3]/ngridcells
  ci = clopper.pearson.ci(power_and_CI[1,3],ngridcells,alpha,CI="two.sided")
  power_and_CI[3,3] = ci$Lower.limit
  power_and_CI[4,3] = ci$Upper.limit
}

saveRDS(error_rates_and_CI,paste(result_folder,"Error_rates_with_CI_Calcium_data.RDS",sep=""))
saveRDS(power_and_CI,paste(result_folder,"Power_estimates_with_CI_Calcium_data.RDS",sep=""))

error_rates_and_CI = round(error_rates_and_CI*1000)/1000
power_and_CI = round(power_and_CI*1000)/1000

t(error_rates_and_CI)
t(power_and_CI)


#### Calculate and save effect sizes ----
# Calcium data

Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Y2.RDS",sep=""))
X1 = readRDS(paste(processed_data_folder,"X1.RDS",sep=""))
X2 = readRDS(paste(processed_data_folder,"X1.RDS",sep=""))
config.20.folds = readRDS(paste(processed_data_folder,"config_20_folds.RDS",sep=""))
config.10.folds = readRDS(paste(processed_data_folder,"config_10_folds.RDS",sep=""))
skips = config.20.folds$cyclic_skip
N = config.20.folds$N
keeps = c((skips+1):(round(N/2)-skips),(round(N/2)+skips+1):(N-skips))
covinds = readRDS(paste(processed_data_folder,"covinds.RDS",sep=""))

TrueResults = readRDS(paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
FakeResults = readRDS(paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))
TrueResults2 = readRDS(paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
FakeResults2 = readRDS(paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))

effect_sizes = list("real"=list(), "fake"=list())

effect_sizes$real$LRT = AddMcFadden(TrueResults, 6, X1, Y1, covinds)
effect_sizes$fake$LRT = AddMcFadden(FakeResults, 6, X1, Y2, covinds)
effect_sizes$real$CS = AddMcFadden(TrueResults, 5, X1, Y1, covinds, keeps=keeps)
effect_sizes$fake$CS = AddMcFadden(FakeResults, 5, X1, Y2, covinds, keeps=keeps)
effect_sizes$real$SR = AddSR(TrueResults2, 3, X1, Y1, covinds, config=config.10.folds)
effect_sizes$fake$SR = AddSR(FakeResults2, 3, X1, Y2, covinds, config=config.10.folds)
effect_sizes$real$mSR = AddSR(TrueResults, 3, X1, Y1, covinds, config=config.20.folds)
effect_sizes$fake$mSR = AddSR(FakeResults, 3, X1, Y2, covinds, config=config.20.folds)
effect_sizes$real$mSRR = AddSRRev(TrueResults, 4, X1, Y1, covinds, config=config.20.folds)
effect_sizes$fake$mSRR = AddSRRev(FakeResults, 4, X1, Y2, covinds, config=config.20.folds)

grid_cells = readRDS(paste(processed_data_folder,"gridcells.RDS",sep=""))
effect_sizes$grid_CS = AddMcFadden(TrueResults[,,grid_cells], 5, X1, Y1, covinds, keeps=keeps)
effect_sizes$nongrid_CS = AddMcFadden(TrueResults[,,setdiff(1:249, grid_cells)], 5, X1, Y1, covinds, keeps=keeps)
effect_sizes$grid_mSRR = AddSRRev(TrueResults[,,grid_cells], 4, X1, Y1, covinds, config=config.20.folds)
effect_sizes$nongrid_mSRR = AddSRRev(TrueResults[,,setdiff(1:249, grid_cells)], 4, X1, Y1, covinds, config=config.20.folds)


# Simulated data

effect_sizes_sim = list("real"=list(), "fake"=list())
X = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_X.RDS",sep=""))
Y1 = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_Y2.RDS",sep=""))
covinds_sim = readRDS(paste(processed_data_folder,"covinds_sim_new.RDS",sep=""))

TrueResultMat = readRDS(paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
FakeResultMat = readRDS(paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))
TrueResultMat2 = readRDS(paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
FakeResultMat2 = readRDS(paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))

effect_sizes_sim$real$LRT = AddMcFadden(TrueResultMat, 6, X, Y1, covinds_sim)
effect_sizes_sim$fake$LRT = AddMcFadden(FakeResultMat, 6, X, Y2, covinds_sim)
effect_sizes_sim$real$CS = AddMcFadden(TrueResultMat, 5, X, Y1, covinds_sim, keeps=keeps)
effect_sizes_sim$fake$CS = AddMcFadden(FakeResultMat, 5, X, Y2, covinds_sim, keeps=keeps)
effect_sizes_sim$real$SR = AddSR(TrueResultMat2, 3, X, Y1, covinds_sim, config=config.10.folds)
effect_sizes_sim$fake$SR = AddSR(FakeResultMat2, 3, X, Y2, covinds_sim, config=config.10.folds)
effect_sizes_sim$real$mSR = AddSR(TrueResultMat, 3, X, Y1, covinds_sim, config=config.20.folds)
effect_sizes_sim$fake$mSR = AddSR(FakeResultMat, 3, X, Y2, covinds_sim, config=config.20.folds)
effect_sizes_sim$real$mSRR = AddSRRev(TrueResultMat, 4, X, Y1, covinds_sim, config=config.20.folds)
effect_sizes_sim$fake$mSRR = AddSRRev(FakeResultMat, 4, X, Y2, covinds_sim, config=config.20.folds)

saveRDS(effect_sizes, paste(result_folder, "effect_sizes.RDS", sep=""))
saveRDS(effect_sizes_sim, paste(result_folder, "effect_sizes_sim.RDS", sep=""))

#
#### Plot effect sizes ----
effect_sizes = readRDS(paste(result_folder, "effect_sizes.RDS", sep=""))
effect_sizes_sim = readRDS(paste(result_folder, "effect_sizes_sim.RDS", sep=""))


Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
meanrates = apply(Y1,2,mean)/(1/7.25)
grid_cells = readRDS(paste(processed_data_folder,"gridcells.RDS",sep=""))
plot_grid_stuff(effect_sizes$grid_CS, effect_sizes$nongrid_CS, T, meanrates1=meanrates[grid_cells], meanrates2=meanrates[setdiff(1:249, grid_cells)])
plot_grid_stuff(effect_sizes$grid_mSRR, effect_sizes$nongrid_mSRR, F)


plot_effect_sizes(effect_sizes, lrt_lims=c(0.2, 0.4, 0.1))
plot_effect_sizes(effect_sizes,T, lrt_lims=c(0.05, 0.15, 0.05))
plot_effect_sizes(effect_sizes_sim, lrt_lims=c(0.03, 0.03, 0.03),SIM=T)
plot_effect_sizes(effect_sizes_sim,T,lrt_lims=c(0.01, 0.01, 0.01))#, lrt_lims=c(0.01, 0.01, 0.01))


#

