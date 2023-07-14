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
library(tidyverse)
library(RColorBrewer)
library(dendextend)
library(rmatio)
library(R.matlab)
library(EnvStats)
library(rhdf5)
library(GenBinomApps)
library(ggvenn)

data_folder = "/Users/fredrine/Documents/CalciumData/Prepped data/" # Folder in which the NAT.mat and NeuronInformation.mat files for the session of data from Zong et. al. (2022) lies
functions_folder = "/Users/fredrine/Documents/CalciumData/Scripts/" # Folder in which the functions_model_selection.R script lies
processed_data_folder = "/Users/fredrine/Documents/CalciumData/Processed_data/" # Folder in which the results from the initial processing will end up (binarized cell activity, covariate matrix with splined versions of the covariates)
result_folder = "/Users/fredrine/Documents/CalciumData/Result_data/" # Folder in which the results from running the forward selection methods will end up

#### Load functions ----
source(paste(functions_folder,"functions_model_selection.R",sep=""))
#### Load data for rat 97045, date 20210317----

# Loads data from the NAT.mat and NeuronInformation.mat files, binarizes the neural activity, and increases the bin size of the tracking data to match the neural data.

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
#### Constructing the design matrix used in the GLM ----

# Creates splined version of the covariates, and filters out some of the neurons that fire very sparsely.

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
  saveRDS(Y1,paste(processed_data_folder,"Y1.RDS",sep=""))
  saveRDS(Y2,paste(processed_data_folder,"Y2.RDS",sep=""))
  saveRDS(X1,paste(processed_data_folder,"X1.RDS",sep=""))
  saveRDS(X2,paste(processed_data_folder,"X2.RDS",sep=""))
  saveRDS(covinds,paste(processed_data_folder,"covinds.RDS",sep=""))
  saveRDS(cellstokeep,paste(processed_data_folder,"cellstokeep.RDS",sep=""))
}
#### Forward selection on simulated data, 20 folds  ----

# Simulates data with properties similar to the calcium data, and performs covarariate selection with various methods, using 20 folds in the cross-validation scheme.

SAVE = FALSE

set.seed(1)
num_trials = 300
N = 12000
sd = 20
bprop = 0.5
alpha = 0.05
chunksperfold = 4
maxstep = 3
scale = 0.25

hdknots_sim = 5
posknots_sim = 2
speknots_sim = 5
knots = c(hdknots_sim,posknots_sim,posknots_sim,speknots_sim)

ncov = 3
covinds_sim = list(`1`=1:(hdknots_sim+1),`2`=(1+hdknots_sim+1):(1+hdknots_sim+(1+posknots_sim)^2),`3`=(2+hdknots_sim+(1+posknots_sim)^2):(2+hdknots_sim+(1+posknots_sim)^2+speknots_sim))

if (TRUE){
  sim.config = config.default()
  sim.config$alpha = alpha
  sim.config$nfolds = 20
  sim.config$chunksize = N/sim.config$nfolds/chunksperfold
  sim.config$family = "binomial"
  sim.config$SILENT = F
  sim.config$SKIP_NEIGHBOR = T
  sim.config$permshuffles = 119
  sim.config$signshuffles = 999
  sim.config$GLMNET = F
  sim.config$lambda = 0.005
  sim.config$glmnetalpha = 0
  sim.config$maxstep = maxstep
  sim.config$cyclic_skip = 75
} # Set configs

methods = c("CV","signrankmax","signrankmax2","cyclicskip")
bettermethodnames =c("CV","mSR","mSRR","CS")
TrueResultMat = array(0,dim=c(length(methods)+1,sim.config$maxstep,num_trials))
FakeResultMat = array(0,dim=c(length(methods)+1,sim.config$maxstep,num_trials))

start = proc.time()[[3]]
for (i in 1:num_trials){
  truedat = make_simulated_data(N,sd,bprop,ncov=ncov,scale=scale)
  fakedat = make_simulated_data(N,sd,ncov=ncov,scale=scale)
  trueY1 = truedat$y
  trueX1 = make_design_matrix(truedat[,2:dim(truedat)[2]],c(),c(),c(2),knots)
  fakeY1 = fakedat$y
  fakeX1 = make_design_matrix(fakedat[,2:dim(fakedat)[2]],c(),c(),c(2),knots)
  
  TrueResultMat[,,i] = forward_selection(Y=trueY1,X=trueX1,cov_inds=covinds_sim,methods=methods,
                                                    config=sim.config)
  FakeResultMat[,,i] = forward_selection(Y=fakeY1,X=fakeX1,cov_inds=covinds_sim,methods=methods,
                                                    config=sim.config)
  
  print_progress(i,num_trials,start)
}


if (SAVE){
  saveRDS(TrueResultMat,paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
  saveRDS(FakeResultMat,paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))
}

#### Plotting results from simulated data, 20 folds ----
TrueResultMat = readRDS(paste(result_folder,"Simulated_data_true_effect_20_folds.RDS",sep=""))
FakeResultMat = readRDS(paste(result_folder,"Simulated_data_no_effect_20_folds.RDS",sep=""))
alpha = 0.05
bettermethodnames = c("CV","mSR","mSRR","CS")

plot_results(TrueResultMat,bettermethodnames,alpha,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
plot_results(FakeResultMat,bettermethodnames,alpha,hline=T,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))

#

#### Forward selection on simulated data, 10 folds  ----

# Simulates data with properties similar to the calcium data, and performs covarariate selection with various methods, using 10 folds in the cross-validation scheme.

SAVE = FALSE

set.seed(1)
num_trials = 300
N = 12000
sd = 20
bprop = 0.5
alpha = 0.05
chunksperfold = 8
maxstep = 3
scale = 0.25

hdknots_sim = 5
posknots_sim = 2
speknots_sim = 5
knots = c(hdknots_sim,posknots_sim,posknots_sim,speknots_sim)

ncov = 3
covinds_sim = list(`1`=1:(hdknots_sim+1),`2`=(1+hdknots_sim+1):(1+hdknots_sim+(1+posknots_sim)^2),`3`=(2+hdknots_sim+(1+posknots_sim)^2):(2+hdknots_sim+(1+posknots_sim)^2+speknots_sim))

if (TRUE){
  sim.config = config.default()
  sim.config$alpha = alpha
  sim.config$nfolds = 10
  sim.config$chunksize = N/sim.config$nfolds/chunksperfold
  sim.config$family = "binomial"
  sim.config$SILENT = F
  sim.config$SKIP_NEIGHBOR = F
  sim.config$permshuffles = 119
  sim.config$signshuffles = 999
  sim.config$GLMNET = F
  sim.config$lambda = 0.005
  sim.config$glmnetalpha = 0
  sim.config$maxstep = maxstep
  sim.config$cyclic_skip = 75
} # Set configs

methods = c("CV","oldHC","oldHCBonf")
bettermethodnames = c("CV","SR","SR_B")
TrueResultMat = array(0,dim=c(length(methods)+1,sim.config$maxstep,num_trials))
FakeResultMat = array(0,dim=c(length(methods)+1,sim.config$maxstep,num_trials))

start = proc.time()[[3]]
for (i in 1:num_trials){
  truedat = make_simulated_data(N,sd,bprop,ncov=ncov,scale=scale)
  fakedat = make_simulated_data(N,sd,ncov=ncov,scale=scale)
  trueY1 = truedat$y
  trueX1 = make_design_matrix(truedat[,2:dim(truedat)[2]],c(),c(),c(2),knots)
  fakeY1 = fakedat$y
  fakeX1 = make_design_matrix(fakedat[,2:dim(fakedat)[2]],c(),c(),c(2),knots)
  
  TrueResultMat[,,i] = forward_selection(Y=trueY1,X=trueX1,cov_inds=covinds_sim,methods=methods,
                                                    config=sim.config)
  FakeResultMat[,,i] = forward_selection(Y=fakeY1,X=fakeX1,cov_inds=covinds_sim,methods=methods,
                                                    config=sim.config)
  
  print_progress(i,num_trials,start)
}

if (SAVE){
  saveRDS(TrueResultMat,paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
  saveRDS(FakeResultMat,paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))
}


#### Plotting results from simulated data, 10 folds ----
TrueResultMat = readRDS(paste(result_folder,"Simulated_data_true_effect_10_folds.RDS",sep=""))
FakeResultMat = readRDS(paste(result_folder,"Simulated_data_no_effect_10_folds.RDS",sep=""))

alpha = 0.05
bettermethodnames = c("CV","SR","SR_B")

plot_results(TrueResultMat,bettermethodnames,alpha,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))
plot_results(FakeResultMat,bettermethodnames,alpha,hline=T,cex=2,nstep=3,covnames=c("B (2D)","A (1D)","C (1D)"))

#
#### Forward selection on calcium data, 20 folds ----

# Performs covarariate selection with various methods on calcium data, using 20 folds in the cross-validation scheme.

#hdknots = 6
#posknots = 4
#speknots = 5
#fam = "binomial" 
#datalength = 12000
#spikethresh = 0.02

SAVE = FALSE

Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Y2.RDS",sep=""))
X1 = readRDS(paste(processed_data_folder,"X1.RDS",sep=""))
X2 = readRDS(paste(processed_data_folder,"X2.RDS",sep=""))
covinds = readRDS(paste(processed_data_folder,"covinds.RDS",sep=""))

numcells = dim(Y1)[2]
methods = c("CV","signrankmax","signrankmax2","cyclicskip")
bettermethodnames = c("CV","mSR","mSRR","CS")
TrueResults = array(0,dim=c(length(methods)+1,3,numcells))
FakeResults = array(0,dim=c(length(methods)+1,3,numcells))
alpha = 0.05
if (TRUE){
  selection.config = config.default()
  selection.config$alpha = 0.05
  selection.config$nfolds = 20
  selection.config$chunksize = 150
  selection.config$family = "binomial"
  selection.config$SILENT = F
  selection.config$SKIP_NEIGHBOR = T
  selection.config$permshuffles = 119
  selection.config$signshuffles = 999
  selection.config$GLMNET = F
  selection.config$lambda = 0.005
  selection.config$maxstep = 3
} # Set configs

start = proc.time()[[3]]
for (i in 1:numcells){
  tryCatch({
    TrueResults[,,i] = forward_selection(Y=Y1[,i],X=X1,cov_inds=covinds,methods=methods,
                                                    config=selection.config)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    FakeResults[,,i] = forward_selection(Y=Y2[,i],X=X1,cov_inds=covinds,methods=methods,
                                                    config=selection.config)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
  print_progress(i,numcells,start)
}

if (SAVE){
  saveRDS(TrueResults,paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
  saveRDS(FakeResults,paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))
}


#### Plotting results from calcium data, 20 folds ----
TrueResults = readRDS(paste(result_folder,"Calcium_data_matched_20_folds.RDS",sep=""))
FakeResults = readRDS(paste(result_folder,"Calcium_data_mismatched_20_folds.RDS",sep=""))

alpha = 0.05
bettermethodnames = c("CV","mSR","mSRR","CS")


plot_results(TrueResults,bettermethodnames,alpha,cex=2,nstep=3)
plot_results(FakeResults,bettermethodnames,alpha,hline = T,cex=2,nstep=3)

#
#### Forward selection on calcium data, 10 folds ----

# Performs covarariate selection with various methods on calcium data, using 10 folds in the cross-validation scheme.

#hdknots = 6
#posknots = 4
#speknots = 5
#fam = "binomial" 
#datalength = 12000
#spikethresh = 0.02

SAVE = FALSE

Y1 = readRDS(paste(processed_data_folder,"Y1.RDS",sep=""))
Y2 = readRDS(paste(processed_data_folder,"Y2.RDS",sep=""))
X1 = readRDS(paste(processed_data_folder,"X1.RDS",sep=""))
X2 = readRDS(paste(processed_data_folder,"X2.RDS",sep=""))
covinds = readRDS(paste(processed_data_folder,"covinds.RDS",sep=""))

numcells = dim(Y1)[2]
methods = c("CV","oldHC","oldHCBonf")
bettermethodnames = c("CV","SR","SR_B")
TrueResults = array(0,dim=c(length(methods)+1,3,numcells))
FakeResults = array(0,dim=c(length(methods)+1,3,numcells))
alpha = 0.05
if (TRUE){
  selection.config = config.default()
  selection.config$alpha = 0.05
  selection.config$nfolds = 10
  selection.config$chunksize = 150
  selection.config$family = "binomial"
  selection.config$SILENT = F
  selection.config$SKIP_NEIGHBOR = F
  selection.config$permshuffles = 119
  selection.config$signshuffles = 999
  selection.config$GLMNET = F
  selection.config$lambda = 0.005
  selection.config$maxstep = 3
} # Set configs

start = proc.time()[[3]]
for (i in 1:numcells){
  tryCatch({
    TrueResults[,,i] = forward_selection(Y=Y1[,i],X=X1,cov_inds=covinds,methods=methods,
                                                    config=selection.config)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  tryCatch({
    FakeResults[,,i] = forward_selection(Y=Y2[,i],X=X1,cov_inds=covinds,methods=methods,
                                                    config=selection.config)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

  print_progress(i,numcells,start)
}

if (SAVE){
  saveRDS(TrueResults,paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
  saveRDS(FakeResults,paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))
}

#### Plotting results from calcium data, 10 folds ----
TrueResults = readRDS(paste(result_folder,"Calcium_data_matched_10_folds.RDS",sep=""))
FakeResults = readRDS(paste(result_folder,"Calcium_data_mismatched_10_folds.RDS",sep=""))

alpha = 0.05
bettermethodnames = c("CV","SR","SR_B")


plot_results(TrueResults,bettermethodnames,alpha,cex=2,nstep=3)
plot_results(FakeResults,bettermethodnames,alpha,hline = T,cex=2,nstep=3)

#
#### Make tuning tables ----

# Summarizes the results from the methods using cyclical shift, signed-rank test, and the modified method using signed-rank.

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

# Calculate error rates of each method (for simulated data, how many times is a covariate included when it shouldn't in the case where there is no connection between the observed covariates and the response,
# and for calcium data when the neural data and behavior is mismatched), with corresponding 95% confidence intervals (Clopper-Pearson).
# Calculate power of each method (for simulated data, how many times is a position included when it should in the case where there is a connection between the position and the response,
# and for calcium data the proportion of grid cells that are classified as being positionally tuned), with corresponding 95% confidence intervals (Clopper-Pearson).

# Simulated data
ntrials = 300
alpha = 0.05

methods = c("CV","mSR","mSRR","CS","SR","SR_B")
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

for (i in 1:4){
  error_rates_and_CI[1,i] = length(which(FakeResultMat1[i+1,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,i] = error_rates_and_CI[1,i]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,i],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,i] = ci$Lower.limit
  error_rates_and_CI[4,i] = ci$Upper.limit
  
  
  power_and_CI[1,i] = length(which(TrueResultMat1[1,1,] == 2 & TrueResultMat1[i+1,1,] <= alpha)) + length(which(TrueResultMat1[1,2,] == 2 & TrueResultMat1[i+1,2,] <= alpha)) + length(which(TrueResultMat1[1,3,] == 2 & TrueResultMat1[i+1,3,] <= alpha))
  # Position (covariate number 2) included means a success
  power_and_CI[2,i] = power_and_CI[1,i]/ntrials
  ci = clopper.pearson.ci(power_and_CI[1,i],ntrials,alpha,CI="two.sided")
  power_and_CI[3,i] = ci$Lower.limit
  power_and_CI[4,i] = ci$Upper.limit
}

for (i in 1:2){
  error_rates_and_CI[1,i+4] = length(which(FakeResultMat2[i+2,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,i+4] = error_rates_and_CI[1,i+4]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,i+4],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,i+4] = ci$Lower.limit
  error_rates_and_CI[4,i+4] = ci$Upper.limit
  
  power_and_CI[1,i+4] = length(which(TrueResultMat2[1,1,] == 2 & TrueResultMat2[i+2,1,] <= alpha)) + length(which(TrueResultMat2[1,2,] == 2 & TrueResultMat2[i+2,2,] <= alpha)) + length(which(TrueResultMat2[1,3,] == 2 & TrueResultMat2[i+2,3,] <= alpha))
  # Position (covariate number 2) included means a success
  power_and_CI[2,i+4] = power_and_CI[1,i+4]/ntrials
  ci = clopper.pearson.ci(power_and_CI[1,i+4],ntrials,alpha,CI="two.sided")
  power_and_CI[3,i+4] = ci$Lower.limit
  power_and_CI[4,i+4] = ci$Upper.limit
}

saveRDS(error_rates_and_CI,paste(result_folder,"Error_rates_with_CI_Simulated_data.RDS",sep=""))
saveRDS(power_and_CI,paste(result_folder,"Power_estimates_with_CI_Simulated_data.RDS",sep=""))

error_rates_and_CI = round(error_rates_and_CI*1000)/1000
power_and_CI = round(power_and_CI*1000)/1000

error_rates_and_CI
power_and_CI

# Calcium data
alpha = 0.05
grid_cells = readRDS(paste(processed_data_folder,"gridcells.RDS",sep=""))
ntrials = 249
ngridcells = length(grid_cells)

methods = c("CV","mSR","mSRR","CS","SR","SR_B")
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

for (i in 1:4){
  error_rates_and_CI[1,i] = length(which(FakeResultMat1[i+1,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,i] = error_rates_and_CI[1,i]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,i],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,i] = ci$Lower.limit
  error_rates_and_CI[4,i] = ci$Upper.limit
  
  
  power_and_CI[1,i] = length(which(TrueResultMat1[1,1,grid_cells] == 2 & TrueResultMat1[i+1,1,grid_cells] <= alpha)) + length(which(TrueResultMat1[1,2,grid_cells] == 2 & TrueResultMat1[i+1,2,grid_cells] <= alpha)) + length(which(TrueResultMat1[1,3,grid_cells] == 2 & TrueResultMat1[i+1,3,grid_cells] <= alpha))
  # Position included for a grid cell is considered a success
  power_and_CI[2,i] = power_and_CI[1,i]/ngridcells
  ci = clopper.pearson.ci(power_and_CI[1,i],ngridcells,alpha,CI="two.sided")
  power_and_CI[3,i] = ci$Lower.limit
  power_and_CI[4,i] = ci$Upper.limit
}

for (i in 1:2){
  error_rates_and_CI[1,i+4] = length(which(FakeResultMat2[i+2,1,] <= alpha))
  # p-value <= alpha in first step when no observed covariate is relevant means you make an error
  error_rates_and_CI[2,i+4] = error_rates_and_CI[1,i+4]/ntrials
  ci = clopper.pearson.ci(error_rates_and_CI[1,i+4],ntrials,alpha,CI="two.sided")
  error_rates_and_CI[3,i+4] = ci$Lower.limit
  error_rates_and_CI[4,i+4] = ci$Upper.limit
  
  power_and_CI[1,i+4] = length(which(TrueResultMat2[1,1,grid_cells] == 2 & TrueResultMat2[i+2,1,grid_cells] <= alpha)) + length(which(TrueResultMat2[1,2,grid_cells] == 2 & TrueResultMat2[i+2,2,grid_cells] <= alpha)) + length(which(TrueResultMat2[1,3,grid_cells] == 2 & TrueResultMat2[i+2,3,grid_cells] <= alpha))
  # Position included for a grid cell is considered a success
  power_and_CI[2,i+4] = power_and_CI[1,i+4]/ngridcells
  ci = clopper.pearson.ci(power_and_CI[1,i+4],ngridcells,alpha,CI="two.sided")
  power_and_CI[3,i+4] = ci$Lower.limit
  power_and_CI[4,i+4] = ci$Upper.limit
}

saveRDS(error_rates_and_CI,paste(result_folder,"Error_rates_with_CI_Calcium_data.RDS",sep=""))
saveRDS(power_and_CI,paste(result_folder,"Power_estimates_with_CI_Calcium_data.RDS",sep=""))

error_rates_and_CI = round(error_rates_and_CI*1000)/1000
power_and_CI = round(power_and_CI*1000)/1000

error_rates_and_CI
power_and_CI
