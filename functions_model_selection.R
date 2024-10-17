#### Data processing functions ----
gaussian_smoother1D = function(dat,sigma,window=NULL,periodic=FALSE,SQUARE=FALSE,na.rm=F){
  if (is.null(window)){
    window = 8*sigma + 1
    mid = 4*sigma+1
  } else {
    mid = round((window-1)/2)+1
  }
  # Gaussian smoothing of a vector
  lendat = length(dat)
  kernel = array(0,dim=window)
  for (i in 1:window){
    if (SQUARE){
      kernel[i] = ((pnorm(0.5-(mid-i),sd=sigma))-(pnorm(-0.5-(mid-i),sd=sigma)))^2#dnorm(i-mid,sd=sigma)#
    } else {
      kernel[i] = dnorm(i-mid,sd=sigma)##pnorm(0.5-(mid-i),sd=sigma)-pnorm(-0.5-(mid-i),sd=sigma)#
    }
  }
  kernel = kernel/sum(kernel)
  newdat = array(0,lendat+2*(mid-1))
  newdat[mid:(lendat+mid-1)] = dat
  if (periodic){
    newdat[1:(mid-1)] = dat[(lendat-mid+2):lendat]
    newdat[(lendat+mid):(lendat+2*mid-2)] = dat[1:(mid-1)]
  } else {
    newdat[1:(mid-1)] = array(dat[1],mid-1)
    newdat[(lendat+mid):(lendat+2*mid-2)] = array(dat[lendat],mid-1)
  }
  result = array(0,lendat)
  for (i in 1:lendat){
    if (na.rm & is.na(newdat[i+mid-1])){
      result[i] = NA
    } else {
      result[i] = sum(kernel*newdat[i:(i+window-1)],na.rm=na.rm)
    }
  }
  return(result)
}

gaussian_smoother_invariant = function(dat,sigma,window=NULL,periodic=FALSE){
  # Scales the kernel so that the sum of squares is one
  # This means that when you filter independent data, the marginal variance is not changed
  if (is.null(window)){
    window = 8*sigma + 1
    mid = 4*sigma+1
  } else {
    mid = round((window-1)/2)+1
  }
  # Gaussian smoothing of a vector
  lendat = length(dat)
  kernel = array(0,dim=window)
  for (i in 1:window){
    kernel[i] = dnorm(i-mid,sd=sigma)##pnorm(0.5-(mid-i),sd=sigma)-pnorm(-0.5-(mid-i),sd=sigma)#
  }
  kernel = kernel/sqrt(sum(kernel^2))
  newdat = array(0,lendat+2*(mid-1))
  newdat[mid:(lendat+mid-1)] = dat
  if (periodic){
    newdat[1:(mid-1)] = dat[(lendat-mid+2):lendat]
    newdat[(lendat+mid):(lendat+2*mid-2)] = dat[1:(mid-1)]
  } else {
    newdat[1:(mid-1)] = array(dat[1],mid-1)
    newdat[(lendat+mid):(lendat+2*mid-2)] = array(dat[lendat],mid-1)
  }
  result = array(0,lendat)
  for (i in 1:lendat){
    result[i] = sum(kernel*newdat[i:(i+window-1)])
  }
  return(result)
}

exponential_smoother1D = function(dat,invrate,power=1,window=NULL,periodic=FALSE){
  if (is.null(window)){
    window = 16*invrate + 1
    mid = 8*invrate+1
  } else {
    mid = round((window-1)/2)+1
  }
  # Exponential smoothing of a vector
  lendat = length(dat)
  kernel = array(0,dim=window)
  #kernel[mid] = dexp(0,rate=1/invrate)
  kernel[mid] = 1/invrate
  for (i in (mid+1):window){
    #kernel[i] = dexp(i-mid,rate=1/invrate)
    kernel[i] = exp(-((i-mid)/invrate)^power)/invrate
    kernel[window-i+1] = kernel[i]
  }
  kernel = kernel/sum(kernel)
  newdat = array(0,lendat+2*(mid-1))
  newdat[mid:(lendat+mid-1)] = dat
  if (periodic){
    newdat[1:(mid-1)] = dat[(lendat-mid+2):lendat]
    newdat[(lendat+mid):(lendat+2*mid-2)] = dat[1:(mid-1)]
  } else {
    newdat[1:(mid-1)] = array(dat[1],mid-1)
    newdat[(lendat+mid):(lendat+2*mid-2)] = array(dat[lendat],mid-1)
  }
  result = array(0,lendat)
  for (i in 1:lendat){
    result[i] = sum(kernel*newdat[i:(i+window-1)])
  }
  return(result)
}

gaussian_smoother2D = function(dat,sigma,window=NULL,periodic=FALSE){
  dim1 = dim(dat)[1]
  dim2 = dim(dat)[2]
  
  if (is.null(window)){
    window = 8*sigma + 1
  } 
  kernel = matrix(0,nrow=window,ncol=window)
  mid = round((window-1)/2)+1
  for (i in 0:(mid-1)){
    for (j in 0:(mid-1)){
      weight = dmvnorm(c(i,j),sigma=diag(sigma,2))#dnorm(sqrt(i^2+j^2),sd=sigma)
      kernel[mid+i,mid+j] = weight
      kernel[mid+i,mid-j] = weight
      kernel[mid-i,mid+j] = weight
      kernel[mid-i,mid-j] = weight
    }
  }
  
  result = matrix(0,nrow=dim1,ncol=dim2)
  padded = matrix(nrow=dim1+window-1,ncol=dim2+window-1)
  padded[mid:(dim1+mid-1),mid:(dim2+mid-1)] = dat
  
  for (i in 1:dim1){
    for (j in 1:dim2){
      if (is.na(dat[i,j])){
        result[i,j] = NA
      } else {
        nonans = which(!is.na(padded[i:(i+window-1),j:(j+window-1)]))
        nonakernel = kernel/sum(kernel[nonans])
        prod = padded[i:(i+window-1),j:(j+window-1)]*nonakernel
        result[i,j] = sum(prod[nonans])
      }
    }
  }
  
  
  #result = t(apply(dat,1,gaussian_smoother1D,sigma=sigma,window=window,periodic=periodic))
  #result = apply(result,2,gaussian_smoother1D,sigma=sigma,window=window,periodic=periodic)
  return(result)
}

unperiodify = function(invec){
  # Makes the input vector no longer periodic in the interval -pi to pi
  vec = invec
  diffs = diff(vec)
  changes = which(abs(diffs)>pi)
  m = length(changes)
  
  prevchange = 0
  for (i in 1:m){
    if (i == m){
      end = length(vec)
    } else {
      end = changes[i+1]
    }
    change = prevchange - sign(diffs[changes[i]])
    vec[(changes[i]+1):end] = vec[(changes[i]+1):end] + 2*pi*change
    prevchange = change
  }
  return(vec)
}

periodify = function(vec){
  # Makes the vector periodic in the interval -pi to pi
  return((vec+pi)%%(2*pi)-pi)
}

fix_periodic = function(df,covname){
  m = dim(df)[1]
  wherena = which(is.na(df[[covname]]))
  wherenotna = which(!is.na(df[[covname]]))
  
  withoutna = df[[covname]][wherenotna]
  withoutnaunwrapped = unperiodify(withoutna)
  withnaunwrapped = array(NA,dim=m)
  withnaunwrapped[wherenotna] = withoutnaunwrapped
  return(periodify(approx(x=(1:m),y=withnaunwrapped,xout=(1:m),rule=2)$y))
}

splinify = function(covariate,min=0,max=1,nint=3,periodic=FALSE,XY=FALSE,perc95both=FALSE,perc95max=FALSE,findmin=FALSE,findmax=FALSE){
  if (XY){
    if (findmin){
      min = c(min(covariate$X),min(covariate$Y))
    }
    if (findmax){
      max = c(max(covariate$X),max(covariate$Y))
    }
    
    X_knots = seq(min[1],max[1],(max[1]-min[1])/(nint[1]+1))[2:(nint[1]+1)]
    X_b_knots = c(min[1],max[1])*1.1
    
    Y_knots = seq(min[2],max[2],(max[2]-min[2])/(nint[2]+1))[2:(nint[2]+1)]
    Y_b_knots = c(min[2],max[2])*1.1
    
    if (periodic[1]){
      X = pbs(covariate$X,knots=X_knots,Boundary.knots=X_b_knots)
    } else {
      X = ns(covariate$X,knots=X_knots,Boundary.knots=X_b_knots)
    }
    if (periodic[2]){
      Y = pbs(covariate$Y,knots=Y_knots,Boundary.knots=Y_b_knots)
    } else {
      Y = ns(covariate$Y,knots=Y_knots,Boundary.knots=Y_b_knots)
    }
    return(kr(X,Y,byrow=TRUE))
  }
  
  if (findmin){
    min = min(covariate)
  }
  if (findmax){
    max = max(covariate)
  }
  
  if (perc95max){
    top = sort(covariate)[round(0.95*length(covariate))]
    knots = seq(min,top,(top-min)/(nint+1))[2:(nint+1)]
  } else if (perc95both){
    sorted = sort(covariate)
    len = length(covariate)
    top = sorted[round(0.975*len)]
    bot = sorted[round(0.25*len)]
    knots = seq(bot,top,(top-bot)/(nint+1))[2:(nint+1)]
  } else {
    knots = seq(min,max,(max-min)/(nint+1))[2:(nint+1)]
  }
  
  b_knots = c(min,max)
  
  if (periodic){
    return(pbs(covariate,knots=knots,Boundary.knots=b_knots))
  } else {
    return(ns(covariate,knots=knots,Boundary.knots=b_knots))
  }
}

make_data = function(spiketimes,new_binsize,covariates=NULL,cov_binsize=1,onlycount=FALSE,periodiccovs = c(),max_T=NULL){
  # Takes spiketimes and covariate data frame, creates dataframe with
  # counts and the desired time bin size
  
  if (!onlycount){
    M1 = dim(covariates)[1]
    max_T = cov_binsize*M1
  } else {
    if (is.null(max_T)){
      max_T = max(spiketimes + 1)
    }
  }
  
  spiketimes = spiketimes[which(spiketimes < max_T)]
  #N = length(spiketimes)
  M2 = ceiling(max_T/new_binsize)
  
  counts = array(0,M2)
  spike_ind = ceiling(spiketimes/new_binsize)
  for (i in spike_ind){
    counts[i] = counts[i] + 1
  }
  
  if (onlycount){
    return(counts)
  }
  
  df = data.frame(counts)
  
  for (j in 1:dim(covariates)[2]){
    if (j %in% periodiccovs){
      maxa = max(covariates[,j],na.rm=T)
      if (maxa > 4){
        if (maxa > 8){
          if (maxa > 200){
            covariates[,j] = covariates[,j]*pi/180 - pi
          } else {
            covariates[,j] = covariates[,j]*pi/180
          }
        } else {
          covariates[,j] = covariates[,j] - pi
        }
      }
      covariates[,j] = fix_periodic(covariates,colnames(covariates)[j])
      nonperiodic = unperiodify(covariates[,j])
      df[[paste("X",j,sep="")]] = periodify(approx(x=(1:M1)*cov_binsize,y=nonperiodic,xout=(1:M2)*new_binsize,rule=2)$y)
    } else {
      df[[paste("X",j,sep="")]] = approx(x=(1:M1)*cov_binsize,y=covariates[,j],xout=(1:M2)*new_binsize,rule=2)$y
    }
  }
  
  colnames(df) = c("counts",colnames(covariates))
  return(df)
}

make_design_matrix = function(covariate_df,positive_ind,periodic_ind,XY_ind,nknots){
  firstfindmin = TRUE
  firstper = FALSE
  firstpos = FALSE
  
  if (length(nknots) == 1){nknots = rep(nknots,dim(covariate_df)[2])}
  
  if (1 %in% positive_ind){firstfindmin = FALSE}
  if (1 %in% periodic_ind){firstper = TRUE}
  if (1 %in% XY_ind){firstpos = TRUE}
  
  if (firstpos){
    df = covariate_df[,1:2]
    colnames(df) = c("X","Y")
    X = splinify(df,nint=nknots[1:2],periodic=rep(firstper,2),XY=TRUE,findmin=firstfindmin,findmax=TRUE)
    i = 3
  } else {
    if (1 %in% positive_ind){
      X = splinify(covariate_df[,1],nint=nknots[1],periodic=firstper,findmin=firstfindmin,findmax=TRUE,perc95max=TRUE)
    } else {
      X = splinify(covariate_df[,1],nint=nknots[1],periodic=firstper,findmin=firstfindmin,findmax=TRUE)
    }
    i = 2
  }
  
  while (i <= dim(covariate_df)[2]){
    if (i %in% XY_ind){
      X = cbind(X,splinify(list("X"=covariate_df[,i],"Y"=covariate_df[,i+1]),nint=rep(nknots[i],2),periodic=rep(FALSE,2),XY=TRUE,findmin=TRUE,findmax=TRUE))
      i = i + 2
    } else {
      findmin = TRUE
      #per = FALSE
      if (i %in% positive_ind){findmin = FALSE}
      if (i %in% periodic_ind){
        X = cbind(X,splinify(covariate_df[,i],min=-pi,max=pi,nint=nknots[i],periodic=TRUE))
      } else {
        if (i %in% positive_ind){
          X = cbind(X,splinify(covariate_df[,i],nint=nknots[i],findmin=findmin,findmax=TRUE,perc95max=TRUE))
        } else {
          X = cbind(X,splinify(covariate_df[,i],nint=nknots[i],findmin=findmin,findmax=TRUE))
        }
      }
      i = i + 1
    }
  }
  return(X)
}

BinMean <- function (vec, every, SUM = FALSE, na.rm = FALSE,MATRIX=FALSE) {
  if (MATRIX){
    m = dim(vec)[2]
    firstone = BinMean(vec[,1],every=every,SUM=SUM,na.rm=na.rm)
    k = length(firstone)
    x = matrix(0,nrow=k,ncol=m)
    x[,1] = firstone
    for (i in 2:m){
      x[,i] = BinMean(vec[,i],every=every,SUM=SUM,na.rm=na.rm)
    }
  } else {
    n <- length(vec)
    if (SUM){
      x <- .colSums(vec, every, n %/% every, na.rm)
      r <- n %% every
      if (r) x <- c(x, sum(vec[(n - r + 1):n], na.rm = na.rm))
    } else {
      x <- .colMeans(vec, every, n %/% every, na.rm)
      r <- n %% every
      if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
    }
  }
  return(x)
}

add_boundaries_to_continuous_signal = function(signal,upper_boundary,lower_boundary,type="mirror"){
  new_signal = signal
  n = length(signal)
  if (lower_boundary > upper_boundary){
    print("Oops, the lower boundary must be smaller than the upper boundary")
    return(signal)
  }
  where_too_large = which(signal > upper_boundary)
  where_too_small = which(signal < lower_boundary)
  
  while (length(where_too_large) > 0 || length(where_too_small) > 0){
    if (type == "mirror"){
      new_signal[where_too_large] = upper_boundary - (new_signal[where_too_large] - upper_boundary)
      new_signal[where_too_small] = lower_boundary + (lower_boundary - new_signal[where_too_small])
    } else {
      new_signal[where_too_large] = upper_boundary
      new_signal[where_too_small] = lower_boundary
    }
    
    where_too_large = which(new_signal > upper_boundary)
    where_too_small = which(new_signal < lower_boundary)
  }
  
  return(new_signal)
}

#### Ratemap functions ----
tuningcurve = function(act1,act2=NULL,act3=NULL,act4=NULL,X,binning=80,binsize,smoothingsd=15,minocc=0.01,xrange1=NULL,PERIODIC=F){
  if (is.null(xrange1)){
    xrange1 = range(X)
  }
  tc1 = array(NA,binning)
  if (!is.null(act2)){
    tc2 = array(NA,binning)
  }
  if (!is.null(act3)){
    tc3 = array(NA,binning)
  }
  if (!is.null(act4)){
    tc4 = array(NA,binning)
  }
  occupancyvec = array(NA,binning)
  Xs = floor((X-xrange1[1])/(xrange1[2]*1.001-xrange1[1])*binning)
  for (j in 0:(binning-1)){
    occupancy = which(Xs==j)
    occupancyvec[j+1] = length(occupancy)*binsize
    if (length(occupancy)*binsize > minocc){
      tc1[j+1] = sum(act1[occupancy])/length(occupancy)/binsize
      if (!is.null(act2)){
        tc2[j+1] = sum(act2[occupancy])/length(occupancy)/binsize
      }
      if (!is.null(act3)){
        tc3[j+1] = sum(act3[occupancy])/length(occupancy)/binsize
      }
      if (!is.null(act4)){
        tc4[j+1] = sum(act4[occupancy])/length(occupancy)/binsize
      }
    }
  }
  pded1 = array(NA,binning+2)
  pded1[2:(binning+1)] = tc1
  if (!is.null(act2)){
    pded2 = array(NA,binning+2)
    pded2[2:(binning+1)] = tc2
  }
  if (!is.null(act3)){
    pded3 = array(NA,binning+2)
    pded3[2:(binning+1)] = tc3
  }
  if (!is.null(act4)){
    pded4 = array(NA,binning+2)
    pded4[2:(binning+1)] = tc4
  }
  for (j in 1:binning){
    if (!is.na(pded1[j+1])){
      numnotna = length(which(!is.na(c(pded1[j],pded1[j+2],pded1[j+1]))))
      if (numnotna == 0){
        tc1[j] = NA
      }
    }
    if (!is.null(act2)){
      if (!is.na(pded2[j+1])){
        numnotna = length(which(!is.na(c(pded2[j],pded2[j+2],pded2[j+1]))))
        if (numnotna == 0){
          tc2[j] = NA
        }
      }
    }
    if (!is.null(act3)){
      if (!is.na(pded3[j+1])){
        numnotna = length(which(!is.na(c(pded3[j],pded3[j+2],pded3[j+1]))))
        if (numnotna == 0){
          tc3[j] = NA
        }
      }
    }
    if (!is.null(act4)){
      if (!is.na(pded4[j+1])){
        numnotna = length(which(!is.na(c(pded4[j],pded4[j+2],pded4[j+1]))))
        if (numnotna == 0){
          tc4[j] = NA
        }
      }
    }
  }
  if (length(smoothingsd) == 1){
    smoothingsd = c(smoothingsd,smoothingsd,smoothingsd,smoothingsd)
  }
  if (smoothingsd[1] > 0){
    tc1 = gaussian_smoother1D(tc1,smoothingsd[1],periodic=PERIODIC)
    occupancy_final = gaussian_smoother1D(occupancyvec,smoothingsd[1],periodic=PERIODIC)
  } else {
    occupancy_final = occupancyvec
  }
  if (smoothingsd[2] > 0){
    if (!is.null(act2)){
      tc2 = gaussian_smoother1D(tc2,smoothingsd[2],periodic=PERIODIC)
    }
  }
  if (smoothingsd[3] > 0){
    if (!is.null(act3)){
      tc3 = gaussian_smoother1D(tc3,smoothingsd[3],periodic=PERIODIC)
    }
  }
  if (smoothingsd[4] > 0){
    if (!is.null(act4)){
      tc4 = gaussian_smoother1D(tc4,smoothingsd[4],periodic=PERIODIC)
    }
  }
  results = list(tc1=tc1)
  if (!is.null(act2)){
    results$tc2 = tc2
  }
  if (!is.null(act3)){
    results$tc3 = tc3
  }
  if (!is.null(act4)){
    results$tc4 = tc4
  }
  results$occupancy = occupancy_final
  return(results)
}

ratemapXY = function(act1,act2=NULL,act3=NULL,act4=NULL,X,Y,binning=80,binsize,smoothingsd=15,minocc=0.01,xrange1=NULL,yrange1=NULL,CORNERVAL=NULL){
  if (is.null(xrange1)){
    xrange1 = range(X)
  }
  if (is.null(yrange1)){
    yrange1 = range(Y)
  }
  rm1 = matrix(NA,nrow=binning,ncol=binning)
  if (!is.null(act2)){
    rm2 = matrix(NA,nrow=binning,ncol=binning)
  }
  if (!is.null(act3)){
    rm3 = matrix(NA,nrow=binning,ncol=binning)
  }
  if (!is.null(act4)){
    rm4 = matrix(NA,nrow=binning,ncol=binning)
  }
  occupancymat = matrix(NA,nrow=binning,ncol=binning)
  Xs = floor((X-xrange1[1])/(xrange1[2]*1.001-xrange1[1])*binning)
  Ys = floor((Y-yrange1[1])/(yrange1[2]*1.001-yrange1[1])*binning)
  for (j in 0:(binning-1)){
    for (k in 0:(binning-1)){
      occupancy = intersect(which(Xs==j),which(Ys==k))
      occupancymat[j+1,k+1] = length(occupancy)*binsize
      if (length(occupancy)*binsize > minocc){
        rm1[j+1,k+1] = sum(act1[occupancy])/length(occupancy)/binsize
        if (!is.null(act2)){
          rm2[j+1,k+1] = sum(act2[occupancy])/length(occupancy)/binsize
        }
        if (!is.null(act3)){
          rm3[j+1,k+1] = sum(act3[occupancy])/length(occupancy)/binsize
        }
        if (!is.null(act4)){
          rm4[j+1,k+1] = sum(act4[occupancy])/length(occupancy)/binsize
        }
      }
    }
  }
  pded1 = matrix(NA,nrow=binning+2,ncol=binning+2)
  pded1[2:(binning+1),2:(binning+1)] = rm1
  if (!is.null(act2)){
    pded2 = matrix(NA,nrow=binning+2,ncol=binning+2)
    pded2[2:(binning+1),2:(binning+1)] = rm2
  }
  if (!is.null(act3)){
    pded3 = matrix(NA,nrow=binning+2,ncol=binning+2)
    pded3[2:(binning+1),2:(binning+1)] = rm3
  }
  if (!is.null(act4)){
    pded4 = matrix(NA,nrow=binning+2,ncol=binning+2)
    pded4[2:(binning+1),2:(binning+1)] = rm4
  }
  for (j in 1:binning){
    for (k in 1:binning){
      if (!is.na(pded1[j+1,k+1])){
        numnotna = length(which(!is.na(c(pded1[j,k+1],
                                         pded1[j,k+2],pded1[j,k],pded1[j+2,k+1],
                                         pded1[j+2,k+2],pded1[j+2,k],pded1[j+1,k],
                                         pded1[j+1,k+2]))))
        if (numnotna < 3){
          rm1[j,k] = NA
        }
      }
      if (!is.null(act2)){
        if (!is.na(pded2[j+1,k+1])){
          numnotna = length(which(!is.na(c(pded2[j,k+1],
                                           pded2[j,k+2],pded2[j,k],pded2[j+2,k+1],
                                           pded2[j+2,k+2],pded2[j+2,k],pded2[j+1,k],
                                           pded2[j+1,k+2]))))
          if (numnotna < 3){
            rm2[j,k] = NA
          }
        }
      }
      if (!is.null(act3)){
        if (!is.na(pded3[j+1,k+1])){
          numnotna = length(which(!is.na(c(pded3[j,k+1],
                                           pded3[j,k+2],pded3[j,k],pded3[j+2,k+1],
                                           pded3[j+2,k+2],pded3[j+2,k],pded3[j+1,k],
                                           pded3[j+1,k+2]))))
          if (numnotna < 3){
            rm3[j,k] = NA
          }
        }
      }
      if (!is.null(act4)){
        if (!is.na(pded4[j+1,k+1])){
          numnotna = length(which(!is.na(c(pded4[j,k+1],
                                           pded4[j,k+2],pded4[j,k],pded4[j+2,k+1],
                                           pded4[j+2,k+2],pded4[j+2,k],pded4[j+1,k],
                                           pded4[j+1,k+2]))))
          if (numnotna < 3){
            rm4[j,k] = NA
          }
        }
      }
    }
  }
  if (length(smoothingsd) == 1){
    smoothingsd = c(smoothingsd,0,0,0)
  }
  if (smoothingsd[1] > 0){
    rm1 = gaussian_smoother2D(rm1,smoothingsd[1])
    occupancy_final = gaussian_smoother2D(occupancymat,smoothingsd[1])
  } else {
    occupancy_final = occupancymat
  }
  if (smoothingsd[2] > 0){
    if (!is.null(act2)){
      rm2 = gaussian_smoother2D(rm2,smoothingsd[2])
    }
  }
  if (smoothingsd[3] > 0){
    if (!is.null(act3)){
      rm3 = gaussian_smoother2D(rm3,smoothingsd[3])
    }
  }
  if (smoothingsd[4] > 0){
    if (!is.null(act4)){
      rm4 = gaussian_smoother2D(rm4,smoothingsd[4])
    }
  }
  if (!is.null(CORNERVAL)){
    rm1[1,1] = CORNERVAL
  }
  results = list(rm1=rm1)
  if (!is.null(act2)){
    if (!is.null(CORNERVAL)){
      rm2[1,1] = CORNERVAL
    }
    results$rm2 = rm2
  }
  if (!is.null(act3)){
    if (!is.null(CORNERVAL)){
      rm3[1,1] = CORNERVAL
    }
    results$rm3 = rm3
  }
  if (!is.null(act4)){
    if (!is.null(CORNERVAL)){
      rm4[1,1] = CORNERVAL
    }
    results$rm4 = rm4
  }
  results$occupancy = occupancy_final
  return(results)
}

plot_rms_from_list_normalised = function(X,Y,neuralmat,cell_list,smoothing_sd=15,binning=80,min_occ=0.01,num_plots=9,mpar=c(1,1),main="",noaxis=F,cex=1,maxrate=NULL,RETURNPLOT=FALSE,cbar=NULL,DONTNORMALIZE=F){
  if (!is.null(mpar)){
    par(mfrow=mpar)
  }
  if (length(cell_list) > num_plots){
    cells = cell_list[floor(seq(1,length(cell_list),length.out = num_plots))]
  } else if (length(cell_list) < num_plots) {
    cells = c(cell_list,rep(cell_list[length(cell_list)],num_plots-length(cell_list)))
  } else {
    cells = cell_list
  }
  for (cell in cells){
    if (!cell %in% cell_list){
      print("WHAT")
    }
  }
  im = matrix(NA,nrow=3*binning+2,ncol=1)
  for (i in 1:round(num_plots/3)){
    rms = ratemapXY(act1=neuralmat[,cells[3*(i-1)+1]],
                    act2=neuralmat[,cells[3*(i-1)+2]],
                    act3=neuralmat[,cells[3*(i-1)+3]],
                    X=X,Y=Y,binning=binning,
                    binsize=1/7.5,smoothingsd=smoothing_sd,minocc=min_occ)
    quantiles = c(quantile(c(rms$rm1),0.95,na.rm=T),quantile(c(rms$rm2),0.95,na.rm=T),quantile(c(rms$rm3),0.95,na.rm=T))
    rms$rm1[which(rms$rm1 > quantiles[1])] = quantiles[1]
    rms$rm2[which(rms$rm2 > quantiles[2])] = quantiles[2]
    rms$rm3[which(rms$rm3 > quantiles[3])] = quantiles[3]
    if (DONTNORMALIZE){
      im = cbind(im,rep(NA,3*binning+2),rbind(rms$rm1,rep(NA,binning),rms$rm2,rep(NA,binning),rms$rm3))
    } else {
      im = cbind(im,rep(NA,3*binning+2),rbind(rms$rm1/quantiles[1],rep(NA,binning),rms$rm2/quantiles[2],rep(NA,binning),rms$rm3/quantiles[3]))
    }
  }
  
  if (RETURNPLOT){
    return(im)
  } else {
    image(im,main=main,col=cbar,axes=!noaxis,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)#,axis.args=list(cex.axis=cex))
  }
}

#### QQ rate plots ----

predict_rates = function(X_train, Y_train, X_test, Y_test, covinds, cells, resmat, alpha, family, ALWAYSPOS=F){
  predicted_rates = matrix(NA, nrow=dim(Y_test)[1], ncol=length(cells))
  for (i in 1:length(cells)){
    cell = cells[i]
    if (ALWAYSPOS){
      significant = c(2)
    } else {
      significant = sort(resmat[1,which(resmat[5,,cell] <= alpha),cell])
    }
    if (length(significant) == 0){
      predicted_rates[, i] = rep(mean(Y_train[, cell]), dim(Y_test)[1])
    } else {
      cols_to_use = c()
      for (s in significant){
        cols_to_use = c(cols_to_use, covinds[[s]])
      }
      X_ = X_train[,cols_to_use]
      X_2 = X_test[,cols_to_use]
      df = data.frame(cbind(Y_train[, cell], X_))
      df_test = data.frame(cbind(Y_test[, cell], X_2))
      colnames(df) = c("y", paste("x", 1:dim(X_)[2], sep=""))
      colnames(df_test) = colnames(df)
      fit = glm(y ~ .,family=family, data=df)
      pred = predict(fit, newdata=df_test, type="response")
      predicted_rates[, i] = pred #fit$fitted.values
    }
  }
  return(predicted_rates)
}

split_cor_average = function(vector1, vector2, num_groups){
  cor_vals = array(NA, num_groups)
  for (i in 1:num_groups){
    indices = seq(i, length(vector1), by = num_groups)
    cor_vals[i] = cor(vector1[indices], vector2[indices])
  }
  return(mean(cor_vals))
}

QQ_plots = function(X_train, Y_train, X_test, Y_test, covinds, cells, resmat, alpha, family, binsize, mfrow, num_groups, ALWAYSPOS=F, QQ=T){
  predicted_rates = predict_rates(X_train, Y_train, X_test, Y_test, covinds, cells, resmat, alpha, family, ALWAYSPOS=ALWAYSPOS)
  predicted_binned = BinMean(predicted_rates, every=binsize, MATRIX=T)*7.25
  rate_binned = BinMean(Y_test[,cells], every=binsize, MATRIX=T)*7.25
  quantiles = seq(0,1,by=0.05)
  old_pty = par("pty")
  r2s = c()
  par(mfrow=mfrow,pty="s")
  for (i in 1:length(cells)){
    if (QQ){
      predicted = quantile(predicted_binned[, i], quantiles)
      truth = quantile(rate_binned[, i], quantiles)
      fit2 = lm(truth~predicted)
      a = fit2$coefficients[1]
      b = fit2$coefficients[2]
      x = c(predicted[1],predicted[length(predicted)])
      y = c(truth[1],truth[length(truth)])
      r2 = cor(sort(predicted_binned[, i]),sort(rate_binned[, i]))
      maxval = max(c(predicted,truth))
      print(maxval)
      plot(predicted, truth, main=paste("Cor, qq:", round(r2*1000)/1000), xlim=c(0,maxval),ylim=c(0,maxval))
      #lines(x, a + b * x,lty=2,col="red")
      lines(c(0,maxval), c(0, maxval),lty=2,col="red")
    } else {
      fit = lm(rate_binned[, i]~predicted_binned[, i])
      a1 = fit$coefficients[1]
      b1 = fit$coefficients[2]
      x1 = range(predicted_binned[, i])
      y1 = range(rate_binned[, i])
      r2 = split_cor_average(predicted_binned[, i], rate_binned[, i], num_groups)
      plot(predicted_binned[, i],rate_binned[, i], main=paste("Cor rates:", round(r2*1000)/1000))
      #lines(x1, a1 + b1 * x1,lty=2,col="red")
      lines(y1, y1,lty=2,col="red")
    }
    r2s = c(r2s, r2)
    #lines(c(truth[1], truth[11]), c(truth[1], truth[11]),lty=2,col="red")
  }
  par(pty=old_pty)
  print(paste("Mean R2:", mean(r2s)))
}

#### Cross-validation functions ----
make_folds = function(n, nfolds, chunksize){
  folds = array(0,n)
  chunksofeach = ceiling(n/nfolds/chunksize)
  if (chunksofeach == 1 & chunksize*nfolds > n){
    toremove = nfolds*chunksize-n
    folds = c(rep(1:(nfolds-toremove),each=chunksize),rep((nfolds-toremove+1):nfolds,each=chunksize-1))
  } else {
    folds = rep(rep(1:nfolds,each=chunksize),chunksofeach)[1:n]
  }
  return(folds)
}

OLD_llh_CV = function(resp,preds,comparemean,NORMAL=FALSE,family=family){
  if (family == "poisson"){
    normal_llh = sum(resp*log(preds)-preds)#-log(factorial(resp))
    if (NORMAL){return(normal_llh)}
    mean_llh = sum(resp*log(comparemean)-comparemean)#-log(factorial(resp))
    diff_llh = normal_llh - mean_llh
    scaled_llh = diff_llh/log(2)#/sum(resp)
  } else if (family == "binomial"){
    normal_llh = sum(resp*log(preds/(1-preds))+log(1-preds))
    if (NORMAL){return(normal_llh)}
    mean_llh = sum(resp*log(comparemean/(1-comparemean))+log(1-comparemean))
    diff_llh = normal_llh - mean_llh
    scaled_llh = diff_llh/log(2)#/sum(resp)
  }
  
  return(scaled_llh)
}

llh_CV = function(resp,preds,family=family){
  if (family == "poisson"){
    llh = sum(resp*log(preds)-preds)#-log(factorial(resp))
  } else if (family == "binomial"){
    llh = sum(resp*log(preds/(1-preds))+log(1-preds))
  }
  return(llh)
}

OLD_CV_folds = function(Y,X,folds,NORMAL=F,ONLY_INTERCEPT=F,family=family,nfolds=10,SKIP_NEIGHBOR=F,indelses=NULL){
  #fam = "poisson"
  n = length(Y)
  
  if (ONLY_INTERCEPT){
    X = rep(1,n)
  }
  
  dat = data.frame(cbind(Y,X))
  colnames(dat)[1] = "Y"
  
  llh_folds = array(0,nfolds)
  
  for (k in 1:nfolds){
    indk = which(folds == k)
    
    if (!is.null(indelses)){
      indelse = indelses[,k]
    } else {
      if (SKIP_NEIGHBOR == TRUE){
        if (k == nfolds){
          indk_plus_1 = which(folds == 1)
        } else {
          indk_plus_1 = which(folds == k+1)
        }
        if (k == 1){
          indk_minus_1 = which(folds == nfolds)
        } else {
          indk_minus_1 = which(folds == k-1)
        }
        indelse = setdiff(1:n,c(indk,indk_plus_1,indk_minus_1))
      } else if (SKIP_NEIGHBOR > 0){
        chunkstarts = c(indk[1],indk[which(diff(indk) > 1)+1])
        chunkends = c(indk[which(diff(indk) > 1)],indk[length(indk)])
        indsafter = c()
        indsbefore = c()
        for (i in 1:SKIP_NEIGHBOR){
          indsafter = c(indsafter,chunkends+i)
          indsbefore = c(indsbefore,chunkstarts-i)
        }
        
        indsafter[which(indsafter > n)] = indsafter[which(indsafter > n)] - n
        indsbefore[which(indsbefore < 1)] = indsbefore[which(indsbefore < 1)] + n
        
        #if (k == 1){
        #  indsbefore = c(indsbefore,(n-SKIP_NEIGHBOR+1):n)
        #}
        #if (k == nfolds){
        #  indsafter = c(indsafter,1:SKIP_NEIGHBOR)
        #}
        
        indelse = setdiff(1:n,c(indk,indsafter,indsbefore))
        if (length(indelse) < length(indk)){
          print("Oops, you're removing too much, you don't want your training set to be smaller than your test set!")
          return(NULL)
        }
      } else {
        indelse = setdiff(1:n,indk)
      }
    }
    test = dat[indk,]
    train = dat[indelse,]
    
    if (ONLY_INTERCEPT){
      preds = rep(mean(train$Y),length(indk))
    } else {
      fitki = glm(Y~.,data=train,family=family)
      preds = predict(fitki,newdata=test,type="response")
      
    }
    resp = test$Y
    llh_folds[k] = llh_CV(resp,preds,mean(test$Y),NORMAL=NORMAL,family=family)
  }
  return(llh_folds)
}

CV_folds = function(Y,X,config,ONLY_INTERCEPT=F,indelses=NULL){
  #fam = "poisson"
  n = length(Y)
  nfolds = config$nfolds
  indice_split = config$indice_split
  family = config$family
  
  if (ONLY_INTERCEPT){
    X = rep(1,n)
  }
  
  dat = data.frame(cbind(Y,X))
  colnames(dat)[1] = "Y"
  
  llh_folds = array(0,nfolds)
  
  for (k in 1:nfolds){
    indk = indice_split[[k]]$indk
    if (!is.null(indelses)){
      indelse = indelses[k]
    } else {
      indelse = indice_split[[k]]$indelse
    }
    
    test = dat[indk,]
    train = dat[indelse,]
    
    if (ONLY_INTERCEPT){
      preds = rep(mean(train$Y),length(indk))
    } else {
      fitki = glm(Y~.,data=train,family=family)
      preds = predict(fitki,newdata=test,type="response")
      
    }
    resp = test$Y
    llh_folds[k] = llh_CV(resp,preds,family=family)
  }
  return(llh_folds)
}

indice_splits = function(n, folds, nfolds, SKIP_NEIGHBOR=F){
  ind_splits = list()
  for (k in 1:nfolds){
    indk = which(folds == k)
    
    if (SKIP_NEIGHBOR == TRUE){
      if (k == nfolds){
        indk_plus_1 = which(folds == 1)
      } else {
        indk_plus_1 = which(folds == k+1)
      }
      if (k == 1){
        indk_minus_1 = which(folds == nfolds)
      } else {
        indk_minus_1 = which(folds == k-1)
      }
      indelse = setdiff(1:n,c(indk,indk_plus_1,indk_minus_1))
    } else if (SKIP_NEIGHBOR > 0){
      chunkstarts = c(indk[1],indk[which(diff(indk) > 1)+1])
      chunkends = c(indk[which(diff(indk) > 1)],indk[length(indk)])
      indsafter = c()
      indsbefore = c()
      for (i in 1:SKIP_NEIGHBOR){
        indsafter = c(indsafter,chunkends+i)
        indsbefore = c(indsbefore,chunkstarts-i)
      }
      
      indsafter[which(indsafter > n)] = indsafter[which(indsafter > n)] - n
      indsbefore[which(indsbefore < 1)] = indsbefore[which(indsbefore < 1)] + n
      
      indelse = setdiff(1:n,c(indk,indsafter,indsbefore))
      if (length(indelse) < length(indk)){
        print("Oops, you're removing too much, you don't want your training set to be smaller than your test set!")
        return(NULL)
      }
    } else {
      indelse = setdiff(1:n,indk)
    }
    ind_splits[[k]] = list("indk"=indk, "indelse"=indelse)
  }
  return(ind_splits)
}

sample_folds = function(n,nfolds,chunksize){
  chunksofeach = ceiling(n/nfolds/chunksize)
  return(rep(sample(rep(1:nfolds,chunksofeach),nfolds*chunksofeach,replace=F),each=chunksize)[1:n])
}

#### Cyclic shift functions ----
cyclic_shift = function(vec,silent=TRUE){
  temp = FALSE
  if (is.null((dim(vec)))){
    n = length(vec)
    new_vec = array(0,dim=n)
  } else if (length(dim(vec))==2){
    n = dim(vec)[1]
    ncols = dim(vec)[2]
    new_vec = matrix(0,nrow=n,ncol=ncols)
  } else if (is.null((dim(dim(vec))))){
    n = length(vec)
    new_vec = array(0,dim=n)
    temp = TRUE
  } else {
    n = dim(vec)[1]
    ncols = dim(vec)[2]
    new_vec = matrix(0,nrow=n,ncol=ncols)
    if (length(dim(vec))==1){
      #print("hmm")
      temp = TRUE
    }
  }
  lag = floor(runif(1,min=0,max=n))
  if (is.null((dim(vec))) || temp){
    new_vec[(lag+1):n] = vec[1:(n-lag)]
    if (lag > 0){
      new_vec[1:lag] = vec[(n-lag+1):n]
    }
  } else {
    for (i in 1:(n-lag)){
      for (j in 1:ncols){
        new_vec[lag+i,j] = vec[i,j]
      }
    }
    #new_vec[(lag+1):n,1:ncols] = vec[1:(n-lag),1:ncols]
    if (lag > 0){
      for (i in 1:lag){
        for (j in 1:ncols){
          new_vec[i,j] = vec[n-lag+i,j]
        }
      }
      #new_vec[1:lag,1:ncols] = vec[(n-lag+1):n,1:ncols]
    }
  }
  if (!silent){
    print(lag)
  }
  return(new_vec)
}

cyclic_shift_ind = function(n,minlag=0,RETURNLAG=FALSE){
  lag = floor(runif(1,min=minlag,max=n-minlag))
  if (lag == 0){
    ind = 1:n
  } else {
    ind = c((lag+1):n,1:lag)
  }
  if (RETURNLAG){
    return(list(lag=lag,ind=ind))
  } else {
    return(ind)
  }
}

#### Forward selection ----
signrank = function(diffs,ONLY_POS=FALSE){
  nonadiffs = diffs[which(!is.na(diffs))]
  abses = abs(nonadiffs)
  ranks = order(order(abses))
  if (ONLY_POS){
    signrank = sum(ranks[which(nonadiffs > 0)])
  } else {
    signrank = sum(ranks*sign(nonadiffs))
  }
  return(signrank)
}

CV_scores = function(Y,X,cov_inds,config,covariates,reverse=NULL,indelses=NULL){
  #Calculate cross validated scores for selected model (covariates)
  if (length(covariates) == 0){
    scores = suppressWarnings(CV_folds(Y,NULL,config,ONLY_INTERCEPT=T,indelses=indelses))
    return(scores)
  }
  
  num_cols = 0
  for (i in covariates){
    num_cols = num_cols + length(cov_inds[[i]])
  }
  X_ = matrix(NA, nrow=length(Y), ncol=num_cols)
  
  k = 1
  for (i in covariates){
    if (! is.null(reverse)){
      if (reverse == i){
        for (j in cov_inds[[i]]){
          X_[, k] = rev(X[, j])
          k = k + 1
        }
      }
    } else {
      for (j in cov_inds[[i]]){
        X_[, k] = X[, j]
        k = k + 1
      }
    }
  }
  scores = suppressWarnings(CV_folds(Y,X_,config,indelses=indelses))
  return(scores)
}

scores_from_list = function(cells, Y, X, covinds, config){
  if (length(cells) == 1){
    scores = array(NA, 8)
  } else {
    scores = matrix(NA, nrow=8, ncol=length(cells))
  }
  
  start = proc.time()[[3]]
  for (i in 1:length(cells)){
    cell = cells[i]
    for (i1 in c(0,1)){
      for (i2 in c(0,1)){
        for (i3 in c(0,1)){
          covs = c()
          row = 1
          if (i1 == 1){
            covs = c(covs, 1)
            row = row + 1
          }
          if (i2 == 1){
            covs = c(covs, 2)
            row = row + 2
          }
          if (i3 == 1){
            covs = c(covs, 3)
            row = row + 4
          }
          if (length(cells) == 1){
            scores[row] = mean(CV_scores(Y,X,covinds,config,covs))
          } else {
            scores[row, i] = mean(CV_scores(Y[,cell],X,covinds,config,covs))
          }
        }
      }
    }
    if (length(cells) > 1){
      print_progress(i, length(cells), start)
    }
  }
  if (length(cells) == 1){
    scores[1:8] = scores[c(1,2,3,5,4,6,7,8)]
  } else {
    scores[1:8,] = scores[c(1,2,3,5,4,6,7,8),]
  }
  return(scores)
}

score_vecs_from_list = function(cells, Y, X, covinds, config){
  if (length(cells) == 1){
    scores = array(NA, c(nfolds, 8))
  } else {
    scores = array(NA, c(nfolds, 8, length(cells)))
  }
  
  start = proc.time()[[3]]
  for (i in 1:length(cells)){
    cell = cells[i]
    for (i1 in c(0,1)){
      for (i2 in c(0,1)){
        for (i3 in c(0,1)){
          covs = c()
          row = 1
          if (i1 == 1){
            covs = c(covs, 1)
            row = row + 1
          }
          if (i2 == 1){
            covs = c(covs, 2)
            row = row + 2
          }
          if (i3 == 1){
            covs = c(covs, 3)
            row = row + 4
          }
          if (length(cells) == 1){
            scores[,row] = CV_scores(Y,X,covinds,config,covs)
          } else {
            scores[,row, i] = CV_scores(Y[,cell],X,covinds,config,covs)
          }
        }
      }
    }
    if (length(cells) > 1){
      print_progress(i, length(cells), start)
    }
  }
  if (length(cells) == 1){
    scores[,1:8] = scores[,c(1,2,3,5,4,6,7,8)]
  } else {
    scores[,1:8,] = scores[,c(1,2,3,5,4,6,7,8),]
  }
  return(scores)
}

one_step = function(Y,X,cov_inds,config,methods,current=c(),indelses=NULL,FAKE=F){
  #cov_inds is a list with length = total number of considered covariates
  #Each element in the list is an array of column-indices from X
  #Current is an array of list-indices
  alpha = config$alpha
  nfolds = config$nfolds
  family = config$family
  permshuffles = config$permshuffles
  signshuffles = config$signshuffles
  cyclic_skip = config$cyclic_skip
  
  
  #Select the one with the largest llh increase
  #Do a model label permutation test for this one
  
  if (length(current) == 0){
    current = list()
    new_candidates = 1:length(cov_inds)
    current_llh_folds = CV_folds(Y,NULL,config,ONLY_INTERCEPT=T,indelses=indelses)
    current_llh = mean(current_llh_folds,na.rm=T)
    new_llh_folds_ = array(0,dim=c(length(new_candidates),nfolds))
    if ("signrankmax2" %in% methods){
      new_llh_folds_2 = array(0,dim=c(length(new_candidates),nfolds))
    }
    for (i in new_candidates){
      X_ = X[,cov_inds[[i]]]
      new_llh_folds_[i,] = CV_folds(Y,X_,config,indelses=indelses)
      if ("signrankmax2" %in% methods){
        X_ = X[rev(1:dim(X)[1]),cov_inds[[i]]]
        new_llh_folds_2[i,] = CV_folds(Y,X_,config,indelses=indelses)
      }
    }
    means = apply(new_llh_folds_,1,mean,na.rm=T)
    selected = new_candidates[which.max(means)]
    new_llh_folds = new_llh_folds_[selected,]
    Tval = means[selected]-current_llh
    signedrank = signrank(new_llh_folds-current_llh_folds)
    
    if ("signrankmax2" %in% methods & length(methods) == 2){
      ranks = array(NA, length(new_candidates))
      for (i in 1:length(new_candidates)){
        ranks[i] = signrank(new_llh_folds_[i,] - new_llh_folds_2[i,])
      }
      selected = new_candidates[which.max(ranks)]
      Tval = ranks[selected]
    }
    
    if ("lrt" %in% methods || "lrtbonf" %in% methods || "cyclic" %in% methods){
      X_ = X[,cov_inds[[selected]]]
      fit = glm(Y~X_,family=family)
      fit0 = glm(Y~1,family=family)
      current$llh = logLik(fit)
      current_llh = current$llh
      logLik0 = logLik(fit0)
      McFadden = 1 - current_llh/logLik0
      lr = current_llh-logLik0
    }
    
    current$llh_folds = new_llh_folds
    current$model = c(selected)
    current$Tval = Tval
    pvals = array(NA,length(methods))
    effect_sizes = array(NA,length(methods))
    
    for (m in 1:length(methods)){
      method = methods[m]
      if (method == "CV"){
        if (Tval > 0){
          pvals[m] = 0
        } else {
          pvals[m] = 1
        }
      } else if (method == "cyclic"){
        effect_sizes[m] = McFadden
        if (FAKE || Tval > 0){#(Tval > 0){
          cyclic_t = lr
          
          # Then shuffle and test
          ts_cyclic = array(0,permshuffles)
          for (i in 1:permshuffles){
            X_ = cyclic_shift(X[,cov_inds[[selected]]])
            fit = glm(Y~X_,family=family)
            ts_cyclic[i] = logLik(fit)-logLik0
            
            if (!FAKE){
              # For speedup:
              pval_sofar = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
              if (pval_sofar > alpha){
                EARLY_BREAK = TRUE
                pval_estimate = min(c(1, pval_sofar / i * permshuffles))
                break
              }
            }
          }
          if (EARLY_BREAK){
            pvals[m] = pval_estimate
          } else {
            pvals[m] = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
          }
        } else {
          pvals[m] = 1
        }
      } else if (method == "lrt"){
        effect_sizes[m] = McFadden
        if (TRUE){#(Tval > 0){
          # Do likelihood ratio test
          pvals[m] = 1-pchisq(2*lr,df=length(cov_inds[[selected]]))
        } else {
          pvals[m] = 1
        }
      } else if (method == "lrtbonf"){
        effect_sizes[m] = McFadden
        if (TRUE){#(Tval > 0){
          # Do likelihood ratio test with Bonferroni correction
          pvals[m] = min(c((1-pchisq(2*lr,df=length(cov_inds[[selected]])))*length(new_candidates),1))
        } else {
          pvals[m] = 1
        }
      } else if (method == "cyclicskip"){
        n = length(Y)
        keep = c((cyclic_skip+1):(round(n/2)-cyclic_skip),(round(n/2)+cyclic_skip+1):(n-cyclic_skip))
        X_ = X[keep,cov_inds[[selected]]]
        Y_ = Y[keep]
        fit = glm(Y_~X_,family=family)
        fit0 = glm(Y_~1,family=family)
        current$skip_llh = logLik(fit)
        logLik0 = logLik(fit0)
        effect_sizes[m] = 1 - current$skip_llh / logLik0
        if (FAKE || Tval > 0){#(Tval > 0){
          cyclic_t = current$skip_llh-logLik0
          
          # Then shuffle and test
          ts_cyclic = array(0,permshuffles)
          EARLY_BREAK = FALSE
          for (i in 1:permshuffles){
            li = cyclic_shift_ind(n,minlag=2*cyclic_skip,RETURNLAG=T)
            keep = setdiff((cyclic_skip+1):(n-cyclic_skip),(n-li$lag-cyclic_skip+1):(n-li$lag+cyclic_skip))
            X_ = X[li$ind[keep], cov_inds[[selected]]]
            Y_ = Y[keep]
            fit0 = glm(Y_~1,family=family)
            logLik0 = logLik(fit0)
            fit = glm(Y_~X_,family=family)
            ts_cyclic[i] = logLik(fit)-logLik0
            
            if (!FAKE){
              # For speedup:
              pval_sofar = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
              if (pval_sofar > alpha){
                EARLY_BREAK = TRUE
                pval_estimate = min(c(1, pval_sofar / i * permshuffles))
                break
              }
            }
          }
          if (EARLY_BREAK){
            pvals[m] = pval_estimate
          } else {
            pvals[m] = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
          }
        } else {
          pvals[m] = 1
        }
      } else if (method == "signrank"){
        effect_sizes[m] = signedrank
        if (TRUE){#(Tval > 0){#
          test_stats = array(0,signshuffles)
          for (i in 1:signshuffles){
            u = round(runif(nfolds))*2-1
            test_stats[i] = signrank(u*(new_llh_folds-current_llh_folds))
          }
          pvals[m] = (length(which(test_stats >= signrank(new_llh_folds-current_llh_folds)))+1)/(signshuffles+1)
        } else {
          pvals[m] = 1
        }
      } else if (method == "signrankmax"){
        effect_sizes[m] = signedrank
        if (TRUE){#(Tval > 0){#
          test_stats = matrix(0,nrow=signshuffles,ncol=length(new_candidates))
          for (i in 1:signshuffles){
            u = round(runif(nfolds))*2-1
            for (j in 1:length(new_candidates)){
              test_stats[i,j] = signrank(u*(new_llh_folds_[j,]-current_llh_folds))
            }
          }
          maxineach = apply(test_stats,1,max)
          pvals[m] = (length(which(maxineach >= signrank(new_llh_folds-current_llh_folds)))+1)/(signshuffles+1)
        } else {
          pvals[m] = 1
        }
      } else if (method == "signrankmax2"){
        signedrankrev = signrank(new_llh_folds_[selected,]-new_llh_folds_2[selected,])
        effect_sizes[m] = signedrankrev
        if (TRUE){#(Tval > 0){#
          test_stats = matrix(0,nrow=signshuffles,ncol=length(new_candidates))
          for (i in 1:signshuffles){
            u = round(runif(nfolds))*2-1
            for (j in 1:length(new_candidates)){
              test_stats[i,j] = signrank(u*(new_llh_folds_[j,]-new_llh_folds_2[j,]))
            }
          }
          maxineach = apply(test_stats,1,max)
          pvals[m] = (length(which(maxineach >= signedrankrev))+1)/(signshuffles+1)
        } else {
          pvals[m] = 1
        }
      } else {
        print("METHOD NOT IMPLEMENTED!!!!!")
        print(method)
      }
    }
    current$pvals = pvals
    current$effect_sizes = effect_sizes
    candidate = selected
  } else {
    current_cols = c()
    for (i in current$model){
      current_cols = c(current_cols,cov_inds[[i]])
    }
    current_llh_folds = current$llh_folds
    new_candidates = setdiff(1:length(cov_inds),current$model)
    new_llh_folds_ = array(0,dim=c(length(new_candidates),nfolds))
    if ("signrankmax2" %in% methods){
      new_llh_folds_2 = array(0,dim=c(length(new_candidates),nfolds))
    }
    j = 1
    for (i in new_candidates){
      X_ = X[,c(current_cols,cov_inds[[i]])]
      new_llh_folds_[j,] = CV_folds(Y,X_,config,indelses=indelses)
      if ("signrankmax2" %in% methods){
        X_ = cbind(X[,current_cols],X[rev(1:dim(X)[1]),cov_inds[[i]]])
        new_llh_folds_2[j,] = CV_folds(Y,X_,config,indelses=indelses)
      }
      j = j + 1
    }
    means = apply(new_llh_folds_,1,mean,na.rm=T)
    selected_ind = which.max(means)
    selected = new_candidates[selected_ind]
    new_llh_folds = new_llh_folds_[selected_ind,]
    
    current_llh = mean(current_llh_folds,na.rm=T)
    Tval = means[selected_ind]-current_llh
    signedrank = signrank(new_llh_folds-current_llh_folds)
    
    if ("signrankmax2" %in% methods & length(methods) == 2){
      ranks = array(NA, length(new_candidates))
      for (i in 1:length(new_candidates)){
        ranks[i] = signrank(new_llh_folds_[i,] - new_llh_folds_2[i,])
      }
      selected_ind = which.max(ranks)
      selected = new_candidates[selected_ind]
      new_llh_folds = new_llh_folds_[selected_ind,]
      Tval = ranks[selected_ind]
    }
    
    if ("lrt" %in% methods || "lrtbonf" %in% methods || "cyclic" %in% methods){
      X_ = X[,c(current_cols,cov_inds[[selected]])]
      fit = glm(Y~X_,family=family)
      new_llh = logLik(fit)
      current_llh = current$llh
      lr = new_llh - current_llh
      McFadden = 1 - new_llh/current_llh
      current$llh = new_llh
    }
    
    current$llh_folds = new_llh_folds
    current$model = c(current$model,selected)
    pvals = array(NA,length(methods))
    effect_sizes = array(NA,length(methods))
    
    for (m in 1:length(methods)){
      if (!is.na(current$pvals[m])){
        if (current$pvals[m] <= alpha){
          method = methods[m]
          if (method == "CV"){
            if (Tval > 0){
              pvals[m] = 0
            } else {
              pvals[m] = 1
            }
          } else if (method == "cyclic"){
            effect_sizes[m] = McFadden
            if (Tval > 0){
              cyclic_t = lr
              
              # Then shuffle and test
              ts_cyclic = array(0,permshuffles)
              EARLY_BREAK = FALSE
              for (i in 1:permshuffles){
                X_ = cbind(X[,current_cols],cyclic_shift(X[,cov_inds[[selected]]]))
                fit = glm(Y~X_,family=family)
                ts_cyclic[i] = logLik(fit)-current_llh
                
                # For speedup:
                pval_sofar = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
                if (pval_sofar > alpha){
                  EARLY_BREAK = TRUE
                  pval_estimate = min(c(1, pval_sofar / i * permshuffles))
                  break
                }
              }
              if (EARLY_BREAK){
                pvals[m] = pval_estimate
              } else {
                pvals[m] = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
              }
            } else {
              pvals[m] = 1
            }
          } else if (method == "lrt"){
            effect_sizes[m] = McFadden
            if (TRUE){#(Tval > 0){
              # Do likelihood ratio test
              pvals[m] = 1-pchisq(2*lr,df=length(cov_inds[[selected]]))
            } else {
              pvals[m] = 1
            }
          } else if (method == "lrtbonf"){
            effect_sizes[m] = McFadden
            if (TRUE){#(Tval > 0){
              # Do likelihood ratio test with Bonferroni correction
              pvals[m] = min(c((1-pchisq(2*lr,df=length(cov_inds[[selected]])))*length(new_candidates),1))
            } else {
              pvals[m] = 1
            }
          } else if (method == "cyclicskip"){
            n = length(Y)
            keep = c((cyclic_skip+1):(round(n/2)-cyclic_skip),(round(n/2)+cyclic_skip+1):(n-cyclic_skip))
            X_ = X[keep,c(current_cols,cov_inds[[selected]])]
            Y_ = Y[keep]
            fit = glm(Y_~X_,family=family)
            new_llh = logLik(fit)
            effect_sizes[m] = 1 - new_llh / current$skip_llh
            if (Tval > 0){#(Tval > 0){
              cyclic_t = new_llh - current$skip_llh
              current$skip_llh = new_llh
              
              # Then shuffle and test
              ts_cyclic = array(0,permshuffles)
              EARLY_BREAK = FALSE
              for (i in 1:permshuffles){
                li = cyclic_shift_ind(n,minlag=2*cyclic_skip,RETURNLAG=T)
                keep = setdiff((cyclic_skip+1):(n-cyclic_skip),(n-li$lag-cyclic_skip+1):(n-li$lag+cyclic_skip))
                Y_ = Y[keep]
                X_ = cbind(X[keep,current_cols],X[li$ind[keep],cov_inds[[selected]]])
                X_0 = X[keep,current_cols]
                Y_ = Y[keep]
                fit = glm(Y_~X_,family=family)
                fit0 = glm(Y_~X_0,family=family)
                ts_cyclic[i] = logLik(fit)-logLik(fit0)
                
                # For speedup:
                pval_sofar = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
                if (pval_sofar > alpha){
                  EARLY_BREAK = TRUE
                  pval_estimate = min(c(1, pval_sofar / i * permshuffles))
                  break
                }
              }
              if (EARLY_BREAK){
                pvals[m] = pval_estimate
              } else {
                pvals[m] = min(c((length(which(ts_cyclic >= cyclic_t))+1)/(permshuffles+1)*length(new_candidates),1)) 
              }
            } else {
              pvals[m] = 1
            }
          } else if (method == "signrank"){
            effect_sizes[m] = signedrank
            if (TRUE){
              test_stats = array(0,signshuffles)
              for (i in 1:signshuffles){
                u = round(runif(nfolds))*2-1
                test_stats[i] = signrank(u*(new_llh_folds-current_llh_folds))
              }
              pvals[m] = (length(which(test_stats >= signedrank))+1)/(signshuffles+1)
            } else {
              pvals[m] = 1
            }
          } else if (method == "signrankmax"){
            effect_sizes[m] = signedrank
            if (TRUE){#(Tval > 0){#(TRUE){#
              test_stats = matrix(0,nrow=signshuffles,ncol=length(new_candidates))
              for (i in 1:signshuffles){
                u = round(runif(nfolds))*2-1
                for (j in 1:length(new_candidates)){
                  test_stats[i,j] = signrank(u*(new_llh_folds_[j,]-current_llh_folds))
                }
              }
              maxineach = apply(test_stats,1,max)
              pvals[m] = (length(which(maxineach >= signedrank))+1)/(signshuffles+1)
            } else {
              pvals[m] = 1
            }
          } else if (method == "signrankmax2"){
            signedrankrev = signrank(new_llh_folds_[selected_ind,]-new_llh_folds_2[selected_ind,])
            effect_sizes[m] = signedrankrev
            if (TRUE){#(Tval > 0){#(TRUE){#
              test_stats = matrix(0,nrow=signshuffles,ncol=length(new_candidates))
              for (i in 1:signshuffles){
                u = round(runif(nfolds))*2-1
                for (j in 1:length(new_candidates)){
                  test_stats[i,j] = signrank(u*(new_llh_folds_[j,]-new_llh_folds_2[j,]))
                }
              }
              maxineach = apply(test_stats,1,max)
              pvals[m] = (length(which(maxineach >= signedrankrev))+1)/(signshuffles+1)
            } else {
              pvals[m] = 1
            }
          }
        }
      }
    }
    candidate = selected
    current$pvals = pvals
    current$effect_sizes = effect_sizes
  }
  if (length(methods)==1 & methods[1] == "CV"){
    return(list("candidate"=candidate,"current"=current,
                CVdiffs1=new_llh_folds_[1,]-current_llh_folds,CVdiffsbest=new_llh_folds-current_llh_folds))
  }
  return(list("candidate"=candidate,"current"=current))
}

config.default = function(){
  config = list()
  config$N = 12000
  config$alpha = 0.05
  config$nfolds = 20
  config$chunksize = 150
  config$family = "binomial"
  config$SILENT = T
  config$SKIP_NEIGHBOR = T
  config$permshuffles = 119
  config$signshuffles = 999
  config$maxstep = 3
  config$cyclic_skip = 75
  return(config)
}

prepare_folds = function(config){
  config$folds = make_folds(config$N, config$nfolds, config$chunksize)
  config$indice_split = indice_splits(config$N, config$folds, config$nfolds, config$SKIP_NEIGHBOR)
  return(config)
}

forward_selection = function(Y,X,cov_inds,methods=c("CV","signrank","cyclic"),config=config.default(),FAKE=F){
  alpha = config$alpha
  SILENT = config$SILENT
  maxstep = config$maxstep
  
  ncov = length(cov_inds)
  nrows1 = length(methods) + 1
  nrows2 = length(methods) * 2 + 1
  results = matrix(NA,nrow=nrows2,ncol=maxstep)
  
  #print("Starting first step!")
  current = c()
  step = one_step(Y,X,cov_inds,config,methods,current=current,FAKE=FAKE)
  
  current = step$current
  results[1, 1] = step$candidate
  results[2:nrows1, 1] = current$pvals
  results[(nrows1 + 1):nrows2, 1] = current$effect_sizes
  
  finished = FALSE
  if (sum(current$pvals <= alpha,na.rm=T) == 0 || ncov == 1 || maxstep == 1){finished = TRUE}
  step_num = 1
  while (!finished){
    step_num = step_num + 1
    if (!SILENT){
      print(paste("Starting step number ",step_num,"!",sep=""))
    }
    step = one_step(Y,X,cov_inds,config,methods,current=current,FAKE=FAKE)
    current = step$current
    results[1, step_num] = step$candidate
    results[2:nrows1, step_num] = current$pvals
    results[(nrows1 + 1):nrows2, step_num] = current$effect_sizes
    
    if (sum(current$pvals <= alpha,na.rm=T) == 0 || step_num == ncov || step_num == maxstep){
      finished = TRUE
    }
  }
  return(results)
}

forward_selection_CV_variants = function(Y,X,cov_inds,fold_nums,chunksizes,neighbor_skips,RANDOMFOLDS=FALSE,config=config.default(),foldses=NULL,indelseses=NULL,RANDOMFOLDSONLAST=FALSE,BUTSIGNEDRANK=FALSE){
  if (TRUE){
    alpha = config$alpha
    family = config$family
  }#use configs
  n = length(Y)
  num_versions = length(fold_nums)
  results = matrix(NA,nrow=num_versions,ncol=3)
  CVdiffs = matrix(NA,nrow=num_versions,ncol=max(fold_nums)*2)
  # Which covariate is selected, and if the mean CV score is better than nullmodel
  
  for (i in 1:num_versions){
    if (i == num_versions & RANDOMFOLDSONLAST){
      RANDOMFOLDS = TRUE
    }
    nfolds = fold_nums[i]
    if (!is.null(foldses)){
      folds = foldses[[i]]
    } else {
      folds = array(0,n)
      chunksize = chunksizes[i]
      if (chunksize == 1 & RANDOMFOLDS){
        folds = sample(1:nfolds,n,replace=T)
      } else {
        chunksofeach = ceiling(n/nfolds/chunksize)
        if (chunksofeach == 1 & chunksize*nfolds > n){
          toremove = nfolds*chunksize-n
          folds = c(rep(1:(nfolds-toremove),each=chunksize),rep((nfolds-toremove+1):nfolds,each=chunksize-1))
        } else {
          folds = rep(rep(1:nfolds,each=chunksize),chunksofeach)[1:n]
        }
      }
    }
    if (BUTSIGNEDRANK){
      if (!is.null(indelseses)){
        step = one_step(Y,X,cov_inds,alpha,folds,nfolds,family,methods=c("oldHC"),current=c(),SKIP_NEIGHBOR=neighbor_skips[i],indelses=indelseses[[i]])
      } else {
        step = one_step(Y,X,cov_inds,alpha,folds,nfolds,family,methods=c("oldHC"),current=c(),SKIP_NEIGHBOR=neighbor_skips[i])
      }
    } else {
      if (!is.null(indelseses)){
        step = one_step(Y,X,cov_inds,alpha,folds,nfolds,family,methods=c("CV"),current=c(),SKIP_NEIGHBOR=neighbor_skips[i],indelses=indelseses[[i]])
      } else {
        step = one_step(Y,X,cov_inds,alpha,folds,nfolds,family,methods=c("CV"),current=c(),SKIP_NEIGHBOR=neighbor_skips[i])
      }
    }
    
    results[i,1] = step$candidate
    results[i,2] = step$current$pvals
    results[i,3] = step$current$Tval
    
    if (!BUTSIGNEDRANK){
      CVdiffs[i,1:nfolds] = step$CVdiffs1
      CVdiffs[i,(max(fold_nums)+1):(max(fold_nums)+nfolds)] = step$CVdiffsbest
    }
  }
  results = cbind(results,fold_nums,chunksizes,neighbor_skips)
  colnames(results) = c("candidate","pval","CV_score","num_folds","chunksize","skip_neighbor")
  if (!BUTSIGNEDRANK){
    results = cbind(results,CVdiffs)
  }
  return(results)
}

plot_results = function(resultmat,methods,alpha=0.05,ncov=3,nstep=NULL,covnames=c("Position","HD","Speed"),pvaldiststep=1,hline=FALSE,nonsig=FALSE,cex=1){
  
  nummethod = length(methods)
  if (is.null(nstep)){
    nstep = ncov
  }
  # First step, significant selection, for each method for each covariate
  # Same for second, third (multiple columns)
  
  # p-value distribution for first step in mismatched
  ncells = dim(resultmat)[3]
  mainframe = data.frame(array(0,dim=nummethod*nstep*ncells))
  colnames(mainframe) = c("Pvalue")
  mainframe$Method = rep(methods,each=nstep*ncells)
  mainframe$Step = array(0,nummethod*nstep*ncells)
  mainframe$Covariate = array(0,nummethod*nstep*ncells)
  mainframe$Significant = array(NA,nummethod*nstep*ncells)
  ind = 0
  for (m in 1:nummethod){
    for (s in 1:nstep){
      for (c in 1:ncells){
        ind = ind + 1
        mainframe$Pvalue[ind] = resultmat[1+m,s,c]
        mainframe$Step[ind] = s
        mainframe$Covariate[ind] = resultmat[1,s,c]
        if (!is.na(mainframe$Pvalue[ind])){
          if (nonsig){
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s-(nummethod-m)*0.1-0.05
            } else {
              mainframe$Significant[ind] = s+(m-1)*0.1+0.05
            }
          } else {
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s+(m-1)*0.1-0.1*(nummethod-1)/2
            }
          }
        } else {
          mainframe$Significant[ind] = NA
        }
      }
    }
  }
  
  if (length(which(mainframe$Pvalue <= alpha))==0){
    SKIPTOP = TRUE
  } else {
    SKIPTOP = FALSE
  }
  
  if (!SKIPTOP){
    breaks = c()
    for (m in 1:nummethod){
      for (s in 1:nstep){
        if (nonsig){
          breaks = c(breaks,s-(nummethod-m)*0.1-0.01,s-(nummethod-m)*0.1-0.09,s+(m-1)*0.1+0.01,s+(m-1)*0.1+0.09)
        } else {
          breaks = c(breaks,s+(m-1)*0.1-0.1*(nummethod-1)/2-0.045,s+(m-1)*0.1-0.1*(nummethod-1)/2+0.045)
        }
      }
    }
    breaks = sort(breaks)
    par(mfrow=c(2,1))
    if (nonsig){
      titl = "Number of sign. and non-sign. tests for each method in each step"
    } else {
      titl = "Number of significant tests for each method in each step"
    }
    hist(mainframe$Significant,xlab="Step", main=titl, col="firebrick",lty="blank",ylim=c(0,ncells),breaks=breaks,freq = T, xaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
    axis(side=1, at=0:4, labels=c("",1,2,3,""),cex.lab=cex,cex=cex, cex.axis=cex)
    hist(mainframe$Significant[which(mainframe$Covariate > 1)], col="royalblue1",lty="blank", add=TRUE,breaks=breaks,freq=T)
    hist(mainframe$Significant[which(mainframe$Covariate > 2)], col="darkcyan",lty="blank", add=TRUE,breaks=breaks,freq=T)
    legend("topright",legend=covnames,col=c("royalblue1","firebrick","darkcyan"),pch=15,cex=cex)
  }
  pbreaksbase = seq(0,1.0001,length.out = 21)
  pbreaks = c()
  methodlab = c("")
  for (m in 1:nummethod){
    pbreaks = c(pbreaks, pbreaksbase+1.1*(m-1))
    methodlab = c(methodlab,methods[m],"")
  }
  if (hline){
    yl = c(0,ncells/20*8)
    hist(mainframe$Pvalue[which(mainframe$Method == methods[1] & mainframe$Step == pvaldiststep)],main="Distribution of p-values in first step for each method",xlab="",breaks=pbreaks,freq=T,xaxt="n",ylim=yl,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
    axis(side=1, at=seq(0.5,0.5+1.1*(nummethod-1),by=1.1),tick=FALSE, labels=methods,cex.lab=cex,cex=cex, cex.axis=cex)
    for (m in 2:nummethod){
      hist(1.1*(m-1)+mainframe$Pvalue[which(mainframe$Method == methods[m] & mainframe$Step == pvaldiststep)],add=T,breaks=pbreaks,freq=T)
    }
    abline(h=ncells/20,lty=2,col="red")
  } else {
    hist(mainframe$Pvalue[which(mainframe$Method == methods[1] & mainframe$Step == pvaldiststep)],main="Distribution of p-values in first step for each method",xlab="",breaks=pbreaks,freq=T,xaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
    axis(side=1, at=seq(0.5,0.5+1.1*(nummethod-1),by=1.1),tick=FALSE, labels=methods,cex.lab=cex,cex=cex, cex.axis=cex)
    for (m in 2:nummethod){
      hist(1.1*(m-1)+mainframe$Pvalue[which(mainframe$Method == methods[m] & mainframe$Step == pvaldiststep)],add=T,breaks=pbreaks,freq=T)
    }
  }
}

print_progress = function(i,n,start,precision=1){
  pad=2*ceiling(log(n)/log(10))+1
  sofar = round((proc.time()[3] - start)*10^precision)/10^precision
  remain = round((sofar/i*(n-i))*10^precision)/10^precision
  output = paste(str_pad(paste(i,"/",n,sep=""),pad,"left"),", Elapsed: ",str_pad(formatC(sofar,format = "f",digits=precision),7+precision,"left"),", remaining: ",str_pad(formatC(remain,format = "f",digits=precision),7+precision,"left"),". ",sep="")
  percent = round(i/n*100)
  halfpercent = ceiling(percent/2)
  output=paste(output,paste("Progress: |",paste(rep("|",halfpercent),collapse=""),paste(rep(" ",50-halfpercent),collapse=""),"| ",percent,"%",sep=""),sep=" ")
  if (i < n){
    message(output,"\r",appendLF=FALSE)
    flush.console()
  } else {
    print(output)
  }
}

history_terms = function(y, num=5, numeach=3){
  n = length(y)
  hist = matrix(0, nrow=n, ncol=num)
  
  for (i in 1:num){
    for (j in 1:numeach){
      hist[((i - 1) * numeach + j + 1) : n, i] = hist[((i - 1) * numeach + j + 1) : n, i] + y[1 : (n - ((i - 1) * numeach + j))]
    }
  }
  return(hist)
}

#### Barplots ----
plot_combined_results_better_bars = function(resultmat,resultmat2,methods,methods2,newlabs,alpha=0.05,ncov=3,nstep=NULL,covnames=c("HD","Pos","Spe"),pvaldiststep=1,hline=FALSE,nonsig=FALSE,cex=1,ylim1=NULL,ylim2=NULL,logscale=FALSE,line=-1,ymax=NULL){
  if (hline){
    par(mfrow=c(2,1))
  } else {
    par(mfrow=c(1,1))
  }
  parlwd = par("lwd")
  par(lwd=cex)
  nummethod = length(methods)
  nummethod2 = length(methods2)-1
  totmethods = c(methods[c(1,5)],methods2[2:(nummethod2+1)], methods[2:4])#c(methods,methods2[2:(nummethod2+1)])
  
  totnummeth = nummethod + nummethod2
  if (is.null(nstep)){
    nstep = ncov
  }
  ncells = dim(resultmat)[3]
  mainframe = data.frame(array(0,dim=totnummeth*nstep*ncells))
  colnames(mainframe) = c("Pvalue")
  mainframe$Method = rep(totmethods,each=nstep*ncells)
  mainframe$Step = array(0,totnummeth*nstep*ncells)
  mainframe$Covariate = array(0,totnummeth*nstep*ncells)
  mainframe$Significant = array(NA,totnummeth*nstep*ncells)
  ind = 0
  m_inds = c(2,6,3,3,4,5)
  for (m in 1:totnummeth){
    mj = m_inds[m]
    if (m == 3){
      resmat = resultmat2
    } else {
      resmat = resultmat
    }
    for (s in 1:nstep){
      for (c in 1:ncells){
        ind = ind + 1
        mainframe$Pvalue[ind] = resmat[mj,s,c]
        mainframe$Step[ind] = s
        mainframe$Covariate[ind] = resmat[1,s,c]
        if (!is.na(mainframe$Pvalue[ind])){
          if (nonsig){
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s-(totnummeth-m)*0.1-0.05
            } else {
              mainframe$Significant[ind] = s+(m-1)*0.1+0.05
            }
          } else {
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s+(m-1)*0.1-0.1*(totnummeth-1)/2
            }
          }
        } else {
          mainframe$Significant[ind] = NA
        }
      }
    }
  }
  if (length(which(mainframe$Pvalue <= alpha))==0){
    SKIPTOP = TRUE
  } else {
    SKIPTOP = FALSE
  }
  
  ###
  
  nummethod1 = length(methods)
  nummethod2 = length(methods2)-1
  combmethod = c(methods[c(1,5)],methods2[2:(nummethod2+1)], methods[2:4])#c(methods,methods2[2:length(methods2)])
  tunings = list()
  
  for (i in 1:nummethod1){
    tunings[[methods[i]]] = c()
    
    for (j in 1:dim(resultmat)[3]){
      if (resultmat[i+1,1,j] <= alpha){
        tuning = covnames[resultmat[1,1,j]]
        if (resultmat[i+1,2,j] <= alpha){
          tuning = paste(tuning,covnames[resultmat[1,2,j]],sep="")
          if (resultmat[i+1,3,j] <= alpha){
            tuning = paste(tuning,covnames[resultmat[1,3,j]],sep="")
          }
        }
      } else {
        tuning = "None"
      }
      if (tuning %in% c("HDPos","PosHD")){tuning = "PosHD"}
      if (tuning %in% c("SpePos","PosSpe")){tuning = "PosSpe"}
      if (tuning %in% c("HDSpe","SpeHD")){tuning = "HDSpe"}
      if (tuning %in% c("HDPosSpe","PosHDSpe","PosSpeHD","SpePosHD","HDSpePos","SpeHDPos")){tuning = "All"}
      if (tuning %in% c("1D-A 2D ","2D 1D-A ")){tuning = "2D 1D-A "}
      if (tuning %in% c("1D-B2D ","2D 1D-B")){tuning = "2D 1D-B"}
      if (tuning %in% c("1D-A 1D-B","1D-B1D-A ")){tuning = "1D-A 1D-B"}
      if (tuning %in% c("1D-A 2D 1D-B","1D-A 1D-B2D ","2D 1D-A 1D-B","2D 1D-B1D-A ","1D-B1D-A 2D ","1D-B2D 1D-A ")){tuning = "All"}
      if (tuning %in% c("B (1D)A (2D)","A (2D)B (1D)")){tuning = "AB"}
      if (tuning %in% c("C (1D)A (2D)","A (2D)C (1D)")){tuning = "AC"}
      if (tuning %in% c("C (1D)B (1D)","B (1D)C (1D)")){tuning = "BC"}
      if (tuning %in% c("B (1D)A (2D)C (1D)","C (1D)B (1D)A (2D)","C (1D)A (2D)B (1D)","B (1D)C (1D)A (2D)","A (2D)B (1D)C (1D)","A (2D)C (1D)B (1D)")){tuning = "All"}
      tunings[[methods[i]]] = c(tunings[[methods[i]]],tuning)
    }
  }
  for (i in 1:nummethod2){
    tunings[[methods2[i+1]]] = c()
    
    for (j in 1:dim(resultmat2)[3]){
      if (resultmat2[i+2,1,j] <= alpha){
        tuning = covnames[resultmat2[1,1,j]]
        if (resultmat2[i+2,2,j] <= alpha){
          tuning = paste(tuning,covnames[resultmat2[1,2,j]],sep="")
          if (resultmat2[i+2,3,j] <= alpha){
            tuning = paste(tuning,covnames[resultmat2[1,3,j]],sep="")
          }
        }
      } else {
        tuning = "None"
      }
      if (tuning %in% c("HDPos","PosHD")){tuning = "PosHD"}
      if (tuning %in% c("SpePos","PosSpe")){tuning = "PosSpe"}
      if (tuning %in% c("HDSpe","SpeHD")){tuning = "HDSpe"}
      if (tuning %in% c("HDPosSpe","PosHDSpe","PosSpeHD","SpePosHD","HDSpePos","SpeHDPos")){tuning = "All"}
      if (tuning %in% c("1D-A 2D ","2D 1D-A ")){tuning = "2D 1D-A "}
      if (tuning %in% c("1D-B2D ","2D 1D-B")){tuning = "2D 1D-B"}
      if (tuning %in% c("1D-A 1D-B","1D-B1D-A ")){tuning = "1D-A 1D-B"}
      if (tuning %in% c("1D-A 2D 1D-B","1D-A 1D-B2D ","2D 1D-A 1D-B","2D 1D-B1D-A ","1D-B1D-A 2D ","1D-B2D 1D-A ")){tuning = "All"}
      if (tuning %in% c("B (1D)A (2D)","A (2D)B (1D)")){tuning = "AB"}
      if (tuning %in% c("C (1D)A (2D)","A (2D)C (1D)")){tuning = "AC"}
      if (tuning %in% c("C (1D)B (1D)","B (1D)C (1D)")){tuning = "BC"}
      if (tuning %in% c("B (1D)A (2D)C (1D)","C (1D)B (1D)A (2D)","C (1D)A (2D)B (1D)","B (1D)C (1D)A (2D)","A (2D)B (1D)C (1D)","A (2D)C (1D)B (1D)")){tuning = "All"}
      tunings[[methods2[i+1]]] = c(tunings[[methods2[i+1]]],tuning)
    }
  }
  
  otherframe = data.frame(array(0,dim=totnummeth*ncells))
  otherframe$method = rep(combmethod,each=ncells)
  otherframe$tuning = rep("",totnummeth*ncells)
  otherframe$val = rep(0,totnummeth*ncells)
  if (covnames[1] == "HD"){
    labels = c("None","Pos","HD","Spe","PosHD","PosSpe","HDSpe","All")
  } else if (covnames[1] == "1D-A ") {
    labels = c("None","2D ","1D-A ","1D-B","2D 1D-A ","2D 1D-B","1D-A 1D-B","All")
  } else {
    labels = c("None","A (2D)","B (1D)","C (1D)","AB","AC","BC","All")
  }
  
  for (j in 1:length(combmethod)){
    m = combmethod[j]
    tuning = tunings[[m]]
    
    otherframe$tuning[((j-1)*ncells+1):(j*ncells)] = tuning
    
    for (k in 1:ncells){
      otherframe$val[(j-1)*ncells+k] = which(labels == otherframe$tuning[(j-1)*ncells+k])/9 + (j-1)*1.1
    }
  }
  
  xl = c(0,1.1*(totnummeth)+0.6)
  if (!SKIPTOP){
    cols = c("#999999","#56B4E9","#D55E00","#F0E442","#CC79A7","#009E73","#E69F00","#000000")
    breaksbase = seq(0,1.0001,length.out = 18)
    for (k in seq(1,18,by=2)){
      breaksbase[k] = breaksbase[k] + 1/25
    }
    breaks = c()
    methodlab = c("")
    for (m in 1:totnummeth){
      breaks = c(breaks, breaksbase+1.1*(m-1))
      methodlab = c(methodlab,totmethods[m],"")
    }
    
    if (logscale){
      header = paste("Distribution of classification for each method, ",ncells," cells in total (log-scale)",sep="")
    } else {
      header = paste("Distribution of classification for each method, ",ncells," cells in total",sep="")
    }
    if (logscale){
      ylim1 = c(0,log(1+ncells)/log(10)*1.1)
      
      h1 = hist(otherframe$val[which(otherframe$tuning == labels[1])],breaks=breaks,plot=F)
      counts = h1$counts
      h1$counts = log(1+h1$counts)/log(10)
      h1$density = h1$counts
      plot(h1,ylim=ylim1,main=header,xlim=xl,xlab="",ylab="",col=cols[1],lty=1, add=FALSE,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
      title(ylab="Frequency",cex=cex,cex.lab=cex,line=1)
      q = length(counts)
      for (j in 1:q){
        if (counts[j] > 0){
          x = h1$mids[j]
          y = log(1+counts[j])/log(10)
          text(labels = paste(counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
        }
      }
      
      for (k in 2:8){
        h2 = hist(otherframe$val[which(otherframe$tuning == labels[k])],breaks=breaks,plot=F)
        counts = h2$counts
        h2$counts = log(1+h2$counts)/log(10)
        h2$density = h2$counts
        plot(h2,ylim=ylim1,col=cols[k],lty=1, add=TRUE,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        
        q = length(counts)
        for (j in 1:q){
          if (counts[j] > 0){
            x = h2$mids[j]
            y = log(1+counts[j])/log(10)
            text(labels = paste(counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
          }
        }
      }
      ylabs = c(0,1,2,4,8,16,32,64,128,256)
      axis(side=2, at=log(ylabs+1)/log(10),tick=TRUE, labels=ylabs,lwd=cex,cex.lab=cex,cex=cex, cex.axis=cex,pos=line)
    } else {
      if (is.null(ymax)){
        ylim1 = c(0,ncells*1.1)
      } else {
        ylim1 = c(0,ymax*1.1)
      }
      h1=hist(otherframe$val[which(otherframe$tuning == labels[1])],xlim=xl,ylim=ylim1,main=header,ylab="",xlab="",col=cols[1],lty=1, add=FALSE,breaks=breaks,freq=T,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
      q = length(h1$counts)
      for (j in 1:q){
        if (h1$counts[j] > 0){
          x = h1$mids[j]
          y = h1$counts[j]
          text(labels = paste(h1$counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
        }
      }
      for (k in 2:8){
        h2=hist(otherframe$val[which(otherframe$tuning == labels[k])],breaks=breaks,col=cols[k],lty=1, add=TRUE,freq=T,xaxt="n",xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        q = length(h2$counts)
        for (j in 1:q){
          if (h2$counts[j] > 0){
            x = h2$mids[j]
            y = h2$counts[j]
            text(labels = paste(h2$counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
          }
        }
      }
      title(ylab="Frequency",cex=cex,cex.lab=cex,line=1)
      ylabs = seq(0,50*round(ncells/50),by=50)
      axis(side=2, at=ylabs,tick=TRUE, labels=ylabs,cex.lab=cex,cex=cex,lwd=cex, cex.axis=cex,pos=line)
    }
    
    ylimthing = par("yaxp")[1:2]
    pos = ylimthing[1] - 0.1*(ylimthing[2]-ylimthing[1])
    axis(side=1, at=seq(0.5,0.5+1.1*(totnummeth-1),by=1.1),tick=FALSE,lwd=cex,pos = pos, labels=newlabs,cex.lab=cex,cex=cex, cex.axis=cex)
    ticks1 = seq(breaksbase[1],breaksbase[1]+1.1*(totnummeth-1),by=1.1)
    ticks2 = seq(1.0001,1.0001+1.1*(totnummeth-1),by=1.1)
    axis(side=1, at=sort(c(ticks1,ticks2)),tick=T, labels=rep("",length(ticks1)+length(ticks2)),cex.lab=cex,cex=cex, cex.axis=cex,pos=0,lwd=cex)
    par(xpd=T)
    #legend("topright",legend=labels,cex=cex,fill=cols)
    legend(totnummeth + 0.7, ylim1[2]*1.1,legend=labels,cex=cex,fill=cols)
    par(xpd=F)
  }
  pbreaksbase = c(-0.000001,seq(0.050001,1.000001,length.out = 20))
  pbreaks = c()
  methodlab = c("")
  for (m in 1:totnummeth){
    pbreaks = c(pbreaks, pbreaksbase+1.1*(m-1))
    methodlab = c(methodlab,totmethods[m],"")
  }
  if (hline){
    if (is.null(ylim2)){
      yl = c(0,3)
    } else {
      yl = ylim2
    }
    hist(mainframe$Pvalue[which(mainframe$Method == totmethods[1] & mainframe$Step == pvaldiststep)],main=paste("Distribution of ",expression(p),"-values in first step for each method",sep=""),xlab="",ylab="",breaks=pbreaks,freq=F,lwd=cex,xaxt="n",yaxt="n",xlim=xl,ylim=yl,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
    title(ylab="Density",cex=cex,cex.lab=cex,line=1)
    ylimthing = par("yaxp")[1:2]
    pos = ylimthing[1] - 0.1*(ylimthing[2]-ylimthing[1])
    axis(side=1, at=seq(0.5,0.5+1.1*(totnummeth-1),by=1.1),pos=pos,tick=FALSE, labels=newlabs,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    ticks1 = seq(0,1.1*(totnummeth-1),by=1.1)
    ticks2 = seq(1.0001,1.0001+1.1*(totnummeth-1),by=1.1)
    axis(side=1, at=sort(c(ticks1,ticks2)),tick=T, labels=rep("",length(ticks1)+length(ticks2)),cex.lab=cex,cex=cex, cex.axis=cex,pos=0,lwd=cex)
    #ylabs = seq(0,20*round(yl[2]/20),by=20)
    ylabs = seq(0,3,by=0.5)
    axis(side=2, at=ylabs,tick=TRUE, labels=ylabs,lwd=cex,pos=line,cex.lab=cex,cex=cex, cex.axis=cex)
    for (m in 2:totnummeth){
      hist(1.1*(m-1)+mainframe$Pvalue[which(mainframe$Method == totmethods[m] & mainframe$Step == pvaldiststep)],add=T,breaks=pbreaks,freq=F,lwd=cex)
    }
    lines(c(0,1.1*(totnummeth-1)+1.0001),c(1,1),lty=1,col="red",lwd=cex*2)
    legend("topright",legend=c(expression(paste(alpha," = 0.05",sep=""))),lty=1,col="red",lwd=cex*2,cex=cex)
  } else {
    if (FALSE){
      hist(mainframe$Pvalue[which(mainframe$Method == totmethods[1] & mainframe$Step == pvaldiststep)],main=paste("Distribution of ",expression(p),"-values in first step for each method",sep=""),xlab="",xlim=xl,breaks=pbreaks,freq=T,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
      ylimthing = par("yaxp")[1:2]
      pos = ylimthing[1] - 0.1*(ylimthing[2]-ylimthing[1])
      axis(side=1, at=seq(0.5,0.5+1.1*(totnummeth-1),by=1.1),pos=pos,tick=FALSE, labels=newlabs,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
      
      ticks1 = seq(0,1.1*(totnummeth-1),by=1.1)
      ticks2 = seq(1.0001,1.0001+1.1*(totnummeth-1),by=1.1)
      axis(side=1, at=sort(c(ticks1,ticks2)),tick=T, labels=rep("",length(ticks1)+length(ticks2)),cex.lab=cex,cex=cex, cex.axis=cex,pos=0,lwd=cex)
      ylabs = seq(0,50*round(ncells/50),by=50)
      axis(side=2, at=ylabs,tick=TRUE, labels=ylabs,line=line,cex.lab=cex,cex=cex, cex.axis=cex)
      for (m in 2:totnummeth){
        hist(1.1*(m-1)+mainframe$Pvalue[which(mainframe$Method == totmethods[m] & mainframe$Step == pvaldiststep)],add=T,breaks=pbreaks,freq=T)
      }
    }
  }
  par(lwd=parlwd)
}

OLDplot_combined_results_better_bars = function(resultmat,resultmat2,methods,methods2,newlabs,alpha=0.05,ncov=3,nstep=NULL,covnames=c("HD","Pos","Spe"),pvaldiststep=1,hline=FALSE,nonsig=FALSE,cex=1,ylim1=NULL,ylim2=NULL,logscale=FALSE,line=-1){
  if (hline){
    par(mfrow=c(2,1))
  } else {
    par(mfrow=c(1,1))
  }
  parlwd = par("lwd")
  par(lwd=cex)
  nummethod = length(methods)
  nummethod2 = length(methods2)-1
  totmethods = c(methods,methods2[2:(nummethod2+1)])
  
  totnummeth = nummethod + nummethod2
  if (is.null(nstep)){
    nstep = ncov
  }
  ncells = dim(resultmat)[3]
  mainframe = data.frame(array(0,dim=totnummeth*nstep*ncells))
  colnames(mainframe) = c("Pvalue")
  mainframe$Method = rep(totmethods,each=nstep*ncells)
  mainframe$Step = array(0,totnummeth*nstep*ncells)
  mainframe$Covariate = array(0,totnummeth*nstep*ncells)
  mainframe$Significant = array(NA,totnummeth*nstep*ncells)
  ind = 0
  for (m in 1:nummethod){
    for (s in 1:nstep){
      for (c in 1:ncells){
        ind = ind + 1
        mainframe$Pvalue[ind] = resultmat[1+m,s,c]
        mainframe$Step[ind] = s
        mainframe$Covariate[ind] = resultmat[1,s,c]
        if (!is.na(mainframe$Pvalue[ind])){
          if (nonsig){
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s-(totnummeth-m)*0.1-0.05
            } else {
              mainframe$Significant[ind] = s+(m-1)*0.1+0.05
            }
          } else {
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s+(m-1)*0.1-0.1*(totnummeth-1)/2
            }
          }
        } else {
          mainframe$Significant[ind] = NA
        }
      }
    }
  }
  
  for (m in (nummethod+1:nummethod2)){
    for (s in 1:nstep){
      for (c in 1:ncells){
        ind = ind + 1
        mainframe$Pvalue[ind] = resultmat2[2+m-nummethod,s,c]
        mainframe$Step[ind] = s
        mainframe$Covariate[ind] = resultmat2[1,s,c]
        if (!is.na(mainframe$Pvalue[ind])){
          if (nonsig){
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s-(totnummeth-m)*0.1-0.05
            } else {
              mainframe$Significant[ind] = s+(m-1)*0.1+0.05
            }
          } else {
            if (mainframe$Pvalue[ind] <= alpha){
              mainframe$Significant[ind] = s+(m-1)*0.1-0.1*(totnummeth-1)/2
            }
          }
        } else {
          mainframe$Significant[ind] = NA
        }
      }
    }
  }
  if (length(which(mainframe$Pvalue <= alpha))==0){
    SKIPTOP = TRUE
  } else {
    SKIPTOP = FALSE
  }
  
  ###
  
  nummethod1 = length(methods)
  nummethod2 = length(methods2)-1
  combmethod = c(methods,methods2[2:length(methods2)])
  tunings = list()
  
  for (i in 1:nummethod1){
    tunings[[methods[i]]] = c()
    
    for (j in 1:dim(resultmat)[3]){
      if (resultmat[i+1,1,j] <= alpha){
        tuning = covnames[resultmat[1,1,j]]
        if (resultmat[i+1,2,j] <= alpha){
          tuning = paste(tuning,covnames[resultmat[1,2,j]],sep="")
          if (resultmat[i+1,3,j] <= alpha){
            tuning = paste(tuning,covnames[resultmat[1,3,j]],sep="")
          }
        }
      } else {
        tuning = "None"
      }
      if (tuning %in% c("HDPos","PosHD")){tuning = "PosHD"}
      if (tuning %in% c("SpePos","PosSpe")){tuning = "PosSpe"}
      if (tuning %in% c("HDSpe","SpeHD")){tuning = "HDSpe"}
      if (tuning %in% c("HDPosSpe","PosHDSpe","PosSpeHD","SpePosHD","HDSpePos","SpeHDPos")){tuning = "All"}
      if (tuning %in% c("1D-A 2D ","2D 1D-A ")){tuning = "2D 1D-A "}
      if (tuning %in% c("1D-B2D ","2D 1D-B")){tuning = "2D 1D-B"}
      if (tuning %in% c("1D-A 1D-B","1D-B1D-A ")){tuning = "1D-A 1D-B"}
      if (tuning %in% c("1D-A 2D 1D-B","1D-A 1D-B2D ","2D 1D-A 1D-B","2D 1D-B1D-A ","1D-B1D-A 2D ","1D-B2D 1D-A ")){tuning = "All"}
      if (tuning %in% c("B (1D)A (2D)","A (2D)B (1D)")){tuning = "AB"}
      if (tuning %in% c("C (1D)A (2D)","A (2D)C (1D)")){tuning = "AC"}
      if (tuning %in% c("C (1D)B (1D)","B (1D)C (1D)")){tuning = "BC"}
      if (tuning %in% c("B (1D)A (2D)C (1D)","C (1D)B (1D)A (2D)","C (1D)A (2D)B (1D)","B (1D)C (1D)A (2D)","A (2D)B (1D)C (1D)","A (2D)C (1D)B (1D)")){tuning = "All"}
      tunings[[methods[i]]] = c(tunings[[methods[i]]],tuning)
    }
  }
  for (i in 1:nummethod2){
    tunings[[methods2[i+1]]] = c()
    
    for (j in 1:dim(resultmat2)[3]){
      if (resultmat2[i+2,1,j] <= alpha){
        tuning = covnames[resultmat2[1,1,j]]
        if (resultmat2[i+2,2,j] <= alpha){
          tuning = paste(tuning,covnames[resultmat2[1,2,j]],sep="")
          if (resultmat2[i+2,3,j] <= alpha){
            tuning = paste(tuning,covnames[resultmat2[1,3,j]],sep="")
          }
        }
      } else {
        tuning = "None"
      }
      if (tuning %in% c("HDPos","PosHD")){tuning = "PosHD"}
      if (tuning %in% c("SpePos","PosSpe")){tuning = "PosSpe"}
      if (tuning %in% c("HDSpe","SpeHD")){tuning = "HDSpe"}
      if (tuning %in% c("HDPosSpe","PosHDSpe","PosSpeHD","SpePosHD","HDSpePos","SpeHDPos")){tuning = "All"}
      if (tuning %in% c("1D-A 2D ","2D 1D-A ")){tuning = "2D 1D-A "}
      if (tuning %in% c("1D-B2D ","2D 1D-B")){tuning = "2D 1D-B"}
      if (tuning %in% c("1D-A 1D-B","1D-B1D-A ")){tuning = "1D-A 1D-B"}
      if (tuning %in% c("1D-A 2D 1D-B","1D-A 1D-B2D ","2D 1D-A 1D-B","2D 1D-B1D-A ","1D-B1D-A 2D ","1D-B2D 1D-A ")){tuning = "All"}
      if (tuning %in% c("B (1D)A (2D)","A (2D)B (1D)")){tuning = "AB"}
      if (tuning %in% c("C (1D)A (2D)","A (2D)C (1D)")){tuning = "AC"}
      if (tuning %in% c("C (1D)B (1D)","B (1D)C (1D)")){tuning = "BC"}
      if (tuning %in% c("B (1D)A (2D)C (1D)","C (1D)B (1D)A (2D)","C (1D)A (2D)B (1D)","B (1D)C (1D)A (2D)","A (2D)B (1D)C (1D)","A (2D)C (1D)B (1D)")){tuning = "All"}
      tunings[[methods2[i+1]]] = c(tunings[[methods2[i+1]]],tuning)
    }
  }
  
  otherframe = data.frame(array(0,dim=totnummeth*ncells))
  otherframe$method = rep(combmethod,each=ncells)
  otherframe$tuning = rep("",totnummeth*ncells)
  otherframe$val = rep(0,totnummeth*ncells)
  if (covnames[1] == "HD"){
    labels = c("None","Pos","HD","Spe","PosHD","PosSpe","HDSpe","All")
  } else if (covnames[1] == "1D-A ") {
    labels = c("None","2D ","1D-A ","1D-B","2D 1D-A ","2D 1D-B","1D-A 1D-B","All")
  } else {
    labels = c("None","A (2D)","B (1D)","C (1D)","AB","AC","BC","All")
  }
  for (j in 1:length(combmethod)){
    m = combmethod[j]
    tuning = tunings[[m]]
    
    otherframe$tuning[((j-1)*ncells+1):(j*ncells)] = tuning
    
    for (k in 1:ncells){
      otherframe$val[(j-1)*ncells+k] = which(labels == otherframe$tuning[(j-1)*ncells+k])/9 + (j-1)*1.1
    }
  }
  
  xl = c(0,1.1*(totnummeth)+0.6)
  if (!SKIPTOP){
    cols = c("#999999","#56B4E9","#D55E00","#F0E442","#CC79A7","#009E73","#E69F00","#000000")
    breaksbase = seq(0,1.0001,length.out = 18)
    for (k in seq(1,18,by=2)){
      breaksbase[k] = breaksbase[k] + 1/25
    }
    breaks = c()
    methodlab = c("")
    for (m in 1:totnummeth){
      breaks = c(breaks, breaksbase+1.1*(m-1))
      methodlab = c(methodlab,totmethods[m],"")
    }
    
    if (logscale){
      header = paste("Distribution of classification for each method, ",ncells," cells in total (log-scale)",sep="")
    } else {
      header = paste("Distribution of classification for each method, ",ncells," cells in total",sep="")
    }
    if (logscale){
      ylim1 = c(0,log(1+ncells)/log(10)*1.1)
      
      h1 = hist(otherframe$val[which(otherframe$tuning == labels[1])],breaks=breaks,plot=F)
      counts = h1$counts
      h1$counts = log(1+h1$counts)/log(10)
      h1$density = h1$counts
      plot(h1,ylim=ylim1,main=header,xlim=xl,xlab="",ylab="",col=cols[1],lty=1, add=FALSE,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
      title(ylab="Frequency",cex=cex,cex.lab=cex,line=1)
      q = length(counts)
      for (j in 1:q){
        if (counts[j] > 0){
          x = h1$mids[j]
          y = log(1+counts[j])/log(10)
          text(labels = paste(counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
        }
      }
      
      for (k in 2:8){
        h2 = hist(otherframe$val[which(otherframe$tuning == labels[k])],breaks=breaks,plot=F)
        counts = h2$counts
        h2$counts = log(1+h2$counts)/log(10)
        h2$density = h2$counts
        plot(h2,ylim=ylim1,col=cols[k],lty=1, add=TRUE,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        
        q = length(counts)
        for (j in 1:q){
          if (counts[j] > 0){
            x = h2$mids[j]
            y = log(1+counts[j])/log(10)
            text(labels = paste(counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
          }
        }
      }
      ylabs = c(0,1,2,4,8,16,32,64,128,256)
      axis(side=2, at=log(ylabs+1)/log(10),tick=TRUE, labels=ylabs,lwd=cex,cex.lab=cex,cex=cex, cex.axis=cex,pos=line)
    } else {
      ylim1 = c(0,ncells*1.1)
      h1=hist(otherframe$val[which(otherframe$tuning == labels[1])],xlim=xl,ylim=ylim1,main=header,ylab="",xlab="",col=cols[1],lty=1, add=FALSE,breaks=breaks,freq=T,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
      q = length(h1$counts)
      for (j in 1:q){
        if (h1$counts[j] > 0){
          x = h1$mids[j]
          y = h1$counts[j]
          text(labels = paste(h1$counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
        }
      }
      for (k in 2:8){
        h2=hist(otherframe$val[which(otherframe$tuning == labels[k])],breaks=breaks,col=cols[k],lty=1, add=TRUE,freq=T,xaxt="n",xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        q = length(h2$counts)
        for (j in 1:q){
          if (h2$counts[j] > 0){
            x = h2$mids[j]
            y = h2$counts[j]
            text(labels = paste(h2$counts[j]),x = x,y = y,pos=3,cex=cex*0.6)
          }
        }
      }
      title(ylab="Frequency",cex=cex,cex.lab=cex,line=1)
      ylabs = seq(0,50*round(ncells/50),by=50)
      axis(side=2, at=ylabs,tick=TRUE, labels=ylabs,cex.lab=cex,cex=cex,lwd=cex, cex.axis=cex,pos=line)
    }
    
    ylimthing = par("yaxp")[1:2]
    pos = ylimthing[1] - 0.1*(ylimthing[2]-ylimthing[1])
    axis(side=1, at=seq(0.5,0.5+1.1*(totnummeth-1),by=1.1),tick=FALSE,lwd=cex,pos = pos, labels=newlabs,cex.lab=cex,cex=cex, cex.axis=cex)
    ticks1 = seq(breaksbase[1],breaksbase[1]+1.1*(totnummeth-1),by=1.1)
    ticks2 = seq(1.0001,1.0001+1.1*(totnummeth-1),by=1.1)
    axis(side=1, at=sort(c(ticks1,ticks2)),tick=T, labels=rep("",length(ticks1)+length(ticks2)),cex.lab=cex,cex=cex, cex.axis=cex,pos=0,lwd=cex)
    par(xpd=T)
    legend("topright",legend=labels,cex=cex,fill=cols)
    par(xpd=F)
  }
  pbreaksbase = seq(0,1.0001,length.out = 21)
  pbreaks = c()
  methodlab = c("")
  for (m in 1:totnummeth){
    pbreaks = c(pbreaks, pbreaksbase+1.1*(m-1))
    methodlab = c(methodlab,totmethods[m],"")
  }
  if (hline){
    if (is.null(ylim2)){
      yl = c(0,3)
    } else {
      yl = ylim2
    }
    hist(mainframe$Pvalue[which(mainframe$Method == totmethods[1] & mainframe$Step == pvaldiststep)],main=paste("Distribution of ",expression(p),"-values in first step for each method",sep=""),xlab="",ylab="",breaks=pbreaks,freq=F,lwd=cex,xaxt="n",yaxt="n",xlim=xl,ylim=yl,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
    title(ylab="Density",cex=cex,cex.lab=cex,line=1)
    ylimthing = par("yaxp")[1:2]
    pos = ylimthing[1] - 0.1*(ylimthing[2]-ylimthing[1])
    axis(side=1, at=seq(0.5,0.5+1.1*(totnummeth-1),by=1.1),pos=pos,tick=FALSE, labels=newlabs,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
    ticks1 = seq(0,1.1*(totnummeth-1),by=1.1)
    ticks2 = seq(1.0001,1.0001+1.1*(totnummeth-1),by=1.1)
    axis(side=1, at=sort(c(ticks1,ticks2)),tick=T, labels=rep("",length(ticks1)+length(ticks2)),cex.lab=cex,cex=cex, cex.axis=cex,pos=0,lwd=cex)
    #ylabs = seq(0,20*round(yl[2]/20),by=20)
    ylabs = seq(0,3,by=0.5)
    axis(side=2, at=ylabs,tick=TRUE, labels=ylabs,lwd=cex,pos=line,cex.lab=cex,cex=cex, cex.axis=cex)
    for (m in 2:totnummeth){
      hist(1.1*(m-1)+mainframe$Pvalue[which(mainframe$Method == totmethods[m] & mainframe$Step == pvaldiststep)],add=T,breaks=pbreaks,freq=F,lwd=cex)
    }
    lines(c(0,1.1*(totnummeth-1)+1.0001),c(1,1),lty=1,col="red",lwd=cex*2)
    legend("topright",legend=c(expression(paste(alpha," = 0.05",sep=""))),lty=1,col="red",lwd=cex*2,cex=cex)
  } else {
    if (FALSE){
      hist(mainframe$Pvalue[which(mainframe$Method == totmethods[1] & mainframe$Step == pvaldiststep)],main=paste("Distribution of ",expression(p),"-values in first step for each method",sep=""),xlab="",xlim=xl,breaks=pbreaks,freq=T,xaxt="n",yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
      ylimthing = par("yaxp")[1:2]
      pos = ylimthing[1] - 0.1*(ylimthing[2]-ylimthing[1])
      axis(side=1, at=seq(0.5,0.5+1.1*(totnummeth-1),by=1.1),pos=pos,tick=FALSE, labels=newlabs,cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
      
      ticks1 = seq(0,1.1*(totnummeth-1),by=1.1)
      ticks2 = seq(1.0001,1.0001+1.1*(totnummeth-1),by=1.1)
      axis(side=1, at=sort(c(ticks1,ticks2)),tick=T, labels=rep("",length(ticks1)+length(ticks2)),cex.lab=cex,cex=cex, cex.axis=cex,pos=0,lwd=cex)
      ylabs = seq(0,50*round(ncells/50),by=50)
      axis(side=2, at=ylabs,tick=TRUE, labels=ylabs,line=line,cex.lab=cex,cex=cex, cex.axis=cex)
      for (m in 2:totnummeth){
        hist(1.1*(m-1)+mainframe$Pvalue[which(mainframe$Method == totmethods[m] & mainframe$Step == pvaldiststep)],add=T,breaks=pbreaks,freq=T)
      }
    }
  }
  par(lwd=parlwd)
}

#### Construct simulated data ----
make_simulated_data = function(N,sd,b_proportion=0,numys=1,ncov=3,scale=0.25){
  #h is hidden, a, b, c are observed
  bx = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  by = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  h = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  a = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  c = make_1D_cov(N,sd,dist="uniform",smoother="exponential",uniscale=2.5)
  
  
  baserate = 0.03
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
  
  maxcomb = max(p_combined)
  if (maxcomb > 1){
    p_combined = p_combined/max(p_combined)
  }
  
  p = p_combined
  
  if (numys > 1){
    ys = matrix(0,nrow=N,ncol=numys)
    for (i in 1:numys){
      ys[,i] = rbinom(N,1,p)
    }
    return(ys)
  }
  y = rbinom(N,1,p)
  dat = data.frame(cbind(y,a,bx,by,c))
  
  if (ncov > 3){
    for (i in 4:ncov){
      d = make_1D_cov(N,sd,dist="uniform",smoother="exponential")
      dat[[paste("d",i,sep="")]] = d
    }
  }
  
  return(dat)
}

contribution = function(norm, decay, x, center_x, y=NULL, center_y=NULL, periodic=F){
  if (is.null(y)){
    diffs = abs(x - center_x)
    if (periodic){
      diffs = sapply(diffs, periodic_dist)
    }
    return (norm * decay_func(diffs^2, decay))
  }
  return(norm * decay_func(((x - center_x)^2 + (y - center_y)^2), decay))
}

normalize_speed = function(speed){
  #q = quantile(speed, 0.99)
  #speed[which(speed > q)] = q
  speed = speed - mean(speed)
  maxabsq = max(c(abs(quantile(speed, 0.05)), quantile(speed, 0.95)))
  speed[which(speed < -maxabsq)] = -maxabsq
  speed[which(speed > maxabsq)] = maxabsq
  speed = speed / maxabsq
  return(speed)
}

decay_func = function(square_diff, decay, EXP=F){
  if (EXP){
    return(exp(-decay * square_diff))
  } else {
    return(1 / (1 + decay * square_diff))
  }
}

periodic_dist = function(d){
  return(min(d, 2 * pi - d))
}

BACKUPmake_simulated_data_from_real_covariates = function(num_cells=300,HIDDEN_ONLY=F){
  
  PLOTTHINGS1 = T
  PLOTTHINGS2 = F
  if (PLOTTHINGS2){
    X = readRDS(paste(processed_data_folder,"Simulated_from_real_cov_X.RDS",sep=""))
    covinds_sim = readRDS(paste(processed_data_folder,"covinds_sim_new.RDS",sep=""))
    folds = make_folds(N, 20, 150)
  }
  
  #HD from second track data is hidden
  track_obs = readRDS(paste(processed_data_folder,"trackdata2.RDS",sep=""))
  track_hid = readRDS(paste(processed_data_folder,"trackdata2_other_date.RDS",sep=""))
  
  N = 12000
  hd_hid = track_hid$HD[1:N]
  x_hid = track_hid$headX[1:N]
  y_hid = track_hid$headY[1:N]
  s_hid = normalize_speed(track_hid$bodyspeed[1:N])
  
  intercept = 4.6
  
  norm_pos = 1.5
  norm_hd = 1
  norm_spe = 1.5
  hidden_boost = 0
  
  decay_pos = 0.003
  decay_hd = 3
  decay_spe = 1
  
  # Draw random centers for the tuning peaks
  centers_hd = runif(num_cells, min=-pi, max=pi)
  centers_x = runif(num_cells, min=-35, max=35)
  centers_y = runif(num_cells, min=-35, max=35)
  centers_s = runif(num_cells, min=-0.75, max=0.75)
  
  if (HIDDEN_ONLY){
    Y = matrix(0,nrow=N,ncol=num_cells)
    for (i in 1:num_cells){
      hidden_hd_contribution = contribution(norm_hd + hidden_boost, decay_hd, hd_hid, centers_hd[i], periodic=T)
      hidden_pos_contribution = contribution(norm_pos + hidden_boost, decay_pos, x_hid, centers_x[i], y_hid, centers_y[i])
      hidden_speed_contribution = contribution(norm_spe + hidden_boost, decay_spe, s_hid, centers_s[i])
      hidden_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - hidden_pos_contribution - hidden_speed_contribution))
      Y[,i] = rbinom(N, 1, hidden_rate)
    }
    return(Y)
  }
  
  hd_obs = track_obs$HD[1:N]
  x_obs = track_obs$headX[1:N]
  y_obs = track_obs$headY[1:N]
  s_obs = normalize_speed(track_obs$bodyspeed[1:N])
  
  if (PLOTTHINGS1){
    hidden_hd_contribution = contribution(norm_hd + hidden_boost, decay_hd, hd_hid, centers_hd[1], periodic=T)
    pos_contribution = contribution(norm_pos, decay_pos, x_obs, centers_x[1], y_obs, centers_y[1])
    hidden_speed_contribution = contribution(norm_spe + hidden_boost, decay_spe, s_hid, centers_s[1])
    combined_pos_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - pos_contribution - hidden_speed_contribution))
    
    hidden_hd_contribution = contribution(norm_hd + hidden_boost, decay_hd, hd_hid, centers_hd[4], periodic=T)
    hidden_pos_contribution = contribution(norm_pos, decay_pos, x_hid, centers_x[4], y_hid, centers_y[4])
    spe_contribution = contribution(norm_spe + hidden_boost, decay_spe, s_obs, centers_s[4])
    combined_spe_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - hidden_pos_contribution - spe_contribution))
    
    hd_contribution = contribution(norm_hd, decay_hd, hd_obs, centers_hd[2], periodic=T)
    hidden_pos_contribution = contribution(norm_pos + hidden_boost, decay_pos, x_hid, centers_x[2], y_hid, centers_y[2])
    hidden_speed_contribution = contribution(norm_spe + hidden_boost, decay_spe, s_hid, centers_s[2])
    combined_hd_rate = 1 / (1 + exp(intercept - hd_contribution - hidden_pos_contribution - hidden_speed_contribution))
    
    hidden_hd_contribution = contribution(norm_hd + hidden_boost, decay_hd, hd_hid, centers_hd[3])
    hidden_pos_contribution = contribution(norm_pos + hidden_boost, decay_pos, x_hid, centers_x[3], y_hid, centers_y[3])
    only_hidden_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - hidden_pos_contribution - hidden_speed_contribution))
    if (PLOTTHINGS2){
      hdknots = 6
      posknots = 4
      speknots = 5
      x_vec = seq(min(track_obs$headY[1:12000]), max(track_obs$headX[1:12000]),length.out = 120)
      
      hd_vec = seq(-pi,pi,length.out = 200)
      xvec = rep(x_vec, 120)
      yvec = rep(x_vec, each=120)
      
      XHD1 = make_design_matrix(data.frame(hd_vec), c(),c(1),c(),c(6))
      XHD2 = make_design_matrix(data.frame(hd_vec), c(),c(1),c(),c(4))
      XHD3 = make_design_matrix(data.frame(hd_vec), c(),c(1),c(),c(3))
      XPOS1 = make_design_matrix(data.frame(cbind(xvec,yvec)), c(),c(),c(1,2),c(4,4))
      XPOS2 = make_design_matrix(data.frame(cbind(xvec,yvec)), c(),c(),c(1,2),c(3,3))
      XPOS3 = make_design_matrix(data.frame(cbind(xvec,yvec)), c(),c(),c(1,2),c(2,2))
      
      X1 = make_design_matrix(track_obs[1:12000,c(3,1,2,9)],c(4),c(1),c(2),c(hdknots,posknots,posknots,speknots))
      X2 = make_design_matrix(track_obs[1:12000,c(3,1,2,9)],c(4),c(1),c(2),c(4,3,3,4))
      X3 = make_design_matrix(track_obs[1:12000,c(3,1,2,9)],c(4),c(1),c(2),c(3,2,2,3))
      
      X1 = make_design_matrix(track_obs[1:12000,c(3,1,2,9)],c(4),c(1),c(2),c(hdknots,posknots,posknots,speknots))
      X2 = make_design_matrix(track_obs[1:12000,c(3,1,2,9)],c(4),c(1),c(2),c(4,3,3,4))
      X3 = make_design_matrix(track_obs[1:12000,c(3,1,2,9)],c(4),c(1),c(2),c(3,2,2,3))
      covinds1 = list(`1`=1:hdknots,`2`=(hdknots+1):(hdknots+(1+posknots)^2),`3`=(1+hdknots+(1+posknots)^2):dim(X1)[2])
      covinds2 = list(`1`=1:4,`2`=(4+1):(4+(1+3)^2),`3`=(1+4+(1+3)^2):dim(X2)[2])
      covinds3 = list(`1`=1:3,`2`=(3+1):(3+(1+2)^2),`3`=(1+3+(1+2)^2):dim(X3)[2])
    }
    hd_rate = 1 / (1 + exp(intercept - hd_contribution))
    pos_rate = 1 / (1 + exp(intercept - pos_contribution))
    y1 = rbinom(N, 1, pos_rate)
    y2 = rbinom(N, 1, combined_pos_rate)
    y3 = rbinom(N, 1, hd_rate)
    y4 = rbinom(N, 1, combined_hd_rate)
    y5 = rbinom(N, 1, only_hidden_rate)
    y6 = rbinom(N, 1, combined_spe_rate)
    if (PLOTTHINGS2){
      df1hd = data.frame(cbind(y4,X1[,covinds1$`1`]))
      df2hd = data.frame(cbind(y4,X2[,covinds2$`1`]))
      df3hd = data.frame(cbind(y4,X3[,covinds3$`1`]))
      ndf1hd = data.frame(cbind(rep(1,200),XHD1))
      ndf2hd = data.frame(cbind(rep(1,200),XHD2))
      ndf3hd = data.frame(cbind(rep(1,200),XHD3))
      df1pos = data.frame(cbind(y2,X1[,covinds1$`2`]))
      df2pos = data.frame(cbind(y2,X2[,covinds2$`2`]))
      df3pos = data.frame(cbind(y2,X3[,covinds3$`2`]))
      ndf1pos = data.frame(cbind(rep(1,14400),XPOS1))
      ndf2pos = data.frame(cbind(rep(1,14400),XPOS2))
      ndf3pos = data.frame(cbind(rep(1,14400),XPOS3))
      
      fit1hd = glm(y4~., data=df1hd,family = "binomial")
      fit2hd = glm(y4~., data=df2hd,family = "binomial")
      fit3hd = glm(y4~., data=df3hd,family = "binomial")
      fit1pos = glm(y2~., data=df1pos,family = "binomial")
      fit2pos = glm(y2~., data=df2pos,family = "binomial")
      fit3pos = glm(y2~., data=df3pos,family = "binomial")
      
      p1hd = predict(fit1hd, newdata=ndf1hd, type="response")*7.25
      p2hd = predict(fit2hd, newdata=ndf2hd, type="response")*7.25
      p3hd = predict(fit3hd, newdata=ndf3hd, type="response")*7.25
      
      p1pos = predict(fit1pos, newdata=ndf1pos, type="response")*7.25
      p2pos = predict(fit2pos, newdata=ndf2pos, type="response")*7.25
      p3pos = predict(fit3pos, newdata=ndf3pos, type="response")*7.25
      rm3 = ratemapXY(p1pos, p2pos, p3pos, X=xvec,Y=yvec, binning=120, smoothingsd = c(0,0,0,0), binsize = 1/7.25)
    }
    
    rm = ratemapXY(pos_rate, pos_contribution,Y1[,5],combined_pos_rate, X=x_obs,Y=y_obs, binning=40, smoothingsd = c(0,0,1,0), binsize = 1/7.25)
    rm2 = ratemapXY(combined_pos_rate, y2, y5, X=x_obs,Y=y_obs, binning=40, smoothingsd = c(0,1,1,0), binsize = 1/7.25)
    
    rm4 = ratemapXY(hd_obs,hd_obs, X=x_hid,Y=y_hid, binning=40, smoothingsd = c(0,1,0,0), binsize = 1/7.25)
    
    par(mfrow=c(4,2),pty="s",mar=c(1,1,1,1))
    image.plot(rm$rm1)
    #image.plot(rm$rm2)
    image.plot(rm$rm3)
    image.plot(rm$rm4)
    #image.plot(rm2$rm1)
    image.plot(rm2$rm2)
    image.plot(rm2$rm3)
    if (PLOTTHINGS2){
      image.plot(rm3$rm1)
      image.plot(rm3$rm2)
      image.plot(rm3$rm3)
    }  
    tc = tuningcurve(y3,y4,cyclic_shift(y4),y5, X=hd_obs,binning=40, binsize=1/7.25, smoothingsd=c(0,0,0,0),PERIODIC = T)
    tc2 = tuningcurve(y6, X=s_obs,binning=40, binsize=1/7.25, smoothingsd=c(0,0,0,0),PERIODIC =F)
    #plot(tc$tc1,type="l",ylim=c(0, 1.05*max(tc$tc1)))
    plot(tc$tc2,type="l",ylim=c(0, 1.05*max(tc$tc2)))
    #plot(tc$tc3,type="l",ylim=c(0, 1.05*max(tc$tc2)))
    plot(tc$tc4,type="l",ylim=c(0, 1.05*max(tc$tc2)))
    plot(tc2$tc1,type="l",ylim=c(0, 1.05*max(tc2$tc1)))
    if (PLOTTHINGS2){
      lines(BinMean(p1hd, every=5),col="red")
      lines(BinMean(p2hd, every=5),col="green")
      lines(BinMean(p3hd, every=5),col="blue")
    }
    #plot(tc$tc3,type="l",ylim=c(0, 1.05*max(tc$tc3)))
    #plot(tc$tc4,type="l",ylim=c(0, 1.05*max(tc$tc4)))
    
    print(paste(mean(y2), mean(y4), mean(y5)))
  } # For testing and plotting
  
  half = round(num_cells / 2)
  Y = matrix(0,nrow=N,ncol=num_cells)
  for (i in 1:half){
    hidden_hd_contribution = contribution(norm_hd + hidden_boost, decay_hd, hd_hid, centers_hd[i], periodic=T)
    pos_contribution = contribution(norm_pos, decay_pos, x_obs, centers_x[i], y_obs, centers_y[i])
    combined_pos_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - pos_contribution))
    Y[,i] = rbinom(N, 1, combined_pos_rate)
  }
  for (i in (half + 1):num_cells){
    hd_contribution = contribution(norm_hd, decay_hd, hd_obs, centers_hd[i], periodic=T)
    hidden_pos_contribution = contribution(norm_pos + hidden_boost, decay_pos, x_hid, centers_x[i], y_hid, centers_y[i])
    combined_hd_rate = 1 / (1 + exp(intercept - hd_contribution - hidden_pos_contribution))
    Y[,i] = rbinom(N, 1, combined_hd_rate)
  }
  return(Y)
}

make_simulated_data_from_real_covariates = function(num_cells=300,HIDDEN_ONLY=F){
  # Simulates data with signal from hd and position,
  # but where one or both are from a different experiment than
  # the one where the covariates are treated as observed.
  
  track_obs = readRDS(paste(processed_data_folder,"trackdata2.RDS",sep=""))
  track_hid = readRDS(paste(processed_data_folder,"trackdata2_other_date.RDS",sep=""))
  
  N = 12000
  hd_hid = track_hid$HD[1:N]
  x_hid = track_hid$headX[1:N]
  y_hid = track_hid$headY[1:N]
  s_hid = normalize_speed(track_hid$bodyspeed[1:N])
  
  intercept = 4.6
  
  norm_pos = 1.5
  norm_hd = 1
  norm_spe = 1.5
  
  decay_pos = 0.003
  decay_hd = 3
  decay_spe = 1
  
  # Draw random centers for the tuning peaks
  centers_hd = runif(num_cells, min=-pi, max=pi)
  centers_x = runif(num_cells, min=-35, max=35)
  centers_y = runif(num_cells, min=-35, max=35)
  centers_s = runif(num_cells, min=-0.75, max=0.75)
  
  if (HIDDEN_ONLY){
    # Simulate cells only depending on unobserved covariates, for checking type I error of methods
    Y = matrix(0,nrow=N,ncol=num_cells)
    for (i in 1:num_cells){
      hidden_hd_contribution = contribution(norm_hd, decay_hd, hd_hid, centers_hd[i], periodic=T)
      hidden_pos_contribution = contribution(norm_pos, decay_pos, x_hid, centers_x[i], y_hid, centers_y[i])
      hidden_spe_contribution = contribution(norm_spe, decay_spe, s_hid, centers_s[i])
      hidden_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - hidden_pos_contribution - hidden_spe_contribution))
      Y[,i] = rbinom(N, 1, hidden_rate)
    }
    return(Y)
  }
  
  hd_obs = track_obs$HD[1:N]
  x_obs = track_obs$headX[1:N]
  y_obs = track_obs$headY[1:N]
  s_obs = normalize_speed(track_obs$bodyspeed[1:N])
  
  third = round(num_cells / 3)
  twothirds = round(num_cells * 2 / 3)
  Y = matrix(0,nrow=N,ncol=num_cells)
  for (i in 1:third){
    # Half of the cells are simulated as position tuned,
    # where the cells also are responsive to the unobserved head direction
    # in the separate experiment
    hidden_hd_contribution = contribution(norm_hd, decay_hd, hd_hid, centers_hd[i], periodic=T)
    pos_contribution = contribution(norm_pos, decay_pos, x_obs, centers_x[i], y_obs, centers_y[i])
    hidden_spe_contribution = contribution(norm_spe, decay_spe, s_hid, centers_s[i])
    combined_pos_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - pos_contribution - hidden_spe_contribution))
    Y[,i] = rbinom(N, 1, combined_pos_rate)
  }
  for (i in (third + 1):twothirds){
    # Vice versa
    hd_contribution = contribution(norm_hd, decay_hd, hd_obs, centers_hd[i], periodic=T)
    hidden_pos_contribution = contribution(norm_pos, decay_pos, x_hid, centers_x[i], y_hid, centers_y[i])
    hidden_spe_contribution = contribution(norm_spe, decay_spe, s_hid, centers_s[i])
    combined_hd_rate = 1 / (1 + exp(intercept - hd_contribution - hidden_pos_contribution - hidden_spe_contribution))
    Y[,i] = rbinom(N, 1, combined_hd_rate)
  }
  for (i in (twothirds + 1):num_cells){
    # Vice versa
    hidden_hd_contribution = contribution(norm_hd, decay_hd, hd_hid, centers_hd[i], periodic=T)
    hidden_pos_contribution = contribution(norm_pos, decay_pos, x_hid, centers_x[i], y_hid, centers_y[i])
    spe_contribution = contribution(norm_spe, decay_spe, s_obs, centers_s[i])
    combined_hd_rate = 1 / (1 + exp(intercept - hidden_hd_contribution - hidden_pos_contribution - spe_contribution))
    Y[,i] = rbinom(N, 1, combined_hd_rate)
  }
  return(Y)
}

make_simpler_simulated_data = function(N,sd,b_proportion=0,numys=1,ncov=3){
  #h is hidden, a, b, c are observed
  h = make_1D_cov(N,sd,dist="uniform",smoother="exponential",CENTER_AND_SCALE=T)
  a = make_1D_cov(N,sd,dist="uniform",smoother="exponential",CENTER_AND_SCALE=T)
  b = make_1D_cov(N,sd,dist="uniform",smoother="exponential",CENTER_AND_SCALE=T)
  c = make_1D_cov(N,sd,dist="uniform",smoother="exponential",CENTER_AND_SCALE=T)
  
  const_term = -3 # rate ≈ 0.05
  total_factor = 10 # maxrate ≈ 0.50 (if relevant things are maximal)
  
  p = 1 / (1 + exp(-(const_term + total_factor * ((1 - b_proportion) * h + b_proportion * b))))
  
  if (numys > 1){
    ys = matrix(0,nrow=N,ncol=numys)
    for (i in 1:numys){
      ys[,i] = rbinom(N,1,p)
    }
    return(ys)
  }
  y = rbinom(N,1,p)
  dat = data.frame(cbind(y,a,b,c))
  
  if (ncov > 3){
    for (i in 4:ncov){
      d = make_1D_cov(N,sd,dist="uniform",smoother="exponential",CENTER_AND_SCALE=T)
      dat[[paste("d",i,sep="")]] = d
    }
  }
  
  return(dat)
}

make_1D_cov = function(N,sd,dist="normal",smoother="exponential",bound=0.3,power=1,uniscale=3,CENTER_AND_SCALE=F){
  if (smoother == "exponential"){
    extra = 8
  } else {
    extra = 4
  }
  if (dist == "normal"){
    cov = rnorm(N+extra*sd)
  } else {
    cov = runif(N+extra*sd)*2*uniscale - uniscale
  }
  if (smoother == "exponential"){
    cov = exponential_smoother1D(cov,sd,power=power)
  } else {
    cov = gaussian_smoother1D(cov,sd)
  }
  cov = cov[(extra/2*sd+1):(extra/2*sd+N)]
  if (CENTER_AND_SCALE){
    cov = cov - median(cov)
    cov = cov / max(abs(cov)) * bound
  }
  cov = add_boundaries_to_continuous_signal(cov,bound,-bound)
  return(cov)
}

#### Generate some autocorrelated things ----

generate_X = function(n,smoothing_kernel="gaussian",smoothing_parameter=5,ncov=3){
  if (smoothing_parameter==0){
    X = matrix(rnorm(n*ncov),nrow=n)
    return(X)
  }
  X_ = matrix(rnorm((n+100)*ncov),nrow=n+100)
  X = matrix(0,nrow=n,ncol=ncov)
  for (i in 1:ncov){
    if (smoothing_kernel == "gaussian"){
      X[,i] = gaussian_smoother_invariant(X_[,i],smoothing_parameter)[51:(n+50)]
    } else if (smoothing_kernel == "exponential"){
      X[,i] = exponential_smoother1D(X_[,i],smoothing_parameter)[51:(n+50)]
    }
  }
  return(X)
}

generate_Y = function(n,X,coef1=1.5,coef2=0.5,base=1){
  Y = rpois(n,base*exp(coef1*X[,1]+coef2*X[,2]))
  return(Y)
}
#### Effect sizes ----
AddMcFadden = function(resmat, m_ind, X, Ymat, cov_inds, keeps=NULL){
  n = dim(resmat)[3]
  reslist = list(
    "hd"=array(NA,n),
    "nohd"=array(NA,n),
    "pos"=array(NA,n),
    "nopos"=array(NA,n),
    "spe"=array(NA,n),
    "nospe"=array(NA,n)
  )
  family="binomial"
  in_names = c("hd", "pos", "spe")
  out_names = c("nohd", "nopos", "nospe")
  num_methods = (dim(resmat)[1] - 1)/2
  
  start = proc.time()[3]
  
  for (i in 1:n){
    Y = Ymat[,i]
    included = resmat[1,which(resmat[m_ind,,i] <= 0.05),i]
    for (j in 1:3){
      if (j %in% included){
        when = which(resmat[1,,i] == j)
        if (FALSE){
          if (when==1){
            X_ = X[,cov_inds[[j]]]
            fit0 = glm(Y~1,family=family)
          } else {
            X0 = X[,cov_inds[[resmat[1,1,i]]]]
            if (when == 3){
              X0 = cbind(X0, X[,cov_inds[[resmat[1,2,i]]]])
            }
            X_ = cbind(X0, X[,cov_inds[[j]]])
            fit0 = glm(Y~X0,family=family)
          }
          fit = glm(Y~X_,family=family)
          mcFadden = 1 - logLik(fit) / logLik(fit0)#(logLik(fit) - length(cov_inds[[j]])) / logLik(fit0)
        } else {
          mcFadden = resmat[num_methods + m_ind, when, i]
        }
        reslist[[in_names[j]]][i] = mcFadden
      } else {
        len_incl = length(included)
        if (resmat[1, len_incl + 1, i] == j){
          mcFadden = resmat[num_methods + m_ind, len_incl + 1, i]
          reslist[[out_names[j]]][i] = mcFadden
        }
        else {
          if (is.null(keeps)){
            keeps = 1:length(Y)
          }
          if (len_incl == 0){
            X_ = X[keeps,cov_inds[[j]]]
            fit0 = glm(Y[keeps]~1,family=family)
          } else {
            X0 = X[keeps,cov_inds[[included[1]]]]
            if (length(included) == 2){
              X0 = cbind(X0, X[keeps,cov_inds[[included[2]]]])
            }
            X_ = cbind(X0, X[keeps,cov_inds[[j]]])
            fit0 = glm(Y[keeps]~X0,family=family)
          }
          fit = glm(Y[keeps]~X_,family=family)
          mcFadden = 1 - logLik(fit) / logLik(fit0)#(logLik(fit) - length(cov_inds[[j]])) / logLik(fit0)
          reslist[[out_names[j]]][i] = mcFadden
        }
      }
    }
    print_progress(i, n, start)
  }
  return(reslist)
}
AddSRRev = function(resmat, m_ind, X, Ymat, cov_inds, config){
  n = dim(resmat)[3]
  reslist = list(
    "hd"=array(NA,n),
    "nohd"=array(NA,n),
    "pos"=array(NA,n),
    "nopos"=array(NA,n),
    "spe"=array(NA,n),
    "nospe"=array(NA,n)
  )
  in_names = c("hd", "pos", "spe")
  out_names = c("nohd", "nopos", "nospe")
  N = config$N
  num_methods = (dim(resmat)[1] - 1)/2
  
  start = proc.time()[3]
  for (i in 1:n){
    Y = Ymat[,i]
    included = resmat[1,which(resmat[m_ind,,i] <= 0.05),i]
    for (j in 1:3){
      if (j %in% included){
        when = which(resmat[1,,i] == j)
        if (FALSE){
          if (when==1){
            X_ = X[,cov_inds[[j]]]
            X0 = X[rev(1:N),cov_inds[[j]]]
          } else {
            X0 = X[,cov_inds[[resmat[1,1,i]]]]
            if (when == 3){
              X0 = cbind(X0, X[,cov_inds[[resmat[1,2,i]]]])
            }
            X_ = cbind(X0, X[,cov_inds[[j]]])
            X0 = cbind(X0, X[rev(1:N),cov_inds[[j]]])
          }
          folds0 = CV_folds(Y,X0,config)
          folds1 = CV_folds(Y,X_,config)
          SR = signrank(folds1 - folds0)
        } else {
          SR = resmat[num_methods + m_ind, when, i]
        }
        reslist[[in_names[j]]][i] = SR
      } else {
        len_incl = length(included)
        if (resmat[1, len_incl + 1, i] == j){
          SR = resmat[num_methods + m_ind, len_incl + 1, i]
        } else {
          if (len_incl == 0){
            X_ = X[,cov_inds[[j]]]
            X0 = X[rev(1:N),cov_inds[[j]]]
          } else {
            X0 = X[,cov_inds[[included[1]]]]
            if (length(included) == 2){
              X0 = cbind(X0, X[,cov_inds[[included[2]]]])
            }
            X_ = cbind(X0, X[,cov_inds[[j]]])
            X0 = cbind(X0, X[rev(1:N),cov_inds[[j]]])
          }
          folds0 = CV_folds(Y,X0,config)
          folds1 = CV_folds(Y,X_,config)
          SR = signrank(folds1 - folds0)
        }
        reslist[[out_names[j]]][i] = SR
      }
    }
    print_progress(i, n, start)
  }
  return(reslist)
}
AddSR = function(resmat, m_ind, X, Ymat, cov_inds, config){
  n = dim(resmat)[3]
  reslist = list(
    "hd"=array(NA,n),
    "nohd"=array(NA,n),
    "pos"=array(NA,n),
    "nopos"=array(NA,n),
    "spe"=array(NA,n),
    "nospe"=array(NA,n)
  )
  in_names = c("hd", "pos", "spe")
  out_names = c("nohd", "nopos", "nospe")
  num_methods = (dim(resmat)[1] - 1)/2
  
  start = proc.time()[3]
  for (i in 1:n){
    Y = Ymat[,i]
    included = resmat[1,which(resmat[m_ind,,i] <= 0.05),i]
    for (j in 1:3){
      if (j %in% included){
        when = which(resmat[1,,i] == j)
        if (FALSE){
          if (when == 1){
            X_ = X[,cov_inds[[j]]]
            folds0 = CV_folds(Y,NULL,config,ONLY_INTERCEPT=T)
          } else {
            X0 = X[,cov_inds[[resmat[1,1,i]]]]
            if (when == 3){
              X0 = cbind(X0, X[,cov_inds[[resmat[1,2,i]]]])
            }
            X_ = cbind(X0, X[,cov_inds[[j]]])
            folds0 = CV_folds(Y,X0,config)
          }
          folds1 = CV_folds(Y,X_,config)
          SR = signrank(folds1 - folds0)
        } else {
          SR = resmat[num_methods + m_ind, when, i]
        }
        reslist[[in_names[j]]][i] = SR
      } else {
        len_incl = length(included)
        if (resmat[1, len_incl + 1, i] == j){
          SR = resmat[num_methods + m_ind, len_incl + 1, i]
        } else {
          if (len_incl == 0){
            X_ = X[,cov_inds[[j]]]
            folds0 = CV_folds(Y,NULL,config,ONLY_INTERCEPT=T)
          } else {
            X0 = X[,cov_inds[[included[1]]]]
            if (length(included) == 2){
              X0 = cbind(X0, X[,cov_inds[[included[2]]]])
            }
            X_ = cbind(X0, X[,cov_inds[[j]]])
            folds0 = CV_folds(Y,X0,config)
          }
          folds1 = CV_folds(Y,X_,config)
          SR = signrank(folds1 - folds0)
        }
        reslist[[out_names[j]]][i] = SR
      }
    }
    print_progress(i, n, start)
  }
  return(reslist)
}

plot_effect_sizes = function(filt_list, FAKE=F, lrt_lims=NULL, meanrates=NULL,SIM=F,cex=NULL){
  if (is.null(cex)){
    cex = par("cex")
  }
  methods = c("LRT", "CS", "SR", "mSR", "mSRR")
  par(mfrow=c(length(methods), 3))
  
  
  xlabs = list(
    "LRT"="McFadden's likelihood-ratio index",
    "CS"="McFadden's likelihood-ratio index",
    "SR"="Rank biserial correlation",
    "mSR"="Rank biserial correlation",
    "mSRR"="Rank biserial correlation"
  )
  
  hdcont = "Effect of HD"
  poscont = "Effect of position"
  specont = "Effect of speed"
  opaq = 0.7
  c_green = rgb(0/255,158/255,115/255,opaq) # True positive
  c_blue = rgb(86/255,180/255,233/255,opaq) # True negative
  c_red = rgb(213/255,94/255,0/255,opaq) # False positive
  c_yellow = rgb(240/255,228/255,66/255,opaq) # False negative
  
  if (FAKE){
    which_filt = filt_list$fake
    leg = c("False positive", "True negative")
    c1 = c_red
    c2 = c_blue
  } else {
    which_filt = filt_list$real
    if (SIM){
      leg = c("False positive", "True positive", "False negative", "True negative")
      c1 = c_red
      c2 = c_green
      c3 = c_yellow
      c4 = c_blue
    } else {
      leg = c("In model", "Not in model")
      c1 = c_green
      c2 = c_red
    }
  }
  breaks_list = list(
    "LRT"=seq(0,1,by=0.025),
    "CS"=seq(0,1,by=0.025),
    "SR"=seq(-1,1,by=0.05),
    "mSR"=seq(-1,1,by=0.05),
    "mSRR"=seq(-1,1,by=0.05)
  )
  
  for (i in 1:length(methods)){
    method = methods[i]
    breaks1 = breaks2 = breaks3 = breaks_list[[method]]
    xlab = xlabs[[method]]
    fac = 1
    if (method %in% c("mSR", "mSRR")){
      fac = 210
    } else if (method == "SR"){
      fac = 55
    }
    if (method %in% names(which_filt)){
      filt = which_filt[[method]]
      if (method %in% c("LRT", "CS")){
        if (!is.null(meanrates)){
          lists = names(filt)
          for (i in 1:length(lists)){
            filt[[lists[i]]] = filt[[lists[i]]] * meanrates
          }
        }
        if (!is.null(lrt_lims)){
          breaks1 = seq(0, lrt_lims[1],by=lrt_lims[1]/40)
          breaks2 = seq(0, lrt_lims[2],by=lrt_lims[2]/40)
          breaks3 = seq(0, lrt_lims[3],by=lrt_lims[3]/40)
        }
      }
      
      if (SIM & !FAKE){
        maxcounthddown = max(c(
          max(hist(filt$hd[101:200]/fac, breaks=breaks1, plot=F)$counts),
          max(hist(filt$nohd[101:200]/fac, breaks=breaks1, plot=F)$counts)
        ))
        maxcounthdup = max(c(
          max(hist(filt$hd[c(1:100, 201:300)]/fac, breaks=breaks1, plot=F)$counts),
          max(hist(filt$nohd[c(1:100, 201:300)]/fac, breaks=breaks1, plot=F)$counts)
        ))
        h1 = hist(filt$hd[c(1:100, 201:300)]/fac,breaks=breaks1, plot=F)
        h2 = hist(filt$hd[101:200]/fac,breaks=breaks1, plot=F)
        h3 = hist(filt$nohd[101:200]/fac,breaks=breaks1, plot=F)
        h4 = hist(filt$nohd[c(1:100, 201:300)]/fac,breaks=breaks1, plot=F)
        
        h2$counts = - h2$counts
        h3$counts = - h3$counts
        
        plot(h1, col=c1, main=hdcont,xlab=xlab,ylab=method,ylim=c(-maxcounthddown,maxcounthdup),yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        axis(2, at=c(0, maxcounthdup), tick=T, labels=c(0, maxcounthdup),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        axis(4, at=c(-maxcounthddown, 0), tick=T, labels=c(maxcounthddown, 0),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        plot(h2, col=c2, add=T)
        plot(h3, col=c3, add=T)
        plot(h4, col=c4, add=T)
        
        maxcountposdown = max(c(
          max(hist(filt$pos[1:100]/fac, breaks=breaks2, plot=F)$counts),
          max(hist(filt$nopos[1:100]/fac, breaks=breaks2, plot=F)$counts)
        ))
        maxcountposup = max(c(
          max(hist(filt$pos[101:300]/fac, breaks=breaks2, plot=F)$counts),
          max(hist(filt$nopos[101:300]/fac, breaks=breaks2, plot=F)$counts)
        ))
        
        h1 = hist(filt$pos[101:300]/fac,breaks=breaks2, plot=F)
        h2 = hist(filt$pos[1:100]/fac,breaks=breaks2, plot=F)
        h3 = hist(filt$nopos[1:100]/fac,breaks=breaks2, plot=F)
        h4 = hist(filt$nopos[101:300]/fac,breaks=breaks2, plot=F)
        
        h2$counts = - h2$counts
        h3$counts = - h3$counts
        
        plot(h1, col=c1, main=poscont,xlab=xlab,ylim=c(-maxcountposdown,maxcountposup),yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        axis(2, at=c(0, maxcountposup), tick=T, labels=c(0, maxcountposup),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        axis(4, at=c(-maxcountposdown, 0), tick=T, labels=c(maxcountposdown, 0),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        plot(h2, col=c2, add=T)
        plot(h3, col=c3, add=T)
        plot(h4, col=c4, add=T)
        
        maxcountspedown = max(c(
          max(hist(filt$spe[201:300]/fac, breaks=breaks3, plot=F)$counts),
          max(hist(filt$nospe[201:300]/fac, breaks=breaks3, plot=F)$counts)
        ))
        maxcountspeup = max(c(
          max(hist(filt$spe[1:200]/fac, breaks=breaks3, plot=F)$counts),
          max(hist(filt$nospe[1:200]/fac, breaks=breaks3, plot=F)$counts)
        ))
        
        h1 = hist(filt$spe[1:200]/fac,breaks=breaks3, plot=F)
        h2 = hist(filt$spe[201:300]/fac,breaks=breaks3, plot=F)
        h3 = hist(filt$nospe[201:300]/fac,breaks=breaks3, plot=F)
        h4 = hist(filt$nospe[1:200]/fac,breaks=breaks3, plot=F)
        
        h2$counts = - h2$counts
        h3$counts = - h3$counts
        
        plot(h1, col=c1, main=specont,xlab=xlab,ylim=c(-maxcountspedown,maxcountspeup),yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        axis(2, at=c(0, maxcountspeup), tick=T, labels=c(0, maxcountspeup),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        axis(4, at=c(-maxcountspedown, 0), tick=T, labels=c(maxcountspedown, 0),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        plot(h2, col=c2, add=T)
        plot(h3, col=c3, add=T)
        plot(h4, col=c4, add=T)
        
        if (i == 1){
          legend("topright", legend=leg,fill=c(c1, c2, c3, c4),cex=cex*0.9)
        }
      } else {
        maxcounthd = max(c(max(hist(filt$hd/fac, breaks=breaks1, plot=F)$counts), max(hist(filt$nohd/fac, breaks=breaks1, plot=F)$counts)))
        hist(filt$hd/fac,breaks=breaks1, main=hdcont,xlab=xlab,ylab=method, col=c1,ylim=c(0,maxcounthd),yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        axis(2, at=c(0, maxcounthd), tick=T, labels=c(0, maxcounthd),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        hist(filt$nohd/fac,breaks=breaks1, main=hdcont,xlab=xlab, col=c2, add=T)
        maxcountpos = max(c(max(hist(filt$pos/fac, breaks=breaks2, plot=F)$counts), max(hist(filt$nopos/fac, breaks=breaks2, plot=F)$counts)))
        hist(filt$pos/fac,breaks=breaks2, main=poscont,xlab=xlab, col=c1,ylim=c(0,maxcountpos),yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        axis(2, at=c(0, maxcountpos), tick=T, labels=c(0, maxcountpos),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        hist(filt$nopos/fac,breaks=breaks2, main=poscont,xlab=xlab, col=c2, add=T)
        maxcountspe = max(c(max(hist(filt$spe/fac, breaks=breaks3, plot=F)$counts), max(hist(filt$nospe/fac, breaks=breaks3, plot=F)$counts)))
        hist(filt$spe/fac,breaks=breaks3, main=specont,xlab=xlab, col=c1,ylim=c(0,maxcountspe),yaxt="n",cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
        axis(2, at=c(0, maxcountspe), tick=T, labels=c(0, maxcountspe),cex.lab=cex,cex=cex, cex.axis=cex,lwd=cex)
        hist(filt$nospe/fac,breaks=breaks3, main=specont,xlab=xlab, col=c2, add=T)
        if (i == 1){
          legend("topright", legend=leg,fill=c(c1, c2),cex=cex*0.9)
        }
      }
    } else {
      print(paste("Method ", method, " not found in list of effect sizes", sep=""))
    }
  }
}

