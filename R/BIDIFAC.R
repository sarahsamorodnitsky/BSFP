BIDIFAC=function(data,rmt=T, sigma=NULL,
                 start=NULL, out=FALSE,
                 eps=1e-3, max.iter=1000, pbar=TRUE, seed=NULL, scale_back = TRUE, ...){
  if (!is.null(seed)){set.seed(seed)}
  # if (!rmt & class(sigma)!="matrix") stop("sigma must be a matrix.")
  if (!rmt & !("matrix" %in% class(sigma))) stop("sigma must be a matrix.")

  fit=data.rearrange(data, rmt, sigma)
  sigma.mat=fit$sigma.mat
  X00=fit$out

  mvec=fit$nrows; nvec=fit$ncols
  p=length(mvec); q=length(nvec)
  rm(fit)

  start.ind.m=c(1, cumsum(mvec)[1:(p-1)]+1)
  end.ind.m=cumsum(mvec)

  start.ind.n=c(1, cumsum(nvec)[1:(q-1)]+1)
  end.ind.n=cumsum(nvec)

  lambda.G=sqrt(sum(mvec))+sqrt(sum(nvec))
  lambda.R=sqrt(mvec)+sqrt(sum(nvec))
  lambda.C=sqrt(sum(mvec))+sqrt(nvec)
  lambda.I=tcrossprod(sqrt(mvec), rep(1, length(nvec)))+
    tcrossprod(rep(1, length(mvec)),sqrt(nvec))

  if (!is.null(start)){
    G00=start[[1]]; R00=start[[2]]
    C00=start[[3]]; I00=start[[4]]
  } else {
    G00= matrix(0, nrow = sum(mvec), ncol = sum(nvec)) # replicate(sum(nvec),rnorm(sum(mvec)))
    R00= matrix(0, nrow = sum(mvec), ncol = sum(nvec)) # replicate(sum(nvec),rnorm(sum(mvec)))
    C00=replicate(sum(nvec),rnorm(sum(mvec)))
    I00=replicate(sum(nvec),rnorm(sum(mvec)))
  }

  G00.nuc=NA; R00.nuc=rep(NA, p)
  C00.nuc=rep(NA, q); I00.nuc=matrix(NA,p,q)

  bool=TRUE
  count=1;crit0=0
  if (pbar) pb = txtProgressBar(min = 0, max=max.iter, initial=0, char="-", style = 3)
  while (bool){
    if (pbar){  setTxtProgressBar(pb, count)  }
    crit0.old = crit0

    #Update G to 0
    fit1=softSVD(X00-R00-C00-I00,lambda.G)
    G00 <- matrix(0, nrow = sum(mvec), ncol = sum(nvec)) # G00=fit1$out;
    G00.nuc=fit1$nuc

    #update R to 0
    for (i in 1:p){
      ind=start.ind.m[i]:end.ind.m[i]
      fit1=softSVD(X00[ind,]-G00[ind,]-C00[ind,]-I00[ind,], lambda.R[i])
      # R00[ind,]=fit1$out;
      R00 <- matrix(0, nrow = sum(mvec), ncol = sum(nvec))
      R00.nuc[i]=fit1$nuc
    }

    for (j in 1:q){
      ind=start.ind.n[j]:end.ind.n[j]
      fit1=softSVD(X00[,ind]-G00[,ind]-R00[,ind]-I00[,ind], lambda.C[j])
      C00[,ind]=fit1$out; C00.nuc[j]=fit1$nuc
    }

    for (i in 1:p){
      for (j in 1:q){
        ind1= start.ind.m[i]:end.ind.m[i]
        ind2=start.ind.n[j]:end.ind.n[j]
        fit1=softSVD(X00[ind1,ind2]-G00[ind1,ind2]-R00[ind1,ind2]-C00[ind1,ind2], lambda.I[i,j])
        I00[ind1,ind2]=fit1$out; I00.nuc[i,j]=fit1$nuc
      }
    }

    crit0 = frob(X00-G00-R00-C00-I00)+
      2*lambda.G*G00.nuc+2*sum(lambda.R*R00.nuc)+
      2*sum(lambda.C*C00.nuc)+2*sum(lambda.I*I00.nuc)

    if (abs(crit0.old-crit0)<eps){ bool=FALSE }
    else if (count==max.iter){ bool=FALSE}
    else{ count = count+1 }
  }

  if (scale_back) {
    S00.mat=G00.mat=R00.mat=C00.mat=I00.mat=data
    for (i in 1:p){
      ind1= start.ind.m[i]:end.ind.m[i]
      for (j in 1:q){
        ind2=start.ind.n[j]:end.ind.n[j]
        G00.mat[[i,j]]=G00[ind1,ind2] *sigma.mat[i,j]
        R00.mat[[i,j]]=R00[ind1,ind2] *sigma.mat[i,j]
        C00.mat[[i,j]]=C00[ind1,ind2] *sigma.mat[i,j]
        I00.mat[[i,j]]=I00[ind1,ind2] *sigma.mat[i,j]
        S00.mat[[i,j]]=G00.mat[[i,j]]+R00.mat[[i,j]]+C00.mat[[i,j]]+I00.mat[[i,j]]
      }
    }
  }

  if (!scale_back) {
    S00.mat=G00.mat=R00.mat=C00.mat=I00.mat=data
    for (i in 1:p){
      ind1= start.ind.m[i]:end.ind.m[i]
      for (j in 1:q){
        ind2=start.ind.n[j]:end.ind.n[j]
        G00.mat[[i,j]]=G00[ind1,ind2] #*sigma.mat[i,j] Remove the scaling back by the error variance
        R00.mat[[i,j]]=R00[ind1,ind2] #*sigma.mat[i,j]
        C00.mat[[i,j]]=C00[ind1,ind2] #*sigma.mat[i,j]
        I00.mat[[i,j]]=I00[ind1,ind2] #*sigma.mat[i,j]
        S00.mat[[i,j]]=G00.mat[[i,j]]+R00.mat[[i,j]]+C00.mat[[i,j]]+I00.mat[[i,j]]
      }
    }
  }

  return(list(X=data, S=S00.mat,
              G=G00.mat, R=R00.mat, C=C00.mat, I=I00.mat,
              sigma.mat=sigma.mat, n.vec=nvec,m.vec=mvec))
}

impute.BIDIFAC=function(data,
                        rmt=T, sigma=NULL,
                        pbar=TRUE,
                        start=NULL, max.iter.impute=100,
                        eps.impute=1e-3, scale_back=TRUE,...){
  dim.data=dim(data)
  p=dim.data[1]; q=dim.data[2]

  dim.list=do.call(cbind,lapply(data, dim))
  mvec=apply(matrix(dim.list[1,],p),1,unique)
  nvec=apply(matrix(dim.list[2,],p),2,unique)
  if (class(mvec)=="list" ) stop("the number of rows do not match")
  if (class(nvec)=="list" ) stop("the number of columns do not match")

  impute.index=matrix(list(), nrow = p, ncol=q)
  if (is.null(sigma)) sigma=matrix(1,p,q)

  for (i in 1:p){
    for (j in 1:q){

      if (any(is.na(data[[i,j]]))) {
        fillmat=fill.matrix(data[[i,j]])
        impute.index[[i,j]]=fillmat$na.ind
        if (rmt) sigma[i,j]=sigma.rmt(fillmat$X.fill)
        data[[i,j]]=fillmat$X.fill/sigma[i,j]
      } else {
        data[[i,j]] <- data[[i,j]]/sigma[i,j]
      }

    }
  }

  bool2=TRUE
  count2=1; impute.vec=0
  if (pbar) pb = txtProgressBar(min = 0, max=max.iter.impute, initial=0, char="-", style = 3)
  while (bool2){
    if (pbar){  setTxtProgressBar(pb, count2)  }
    impute.vec.old=impute.vec
    fit=BIDIFAC(data, rmt=F, sigma=matrix(1,p,q),start=start, pbar = F)

    start=list(
      data.rearrange(fit$G)$out, data.rearrange(fit$R)$out,
      data.rearrange(fit$C)$out, data.rearrange(fit$I)$out)

    impute.vec=NULL
    for (i in 1:p){
      for (j in 1:q){
        imp=fit$S[[i,j]][impute.index[[i,j]]]
        data[[i,j]][impute.index[[i,j]]]=imp
        impute.vec=c(impute.vec,imp)
      }
    }

    crit2=frob(impute.vec.old-impute.vec)/length(impute.vec)
    if (crit2<eps.impute){ bool2=FALSE }
    else if (count2==max.iter.impute){
      bool2=FALSE
      cat("\n The algorithm did not converge. \n")
      cat("Try to increase max.iter.impute or relax eps/eps.impute.")
    }
    else{ count2=count2+1 }
  }

  if (scale_back) {
    for (i in 1:p){
      for (j in 1:q){
        fit$X[[i,j]]=fit$X[[i,j]]*sigma[i,j]
        fit$S[[i,j]]=fit$S[[i,j]]*sigma[i,j]
        fit$G[[i,j]]=fit$G[[i,j]]*sigma[i,j]
        fit$R[[i,j]]=fit$R[[i,j]]*sigma[i,j]
        fit$C[[i,j]]=fit$C[[i,j]]*sigma[i,j]
        fit$I[[i,j]]=fit$I[[i,j]]*sigma[i,j]
      }
    }
  }

  if (!scale_back) {
    for (i in 1:p){
      for (j in 1:q){
        fit$X[[i,j]]=fit$X[[i,j]]#*sigma[i,j]
        fit$S[[i,j]]=fit$S[[i,j]]#*sigma[i,j]
        fit$G[[i,j]]=fit$G[[i,j]]#*sigma[i,j]
        fit$R[[i,j]]=fit$R[[i,j]]#*sigma[i,j]
        fit$C[[i,j]]=fit$C[[i,j]]#*sigma[i,j]
        fit$I[[i,j]]=fit$I[[i,j]]#*sigma[i,j]
      }
    }
  }

  fit$sigma.mat=sigma

  return(fit)
}

data.rearrange=function(data,rmt=F,sigma=NULL){
  out=NULL
  p=nrow(data)
  q=ncol(data)

  m.vec=rep(NA,p)
  n.vec= ncol(data[[1,1]]) # do.call(c, lapply(data[1,], ncol))

  if (is.null(sigma)) sigma=matrix(1,p,q)

  for (i in 1:p){
    dimm=do.call(cbind, lapply(data[i,],dim))
    m1=unique(dimm[1,])
    if (length(m1)==1 ){m.vec[i]=m1 }
    else{ stop("the number of rows do not match.") }
    if (!all(dimm[2,], n.vec)){ stop("the number of columns do not match")}

    for (j in 1:q){
      if (rmt) sigma[i,j]=sigma.rmt(data[[i,j]])
      data[[i,j]]=data[[i,j]]/sigma[i,j]
    }

    out=rbind(out,do.call(cbind,data[i,]))
  }

  return(list(out=out, nrows=m.vec, ncols=n.vec, sigma.mat=sigma))
}

frob <- function(X){ sum(X^2,na.rm=T) }

diag2 <- function(x) {
  if (length(x)==1) return(as.matrix(x))
  else if (length(x)>1) return(diag(x))
}

sample2 <- function(x) {
  if (length(x)==1) return(x)
  else if (length(x)>1) return(sample(x))
}

sigma.rmt=function(X){ estim_sigma(X,method="MAD") }

softSVD=function(X, lambda){
  svdX=svd(X)
  nuc=pmax(svdX$d-lambda,0)
  out=tcrossprod(svdX$u, tcrossprod( svdX$v,diag(nuc) ))
  return(list(out=out, nuc=sum(nuc)))
}

fill.matrix=function(X){
  na.ind.original=na.ind=which(is.na(X),arr.ind = T)
  bool=T
  while (bool){
    impute.X=rep(NA, nrow(na.ind))
    for (j in 1:nrow(na.ind)){
      colmean=mean(X[,na.ind[j,2]], na.rm=T)
      rowmean=mean(X[na.ind[j,1],], na.rm=T)
      impute.X[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X[na.ind]=impute.X

    na.ind=which(is.na(X),arr.ind = T)
    if (length(na.ind)==0) bool=F
  }

  return(list(X.fill=X, na.ind=na.ind.original))
}

estim_sigma <- function (X, k = NA, method = c("LN", "MAD"), center = "TRUE") {
  method <- match.arg(method, c("LN", "MAD", "ln", "mad", "Ln",
                                "Mad"), several.ok = T)[1]
  method <- tolower(method)
  if (inherits(X, "data.frame")) {
    X <- as.matrix(X)
  }
  if (sum(sapply(X, is.numeric)) < ncol(X)) {
    stop("all the variables must be numeric")
  }
  if (center == "TRUE") {
    X <- scale(X, scale = F)
  }
  n = nrow(X)
  p = ncol(X)
  svdX = svd(X)
  if (method == "ln" & is.na(k)) {
    warning("Since you did not specify k, k was estimated using the FactoMineR estim_ncp function")
    k <- estim_ncp(X, scale = F)$ncp
    print(paste("k = ", k))
  }
  if (center == "TRUE") {
    N <- (n - 1)
  }
  else {
    N <- n
  }
  if ((k >= min(N, p)) & (method == "ln")) {
    stop("the number k specified has to be smaller than the minimum of the number of rows or columns")
  }
  if (method == "ln") {
    if (k == 0) {
      sigma = sqrt(sum(svdX$d^2)/(N * p))
    }
    else {
      sigma <- sqrt(sum(svdX$d[-c(1:k)]^2)/(N * p - N *
                                              k - p * k + k^2))
    }
  }
  else {
    beta <- min(n, p)/max(n, p)
    lambdastar <- sqrt(2 * (beta + 1) + 8 * beta/((beta +
                                                     1 + (sqrt(beta^2 + 14 * beta + 1)))))
    wbstar <- 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta +
      1.43
    sigma <- median(svdX$d)/(sqrt(max(n, p)) * (lambdastar/wbstar))
  }
  return(sigma)
}
