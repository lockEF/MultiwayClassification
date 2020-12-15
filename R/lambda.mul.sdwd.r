#' Penalty parameters selection for the Multiway Sparse DWD
#'
#' Conduct a k-fold cross-validation for \code{\link{mul.sdwd}} and returns the optimal pair of L1 and L2 parameters.
#'
#' @param X A multiway array with dimensions \eqn{N \times P_1 \times ... \times P_K}. The first dimension (N) give the cases (e.g.subjects) to be classified.
#' @param y A vector of length N with class label (-1 or 1) for each case.
#' @param R Assumed rank of the coefficient array.
#' @param lambda1.vec A vector of L1 candidates. Default is c(0, 1e-4, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25).
#' @param lambda2.vec A vector of L2 candidates. Default is c(0.25, 0.50, 0.75, 1.00, 3, 5).
#' @param nfolds The number of folds. Default value is 10.
#' @param convThresh The algorithm stops when the distance between B_new and B_old is less than convThresh. Default is 1e-7.
#' @param nmax Restrics how many iterations are allowed. Default is 500.
#' @return A list of components
#'
#' \code{par} Optimal pair of L1 and L2 and the maximum t test statistic.
#'
#' \code{tstats} Test statistics for all L1 and L2 condidates.
#'
#' @export
#'
#' @examples
#' ## Load gene expression time course data (?IFNB_Data for more info)
#' data(IFNB_Data)
#'
#'## Select penalty parameters by cross-validation
#' lambda.msdwd1 <- lambda.mul.sdwd(DataArray,y=Class,R=1,lambda1.vec=c(0, 1e-4, 0.001, 0.005, 0.01),
#' lambda2.vec=c(0.50, 1.00, 3, 5), nfolds = 10,convThresh=10^(-5),nmax=500)
#' lambda.msdwd1$par
#'
#' @import sdwd
#' @import stats
#' @import roxygen2
#'


lambda.mul.sdwd <- function(X,y,R=1,lambda1.vec=c(0, 1e-4, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25), lambda2.vec=c(0.25, 0.50, 0.75, 1.00, 3, 5), nfolds = 10,convThresh=10^(-7),nmax=500) {
  if (R==1) {
    lambda <- lambda.select.rank1(X,y,R,lambda1.vec,lambda2.vec,nfolds,convThresh=convThresh,nmax=nmax)
  } else {
    lambda <- lambda.select.rankR(X,y,R,lambda1.vec,lambda2.vec,nfolds,thresh.outer=convThresh,nmax=nmax)
  }
  return(lambda)
}

################################################################
####### Functions for Rank 1 (R=1) Multiway SDWD model #########
################################################################

# Select L1 and L2 using cross validation for rank 1 model
lambda.select.rank1 <- function(X,y,R=1,lambda1.vec, lambda2.vec, nfolds = 10,convThresh=1e-7,nmax=500) {
  res <- mapply(function(iter) lambda1.select.rank1(iter, X,y,R,lambda1.vec, lambda2.vec, nfolds,convThresh,nmax),iter=1:length(lambda2.vec))
  nd <- length(res)
  out.l2 <- do.call(rbind,res[seq(1,nd,2)])
  out.all <- do.call(rbind,res[seq(2,nd,2)])
  par <- out.l2[which.max(out.l2[,3]),]
  return(list('par'=par, 'tstats' = out.all))
}


# Select L1 for fixed L2. This function is used in lambda1.select.rank1
lambda1.select.rank1 <- function(iter=1, X,y,R=1,lambda1.vec=c(0, 1e-4, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),
                                 lambda2.vec=c(0.25, 0.50, 0.75, 1.00, 3, 5), nfolds = 10, convThresh=10^(-7),nmax=500) {
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  r1 <- length(lambda1.vec)
  scores.cv<- matrix(0,nrow=N, ncol=r1)
  foldid = sample(rep(seq(nfolds), length = N))
  for (i in seq(nfolds)) {
    which = foldid ==i
    Xmati <- Xmat[!which,]
    Xi <- array(Xmati,dim=c(dim(Xmati)[1],P))
    Yi <- y[!which]
    t.stats <- matrix(NA, 1, r1)
    res = try(msdwd.rank1.l1(Xi, Yi, R,lambda1.vec,lambda2.vec[iter],convThresh,nmax=nmax))
    for(ii in which(which==TRUE)) {
      for(m in 1:r1) {
        scores.cv[ii,m] <- as.vector(rowSums(Xmat[ii,]%*%t(res$beta.list[[m]])))
      }
    }
  }

  for(i in 1:r1){
    ttest <- try(t.test(scores.cv[,i]~y))
    if(class(ttest)=="try-error") {
      t.stats[1,i] <- 0
    } else {
      t.stats[1,i] <- abs(t.test(scores.cv[,i]~y)$statistic)
    }
  }
  ind.min <- which(t.stats == max(t.stats, na.rm = TRUE), arr.ind = TRUE)
  par <- c('lambda1'=lambda1.vec[ind.min[2]],'lambda2'=lambda2.vec[iter], 'tstats'=max(t.stats, na.rm = TRUE))
  out <- as.data.frame(matrix(0,length(lambda1.vec),3))
  out[,1] <- lambda1.vec;out[,2] <- lambda2.vec[iter];out[,3] <- t(t.stats)
  names(out) <- c("lambda1","lambda2","tstats")
  return(list(par=par, t.stats=out) )
}

## Run msdwd for a series of L1 parameters. The results are used in lambda1.select.rank1.
msdwd.rank1.l1 <- function(X,y,R=1,lambda1.vec=c(0,0.01),lambda2=1,convThresh=10^(-7),nmax=500,num.init=5,nite=20,eps = 1e-8){
  n.l1 <- length(lambda1.vec)
  res=list()
  obj = c()
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  y = as.factor(y)
  Xarrays = list()
  for(l in 1:L){
    Xarrays[[l]] <- array(dim=c(N,P[l],prod(P[-l])))
    perm = c(l,c(1:L)[-l])
    for(i in 1:N){
      X_slice = array(Xmat[i,],dim=c(P))
      X_slice_perm = aperm(X_slice,perm)
      Xarrays[[l]][i,,] = array(X_slice_perm,dim=c(P[l],prod(P[-l])))
    }
  }

  # Initial step: start with multiple initial values
  for (i in 1:num.init) {
    res[[i]] <-initial.step.rank1(X,y,R,nite)
    obj[i] <- res[[i]]$obj[length(res[[i]]$obj)]
  }

  ii <- which.min(obj)
  U <- res[[ii]]$U
  j=0;
  conv=convThresh+1
  Bmat = res[[ii]]$B
  beta <- Bmat
  s1 <- s2 <- NA;obj.value1 <- NA
  l1 <- 0;l2 <- 1;

  # Step 1 continue with the selected path
  while(conv>convThresh) {
    j=j+1
    Bpre = beta
    for(l in 1:L){
      ###Matrix B
      B_red = matrix(nrow=prod(P[-l]),ncol=R)
      for(r in 1:R){
        Ured_r <- lapply(U[-l], function(x) x[,r])
        B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l]))
        s1[r] <- do.call(prod,lapply(U[-l],function(x) sum(abs(x[,r]))))
        s2[r] <- do.call(prod,lapply(U[-l], function(x) sum(x[,r]^2)))
      }

      X_red = matrix(nrow=N,ncol=P[l]*r)
      for(i in 1:N){X_red[i,] = as.vector(Xarrays[[l]][i,,] %*% B_red)}

      pf1weights = pf2weights = c()
      for(r in 1:R) pf2weights = c(pf2weights,rep(s2[r],P[l]))
      for(r in 1:R) pf1weights = c(pf1weights,rep(s1[r],P[l]))
      fit=sdwd(X_red,y,lambda=l1, lambda2 =l2,pf=pf1weights,pf2=pf2weights, standardize = FALSE, eps=eps, strong = FALSE)
      int=fit$b0
      Ulvec <- as.vector(fit$beta)
      U[[l]] = matrix(Ulvec,nrow=P[l])
      if(any(colSums(abs(U[[l]]))==0) | is.na(sum(U[[l]]))) {
        warning('All predictors are estimated to be zero; try a smaller lambda 1 or lower Rank')
        U[[l]][is.na(U[[l]])]=0
        int=0
        return(list(U.list=U.list,int.mat=int.mat, beta.list=beta.list))
      }
    }
    Us <- Uc <- matrix(ncol=L,nrow = R)
    score <- NA
    if(all(U[[l]]==0) | all(is.na(U[[l]]))) {
      break
    } else {
      for(r in 1:R) {
        Us[r,] <- sapply(U, function(x) sum(x[,r]^2))
        score[r] <- prod(Us[r,])^(1/L)
        Uc[r,] <- sqrt(score[r]/Us[r,])
      }

      for (l in 1:L) {
        U[[l]] <- matrix(t(apply(U[[l]], 1, function(x) Uc[,l]*x)),ncol=R)
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }
      beta <- Bmat
      obj.value1[j] <- obj.fun(X,y,R,beta,b0=int,U,l1, l2)
      conv=dist(rbind(as.vector(Bpre),as.vector(beta)))
    }

  }

  U.list <- beta.list <- vector("list",n.l1)
  int.mat <- matrix(NA, n.l1, 1)
  # Step 2 Run the algorithm with fixed penalty parameters
  for(kk in 1:n.l1) {
    j=0
    conv=convThresh+1
    obj.value2 <- NA
    s1<- s2 <- NA
    if(sum(sapply(U,function(x) colSums(x))==0,na.rm = T)>0) {
      for(kk1 in kk:n.l1) {
        for(l in 1:L) {
          u <- rep(0,P[l]*R)
          U[[l]] = matrix(u,ncol=R)
        }
        U.list[[kk1]] <- U
        int.mat[kk1,1] <- 0
        beta.list[[kk1]] <- matrix(0,R,prod(P))

      }
      break
    } else {
      while(conv>convThresh){
        j=j+1
        Bpre = beta
        for(l in 1:L){
          ###Matrix B
          B_red = matrix(nrow=prod(P[-l]),ncol=R)
          for(r in 1:R){
            Ured_r <- lapply(U[-l], function(x) x[,r])
            B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l]))
            s1[r] <- do.call(prod,lapply(U[-l],function(x) sum(abs(x[,r]))))
            s2[r] <- do.call(prod,lapply(U[-l], function(x) sum(x[,r]^2)))
          }

          X_red = matrix(nrow=N,ncol=P[l]*r)
          for(i in 1:N){X_red[i,] = as.vector(Xarrays[[l]][i,,] %*% B_red)}
          pf1weights = pf2weights = c()
          for(r in 1:R) pf2weights = c(pf2weights,rep(s2[r],P[l]))
          for(r in 1:R) pf1weights = c(pf1weights,rep(s1[r],P[l]))
          fit=sdwd(X_red,y,lambda=lambda1.vec[kk], lambda2 =lambda2,pf=pf1weights,pf2=pf2weights, standardize = FALSE, eps=eps, strong = FALSE)
          int=fit$b0
          Ulvec <- as.vector(fit$beta)
          U[[l]] = matrix(Ulvec,nrow=P[l])
          if(any(colSums(abs(U[[l]]))==0) | is.na(sum(U[[l]]))) {
            warning('All predictors are estimated to be zero; try a smaller lambda 1 or lower Rank')
            U[[l]][is.na(U[[l]])]=0
            int=0
            return(list(U.list=U.list,int.mat=int.mat, beta.list=beta.list))
          }
        }

        Us <- Uc <- matrix(ncol=L,nrow = R)
        score <- NA
        if(sum(sapply(U,function(x) colSums(x))==0,na.rm = T)>0){
          Bmat = matrix(nrow=R,ncol=prod(P))
          for(r in 1:R) {
            Ur = lapply(U, function(x) x[,r])
            Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
            Bmat[r,] = array(Br,dim=c(1,prod(P)))
          }
          beta <- Bmat
          break} else {
            for(r in 1:R) {
              Us[r,] <- sapply(U, function(x) sum(x[,r]^2))
              score[r] <- prod(Us[r,])^(1/L)
              Uc[r,] <- sqrt(score[r]/Us[r,])
            }

            for (l in 1:L) {
              U[[l]] <- matrix(t(apply(U[[l]], 1, function(x) Uc[,l]*x)),ncol=R)
            }

            Bmat = matrix(nrow=R,ncol=prod(P))
            for(r in 1:R) {
              Ur = lapply(U, function(x) x[,r])
              Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
              Bmat[r,] = array(Br,dim=c(1,prod(P)))
            }
            beta <- Bmat
            conv=dist(rbind(as.vector(Bpre),as.vector(beta)))

            if(j > nmax){
              warning('The number of iterations achieved maximum')
              break
            }
          }
      }
      U.list[[kk]] <- U
      int.mat[kk,1] <- int
      beta.list[[kk]] <- beta
    }
  }
  return(list(U.list=U.list,int.mat=int.mat, beta.list=beta.list))
}


################################################################
####### Functions for Rank R (R>1) Multiway SDWD model #########
################################################################


####### Select L1 and L2 using cross validation for rank R model  #######
lambda.select.rankR <- function(X,y,R=2,lambda1.vec=c(0, 1e-4, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),
                                lambda2.vec=c(0.25, 0.50, 0.75, 1.00, 3, 5), nfolds = 10,thresh.inner=10^(-5),thresh.outer=10^(-5),nmax=500,num.init=5,nite=10) {
  res <- mapply(function(iter) lambda1.select.rankR(iter, X,y,R,lambda1.vec, lambda2.vec, nfolds,thresh.inner,thresh.outer,nmax,num.init,nite),iter=1:length(lambda2.vec))
  nd <- length(res)
  out.l2 <- do.call(rbind,res[seq(1,nd,2)])
  out.all <- do.call(rbind,res[seq(2,nd,2)])
  par <- out.l2[which.max(out.l2[,3]),]
  return(list('par'=par, 'tstats' = out.all))
}

####### Select L1 for fixed L2. This function is used in lambda.select.rankR
lambda1.select.rankR <- function(iter=1, X,y,R=2,lambda1.vec, lambda2.vec, nfolds = 10,thresh.inner=10^(-5),thresh.outer=10^(-5),nmax=500,num.init=5,nite=10) {
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  #y = as.factor(y)
  Xmat = array(X,dim=c(N,prod(P)))
  r1 <- length(lambda1.vec)
  scores.cv<- matrix(0,nrow=N, ncol=r1)
  foldid = sample(rep(seq(nfolds), length = N))
  for (i in seq(nfolds)) {
    which = foldid ==i
    Xmati <- Xmat[!which,]
    Xi <- array(Xmati,dim=c(dim(Xmati)[1],P))
    Yi <- y[!which]
    t.stats <- matrix(NA, 1, r1)
    res = msdwd.rankR.l1(Xi, Yi, R,lambda1.vec,lambda2.vec[iter],thresh.inner,thresh.outer,nmax,num.init,nite)
    for(ii in which(which==TRUE)) {
      for(m in 1:r1) {
        scores.cv[ii,m] <- sum(Xmat[ii,]%*%res$beta.list[[m]])
      }
    }
  }

  for(i in 1:r1){
    ttest <- try(t.test(scores.cv[,i]~y))
    if(class(ttest)=="try-error" | is.na(ttest$statistic)) {
      t.stats[1,i] <- 0
    } else {
      t.stats[1,i] <- abs(ttest$statistic)
    }
  }
  ind.min <- which(t.stats == max(t.stats, na.rm = TRUE), arr.ind = TRUE)
  par <- c('lambda1'=lambda1.vec[ind.min[2]],'lambda2'=lambda2.vec[iter], 'tstats'=max(t.stats, na.rm = TRUE))
  out <- as.data.frame(matrix(0,length(lambda1.vec),3))
  out[,1] <- lambda1.vec;out[,2] <- lambda2.vec[iter];out[,3] <- t(t.stats)
  names(out) <- c("lambda1","lambda2","tstats")
  return(list('par.L1'=par, 'out.L1'=out))
}

####### Multiway sparse DWD model for a series of L1 parameters #######
msdwd.rankR.l1 <- function(X,y,R=2,lambda1.vec=c(0,0.01),lambda2=1,thresh.inner=10^(-5),thresh.outer=10^(-5),nmax=500,num.init=5,nite=10){
  n.l1 <- length(lambda1.vec)
  l1=0;l2=1
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  #y = as.factor(y)
  res.init <- list()
  obj = c()
  # Initial step: start with multiple initial values
  for (i in 1:num.init) {
    res.init[[i]] <-initial.step.rankR(X,y,R,nite=nite)
    obj[i] <- res.init[[i]]$obj.vec[length(res.init[[i]]$obj.vec)]
  }

  ii <- which.min(obj)
  U <- res.init[[ii]]$U
  beta0=res.init[[ii]]$beta0
  Bmat = res.init[[ii]]$Bmat
  beta = res.init[[ii]]$beta
  jj=0; ite=0
  conv.outer=1
  u=y*(beta0+Xmat%*%beta)
  obj.value1 <- obj.value2 <- NA
  s1 <- s2 <- NA
  int.vec=c()

  Xarrays = list()
  for(l in 1:L){
    Xarrays[[l]] <- array(dim=c(N,P[l],prod(P[-l])))
    perm = c(l,c(1:L)[-l])
    for(i in 1:N){
      X_slice = array(Xmat[i,],dim=c(P))
      X_slice_perm = aperm(X_slice,perm)
      Xarrays[[l]][i,,] = array(X_slice_perm,dim=c(P[l],prod(P[-l])))
    }
  }

  # Step 1 continue with the selected path
  while(conv.outer>thresh.outer){
    beta.prev=beta
    jj=jj+1
    for(l in 1:L){
      ite=ite+1
      ###Matrix B
      B_red = matrix(nrow=prod(P[-l]),ncol=R)
      for(r in 1:R){
        Ured_r <- lapply(U[-l], function(x) x[,r])
        B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l]))
        s1[r] <- do.call(prod,lapply(U[-l],function(x) sum(abs(x[,r]))))
      }
      s2mat = crossprod(B_red)
      X_red = array(dim=c(N,P[l],r))
      for(i in 1:N){X_red[i,,] = Xarrays[[l]][i,,] %*% B_red}

      conv.inner = 1
      while(conv.inner> thresh.inner){
        beta0.prev=beta0
        Uprev = U[[l]]
        #update coefficients
        for(i in 1:r){ for(j in 1:P[l]){
          vprime.u = vprime(u)
          z=4*U[[l]][j,i]-mean(vprime.u * X_red[,j,i]*y)-l2*sum(s2mat[i,-i]*U[[l]][j,-i])
          temp=sign(z)*max(abs(z)-s1[i]*l1,0)/(4+l2*s2mat[i,i])
          u = u+y*(temp-U[[l]][j,i])*X_red[,j,i]
          U[[l]][j,i]=temp
        }}
        #update intercept
        vprime.u = vprime(u)
        temp = beta0-mean(vprime.u*y)/4
        u=u+y*(temp-beta0)
        beta0=temp
        conv.inner = sum((U[[l]]-Uprev)^2)+(beta0-beta0.prev)^2
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }
      beta <- colSums(Bmat)
      if(any(colSums(abs(U[[l]]))==0) | is.na(sum(U[[l]]))) {
        warning('All predictors are estimated to be zero for one component; try a smaller lambda1 or lower Rank')
        U[[l]][is.na(U[[l]])]=0
        beta0=0
        return(list('beta'=beta,'Bmat'=Bmat,'U'=U, 'beta0'=beta0,'obj.vec'=obj.value2))
      }
      Us <- Uc <- matrix(ncol=L,nrow = R)
      score <- NA
      for(r in 1:R) {
        Us[r,] <- sapply(U, function(x) sum(x[,r]^2))
        score[r] <- prod(Us[r,])^(1/L)
        Uc[r,] <- sqrt(score[r]/Us[r,])
      }

      for (l in 1:L) {
        U[[l]] <- matrix(t(apply(U[[l]], 1, function(x) Uc[,l]*x)),ncol=R)
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }

      beta <- colSums(Bmat)
      conv.outer = sum((beta-beta.prev)^2)+(beta0-beta0.prev)^2
      # objective fuction for each dimension
      obj.value1[ite] <- obj.fun1(Xmat,y,beta,Bmat,beta0=beta0,l1, l2)
      obj.value2[ite] <- obj.fun2(Xmat,y,beta,Bmat,beta0=beta0,l1, l2)
    }
  }

  U.list <- beta.list <- vector("list",n.l1)
  int.mat <- matrix(NA, n.l1, 1)

  # Step 2 Run the algorithm with fixed penalty parameters
  for(kk in 1:n.l1) {
    t1 <- mdwd.t1(U,beta,beta0,Bmat,X,y,R,lambda1=lambda1.vec[kk],lambda2)
    U=t1$U;beta=t1$beta;Bmat=t1$Bmat;beta0=t1$beta0
    U.list[[kk]] <- U
    int.mat[kk,1] <- beta0
    beta.list[[kk]] <- beta

    if(sum(sapply(U,function(x) colSums(x))==0,na.rm = T)>0) {
      for(kk1 in kk:n.l1) {
        for(l in 1:L) {
          u <- rep(0,P[l]*R)
          U[[l]] = matrix(u,ncol=R)
        }
        U.list[[kk1]] <- U
        int.mat[kk1,1] <- 0
        beta.list[[kk1]] <- rep(0,prod(P))

      }
      return(list(U.list=U.list,int.mat=int.mat, beta.list=beta.list))
    }

  }
  return(list(U.list=U.list,int.mat=int.mat, beta.list=beta.list))
}

mdwd.t1 <- function(U,beta,beta0,Bmat,X,y,R,lambda1,lambda2,thresh.inner=10^(-5),thresh.outer=10^(-5),nmax=500,num.int=5,nite=10) {
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  #y = as.factor(y)
  u=y*(beta0+Xmat%*%beta)
  obj.value1 <- obj.value2 <- NA
  s1 <- s2 <- NA
  int.vec=c()

  Xarrays = list()
  for(l in 1:L){
    Xarrays[[l]] <- array(dim=c(N,P[l],prod(P[-l])))
    perm = c(l,c(1:L)[-l])
    for(i in 1:N){
      X_slice = array(Xmat[i,],dim=c(P))
      X_slice_perm = aperm(X_slice,perm)
      Xarrays[[l]][i,,] = array(X_slice_perm,dim=c(P[l],prod(P[-l])))
    }
  }
  jj=0; ite=0;conv.outer=1;
  while(jj < nmax & conv.outer>thresh.outer){
    beta.prev=beta
    jj=jj+1
    for(l in 1:L){
      ite=ite+1
      ###Matrix B
      B_red = matrix(nrow=prod(P[-l]),ncol=R)
      for(r in 1:R){
        Ured_r <- lapply(U[-l], function(x) x[,r])
        B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l]))
        s1[r] <- do.call(prod,lapply(U[-l],function(x) sum(abs(x[,r]))))
      }
      s2mat = crossprod(B_red)
      X_red = array(dim=c(N,P[l],r))
      for(i in 1:N){X_red[i,,] = Xarrays[[l]][i,,] %*% B_red}

      conv.inner = 1
      while(conv.inner> thresh.inner){
        beta0.prev=beta0
        Uprev = U[[l]]
        #update coefficients
        for(i in 1:r){ for(j in 1:P[l]){
          vprime.u = vprime(u)
          z=4*U[[l]][j,i]-mean(vprime.u * X_red[,j,i]*y)-lambda2*sum(s2mat[i,-i]*U[[l]][j,-i])
          temp=sign(z)*max(abs(z)-s1[i]*lambda1,0)/(4+lambda2*s2mat[i,i])
          u = u+y*(temp-U[[l]][j,i])*X_red[,j,i]
          U[[l]][j,i]=temp
        }
        }

        #update intercept
        vprime.u = vprime(u)
        temp = beta0-mean(vprime.u*y)/4
        u=u+y*(temp-beta0)
        beta0=temp
        conv.inner = sum((U[[l]]-Uprev)^2)+(beta0-beta0.prev)^2
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }
      beta <- colSums(Bmat)
      if(any(colSums(abs(U[[l]]))==0) | is.na(sum(U[[l]])) ) {
        #warning('All predictors are estimated to be zero for one component; try a smaller lambda1 or lower Rank')
        U[[l]][is.na(U[[l]])]=0
        return(list('beta'=beta,'Bmat'=Bmat,'U'=U, 'beta0'=beta0,'obj.vec'=obj.value2))
      }
      Us <- Uc <- matrix(ncol=L,nrow = R)
      score <- NA
      for(r in 1:R) {
        Us[r,] <- sapply(U, function(x) sum(x[,r]^2))
        score[r] <- prod(Us[r,])^(1/L)
        Uc[r,] <- sqrt(score[r]/Us[r,])
      }

      for (l in 1:L) {
        U[[l]] <- matrix(t(apply(U[[l]], 1, function(x) Uc[,l]*x)),ncol=R)
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }

      beta <- colSums(Bmat)
      conv.outer = sum((beta-beta.prev)^2)+(beta0-beta0.prev)^2
      # objective fuction for each dimension
      obj.value1[ite] <- obj.fun1(Xmat,y,beta,Bmat,beta0=beta0,lambda1, lambda2)
      obj.value2[ite] <- obj.fun2(Xmat,y,beta,Bmat,beta0=beta0,lambda1, lambda2)
    }
  }
  return(list('beta'=beta,'Bmat'=Bmat,'U'=U, 'beta0'=beta0,'obj.vec'=obj.value2))
}
