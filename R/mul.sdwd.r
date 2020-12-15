#' Multiway Sparse Distance Weighted Discrimination
#'
#' Optimize sparse DWD objective for a multiway dataset with any dimension. L1 and L2 parameters are imposed to enforce sparsity in the model.
#'
#' @param X A multiway array with dimensions \eqn{N \times P_1 \times ... \times P_K}. The first dimension (N) give the cases (e.g.subjects) to be classified.
#' @param y A vector of length N with class label (-1 or 1) for each case.
#' @param R Assumed rank of the coefficient array.
#' @param lambda1 The L1 tuning paprameter lambda1 that enforce sparsity. lambda1 and lambda2 can be determined by the function \code{\link{lambda.mul.sdwd}}.
#' @param lambda2 The L2 tuning paprameter lambda2.
#' @param convThresh The algorithm stops when the distance between B_new and B_old is less than convThresh. Default is 1e-7.
#' @param nmax Restrics how many iterations are allowed. Default is 500.
#' @return A list of components
#'
#' \code{beta} Vector of coefficients with length \eqn{P_1 \times ... \times P_K}.
#'
#' \code{int} Intercept.
#'
#' \code{U} A list of K matrices. Each matrix (\eqn{P_k \times R}) corresponds to the estimated weights for \eqn{k^th} dimension.
#'
#' If R=1, U is a list of K vectors.
#'
#'
#' @export
#'
#' @examples
#' ## Load gene expression time course data (?IFNB_Data for more info)
#' data(IFNB_Data)
#'
#' ## Run Multiway SDWD
#'res.msdwd1 <- mul.sdwd(DataArray,y=Class,R=1,lambda1=0,lambda2=1,convThresh = 1e-5, nmax=500)
#'
#' ##Compute projection onto the classification direction for each individual:
#' scores <- c()
#' for(i in 1:length(Class)) scores[i] = sum(as.vector(DataArray[i,,])*res.msdwd1$beta)+res.msdwd1$int
#' plot(scores, col=Class)
#'
#' @import sdwd
#' @import stats
#' @import roxygen2
#'
#'
#'

mul.sdwd <- function(X,y,R=1,lambda1=0,lambda2=1,convThresh = 1e-5, nmax=500) {
  if (R==1) {
    res <- msdwd.rank1(X,y,R,lambda1,lambda2,convThresh=convThresh,nmax=nmax)
  } else {
    res <- msdwd.rankR(X,y,R,lambda1,lambda2,thresh.outer=convThresh,nmax=nmax)
  }
  return(res)
}





################################################################
####### Functions for Rank 1 (R=1) Multiway SDWD model #########
################################################################
msdwd.rank1 <- function(X,y,R=1,lambda1=0,lambda2=1,convThresh=10^(-7),eps = 1e-8, nmax=500,num.init=5,nite=20){
  #This function implemented the algorithm for the rank 1 multiway sparse DWD model.
  #It is based on the objective function with separable penalty term (P^U(B)).
  #The rank is restricted to be 1 (R=1). Note that the funcion works for R>1, but the solution is often shrunk to a lower rank. Thus msdwd.rankR is used when R>1.  More details can be found in the paper.
  #X: multiway array
  #y: binary variable
  #R: specify the rank used in the model. R is restrited to be 1.
  # lambda1 and lambda2 are penalty parameters which can be determined by lambda.select.rank1.
  #nmax: Restrics how many iterations are allowed. Default is 500.
  #num.init: The number of initial values to start with. Default is 5.
  #nite: The number of iterations considered for initial function. Default is 20.
  #eps: The convergence thresholds for sdwd function at each dimension. Default is 1e-8
  #convThresh: The algorithm stops when the distance between B_new and B_old is less than convThresh. Default is 1e-4.

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
  if (lambda1>0){
    print('initializing without sparsity')
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
        if(all(U[[l]]==0)) {
          #U <- lapply(U, function(x) x=matrix(0,ncol=R))
          warning('All predictors are estimated to be zero; try a smaller lambda 1')
          break
        }
      }
      Us <- Uc <- matrix(ncol=L,nrow = R)
      score <- NA
      if(all(U[[l]]==0)) {
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
    print('proceeding with sparsity')
  }


  # Step 2 Run the algorithm with fixed penalty parameters
  j=0
  conv=convThresh+1
  obj.value2 <- NA
  s1<- s2 <- NA
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
      fit=sdwd(X_red,y,lambda=lambda1, lambda2 =lambda2,pf=pf1weights,pf2=pf2weights, standardize = FALSE, eps=eps, strong = FALSE)
      int=fit$b0
      Ulvec <- as.vector(fit$beta)
      U[[l]] = matrix(Ulvec,nrow=P[l])
      if(all(U[[l]]==0)) {
        warning('All predictors are estimated to be zero; try a smaller lambda 1 or lower Rank')
        break
      }
    }
    Us <- Uc <- matrix(ncol=L,nrow = R)
    score <- NA
    if(all(U[[l]]==0)) {
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
      obj.value2[j] <- obj.fun(X,y,R,beta,b0=int,U,lambda1, lambda2)
      conv=dist(rbind(as.vector(Bpre),as.vector(beta)))

      if(j > nmax){
        warning('The number of iterations achieved maximum')
        break
      }
    }
  }
  return(list('beta'=beta, 'int'=int,'U'=U))
}


## Initial step for rank 1 model
initial.step.rank1 <- function(X,y,R=1,nite=20){
  l1=0;l2=1
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  y = as.factor(y)
  U = list()
  for(l in 1:L) {
    u0 <- runif(P[l]*R)
    u <- u0/sqrt(sum(u0^2))
    U[[l]] = matrix(u,ncol=R)
  }
  Bmat = matrix(nrow=R,ncol=prod(P))
  for(r in 1:R){
    Ur = lapply(U, function(x) x[,r])
    Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
    Bmat[r,] = array(Br,dim=c(1,prod(P)))
  }
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
  j=0
  obj.value1 <- obj.value2 <- NA
  beta=Bmat
  s1 <- s2 <- NA
  ite=0
  conv=1
  while(j < nite & conv>1e-4){
    j=j+1
    Bpre = beta
    for(l in 1:L){
      ite=ite+1
      ### Matrix B
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
      fit=sdwd(X_red,y,lambda=l1, lambda2 =l2,pf=pf1weights,pf2=pf2weights, standardize = FALSE, strong = FALSE)
      int=fit$b0
      Ulvec <- as.vector(fit$beta)
      U[[l]] = matrix(Ulvec,nrow=P[l])

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }

      beta <- Bmat
      # objective fuction for each dimension
      obj.value1[ite] <- obj.fun(X,y,R,beta,b0=int,U,l1, l2)
      #cat("step 0",ite,obj.value1[ite],'\n')
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
    #beta <- colSums(Bmat)
    beta <- Bmat
    obj.value2[j] <- obj.fun(X,y,R,beta,b0=int,U,l1, l2)
    conv=dist(rbind(as.vector(Bpre),as.vector(beta)))
  }
  #cat("step 0",j,obj.value1[ite],obj.value2[j],'\n')
  return(list(U=U,int=int,B=Bmat,obj=obj.value1,obj2=obj.value2))
}

# Objective function for rank 1 model
obj.fun<- function(X,y,R=1,Bmat,b0,U,lambda1=0, lambda2=1) {
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  y <- as.numeric(levels(y))[y]
  Xmat <- array(X,dim=c(N,prod(P)))
  s1 <- s2 <- matrix(ncol=R,nrow=L);s11 <- s22 <- NA
  vt <- NA
  for(i in 1:N) {
    t <- y[i]*(b0+sum(Xmat[i,]%*%t(Bmat)))
    if(t <= 1/2) {
      vt[i] <- 1-t
    } else {
      vt[i] <- 1/(4*t)
    }
  }

  for(r in 1:R) {
    for(l in 1:L) {
      s1[l,r]<- sum(abs(U[[l]][,r]))
      s2[l,r] <- sum(U[[l]][,r]^2)
    }
    s11[r] <- prod(s1[,r])
    s22[r] <- prod(s2[,r])
  }
  penalty <- lambda1*sum(s11) + lambda2/2*sum(s22)
  obj <- 1/N * sum(vt) + penalty
  return(obj)
}



################################################################
####### Functions for Rank R (R>1) Multiway SDWD model #########
################################################################

msdwd.rankR=function(X,y,R=2,lambda1=0,lambda2=1,thresh.inner=10^(-5),thresh.outer=10^(-5),nmax=500,num.init=5,nite=10){
  #This function implemented the algorithm for the rank R multiway sparse DWD model.
  #It is based on the objective function with non-separable penalty term (P^(B)).
  #X: multiway array
  #y: binary variable
  #R: specify the rank used in the model.
  # lambda1 and lambda2 are penalty parameters which can be determined by lambda.select.rankR.
  #thresh.outer: The algorithm stops when the distance between B_new and B_old is less than thresh.outer. Default is 1e-5.
  #thresh.inner: This value is used for convergence of the inner algorithm that estimates weights for each dimension. (See Algorithm 2 in the paper )
  #nmax: Restrics how many iterations are allowed. Default is 500.
  #num.init: The number of initial values to start with. Default is 5.
  #nite: The number of iterations considered for initial function. Default is 20.

  l1=0;l2=1
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  #y = as.factor(y)
  res.init <- list()
  obj = c()
  # Initial step: start with multiple initial values
  print('initializing without sparsity')
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
  if(lambda1>0){
    print('initializing without sparsity')
    while(conv.outer>thresh.outer){
      beta.prev=beta
      #  beta0.prev=beta0
      jj=jj+1
      #  Bpre = beta
      for(l in 1:L){
        ite=ite+1
        ###Matrix B
        B_red = matrix(nrow=prod(P[-l]),ncol=R)
        for(r in 1:R){
          Ured_r <- lapply(U[-l], function(x) x[,r])
          B_red[,r] = as.vector(array(apply(expand.grid(Ured_r), 1, prod), dim=P[-l]))
          s1[r] <- do.call(prod,lapply(U[-l],function(x) sum(abs(x[,r]))))
          #        s2[r] <- do.call(prod,lapply(U[-l], function(x) sum(x[,r]^2)))
        }
        s2mat = crossprod(B_red)
        X_red = array(dim=c(N,P[l],R))
        for(i in 1:N){X_red[i,,] = Xarrays[[l]][i,,] %*% B_red}

        conv.inner = 1
        while(conv.inner> thresh.inner){
          beta0.prev=beta0
          Uprev = U[[l]]
          #update coefficients
          for(i in 1:R){ for(j in 1:P[l]){
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
          #  print(conv.inner)
        }

        Bmat = matrix(nrow=R,ncol=prod(P))
        for(r in 1:R) {
          Ur = lapply(U, function(x) x[,r])
          Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
          Bmat[r,] = array(Br,dim=c(1,prod(P)))
        }
        beta <- colSums(Bmat)
        if(any(colSums(abs(U[[l]]))==0)) {
          warning('All predictors are estimated to be zero for one component; try a smaller lambda1 or lower Rank')
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
    print('proceeding with sparsity')
  }

  jj=0;conv.outer=1;ite=0
  # Step 2 Run the algorithm with fixed penalty parameters
  while(jj < nmax & conv.outer>thresh.outer){
    beta.prev=beta
    #  beta0.prev=beta0
    jj=jj+1
    #  Bpre = beta
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
        }}
        #update intercept
        vprime.u = vprime(u)
        temp = beta0-mean(vprime.u*y)/4
        u=u+y*(temp-beta0)
        beta0=temp
        conv.inner = sum((U[[l]]-Uprev)^2)+(beta0-beta0.prev)^2
        #  print(conv.inner)
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }
      beta <- colSums(Bmat)
      if(any(colSums(abs(U[[l]]))==0)) {
        warning('All predictors are estimated to be zero for one component; try a smaller lambda1 or lower Rank')
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
  return(list('beta'=beta, 'int'=beta0,'U'=U))
}



####### Intial step for rank R model  #########
initial.step.rankR=function(X,y,R=1,thresh.inner=10^(-5),thresh.outer=10^(-5),nite=10){
  l1 = 0; l2 = 1
  L = length(dim(X))-1
  N = dim(X)[1]
  P = dim(X)[2:(L+1)]
  Xmat = array(X,dim=c(N,prod(P)))
  #y = as.factor(y)
  U = list()
  for(l in 1:L) {
    u0 <- runif(P[l]*R)
    u <- u0/sqrt(sum(u0^2))
    U[[l]] = matrix(u,ncol=R)
  }
  beta0=0
  Bmat = matrix(nrow=R,ncol=prod(P))
  for(r in 1:R){
    Ur = lapply(U, function(x) x[,r])
    Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
    Bmat[r,] = array(Br,dim=c(1,prod(P)))
  }
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
  jj=0
  obj.value1 <- obj.value2 <- NA
  beta=colSums(Bmat)
  s1 <- s2 <- NA
  ite=0
  conv.outer=1
  int.vec=c()
  u=y*(beta0+Xmat%*%beta)
  while(jj < nite & conv.outer>thresh.outer){
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
        #  print(conv.inner)
      }

      Bmat = matrix(nrow=R,ncol=prod(P))
      for(r in 1:R) {
        Ur = lapply(U, function(x) x[,r])
        Br = array(apply(expand.grid(Ur), 1, prod), dim=P)
        Bmat[r,] = array(Br,dim=c(1,prod(P)))
      }
      beta <- colSums(Bmat)
      if(any(colSums(abs(U[[l]]))==0)) {
        warning('All predictors are estimated to be zero for one component; try a smaller lambda1 or lower Rank')
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
  return(list('beta'=beta,'Bmat'=Bmat,'U'=U, 'beta0'=beta0,'obj.vec'=obj.value2))
}

###### Other internal functions #######
vprime=function(u){
  vpu = rep(-1,length(u))
  vpu[u>0.5] = -1/(4*u[u>0.5]^2)
  return(vpu)
}

vfunc =function(u){
  vu = 1-u
  vu[u>0.5]=1/(4*u[u>0.5])
  return(vu)
}


obj.fun1<- function(Xm,y,beta,Bmat,beta0,lambda1,lambda2=1) {
  return(mean(vfunc(y*(beta0+Xm%*%beta)))+lambda1*sum(abs(Bmat))+lambda2/2*sum(Bmat^2))
}

obj.fun2<- function(Xm,y,beta,Bmat,beta0,lambda1,lambda2=1) {
  return(mean(vfunc(y*(beta0+Xm%*%beta)))+lambda1*sum(abs(Bmat))+lambda2/2*sum(beta^2))
}



