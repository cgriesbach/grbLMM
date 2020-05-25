library(lme4)
library(parallel)
library(Matrix)

grbLMM = function(y, X, Z, id, m.stop = 500, ny = .1, cv.dat = NULL, aic = FALSE){
  getV = function(x){solve(x/sigma2 + Qi)}
  getQ = function(x, y){x + tcrossprod(as.numeric(y))}
  cfc = function(x)(length(unique(x)))
  getZ2D = function(x){chol2inv(chol(x + D))}
  
  ### basic definitions
  id = as.numeric(factor(id, levels = unique(id)))
  N = length(id)
  n = length(unique(id))
  id.t = as.numeric(table(id))
  first = rep(FALSE, N)
  for(i in 1:N){
    first[which.max(id==i)] = TRUE
  }
  p = ncol(X)
  q = ncol(Z)
  
  ### extract cluster-constant covariates
  ccc = rep(FALSE, p)
  for(r in 1:p){
    if(Reduce('+', lapply(split(X[,r], id), cfc)) == n){ccc[r] = TRUE}
  }
  Xcc = cbind(1, X[first, ccc])
  Xsw = X[,ccc]
  Xcor = list()
  Xcor[[1]] = Xcc%*%solve(crossprod(Xcc))%*%t(Xcc)
  
  if(q>1){
    for(s in 2:q){
      x = matrix(rep(1, n), n, 1)
      Xcor[[s]] = x%*%solve(crossprod(x))%*%t(x)
    }
  }
  
  p1 = rep(seq(1, q*n, q), q) + rep(0:(q-1), each = n)
  p2 = rep(seq(1, q*n, n), n) + rep(0:(n-1), each = q)
  P1 = sparseMatrix(seq_along(p1), p1)
  P2 = sparseMatrix(seq_along(p2), p2)
  
  Xcor = bdiag(Xcor)
  Xcor = P2%*%(diag(n*q) - Xcor)%*%P1
  
  ### construct random effects design matrix
  Z0 = Z
  Z = Z2 = QQ = list()
  for(i in 1:n){
    if(id.t[i] == 1){
      Z[[i]] = matrix(Z0[id==i,], 1, q)
    }else{
      Z[[i]] = Z0[id==i,]
    }
    Z2[[i]] = crossprod(Z[[i]])
  }
  bZ = Z
  Z = bdiag(Z)
  Z22 = crossprod(Z)
  
  ##### set starting values
  beta = rep(0, p)
  
  ### offset model for intercept and random structure
  if(q==1){
    offset = lme(y ~ 1, random = ~ 1 | id, control = lmeControl(opt = "optim", singular.ok = TRUE, returnObject = TRUE))
  }else{
    offset = lme(y ~ 1, random = ~ Z0[,-1] | id, control = lmeControl(opt = "optim", singular.ok = TRUE, returnObject = TRUE))
  }
  
  ### extract starting values
  int = offset$coefficients$fixed[1]
  gamma = Xcor%*%as.numeric(t(offset$coefficients$random$id))
  sigma2 = offset$sigma^2
  Q = getVarCov(offset)
  
  ### construct initial hat matrix
  C = cbind(1, Z)
  G = kronecker(diag(n), Q)
  B = as.matrix(bdiag(diag(1), solve(G)))
  S = list(C%*%bdiag(1,Xcor)%*%solve(crossprod(C) + sigma2*B)%*%t(C))
  
  ### in case of cv
  clcv = NA
  if(!is.null(cv.dat)){
    idcv = as.numeric(factor(cv.dat$idcv), levels = unique(cv.dat$idcv))
    Ncv = length(idcv)
    ncv = length(unique(idcv))
    idcv.t = as.numeric(table(idcv))
    
    Zcv = list()
    for(i in 1:ncv){
      if(idcv.t[i] == 1){
        z = cv.dat$Zcv[idcv==i,]
        Zcv[[i]] = t(as.matrix(z, 1, q))
      }else{
        Zcv[[i]] = cv.dat$Zcv[idcv==i,]
      }
    }
    
    Zcv = as.matrix(bdiag(Zcv))
  }
  
  ### prepare baselearners
  BL = HM = list()
  for(r in 1:p){
    x = cbind(1, X[,r])
    BL[[r]] = solve(crossprod(x))%*%t(x)
    HM[[r]] = ny*x%*%BL[[r]]
  }
  
  ### define storing matrices/vectors
  INT = int
  BETA = beta
  GAMMA = gamma
  SIGMA2 = sigma2
  CLCV = l = AIC = AICc = c()
  
  for(m in 1:m.stop){
    ###############################################################
    #### S1 #######################################################
    ###############################################################
    eta = as.vector(int + X%*%beta + Z%*%gamma)
    u = y - eta

    fits = matrix(0, 3, p)
    for(r in 1:p){
      fit = BL[[r]]%*%u
      fits[1,r] = fit[1]
      fits[2,r] = fit[2]
      fits[3,r] = sum((u - cbind(1, X[,r])%*%fit)^2)
    }

    best = which.min(fits[3,])
    int = int + ny*fits[1,best]
    beta[best] = beta[best] + ny*fits[2,best]


    ###############################################################
    #### S2 #######################################################
    ###############################################################
    eta = as.vector(int + X%*%beta + Z%*%gamma)
    u = y - eta

    D = solve(Q/sigma2)
    Z2D = lapply(Z2, getZ2D)
    rfit = ny*Xcor%*%bdiag(Z2D)%*%t(Z)
    gamma = gamma + rfit%*%u

    ###############################################################
    #### S3 #######################################################
    ###############################################################
    eta = as.vector(int + X%*%beta + Z%*%gamma)
    
    Qi = solve(Q)
    sigma2 = var(y - eta)

    V = lapply(Z2, getV)
    V = mapply(getQ, V, split(gamma, rep(1:n, each = q)), SIMPLIFY = FALSE)
    Q = Reduce("+", V)/n

    ### predictive risk computation
    if(!is.null(cv.dat)){
      Qcv = diag(Ncv) + Zcv%*%kronecker(diag(ncv), Q/sigma2)%*%t(Zcv)
      clcv = mean((cv.dat$ycv - int - cv.dat$Xcv%*%beta)^2)
    }

    if(aic){
      Sma = HM[[best]]
      Smb = Z%*%rfit
      S[[m+1]] = diag(N) - (diag(N) - Smb)%*%(diag(N) - Sma)%*%(diag(N) - S[[m]])
      df = sum(diag(S[[m+1]]))
      s2 = mean((y - eta)^2)
      AICc[m] = log(s2) + (1 + df/N) / (1 - (df + 2)/N)
    }

    INT = c(INT, int)
    BETA = rbind(BETA, beta)
    GAMMA = cbind(GAMMA, gamma)
    SIGMA2 = c(SIGMA2, sigma2)
    QQ[[m]] = Q
    CLCV = c(CLCV, clcv)

    if(m%%10 == 0){print(m)}
  }
  
  structure(list(int = int, beta = beta, sigma2 = sigma2, gamma = gamma, Q = Q, AICc = AICc,
                 INT = INT, BETA = BETA, SIGMA2 = SIGMA2, GAMMA = t(GAMMA), QQ = QQ, CLCV = CLCV,
                 S = S, eta = eta))
}

cv.grbLMM = function(k, y, X, Z, id, m.stop = 500, ny = .1, cores = 1){
  id = as.numeric(factor(id, levels = unique(id)))
  id.t = as.numeric(table(id))
  n = length(id.t)
  
  ### create an equally sized partition
  sets = sample(cut(1:n, k, labels = F))
  sets = rep(sets, id.t)
  
  ### function for each fold
  k.fold = function(k){
    y_train = y[sets != k]
    X_train = X[sets != k,]
    Z_train = as.matrix(Z[sets != k,])
    id_train = id[sets != k]
    
    cv.dat = list('ycv' = y[sets == k],
                  'Xcv' = X[sets == k,],
                  'Zcv' = as.matrix(Z[sets == k,]),
                  'idcv' = id[sets == k])
    
    model = grbLMM(y_train, X_train, Z_train, id_train, m.stop = m.stop, ny = .1, cv.dat = cv.dat)
    
    return(model$CLCV)
  }
  
  ### execute model on all folds
  if(cores==1){
    cv.ls = lapply(1:k, k.fold)
  }else{
    cv.ls = mclapply(1:k, k.fold, mc.cores = cores)
  }
  
  ### find best performing m.stop averaged over all folds
  cv.MAT = c()
  for(i in 1:k){
    cv.MAT = rbind(cv.MAT, cv.ls[[i]])
  }
  
  pred.risk = colMeans(cv.MAT)
  m.opt = which.min(pred.risk)
  
  model = grbLMM(y, X, Z, id, m.stop = m.opt)
  model$m.opt = m.opt
  model$pred = pred.risk
  model$coef = c(model$int, model$beta)
  model$folds = cv.MAT
  
  return(model)
}
