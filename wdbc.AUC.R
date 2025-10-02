library(pROC)
library(rocbc)
library(dplyr)
library(MASS)
library(emplik)
#library(abind)

#Box-Cox Transformation and Standardization
bc.scale = function(x) {
  if (any(is.na(x))) stop("Input contains NA values.")
  if (any(x <= 0)) stop("Input not positive.")
  bc <- MASS::boxcox(x ~ 1, plotit = FALSE)
  # The optimal lambda is the value on the x-axis corresponding to the peak
  # of the log-likelihood curve on the y-axis.
  lambda_optimal <- bc$x[which.max(bc$y)]
  
  # --- Apply Transformation ---
  # A lambda of 0 corresponds to a natural log transformation.
  if (lambda_optimal == 0) {
    transformed_x <- log(x)
  } else {
    transformed_x <- (x^lambda_optimal - 1) / lambda_optimal
  }
  return(scale(transformed_x))
}

load('wdbc.new')
#Apply Box-Cox transformation and then scale the data
wdbc.bc = as.data.frame(apply(wdbc[,3:32],2,bc.scale)) 
wdbc.bc$D = (wdbc$V2 == 'M')

#Find True AUC and select biomarkers
wdbc_auc = NULL
for (ii in 1:30) {
  auc = auc(roc(response = wdbc.bc$D, predictor = wdbc.bc[,ii]))
  wdbc_auc = c(wdbc_auc,auc)
}

#Necessary parameters
data.vb = wdbc.bc
#colSums(is.na(data.vb))
true_auc = wdbc_auc
n = nrow(data.vb)
precision = 10 ^ (-3)
alpha = 0.05
za = qnorm(0.975)
B = 1000
search_start = 0
search_end = 1
search_step = 0.01
tol = 1e-5

auc_vb_ci = function(bmk,
                     k = 1, a = 6,
                     misv = FALSE) {
  D = data.vb$D
  T = as.numeric(data.vb[, bmk])
  Tc = as.numeric(wdbc[, bmk + 2])
  A = as.numeric(data.vb[, a]) #Covariate
  C = A
  #generate V
  logitp = T + C + k
  pv = 1 / (exp(-logitp) + 1)
  u = runif(n)
  V = (pv > u) + 0
  #Apply rocbc & pROC package to compute CI ignoring vb
  Tv = Tc[V == 1]
  Dv = D[V == 1]
  BC_ci = rocboxcox(Tv, Dv, alpha = 0.05, plots = 'on')$AUCCI
  DeLong_ci = as.numeric(ci.auc(roc(Dv, Tv)))[c(1, 3)]
  get_weights = function(T, D, A, V, n, misv = FALSE) {
    #Fit Disease model
    disease = glm(D[V == 1] ~ T[V == 1] + A[V == 1],
                  family = binomial(link = 'probit'),
                  control = list(maxit = 100))
    beta = disease$coefficients
    rho = drop(pnorm(cbind(rep(1, n), T, A) %*% beta))
    rho[rho > 1 - precision] = 1 - precision
    rho[rho < precision] = precision
    if (misv == FALSE) {
      # Correct V model
      Verif = glm(V ~ T + C,
                  family = binomial(link = 'logit'),
                  control = list(maxit = 100))
      betaV = Verif$coefficients
      pi_hat = (1 + exp(-drop(cbind(rep(
        1, n
      ), T, C) %*% betaV))) ^ (-1)
    }
    else {
      #Wrong V model
      Verif = glm(V ~ T,
                  family = binomial(link = 'logit'),
                  control = list(maxit = 100))
      betaV = Verif$coefficients
      pi_hat = (1 + exp(-drop(cbind(rep(
        1, n
      ), T) %*% betaV))) ^ (-1)
    }
    pi_hat[pi_hat > 1 - precision] = 1 - precision
    pi_hat[pi_hat < precision] = precision
    #weights
    wd.fi = rho
    wd.msi = V * D + (1 - V) * rho
    wd.ipw = V * D / pi_hat
    wd.spe = (V * D - (V - pi_hat) * rho) / pi_hat
    wn.fi = 1 - rho
    wn.msi = V * (1 - D) + (1 - V) * (1 - rho)
    wn.ipw = V * (1 - D) / pi_hat
    wn.spe = (V * (1 - D) - (V - pi_hat) * (1 - rho)) / pi_hat
    return(list(
      'wd' = rbind(wd.fi, wd.msi, wd.ipw, wd.spe),
      'wn' = rbind(wn.fi, wn.msi, wn.ipw, wn.spe)
    ))
  }
  
  BSP = function(dummy) {
    boots = sample(1:n, n, replace = TRUE)
    Tb = T[boots]
    Db = D[boots]
    Ab = A[boots]
    Vb = V[boots]
    VD = sum(Vb * Db)
    VH = sum(Vb * (1 - Db))
    if (VD >= 1 && VH >= 1)
    {
      #get_weights(Tb, Db, Ab, Vb, n, misd = misd, misv = misv)
      wtsb = get_weights(Tb, Db, Ab, Vb, n, misv = misv)
      #calculate delta
      TT = matrix(rep(Tb, n), n)
      Ib = (TT >= t(TT)) + 0
      ghatb = (wtsb$wd %*% Ib) / apply(wtsb$wd, 1, sum)
      deltab = (ghatb %*% t(wtsb$wn)) / apply(wtsb$wn, 1, sum)
      #deltab gives all 16 combs, we only need 4 diagonal elements
      HEL_Ub = drop(wtsb$wn %*% Ib / apply(wtsb$wn, 1, sum))
      HEL_Lb = wtsb$wd * (1 - HEL_Ub - diag(deltahat)) / apply(wtsb$wd, 1, sum)
      HEL_rb = apply(HEL_Lb, 1, function(row) {
        tryCatch(
          el.test(row, mu = 0)$'-2LLR',
          error = function(e)
            NA
        )
      })
      return(c(diag(deltab), unname(HEL_rb)))
    }
  }
  
  HEL_LLR = function(delta) {
    #calculate delta and U for HEL
    #HEL_U = drop(wts$wn[i, ] %*% It / sum(wts$wn[i, ]))
    HEL_L = wts$wd[i, ] * (1 - HEL_U[i, ] - delta) / sum(wts$wd[i, ])
    HEL_r = el.test(HEL_L, mu = 0)$'-2LLR'
    return(HEL_r)
  }
  
  HEL1_LLR = function(delta) {
    HEL_LLR(delta) - BSQ[i]
  }
  
  HEL2_LLR = function(delta) {
    HEL_LLR(delta) * (7 / 9) ^ 3 / BSMedian[i] - qchisq(1 - alpha, 1)
  }
  
  Find_zero = function(f) {
    est_zero = which(sapply(seq(
      search_start, search_end, search_step
    ), f) < 0)
    f0 = f(0)
    f1 = f(1)
    if (f(0) < 0) {
      LB = 0
    } else {
      LB = uniroot(
        f,
        lower = search_start + search_step * (est_zero[1] - 2),
        upper = search_start + search_step * (est_zero[1] - 1),
        tol = tol
      )$root
    }
    if (f(1) < 0) {
      UB = 1
    } else {
      UB = uniroot(
        f,
        lower = search_start + search_step * (max(est_zero) - 1),
        upper = search_start + search_step * max(est_zero),
        tol = tol
      )$root
    }
    # length = UB - LB
    # #p.est = (UB + LB) / 2
    # count = (f(true_auc) <= 0)
    return(c(LB, UB))
  }
  
  wts = get_weights(T, D, A, V, n, misv = misv)
  TT = matrix(rep(T, n), n)
  It = (TT >= t(TT)) + 0
  ghat = (wts$wd %*% It) / apply(wts$wd, 1, sum)
  deltahat = (ghat %*% t(wts$wn)) / apply(wts$wn, 1, sum)
  HEL_U = (wts$wn %*% It) / apply(wts$wn, 1, sum)
  BSP.est = do.call(cbind, Filter(Negate(is.null), lapply(1:B, FUN = BSP)))
  #BSP.est = do.call(cbind, BSP.est[!sapply(BSP.est, is.null)])
  BSP.delta = BSP.est[1:4, ]
  BSP.r = BSP.est[5:8, ]
  BSP.delta[BSP.delta < 0] = 0
  BSP.delta[BSP.delta > 1] = 1
  BPL = apply(BSP.delta, 1, quantile, probs = 0.025, na.rm = TRUE)
  BPU = apply(BSP.delta, 1, quantile, probs = 0.975, na.rm = TRUE)
  BCmean = apply(BSP.delta, 1, mean, na.rm = TRUE)
  BCvar = apply(BSP.delta, 1, var, na.rm = TRUE)
  BCL = BCmean - za * sqrt(BCvar)
  BCU = BCmean + za * sqrt(BCvar)
  BCL = (BCL > 1) + BCL * between(BCL, 0, 1)
  BCU = (BCU > 1) + BCU * between(BCU, 0, 1)
  # LBP = BPU - BPL
  # LBC = BCU - BCL
  # CBP = mapply(between, true_auc, BPL, BPU)
  # CBC = mapply(between, true_auc, BCL, BCU)
  BSQ = apply(BSP.r, 1, quantile, probs = 0.95, na.rm = TRUE)
  BSMedian = apply(BSP.r, 1, median, na.rm = TRUE)
  
  HEL1_CI = NULL
  for (i in 1:4) {
    HEL1_CI = rbind(HEL1_CI, Find_zero(HEL1_LLR))
  }
  HEL2_CI = NULL
  for (i in 1:4) {
    HEL2_CI = rbind(HEL2_CI, Find_zero(HEL2_LLR))
  }
  data_info = rbind(c(true_auc[bmk,1], sum(V) / n), 
                    c(sum(V * D), sum(V * (1 - D))), 
                    DeLong_ci, 
                    BC_ci)
  CIs = cbind(data_info, cbind(BCL, BCU, BPL, BPU), HEL1_CI, HEL2_CI)
  colnames(CIs) = c("True_AUC","Verified%","BCL","BCU","BPL","BPU","HELIL","HELIU","HELIIL","HELIIU")
  rownames(CIs) = c("True_AUC", "Verified D", "DeLong_ci", "BC_ci")
  CI_comp = cbind(
    c(true_auc[bmk,1],sum(V * D),
      between(true_auc[bmk,1], DeLong_ci[1], DeLong_ci[2]),
      between(true_auc[bmk,1], BC_ci[1], BC_ci[2])),
    c(sum(V) / n, sum(V * (1 - D)), 
      DeLong_ci[2] - DeLong_ci[1], BC_ci[2] - BC_ci[1]),
    between(rep(true_auc[bmk,1],4), BCL, BCU),
    BCU - BCL,
    between(rep(true_auc[bmk,1],4), BPL, BPU),
    BPU - BPL,
    between(rep(true_auc[bmk,1],4), HEL1_CI[ , 1], HEL1_CI[ , 2]),
    HEL1_CI[ , 2] - HEL1_CI[ , 1],
    between(rep(true_auc[bmk,1],4), HEL2_CI[ , 1], HEL2_CI[ , 2]),
    HEL2_CI[ , 2] - HEL2_CI[ , 1])
  colnames(CI_comp) = c("Cover","Width","BC","BC","BP","BP","HELI","HELI","HELII","HELII")
  return(list(CIs = CIs, CI_comp = CI_comp))
}

set.seed(123)

auc_vb_ci(bmk = 2, k = 1, a = 6, misv = F)

