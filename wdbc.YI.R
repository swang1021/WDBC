library(rocbc)
library(ThresholdROC)
library(dplyr)

###Load subset of WDBC dataset with prevalence = 0.2
load('wdbc.new')
n = nrow(wdbc)
sum(wdbc$V2 == 'M') / n

###Find true Youden index using complete data
jc = matrix(0, 32, 2)
for (i in 3:32) {
  T = wdbc[, i]
  D = (wdbc$V2 == 'M') + 0
  Tsort = sort(T)
  ctest = Tsort[1:n - 1] + (Tsort[2:n] - Tsort[1:n - 1]) / 2
  s = 0
  for (c in ctest) {
    t_temp_1 = (T >= c) + 0
    t_temp_2 = (T < c) + 0
    sen = drop(t_temp_1 %*% D) / 86
    spe = drop(t_temp_2 %*% (1 - D)) / 344
    s = cbind(s, sen + spe - 1)
  }
  j = max(s)
  cp = ctest[which.max(s) - 1]
  jc[i, ] = c(j, cp) #jc contains true J and co for all biomarkers
}

#Necessary parameters
precision = 10 ^ (-3)
alpha = 0.05
za = qnorm(0.975)
iterations = 2000

#Main function to compute CIs
wdbc.ci = function(seed = 123, 
                   k = 1, 
                   bio, 
                   a = 'SPE', 
                   misv = FALSE) {
  set.seed(seed)
  D = (wdbc$V2 == 'M') + 0
  Tc = as.numeric(wdbc[, bio])
  T = as.numeric(scale(wdbc[, bio]))
  A = as.numeric(scale(wdbc[, 8]))
  #generate V
  logitp = T + A + k
  pv = 1 / (exp(-logitp) + 1)
  u = runif(n)
  V = (pv > u) + 0
  sum(V * D) > 5
  sum(V * (1 - D)) > 5
  #Apply rocbc & ThresholdROC package to compute CI ignoring missing data
  Tv = Tc[V == 1]
  Dv = D[V == 1]
  rocboxcox(Tv, Dv, alpha = 0.05, plots = 'on')$JCI
  Tvsort = sort(Tv)
  cvtest = Tvsort[1:sum(V) - 1] + (Tvsort[2:sum(V)] - Tvsort[1:sum(V) - 1]) / 2
  s = 0
  for (c in cvtest) {
    t_temp_1 = (Tv >= c) + 0
    t_temp_2 = (Tv < c) + 0
    sen = drop(t_temp_1 %*% Dv) / sum(Dv)
    spe = drop(t_temp_2 %*% (1 - Dv)) / sum(1 - Dv)
    s = cbind(s, sen + spe - 1)
  }
  jv = max(s)
  cpv = cvtest[which.max(s) - 1]
  tab = matrix(c(sum((Tv >= cpv) * Dv), sum((Tv >= cpv) * (1 - Dv)), sum((Tv <
                                                                            cpv) * Dv), sum((Tv < cpv) * (1 - Dv))), ncol = 2, byrow = TRUE)
  diagnostic(tab)[8, ]
  #Fit Disease Model
  disease = glm(D[V == 1] ~ T[V == 1] + A[V == 1],
                family = binomial(link = 'probit'),
                control = list(maxit = 100))
  beta = disease$coefficients
  rho = drop(pnorm(cbind(rep(1, n), T, A) %*% beta))
  rho[rho > 1 - precision] = 1 - precision
  rho[rho < precision] = precision
  #Fit Verification Model
  if (misv == FALSE) {
    #Correct V model
    Verif = glm(V ~ T + A,
                family = binomial(link = 'logit'),
                control = list(maxit = 100))
    betaV = Verif$coefficients
    pi_hat = (1 + exp(-drop(cbind(rep(
      1, n
    ), T, A) %*% betaV))) ^ (-1)
    pi_hat[pi_hat > 1 - precision] = 1 - precision
    pi_hat[pi_hat < precision] = precision
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
    pi_hat[pi_hat > 1 - precision] = 1 - precision
    pi_hat[pi_hat < precision] = precision
  }
  #Compute weights
  if (a == 'FI') {
    wd = rho
    wn = 1 - rho
  } else if (a == 'MSI') {
    wd = V * D + (1 - V) * rho
    wn = V * (1 - D) + (1 - V) * (1 - rho)
  } else if (a == 'IPW') {
    wd = V * D / pi_hat
    wn = V * (1 - D) / pi_hat
  } else if (a == 'SPE') {
    wd = (V * D - (V - pi_hat) * rho) / pi_hat
    wn = (V * (1 - D) - (V - pi_hat) * (1 - rho)) / pi_hat
  }
  sum.wd = sum(wd)
  sum.wn = sum(wn)
  #Search for j & cutoff point
  Tsort = sort(T)
  ctest = Tsort[1:n - 1] + (Tsort[2:n] - Tsort[1:n - 1]) / 2
  jm = matrix(0, length(ctest), 3)
  for (k in 1:length(ctest)) {
    c = ctest[k]
    t_temp_1 = (T >= c) + 0
    t_temp_2 = (T < c) + 0
    # Spe
    spe = drop(t_temp_2 %*% wn) / sum.wn
    spe = (spe > 1) + spe * between(spe, 0, 1)
    # Sen
    sen = drop(t_temp_1 %*% wd) / sum.wd
    sen = (sen > 1) + sen * between(sen, 0, 1)
    jm[k, ] = c(spe, sen, sen + spe - 1)
  }
  j = max(jm[, 3])
  cp = ctest[which.max(jm[, 3]) - 1]
  #Find Fhat & Ghat 
  Fhat = jm[which.max(jm[, 3]) - 1, 1]
  Ghat = 1 - jm[which.max(jm[, 3]) - 1, 2]
  Wald.var = Ghat * (1 - Ghat) / sum.wd + Fhat * (1 - Fhat) / sum.wn
  Wald.lb = j - za * sqrt(Wald.var)
  Wald.ub = min(1, j + za * sqrt(Wald.var))
  Wald.l = Wald.ub - Wald.lb
  #MOVER methods
  G_temp = (drop((T < cp) %*% wd) + 0.5 * za ^ 2) / (sum.wd + za ^ 2)
  Gtilde = (G_temp > 1) + G_temp * between(G_temp, 0, 1)
  F_temp = (drop((T < cp) %*% wn) + 0.5 * za ^ 2) / (sum.wn + za ^ 2)
  Ftilde = (F_temp > 1) + F_temp * between(F_temp, 0, 1)
  LGac = Gtilde - za * sqrt(Gtilde * (1 - Gtilde) / (sum.wd + za ^ 2))
  UGac = Gtilde + za * sqrt(Gtilde * (1 - Gtilde) / (sum.wd + za ^ 2))
  LFac = Ftilde - za * sqrt(Ftilde * (1 - Ftilde) / (sum.wn + za ^ 2))
  UFac = Ftilde + za * sqrt(Ftilde * (1 - Ftilde) / (sum.wn + za ^ 2))
  LGws = (Ghat + 0.5 * za ^ 2 / sum.wd - za * sqrt(Ghat * (1 - Ghat) / sum.wd +
                                                     0.25 * za ^ 2 / sum.wd ^ 2)) / (1 + za ^ 2 / sum.wd)
  UGws = (Ghat + 0.5 * za ^ 2 / sum.wd + za * sqrt(Ghat * (1 - Ghat) / sum.wd +
                                                     0.25 * za ^ 2 / sum.wd ^ 2)) / (1 + za ^ 2 / sum.wd)
  LFws = (Fhat + 0.5 * za ^ 2 / sum.wn - za * sqrt(Fhat * (1 - Fhat) / sum.wn +
                                                     0.25 * za ^ 2 / sum.wn ^ 2)) / (1 + za ^ 2 / sum.wn)
  UFws = (Fhat + 0.5 * za ^ 2 / sum.wn + za * sqrt(Fhat * (1 - Fhat) / sum.wn +
                                                     0.25 * za ^ 2 / sum.wn ^ 2)) / (1 + za ^ 2 / sum.wn)
  Mover.ac.lb = j - sqrt((Fhat - LFac) ^ 2 + (UGac - Ghat) ^ 2)
  Mover.ac.ub = min(1, j + sqrt((UFac - Fhat) ^ 2 + (Ghat - LGac) ^ 2))
  Mover.ws.lb = j - sqrt((Fhat - LFws) ^ 2 + (UGws - Ghat) ^ 2)
  Mover.ws.ub = min(1, j + sqrt((UFws - Fhat) ^ 2 + (Ghat - LGws) ^ 2))
  Mover.ac.l = Mover.ac.ub - Mover.ac.lb
  Mover.ws.l = Mover.ws.ub - Mover.ws.lb
  jtilde = Ftilde - Gtilde
  #Generate bootstrap copies of Jtilde
  BSP = function(dummy) {
    boots = sample(1:n, n, replace = TRUE)
    Tb = T[boots]
    Db = D[boots]
    Ab = A[boots]
    Vb = V[boots]
    if (sum(Vb * Db) >= 1 &&
        sum(Vb * (1 - Db)) >= 1) {
      #Correct D model
      disease = glm(
        Db[Vb == 1] ~ Tb[Vb == 1] + Ab[Vb == 1],
        family = binomial(link = 'probit'),
        control = list(maxit = 100)
      )
      beta = disease$coefficients
      rho = drop(pnorm(cbind(rep(1, n), Tb, Ab) %*% beta))
      rho[rho > 1 - precision] = 1 - precision
      rho[rho < precision] = precision
      if (misv == FALSE) {
        #Correct V model
        Verif = glm(Vb ~ Tb + Ab,
                    family = binomial(link = 'logit'),
                    control = list(maxit = 100))
        betaV = Verif$coefficients
        pi_hat = (1 + exp(-drop(cbind(
          rep(1, n), Tb, Ab
        ) %*% betaV))) ^ (-1)
        pi_hat[pi_hat > 1 - precision] = 1 - precision
        pi_hat[pi_hat < precision] = precision
      }
      else{
        #Wrong V model
        Verif = glm(Vb ~ Tb,
                    family = binomial(link = 'logit'),
                    control = list(maxit = 100))
        betaV = Verif$coefficients
        pi_hat = (1 + exp(-drop(cbind(
          rep(1, n), Tb
        ) %*% betaV))) ^ (-1)
        pi_hat[pi_hat > 1 - precision] = 1 - precision
        pi_hat[pi_hat < precision] = precision
      }
      #weights
      if (a == 'FI') {
        bwd = rho
        bwn = 1 - rho
      } else if (a == 'MSI') {
        bwd = Vb * Db + (1 - Vb) * rho
        bwn = Vb * (1 - Db) + (1 - Vb) * (1 - rho)
      } else if (a == 'IPW') {
        bwd = Vb * Db / pi_hat
        bwn = Vb * (1 - Db) / pi_hat
      } else if (a == 'SPE') {
        bwd = (Vb * Db - (Vb - pi_hat) * rho) / pi_hat
        bwn = (Vb * (1 - Db) - (Vb - pi_hat) * (1 - rho)) / pi_hat
      }
      sum.bwd = sum(bwd)
      sum.bwn = sum(bwn)
      Gtildeb = (drop((Tb < cp) %*% bwd) + 0.5 * za ^ 2) / (sum.bwd + za ^ 2)
      Gtildeb = (Gtildeb > 1) + Gtildeb * between(Gtildeb, 0, 1)
      Ftildeb = (drop((Tb < cp) %*% bwn) + 0.5 * za ^ 2) / (sum.bwn + za ^ 2)
      Ftildeb = (Ftildeb > 1) + Ftildeb * between(Ftildeb, 0, 1)
      return(Ftildeb - Gtildeb)
    }
  }
  #Bootstrap intervals
  BSP.est = sapply(1:iterations, FUN = BSP)
  BSP.est = unlist(BSP.est[lengths(BSP.est) != 0]) #remove null elements
  BPL = quantile(BSP.est,
                 probs = 0.025,
                 na.rm = TRUE,
                 names = FALSE)
  BPU = quantile(BSP.est,
                 probs = 0.975,
                 na.rm = TRUE,
                 names = FALSE)
  BCmean = mean(BSP.est, na.rm = TRUE)
  BCvar = var(BSP.est, na.rm = TRUE)
  BCL1 = jtilde - za * sqrt(BCvar)
  BCU1 = min(1, jtilde + za * sqrt(BCvar))
  BCL2 = BCmean - za * sqrt(BCvar)
  BCU2 = min(1, BCmean + za * sqrt(BCvar))
  LBP = BPU - BPL
  LBC1 = BCU1 - BCL1
  LBC2 = BCU2 - BCL2
  return(
    list(
      'True J' = jc[bio,1],
      'Verified' = c(
        '#Diseased' = sum(V * D),
        '%Diseased' = sum(V * D) / sum(D),
        '#Healthy' = sum(V * (1 - D)),
        '%Healthy' = sum(V * (1 - D)) / sum(1 - D),
        '#Verified' = sum(V),
        '%Verified' = sum(V) / n
      ),
      'rocbc.JCI' = rocboxcox(Tv, Dv, alpha = 0.05, plots = 'on')$JCI,
      'TRoc.JCI' = diagnostic(tab)[8, ],
      'BCI' = c('Lower' = BCL1, 'Upper' = BCU1),
      'BCII' = c('Lower' = BCL2, 'Upper' = BCU2),
      'BPac' = c('Lower' = BPL, 'Upper' = BPU),
      'Wald' = c('Lower' = Wald.lb, 'Upper' = Wald.ub),
      'Mover.ac' = c('Lower' = Mover.ac.lb, 'Upper' = Mover.ac.ub),
      'Mover.ws' = c('Lower' = Mover.ws.lb, 'Upper' = Mover.ws.ub)
    )
  )
}
