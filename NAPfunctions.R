
# required packages
if(!('doParallel' %in% rownames(installed.packages()))) install.packages('doParallel')

library(doParallel)

#### helper function ####
mycombine.fixed = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = vector('list', 2)
  for(k in 1:length.list.combined){
    
    list.out[[1]] = rbind(list.out[[1]], list.combined[[k]][[1]])
    list.out[[2]] = rbind(list.out[[2]], list.combined[[k]][[2]])
  }
  
  return(list.out)
}

mycombine.seq.onesample = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = vector('list', 3)
  for(k in 1:length.list.combined){
    
    list.out[[1]] = rbind(list.out[[1]], list.combined[[k]][[1]])
    list.out[[2]] = rbind(list.out[[2]], list.combined[[k]][[2]])
    list.out[[3]] = rbind(list.out[[3]], list.combined[[k]][[3]])
  }
  
  return(list.out)
}

mycombine.seq.twosample = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = vector('list', 4)
  for(k in 1:length.list.combined){
    
    list.out[[1]] = rbind(list.out[[1]], list.combined[[k]][[1]])
    list.out[[2]] = rbind(list.out[[2]], list.combined[[k]][[2]])
    list.out[[3]] = rbind(list.out[[3]], list.combined[[k]][[3]])
    list.out[[4]] = rbind(list.out[[4]], list.combined[[k]][[4]])
  }
  
  return(list.out)
}

#### Bayes factors ####
#### NAP ####
#### one-sample z ####
NAPBF_oneZ = function(obs, nObs, mean.obs, test.statistic,
                      tau.NAP = 0.3/sqrt(2),
                      sigma0 = 1){
  
  if(!missing(obs)){
    
    nObs = length(obs)
    mean.obs = mean(obs)
    
    # constant terms in BF
    r_n = (nObs*(tau.NAP^2))/(nObs*(tau.NAP^2) + 1)
    test.statistic = (sqrt(nObs)*mean.obs)/sigma0
    W = (r_n*(test.statistic^2))/2
    
    return(((nObs*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
  }else if((!missing(nObs))&&(!missing(mean.obs))){
    
    # constant terms in BF
    r_n = (nObs*(tau.NAP^2))/(nObs*(tau.NAP^2) + 1)
    test.statistic = (sqrt(nObs)*mean.obs)/sigma0
    W = (r_n*(test.statistic^2))/2
    
    return(((nObs*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
  }else if((!missing(nObs))&&(!missing(test.statistic))){
    
    # constant terms in BF
    r_n = (nObs*(tau.NAP^2))/(nObs*(tau.NAP^2) + 1)
    W = (r_n*(test.statistic^2))/2
    
    return(((nObs*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
  }
}

#### one-sample t ####
NAPBF_oneT = function(obs, nObs, mean.obs, sd.obs, test.statistic,
                      tau.NAP = 0.3/sqrt(2)){
  
  if(!missing(obs)){
    
    nObs = length(obs)
    mean.obs = mean(obs)
    sd.obs = sd(obs)
    test.statistic = (sqrt(nObs)*mean.obs)/sd.obs
    
    # constant terms in BF
    r_n = (nObs*(tau.NAP^2))/(nObs*(tau.NAP^2) + 1)
    q_n = (r_n*nObs)/(nObs-1)
    G = 1 + (test.statistic^2)/(nObs-1)
    H = 1 + ((1-r_n)*(test.statistic^2))/(nObs-1)
    
    return(((nObs*(tau.NAP^2) + 1)^(-3/2))*
             ((G/H)^(nObs/2))*(1 + (q_n*(test.statistic^2))/H))
    
  }else if((!missing(nObs))&&(!missing(mean.obs))&&(!missing(sd.obs))){
    
    test.statistic = (sqrt(nObs)*mean.obs)/sd.obs
    
    # constant terms in BF
    r_n = (nObs*(tau.NAP^2))/(nObs*(tau.NAP^2) + 1)
    q_n = (r_n*nObs)/(nObs-1)
    G = 1 + (test.statistic^2)/(nObs-1)
    H = 1 + ((1-r_n)*(test.statistic^2))/(nObs-1)
    
    return(((nObs*(tau.NAP^2) + 1)^(-3/2))*
             ((G/H)^(nObs/2))*(1 + (q_n*(test.statistic^2))/H))
    
  }else if((!missing(nObs))&&(!missing(test.statistic))){
    
    # constant terms in BF
    r_n = (nObs*(tau.NAP^2))/(nObs*(tau.NAP^2) + 1)
    q_n = (r_n*nObs)/(nObs-1)
    G = 1 + (test.statistic^2)/(nObs-1)
    H = 1 + ((1-r_n)*(test.statistic^2))/(nObs-1)
    
    return(((nObs*(tau.NAP^2) + 1)^(-3/2))*
             ((G/H)^(nObs/2))*(1 + (q_n*(test.statistic^2))/H))
  }
}

#### two-sample z ####
NAPBF_twoZ = function(obs1, obs2, n1Obs, n2Obs, 
                      mean.obs1, mean.obs2, test.statistic,
                      tau.NAP = 0.3/sqrt(2),
                      sigma0 = 1){
  
  if((!missing(obs1))&&(!missing(obs2))){
    
    n1Obs = length(obs1)
    n2Obs = length(obs2)
    mean.obs1 = mean(obs1)
    mean.obs2 = mean(obs2)
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/sigma0
    W = (r_n*(test.statistic^2))/2
    
    return(((m_n*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(mean.obs1))&&(!missing(mean.obs2))){
    
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/sigma0
    W = (r_n*(test.statistic^2))/2
    
    return(((m_n*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(test.statistic))){
    
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    
    # BF at step n for those not reached decision
    W = (r_n*(test.statistic^2))/2
    
    return(((m_n*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
  }
}

#### two-sample t ####
NAPBF_twoT = function(obs1, obs2, n1Obs, n2Obs, 
                      mean.obs1, mean.obs2, sd.obs1, sd.obs2, test.statistic,
                      tau.NAP = 0.3/sqrt(2)){
  
  if((!missing(obs1))&&(!missing(obs2))){
    
    n1Obs = length(obs1)
    n2Obs = length(obs2)
    mean.obs1 = mean(obs1)
    mean.obs2 = mean(obs2)
    sd.obs1 = sd(obs1)
    sd.obs2 = sd(obs2)
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    q_n = (r_n*(n1Obs+n2Obs-1))/(n1Obs+n2Obs-2)
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt(((n1Obs-1)*(sd.obs1^2) + (n2Obs-1)*(sd.obs2^2))/
                       (n1Obs+n2Obs-2))
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/pooleds_n
    G = 1 + (test.statistic^2)/(n1Obs+n2Obs-2)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n1Obs+n2Obs-2)
    
    return(((m_n*(tau.NAP^2) + 1)^(-3/2))*
             ((G/H)^((n1Obs+n2Obs-1)/2))*
             (1 + (q_n*(test.statistic^2))/H))
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(mean.obs1))&&(!missing(mean.obs2))&&
           (!missing(sd.obs1))&&(!missing(sd.obs2))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    q_n = (r_n*(n1Obs+n2Obs-1))/(n1Obs+n2Obs-2)
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt(((n1Obs-1)*(sd.obs1^2) + (n2Obs-1)*(sd.obs2^2))/
                       (n1Obs+n2Obs-2))
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/pooleds_n
    G = 1 + (test.statistic^2)/(n1Obs+n2Obs-2)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n1Obs+n2Obs-2)
    
    return(((m_n*(tau.NAP^2) + 1)^(-3/2))*
             ((G/H)^((n1Obs+n2Obs-1)/2))*
             (1 + (q_n*(test.statistic^2))/H))
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(test.statistic))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    q_n = (r_n*(n1Obs+n2Obs-1))/(n1Obs+n2Obs-2)
    
    # BF at step n for those not reached decision
    G = 1 + (test.statistic^2)/(n1Obs+n2Obs-2)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n1Obs+n2Obs-2)
    
    return(((m_n*(tau.NAP^2) + 1)^(-3/2))*
             ((G/H)^((n1Obs+n2Obs-1)/2))*
             (1 + (q_n*(test.statistic^2))/H))
    
  }
}


#### Composite alternative and Hajnal's ratio ####
#### one-sample z ####
HajnalBF_oneZ = function(obs, nObs, mean.obs, test.statistic, 
                         es1 = 0.3, sigma0 = 1){
  
  if(!missing(obs)){
    
    nObs = length(obs)
    mean.obs = mean(obs)
    
    test.statistic = (sqrt(nObs)*mean.obs)/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = sqrt(nObs)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(nObs)*abs(es1))/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(nObs))&&(!missing(mean.obs))){
    
    test.statistic = (sqrt(nObs)*mean.obs)/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = sqrt(nObs)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(nObs)*abs(es1))/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(nObs))&&(!missing(test.statistic))){
    
    return((dnorm(x = test.statistic,
                  mean = sqrt(nObs)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(nObs)*abs(es1))/
              dnorm(x = test.statistic))/2)
  }
}

#### one-sample t ####
HajnalBF_oneT = function(obs, nObs, mean.obs, sd.obs, test.statistic, 
                         es1 = 0.3){
  
  if(!missing(obs)){
    
    nObs = length(obs)
    mean.obs = mean(obs)
    sd.obs = sd(obs)
    test.statistic = (sqrt(nObs)*mean.obs)/sd.obs
    
    return(df(x = test.statistic^2, df1 = 1, df2 = nObs-1,
              ncp = nObs*(es1^2))/
             df(x = test.statistic^2, df1 = 1, df2 = nObs-1))
    
  }else if((!missing(nObs))&&(!missing(mean.obs))&&(!missing(sd.obs))){
    
    test.statistic = (sqrt(nObs)*mean.obs)/sd.obs
    
    return(df(x = test.statistic^2, df1 = 1, df2 = nObs-1,
              ncp = nObs*(es1^2))/
             df(x = test.statistic^2, df1 = 1, df2 = nObs-1))
    
  }else if((!missing(nObs))&&(!missing(test.statistic))){
    
    return(df(x = test.statistic^2, df1 = 1, df2 = nObs-1,
              ncp = nObs*(es1^2))/
             df(x = test.statistic^2, df1 = 1, df2 = nObs-1))
  }
}

#### two-sample z ####
HajnalBF_twoZ = function(obs1, obs2, n1Obs, n2Obs, 
                         mean.obs1, mean.obs2, test.statistic,
                         es1 = 0.3, sigma0 = 1){
  
  if((!missing(obs1))&&(!missing(obs2))){
    
    n1Obs = length(obs1)
    n2Obs = length(obs2)
    mean.obs1 = mean(obs1)
    mean.obs2 = mean(obs2)
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(mean.obs1))&&(!missing(mean.obs2))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(test.statistic))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    return((dnorm(x = test.statistic,
                  mean = sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic))/2)
  }
}

#### two-sample t ####
HajnalBF_twoT = function(obs1, obs2, n1Obs, n2Obs, 
                         mean.obs1, mean.obs2, sd.obs1, sd.obs2, test.statistic,
                         es1 = 0.3){
  
  if((!missing(obs1))&&(!missing(obs2))){
    
    n1Obs = length(obs1)
    n2Obs = length(obs2)
    mean.obs1 = mean(obs1)
    mean.obs2 = mean(obs2)
    sd.obs1 = sd(obs1)
    sd.obs2 = sd(obs2)
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt(((n1Obs-1)*(sd.obs1^2) + (n2Obs-1)*(sd.obs2^2))/
                       (n1Obs+n2Obs-2))
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/pooleds_n
    
    return(df(x = test.statistic^2, df1 = 1, df2 = n1Obs+n2Obs-2,
              ncp = m_n*(es1^2))/
             df(x = test.statistic^2, df1 = 1, df2 = n1Obs+n2Obs-2))
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(mean.obs1))&&(!missing(mean.obs2))&&
           (!missing(sd.obs1))&&(!missing(sd.obs2))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt(((n1Obs-1)*(sd.obs1^2) + (n2Obs-1)*(sd.obs2^2))/
                       (n1Obs+n2Obs-2))
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/pooleds_n
    
    return(df(x = test.statistic^2, df1 = 1, df2 = n1Obs+n2Obs-2,
              ncp = m_n*(es1^2))/
             df(x = test.statistic^2, df1 = 1, df2 = n1Obs+n2Obs-2))
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(test.statistic))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    return(df(x = test.statistic^2, df1 = 1, df2 = n1Obs+n2Obs-2,
              ncp = m_n*(es1^2))/
             df(x = test.statistic^2, df1 = 1, df2 = n1Obs+n2Obs-2))
  }
}


#### fixed design NAP for a given n at multiple effect sizes ####
#### one-sample z ####
fixedNAP.oneZ_n = function(es = c(0, .2, .3, .5), n.fixed = 20, 
                           tau.NAP = 0.3/sqrt(2), sigma0 = 1,
                           nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch.size.increment = function(narg){100}
  
  nmin = min(100, n.fixed)
  nmax = n.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumsum_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:nmin,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen*sigma0, sigma0)
                       })
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:n.increment,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen*sigma0, sigma0)
                       })
      }
      
      # sum of observations until step n
      cumsum_n = cumsum_n + rowSums(obs_n)
      
      seq.converged_es = (n==nmax)
    }
    
    # constant terms in BF
    r_n = (n*(tau.NAP^2))/(n*(tau.NAP^2) + 1)
    
    # xbar and S until step n
    xbar_n = cumsum_n/n
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n)*xbar_n)/sigma0
    W = (r_n*(test.statistic^2))/2
    BF_n = ((n*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W)
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### one-sample t ####
fixedNAP.oneT_n = function(es = c(0, .2, .3, .5), n.fixed = 20, 
                           tau.NAP = 0.3/sqrt(2),
                           nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch.size.increment = function(narg){100}
  
  nmin = min(100, n.fixed)
  nmax = n.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumSS_n = cumsum_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:nmin,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen, 1)
                       })
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:n.increment,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen, 1)
                       })
      }
      
      # sum of observations until step n
      cumsum_n = cumsum_n + rowSums(obs_n)
      
      # sum of squares of observations until step n
      cumSS_n = cumSS_n + rowSums(obs_n^2)
      
      seq.converged_es = (n==nmax)
    }
    
    # constant terms in BF
    r_n = (n*(tau.NAP^2))/(n*(tau.NAP^2) + 1)
    q_n = (r_n*n)/(n-1)
    
    # xbar and S until step n
    xbar_n = cumsum_n/n
    s_n = sqrt((cumSS_n - (n*(xbar_n^2)))/(n-1))
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n)*xbar_n)/s_n
    G = 1 + (test.statistic^2)/(n-1)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n-1)
    BF_n = ((n*(tau.NAP^2) + 1)^(-3/2))*
      ((G/H)^(n/2))*(1 + (q_n*(test.statistic^2))/H)
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### two-sample z ####
fixedNAP.twoZ_n = function(es = c(0, .2, .3, .5), n1.fixed = 20, n2.fixed = 20,
                           tau.NAP = 0.3/sqrt(2), sigma0 = 1,
                           nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch1.size.increment = batch2.size.increment = function(narg){100}
  
  n1min = min(100, n1.fixed)
  n2min = min(100, n2.fixed)
  
  n1max = n1.fixed
  n2max = n2.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumsum1_n = cumsum2_n = BF_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        obs1_n = mapply(X = 1:n1,
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es.gen/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        obs2_n = mapply(X = 1:n2,
                        FUN = function(X){
                          
                          rnorm(nReplicate, es.gen/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(n1.increment>0){
          
          obs1_n = mapply(X = 1:n1.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, sigma0)
                          })
          
          # sum of observations until step n
          cumsum1_n = cumsum1_n + rowSums(obs1_n)
        }
        
        if(n2.increment>0){
          
          obs2_n = mapply(X = 1:n2.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, sigma0)
                          })
          
          # sum of observations until step n
          cumsum2_n = cumsum2_n + rowSums(obs2_n)
        }
      }
      
      seq.converged_es = (n1==n1max)&&(n2==n2max)
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    
    # xbar and S until step n
    xbar1_n = cumsum1_n/n1
    xbar2_n = cumsum2_n/n2
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/sigma0
    W = (r_n*(test.statistic^2))/2
    BF_n = ((m_n*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W)
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### two-sample t ####
fixedNAP.twoT_n = function(es = c(0, .2, .3, .5), n1.fixed = 20, n2.fixed = 20,
                           tau.NAP = 0.3/sqrt(2),
                           nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch1.size.increment = batch2.size.increment = function(narg){100}
  
  n1min = min(100, n1.fixed)
  n2min = min(100, n2.fixed)
  
  n1max = n1.fixed
  n2max = n2.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumSS1_n = cumSS2_n = cumsum1_n = cumsum2_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        obs1_n = mapply(X = 1:n1,
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es.gen/2, 1)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        # sum of squares of observations until step n
        cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
        
        obs2_n = mapply(X = 1:n2,
                        FUN = function(X){
                          
                          rnorm(nReplicate, es.gen/2, 1)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
        # sum of squares of observations until step n
        cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(n1.increment>0){
          
          obs1_n = mapply(X = 1:n1.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, 1)
                          })
          
          # sum of observations until step n
          cumsum1_n = cumsum1_n + rowSums(obs1_n)
          
          # sum of squares of observations until step n
          cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
        }
        
        if(n2.increment>0){
          
          obs2_n = mapply(X = 1:n2.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, 1)
                          })
          
          # sum of observations until step n
          cumsum2_n = cumsum2_n + rowSums(obs2_n)
          
          # sum of squares of observations until step n
          cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
        }
      }
      
      seq.converged_es = (n1==n1max)&&(n2==n2max)
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    q_n = (r_n*(n1+n2-1))/(n1+n2-2)
    
    # xbar and S until step n
    xbar1_n = cumsum1_n/n1
    xbar2_n = cumsum2_n/n2
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt((cumSS1_n - n1*(xbar1_n^2) +
                        cumSS2_n - n2*(xbar2_n^2))/(n1+n2-2))
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/pooleds_n
    G = 1 + (test.statistic^2)/(n1+n2-2)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n1+n2-2)
    BF_n = ((m_n*(tau.NAP^2) + 1)^(-3/2))*
      ((G/H)^((n1+n2-1)/2))*(1 + (q_n*(test.statistic^2))/H)
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}


#### fixed design with composite alternative for a given n at multiple effect sizes ####
#### one-sample z ####
fixedHajnal.oneZ_n = function(es1 = 0.3, es = c(0, .2, .3, .5), n.fixed = 20, 
                              sigma0 = 1,
                              nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch.size.increment = function(narg){100}
  
  nmin = min(100, n.fixed)
  nmax = n.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumsum_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:nmin,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen, sigma0)
                       })
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:n.increment,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen, sigma0)
                       })
      }
      
      # sum of observations until step n
      cumsum_n = cumsum_n + rowSums(obs_n)
      
      seq.converged_es = (n==nmax)
    }
    
    # xbar and S until step n
    xbar_n = cumsum_n/n
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n)*xbar_n)/sigma0
    BF_n = (dnorm(x = test.statistic,
                  mean = sqrt(n)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(n)*abs(es1))/
              dnorm(x = test.statistic))/2
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### one-sample t ####
fixedHajnal.oneT_n = function(es1 = 0.3, es = c(0, .2, .3, .5), n.fixed = 20, 
                              nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch.size.increment = function(narg){100}
  
  nmin = min(100, n.fixed)
  nmax = n.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumSS_n = cumsum_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:nmin,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen, 1)
                       })
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        obs_n = mapply(X = 1:n.increment,
                       FUN = function(X){
                         
                         rnorm(nReplicate, es.gen, 1)
                       })
      }
      
      # sum of observations until step n
      cumsum_n = cumsum_n + rowSums(obs_n)
      
      # sum of squares of observations until step n
      cumSS_n = cumSS_n + rowSums(obs_n^2)
      
      seq.converged_es = (n==nmax)
    }
    
    # xbar and S until step n
    xbar_n = cumsum_n/n
    s_n = sqrt((cumSS_n - (n*(xbar_n^2)))/(n-1))
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n)*xbar_n)/s_n
    BF_n = df(x = test.statistic^2, df1 = 1, df2 = n-1,
              ncp = n*(es1^2))/
      df(x = test.statistic^2, df1 = 1, df2 = n-1)
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### two-sample z ####
fixedHajnal.twoZ_n = function(es1 = 0.3, es = c(0, .2, .3, .5), n1.fixed = 20, n2.fixed = 20,
                              sigma0 = 1, nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch1.size.increment = batch2.size.increment = function(narg){100}
  
  n1min = min(100, n1.fixed)
  n2min = min(100, n2.fixed)
  
  n1max = n1.fixed
  n2max = n2.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumsum1_n = cumsum2_n = BF_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        obs1_n = mapply(X = 1:n1,
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es.gen/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        obs2_n = mapply(X = 1:n2,
                        FUN = function(X){
                          
                          rnorm(nReplicate, es.gen/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(n1.increment>0){
          
          obs1_n = mapply(X = 1:n1.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, sigma0)
                          })
          
          # sum of observations until step n
          cumsum1_n = cumsum1_n + rowSums(obs1_n)
        }
        
        if(n2.increment>0){
          
          obs2_n = mapply(X = 1:n2.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, sigma0)
                          })
          
          # sum of observations until step n
          cumsum2_n = cumsum2_n + rowSums(obs2_n)
        }
      }
      
      seq.converged_es = (n1==n1max)&&(n2==n2max)
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    
    # xbar and S until step n
    xbar1_n = cumsum1_n/n1
    xbar2_n = cumsum2_n/n2
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/sigma0
    BF_n = (dnorm(x = test.statistic,
                  mean = sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = -sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic))/2
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### two-sample t ####
fixedHajnal.twoT_n = function(es1 = 0.3, es = c(0, .2, .3, .5), n1.fixed = 20, n2.fixed = 20,
                              nReplicate = 5e+4, nCore){
  
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  batch1.size.increment = batch2.size.increment = function(narg){100}
  
  n1min = min(100, n1.fixed)
  n2min = min(100, n2.fixed)
  
  n1max = n1.fixed
  n2max = n2.fixed
  
  doParallel::registerDoParallel(cores = nCore)
  out.combined = foreach(es.gen = es, .combine = 'mycombine.fixed', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF, posterior prob that H0 is true and
    ## posterior mean of effect size
    # required storages
    cumSS1_n = cumSS2_n = cumsum1_n = cumsum2_n = numeric(nReplicate)
    
    seq.step = 0
    seq.converged_es = F
    while(!seq.converged_es){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        obs1_n = mapply(X = 1:n1,
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es.gen/2, 1)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        # sum of squares of observations until step n
        cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
        
        obs2_n = mapply(X = 1:n2,
                        FUN = function(X){
                          
                          rnorm(nReplicate, es.gen/2, 1)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
        # sum of squares of observations until step n
        cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(n1.increment>0){
          
          obs1_n = mapply(X = 1:n1.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, 1)
                          })
          
          # sum of observations until step n
          cumsum1_n = cumsum1_n + rowSums(obs1_n)
          
          # sum of squares of observations until step n
          cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
        }
        
        if(n2.increment>0){
          
          obs2_n = mapply(X = 1:n2.increment,
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, 1)
                          })
          
          # sum of observations until step n
          cumsum2_n = cumsum2_n + rowSums(obs2_n)
          
          # sum of squares of observations until step n
          cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
        }
      }
      
      seq.converged_es = (n1==n1max)&&(n2==n2max)
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    
    # xbar and S until step n
    xbar1_n = cumsum1_n/n1
    xbar2_n = cumsum2_n/n2
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt((cumSS1_n - n1*(xbar1_n^2) +
                        cumSS2_n - n2*(xbar2_n^2))/(n1+n2-2))
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/pooleds_n
    BF_n = df(x = test.statistic^2, df1 = 1, df2 = n1+n2-2,
              ncp = m_n*(es1^2))/
      df(x = test.statistic^2, df1 = 1, df2 = n1+n2-2)
    
    list(c(es.gen, mean(log(BF_n))), BF_n)
  }
  
  names(out.combined) = c('summary', 'BF')
  
  if(length(es)==1){
    
    out.combined$summary = as.data.frame(matrix(out.combined$summary,
                                                nrow = 1, ncol = 2, byrow = T))
    
    out.combined$BF = matrix(data = out.combined$BF, nrow = 1,
                             ncol = nReplicate, byrow = T)
  }
  colnames(out.combined$summary) = c('effect.size', 'avg.logBF')
  
  return(out.combined)
}

#### fixed design NAP at a given effect size for varied n ####
#### one-sample z ####
fixedNAP.oneZ_es = function(es = 0, nmin = 20, nmax = 5000,
                            tau.NAP = 0.3/sqrt(2), sigma0 = 1,
                            batch.size.increment, nReplicate = 5e+4){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n.seq = nmin
  while(n.seq[length(n.seq)]<nmax){
    
    n.seq = c(n.seq,
              n.seq[length(n.seq)] + 
                min(batch.size.increment(n.seq[length(n.seq)]),
                    nmax-n.seq[length(n.seq)]))
  }
  nStep = length(n.seq)
  
  ## simulating data and calculating the BF
  # required storages
  cumsum_n = numeric(nReplicate)
  BF = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:n.seq[seq.step],
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, sigma0)
                     })
      
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:(n.seq[seq.step]-n.seq[seq.step-1]),
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, sigma0)
                     })
    }
    
    # constant terms in BF
    r_n = (n.seq[seq.step]*(tau.NAP^2))/(n.seq[seq.step]*(tau.NAP^2) + 1)
    
    # sum of observations until step n
    cumsum_n = cumsum_n + rowSums(obs_n)
    
    # xbar and S until step n
    xbar_n = cumsum_n/n.seq[seq.step]
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n.seq[seq.step])*xbar_n)/sigma0
    W = (r_n*(test.statistic^2))/2
    BF = cbind(BF,
               ((n.seq[seq.step]*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n' = n.seq, 
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}

#### one-sample t ####
fixedNAP.oneT_es = function(es = 0, nmin = 20, nmax = 5000,
                            tau.NAP = 0.3/sqrt(2),
                            batch.size.increment, nReplicate = 5e+4){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n.seq = nmin
  while(n.seq[length(n.seq)]<nmax){
    
    n.seq = c(n.seq,
              n.seq[length(n.seq)] + 
                min(batch.size.increment(n.seq[length(n.seq)]),
                    nmax-n.seq[length(n.seq)]))
  }
  nStep = length(n.seq)
  
  ## simulating data and calculating the BF
  # required storages
  cumSS_n = cumsum_n = numeric(nReplicate)
  BF = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:n.seq[seq.step],
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, 1)
                     })
      
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:(n.seq[seq.step]-n.seq[seq.step-1]),
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, 1)
                     })
    }
    
    # constant terms in BF
    r_n = (n.seq[seq.step]*(tau.NAP^2))/(n.seq[seq.step]*(tau.NAP^2) + 1)
    q_n = (r_n*n.seq[seq.step])/(n.seq[seq.step]-1)
    
    # sum of observations until step n
    cumsum_n = cumsum_n + rowSums(obs_n)
    
    # sum of squares of observations until step n
    cumSS_n = cumSS_n + rowSums(obs_n^2)
    
    # xbar and S until step n
    xbar_n = cumsum_n/n.seq[seq.step]
    s_n = sqrt((cumSS_n - (n.seq[seq.step]*(xbar_n^2)))/
                 (n.seq[seq.step]-1))
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n.seq[seq.step])*xbar_n)/s_n
    G = 1 + (test.statistic^2)/(n.seq[seq.step]-1)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n.seq[seq.step]-1)
    BF = cbind(BF,
               ((n.seq[seq.step]*(tau.NAP^2) + 1)^(-3/2))*
                 ((G/H)^(n.seq[seq.step]/2))*(1 + (q_n*(test.statistic^2))/H))
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n' = n.seq, 
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}

#### two-sample z ####
fixedNAP.twoZ_es = function(es = 0, n1min = 20, n2min = 20,
                            n1max = 5000, n2max = 5000,
                            tau.NAP = 0.3/sqrt(2), sigma0 = 1,
                            batch1.size.increment, batch2.size.increment, nReplicate = 5e+4){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n1.seq = n1min
  n2.seq = n2min
  while(any(n1.seq[length(n1.seq)]<n1max, n2.seq[length(n2.seq)]<n2max)){
    
    n1.seq = c(n1.seq,
               n1.seq[length(n1.seq)] + 
                 min(batch1.size.increment(n1.seq[length(n1.seq)]),
                     n1max-n1.seq[length(n1.seq)]))
    n2.seq = c(n2.seq,
               n2.seq[length(n2.seq)] + 
                 min(batch2.size.increment(n2.seq[length(n2.seq)]),
                     n2max-n2.seq[length(n2.seq)]))
  }
  nStep = length(n1.seq)
  
  ## simulating data and calculating the BF, posterior prob that H0 is true and
  ## posterior mean of effect size
  # required storages
  cumsum1_n = cumsum2_n = numeric(nReplicate)
  BF = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs1_n = mapply(X = 1:n1.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, -es/2, sigma0)
                      })
      
      # sum of observations until step n
      cumsum1_n = cumsum1_n + rowSums(obs1_n)
      
      # xbar and S until step n
      xbar1_n = cumsum1_n/n1.seq[seq.step]
      
      obs2_n = mapply(X = 1:n2.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, es/2, sigma0)
                      })
      
      # sum of observations until step n
      cumsum2_n = cumsum2_n + rowSums(obs2_n)
      
      # xbar and S until step n
      xbar2_n = cumsum2_n/n2.seq[seq.step]
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      if(n1.seq[seq.step]>n1.seq[seq.step-1]){
        
        obs1_n = mapply(X = 1:(n1.seq[seq.step]-n1.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        # xbar and S until step n
        xbar1_n = cumsum1_n/n1.seq[seq.step]
      }
      
      if(n2.seq[seq.step]>n2.seq[seq.step-1]){
        
        obs2_n = mapply(X = 1:(n2.seq[seq.step]-n2.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, es/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
        # xbar and S until step n
        xbar2_n = cumsum2_n/n2.seq[seq.step]
      }
    }
    
    # constant terms in BF
    m_n = (n1.seq[seq.step]*n2.seq[seq.step])/(n1.seq[seq.step]+n2.seq[seq.step])
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/sigma0
    W = (r_n*(test.statistic^2))/2
    BF = cbind(BF, ((m_n*(tau.NAP^2) + 1)^(-3/2))*(1 + 2*W)*exp(W))
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n1' = n1.seq, 'n2' = n2.seq,
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}

#### two-sample t ####
fixedNAP.twoT_es = function(es = 0, n1min = 20, n2min = 20,
                            n1max = 5000, n2max = 5000,
                            tau.NAP = 0.3/sqrt(2),
                            batch1.size.increment, batch2.size.increment, nReplicate = 5e+4){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n1.seq = n1min
  n2.seq = n2min
  while(any(n1.seq[length(n1.seq)]<n1max, n2.seq[length(n2.seq)]<n2max)){
    
    n1.seq = c(n1.seq,
               n1.seq[length(n1.seq)] + 
                 min(batch1.size.increment(n1.seq[length(n1.seq)]),
                     n1max-n1.seq[length(n1.seq)]))
    n2.seq = c(n2.seq,
               n2.seq[length(n2.seq)] + 
                 min(batch2.size.increment(n2.seq[length(n2.seq)]),
                     n2max-n2.seq[length(n2.seq)]))
  }
  nStep = length(n1.seq)
  
  ## simulating data and calculating the BF, posterior prob that H0 is true and
  ## posterior mean of effect size
  # required storages
  cumSS1_n = cumSS2_n = cumsum1_n = cumsum2_n = numeric(nReplicate)
  BF = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs1_n = mapply(X = 1:n1.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, -es/2, 1)
                      })
      
      # sum of observations until step n
      cumsum1_n = cumsum1_n + rowSums(obs1_n)
      
      # sum of squares of observations until step n
      cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
      
      # xbar and S until step n
      xbar1_n = cumsum1_n/n1.seq[seq.step]
      
      obs2_n = mapply(X = 1:n2.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, es/2, 1)
                      })
      
      # sum of observations until step n
      cumsum2_n = cumsum2_n + rowSums(obs2_n)
      
      # sum of squares of observations until step n
      cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
      
      # xbar and S until step n
      xbar2_n = cumsum2_n/n2.seq[seq.step]
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      if(n1.seq[seq.step]>n1.seq[seq.step-1]){
        
        obs1_n = mapply(X = 1:(n1.seq[seq.step]-n1.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es/2, 1)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        # sum of squares of observations until step n
        cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
        
        # xbar and S until step n
        xbar1_n = cumsum1_n/n1.seq[seq.step]
      }
      
      if(n2.seq[seq.step]>n2.seq[seq.step-1]){
        
        obs2_n = mapply(X = 1:(n2.seq[seq.step]-n2.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, es/2, 1)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
        # sum of squares of observations until step n
        cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
        
        # xbar and S until step n
        xbar2_n = cumsum2_n/n2.seq[seq.step]
      }
    }
    
    # constant terms in BF
    m_n = (n1.seq[seq.step]*n2.seq[seq.step])/(n1.seq[seq.step]+n2.seq[seq.step])
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    q_n = (r_n*(n1.seq[seq.step]+n2.seq[seq.step]-1))/(n1.seq[seq.step]+n2.seq[seq.step]-2)
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt((cumSS1_n - n1.seq[seq.step]*(xbar1_n^2) +
                        cumSS2_n - n2.seq[seq.step]*(xbar2_n^2))/
                       (n1.seq[seq.step]+n2.seq[seq.step]-2))
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/pooleds_n
    G = 1 + (test.statistic^2)/(n1.seq[seq.step]+n2.seq[seq.step]-2)
    H = 1 + ((1-r_n)*(test.statistic^2))/(n1.seq[seq.step]+n2.seq[seq.step]-2)
    BF = cbind(BF, 
               ((m_n*(tau.NAP^2) + 1)^(-3/2))*
                 ((G/H)^((n1.seq[seq.step]+n2.seq[seq.step]-1)/2))*
                 (1 + (q_n*(test.statistic^2))/H))
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n1' = n1.seq, 'n2' = n2.seq,
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}


#### fixed design with composite alternative at a given effect size for varied n ####
#### one-sample z ####
fixedHajnal.oneZ_es = function(es = 0, es1 = 0.3, nmin = 20, nmax = 5000,
                               sigma0 = 1, batch.size.increment, nReplicate = 5e+4){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n.seq = nmin
  while(n.seq[length(n.seq)]<nmax){
    
    n.seq = c(n.seq,
              n.seq[length(n.seq)] + 
                min(batch.size.increment(n.seq[length(n.seq)]),
                    nmax-n.seq[length(n.seq)]))
  }
  nStep = length(n.seq)
  
  ## simulating data and calculating the BF
  # required storages
  cumsum_n = numeric(nReplicate)
  BF = cohend = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:n.seq[seq.step],
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, sigma0)
                     })
      
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:(n.seq[seq.step]-n.seq[seq.step-1]),
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, sigma0)
                     })
    }
    
    # sum of observations until step n
    cumsum_n = cumsum_n + rowSums(obs_n)
    
    # xbar and S until step n
    xbar_n = cumsum_n/n.seq[seq.step]
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n.seq[seq.step])*xbar_n)/sigma0
    BF = cbind(BF,
               (dnorm(x = test.statistic,
                      mean = sqrt(n.seq[seq.step])*abs(es1))/
                  dnorm(x = test.statistic) +
                  dnorm(x = test.statistic,
                        mean = -sqrt(n.seq[seq.step])*abs(es1))/
                  dnorm(x = test.statistic))/2)
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n' = n.seq, 
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}

#### one-sample t ####
fixedHajnal.oneT_es = function(es = 0, es1 = 0.3, nmin = 20, nmax = 5000,
                               batch.size.increment, nReplicate = 5e+4){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n.seq = nmin
  while(n.seq[length(n.seq)]<nmax){
    
    n.seq = c(n.seq,
              n.seq[length(n.seq)] + 
                min(batch.size.increment(n.seq[length(n.seq)]),
                    nmax-n.seq[length(n.seq)]))
  }
  nStep = length(n.seq)
  
  ## simulating data and calculating the BF
  # required storages
  cumSS_n = cumsum_n = numeric(nReplicate)
  BF = unbiased.cohend = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:n.seq[seq.step],
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, 1)
                     })
      
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      obs_n = mapply(X = 1:(n.seq[seq.step]-n.seq[seq.step-1]),
                     FUN = function(X){
                       
                       rnorm(nReplicate, es, 1)
                     })
    }
    
    # sum of observations until step n
    cumsum_n = cumsum_n + rowSums(obs_n)
    
    # sum of squares of observations until step n
    cumSS_n = cumSS_n + rowSums(obs_n^2)
    
    # xbar and S until step n
    xbar_n = cumsum_n/n.seq[seq.step]
    s_n = sqrt((cumSS_n - (n.seq[seq.step]*(xbar_n^2)))/
                 (n.seq[seq.step]-1))
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(n.seq[seq.step])*xbar_n)/s_n
    BF = cbind(BF,
               df(x = test.statistic^2, df1 = 1, df2 = n.seq[seq.step]-1,
                  ncp = n.seq[seq.step]*(es1^2))/
                 df(x = test.statistic^2, df1 = 1, df2 = n.seq[seq.step]-1))
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n' = n.seq, 
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}

#### two-sample z ####
fixedHajnal.twoZ_es = function(es = 0, es1 = 0.3, n1min = 20, n2min = 20,
                               n1max = 5000, n2max = 5000, sigma0 = 1,
                               batch1.size.increment, batch2.size.increment, nReplicate = 5e+4){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n1.seq = n1min
  n2.seq = n2min
  while(any(n1.seq[length(n1.seq)]<n1max, n2.seq[length(n2.seq)]<n2max)){
    
    n1.seq = c(n1.seq,
               n1.seq[length(n1.seq)] + 
                 min(batch1.size.increment(n1.seq[length(n1.seq)]),
                     n1max-n1.seq[length(n1.seq)]))
    n2.seq = c(n2.seq,
               n2.seq[length(n2.seq)] + 
                 min(batch2.size.increment(n2.seq[length(n2.seq)]),
                     n2max-n2.seq[length(n2.seq)]))
  }
  nStep = length(n1.seq)
  
  ## simulating data and calculating the BF, posterior prob that H0 is true and
  ## posterior mean of effect size
  # required storages
  cumsum1_n = cumsum2_n = numeric(nReplicate)
  BF = cohend = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs1_n = mapply(X = 1:n1.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, -es/2, sigma0)
                      })
      
      # sum of observations until step n
      cumsum1_n = cumsum1_n + rowSums(obs1_n)
      
      # xbar and S until step n
      xbar1_n = cumsum1_n/n1.seq[seq.step]
      
      obs2_n = mapply(X = 1:n2.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, es/2, sigma0)
                      })
      
      # sum of observations until step n
      cumsum2_n = cumsum2_n + rowSums(obs2_n)
      
      # xbar and S until step n
      xbar2_n = cumsum2_n/n2.seq[seq.step]
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      if(n1.seq[seq.step]>n1.seq[seq.step-1]){
        
        obs1_n = mapply(X = 1:(n1.seq[seq.step]-n1.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        # xbar and S until step n
        xbar1_n = cumsum1_n/n1.seq[seq.step]
      }
      
      if(n2.seq[seq.step]>n2.seq[seq.step-1]){
        
        obs2_n = mapply(X = 1:(n2.seq[seq.step]-n2.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, es/2, sigma0)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
        # xbar and S until step n
        xbar2_n = cumsum2_n/n2.seq[seq.step]
      }
    }
    
    # constant terms in BF
    m_n = (n1.seq[seq.step]*n2.seq[seq.step])/(n1.seq[seq.step]+n2.seq[seq.step])
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/sigma0
    BF = cbind(BF, 
               (dnorm(x = test.statistic,
                      mean = sqrt(m_n)*abs(es1))/
                  dnorm(x = test.statistic) +
                  dnorm(x = test.statistic,
                        mean = -sqrt(m_n)*abs(es1))/
                  dnorm(x = test.statistic))/2)
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n1' = n1.seq, 'n2' = n2.seq,
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}

#### two-sample t ####
fixedHajnal.twoT_es = function(es = 0, es1 = 0.3, n1min = 20, n2min = 20,
                               n1max = 5000, n2max = 5000,
                               batch1.size.increment, batch2.size.increment, nReplicate = 5e+4){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
  }
  
  # seq of sample size
  n1.seq = n1min
  n2.seq = n2min
  while(any(n1.seq[length(n1.seq)]<n1max, n2.seq[length(n2.seq)]<n2max)){
    
    n1.seq = c(n1.seq,
               n1.seq[length(n1.seq)] + 
                 min(batch1.size.increment(n1.seq[length(n1.seq)]),
                     n1max-n1.seq[length(n1.seq)]))
    n2.seq = c(n2.seq,
               n2.seq[length(n2.seq)] + 
                 min(batch2.size.increment(n2.seq[length(n2.seq)]),
                     n2max-n2.seq[length(n2.seq)]))
  }
  nStep = length(n1.seq)
  
  ## simulating data and calculating the BF, posterior prob that H0 is true and
  ## posterior mean of effect size
  # required storages
  cumSS1_n = cumSS2_n = cumsum1_n = cumsum2_n = numeric(nReplicate)
  BF = unbiased.cohend = NULL
  
  pb = txtProgressBar(min = 1, max = nStep, style = 3)
  for(seq.step in 1:nStep){
    
    # sample size used at this step
    if(seq.step==1){
      
      # simulating obs at this step
      set.seed(seq.step)
      obs1_n = mapply(X = 1:n1.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, -es/2, 1)
                      })
      
      # sum of observations until step n
      cumsum1_n = cumsum1_n + rowSums(obs1_n)
      
      # sum of squares of observations until step n
      cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
      
      # xbar and S until step n
      xbar1_n = cumsum1_n/n1.seq[seq.step]
      
      obs2_n = mapply(X = 1:n2.seq[seq.step],
                      FUN = function(X){
                        
                        rnorm(nReplicate, es/2, 1)
                      })
      
      # sum of observations until step n
      cumsum2_n = cumsum2_n + rowSums(obs2_n)
      
      # sum of squares of observations until step n
      cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
      
      # xbar and S until step n
      xbar2_n = cumsum2_n/n2.seq[seq.step]
      
    }else{
      
      # simulating obs at this step
      set.seed(seq.step)
      if(n1.seq[seq.step]>n1.seq[seq.step-1]){
        
        obs1_n = mapply(X = 1:(n1.seq[seq.step]-n1.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, -es/2, 1)
                        })
        
        # sum of observations until step n
        cumsum1_n = cumsum1_n + rowSums(obs1_n)
        
        # sum of squares of observations until step n
        cumSS1_n = cumSS1_n + rowSums(obs1_n^2)
        
        # xbar and S until step n
        xbar1_n = cumsum1_n/n1.seq[seq.step]
      }
      
      if(n2.seq[seq.step]>n2.seq[seq.step-1]){
        
        obs2_n = mapply(X = 1:(n2.seq[seq.step]-n2.seq[seq.step-1]),
                        FUN = function(X){
                          
                          rnorm(nReplicate, es/2, 1)
                        })
        
        # sum of observations until step n
        cumsum2_n = cumsum2_n + rowSums(obs2_n)
        
        # sum of squares of observations until step n
        cumSS2_n = cumSS2_n + rowSums(obs2_n^2)
        
        # xbar and S until step n
        xbar2_n = cumsum2_n/n2.seq[seq.step]
      }
    }
    
    # constant terms in BF
    m_n = (n1.seq[seq.step]*n2.seq[seq.step])/(n1.seq[seq.step]+n2.seq[seq.step])
    
    # BF at step n for those not reached decision
    pooleds_n = sqrt((cumSS1_n - n1.seq[seq.step]*(xbar1_n^2) +
                        cumSS2_n - n2.seq[seq.step]*(xbar2_n^2))/
                       (n1.seq[seq.step]+n2.seq[seq.step]-2))
    test.statistic = (sqrt(m_n)*(xbar2_n - xbar1_n))/pooleds_n
    BF = cbind(BF, 
               df(x = test.statistic^2, df1 = 1, df2 = n1.seq[seq.step]+n2.seq[seq.step]-2,
                  ncp = m_n*(es1^2))/
                 df(x = test.statistic^2, df1 = 1, df2 = n1.seq[seq.step]+n2.seq[seq.step]-2))
    
    setTxtProgressBar(pb, seq.step)
  }
  
  return(list('summary' = data.frame('n1' = n1.seq, 'n2' = n2.seq,
                                     'avg.logBF' = colMeans(log(BF))),
              'BF' = BF))
}


#### SBF + NAP ####
#### one-sample z ####
SBFNAP_oneZ = function(es = c(0, .2, .3, .5), nmin = 1, nmax = 5000,
                       tau.NAP = 0.3/sqrt(2), sigma0 = 1,
                       RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                       batch.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
    # batch.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.onesample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumSS.not.reached.decision = 
      cumsum.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N = rep(nmax, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:nmin,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:nmin, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = nmin,
                   byrow = T)
        }
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:n.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:n.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n.increment,
                   byrow = T)
        }
      }
      
      # constant terms in BF
      r_n = (n*(tau.NAP^2))/(n*(tau.NAP^2) + 1)
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + rowSums(obs.not.reached.decision)
      
      # xbar and S until step n
      xbar.not.reached.decision = cumsum.not.reached.decision/n
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(n)*xbar.not.reached.decision)/sigma0
      W.not.reached.decision = (r_n*(test.statistic.not.reached.decision^2))/2
      BF.not.reached.decision = ((n*(tau.NAP^2) + 1)^(-3/2))*
        (1 + 2*W.not.reached.decision)*exp(W.not.reached.decision)
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N[not.reached.decision[reached.decision_n]] = n
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum.not.reached.decision = cumsum.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n==nmax)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N), mean(log(BF))),
         N, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 6, byrow = T)
    
    out.OCandASN$N = matrix(data = out.OCandASN$N, nrow = nReplicate,
                            ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN', 'avg.logBF')
  rownames(out.OCandASN$N) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### one-sample t ####
SBFNAP_oneT = function(es = c(0, .2, .3, .5), nmin = 2, nmax = 5000,
                       tau.NAP = 0.3/sqrt(2),
                       RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                       batch.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
    # batch.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.onesample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumSS.not.reached.decision = 
      cumsum.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N = rep(nmax, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:nmin,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:nmin, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = nmin,
                   byrow = T)
        }
        
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:n.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:n.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n.increment,
                   byrow = T)
        }
      }
      
      # constant terms in BF
      r_n = (n*(tau.NAP^2))/(n*(tau.NAP^2) + 1)
      q_n = (r_n*n)/(n-1)
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + rowSums(obs.not.reached.decision)
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + rowSums(obs.not.reached.decision^2)
      
      # xbar and S until step n
      xbar.not.reached.decision = cumsum.not.reached.decision/n
      s.not.reached.decision = sqrt((cumSS.not.reached.decision - 
                                       n*(xbar.not.reached.decision^2))/(n-1))
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(n)*xbar.not.reached.decision)/s.not.reached.decision
      G.not.reached.decision = 1 + (test.statistic.not.reached.decision^2)/(n-1)
      H.not.reached.decision = 1 + ((1-r_n)*(test.statistic.not.reached.decision^2))/(n-1)
      BF.not.reached.decision = ((n*(tau.NAP^2) + 1)^(-3/2))*
        ((G.not.reached.decision/H.not.reached.decision)^(n/2))*
        (1 + (q_n*(test.statistic.not.reached.decision^2))/
           H.not.reached.decision)
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N[not.reached.decision[reached.decision_n]] = n
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum.not.reached.decision = cumsum.not.reached.decision[-reached.decision_n]
        cumSS.not.reached.decision = cumSS.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n==nmax)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N), mean(log(BF))),
         N, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 6, byrow = T)
    
    out.OCandASN$N = matrix(data = out.OCandASN$N, nrow = nReplicate,
                            ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN', 'avg.logBF')
  rownames(out.OCandASN$N) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### two-sample z ####
SBFNAP_twoZ = function(es = c(0, .2, .3, .5), n1min = 1, n2min = 1,
                       n1max = 5000, n2max = 5000,
                       tau.NAP = 0.3/sqrt(2), sigma0 = 1,
                       RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                       batch1.size.increment, batch2.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
    # batch1.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
    # batch2.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.twosample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumsum1.not.reached.decision = cumsum2.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N1 = rep(n1max, nReplicate)
    N2 = rep(n2max, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1min,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2min,
                   byrow = T)
        }
        
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1.increment,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2.increment,
                   byrow = T)
        }
      }
      
      # constant terms in BF
      m_n = (n1*n2)/(n1+n2)
      r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + rowSums(obs1.not.reached.decision)
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + rowSums(obs2.not.reached.decision)
      
      # xbar and S until step n
      xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
      xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(m_n)*(xbar2.not.reached.decision - xbar1.not.reached.decision))/sigma0
      W.not.reached.decision = (r_n*(test.statistic.not.reached.decision^2))/2
      BF.not.reached.decision = ((m_n*(tau.NAP^2) + 1)^(-3/2))*
        (1 + 2*W.not.reached.decision)*exp(W.not.reached.decision)
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N1[not.reached.decision[reached.decision_n]] = n1
        N2[not.reached.decision[reached.decision_n]] = n2
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum1.not.reached.decision = cumsum1.not.reached.decision[-reached.decision_n]
        cumsum2.not.reached.decision = cumsum2.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n1==n1max)||(n2==n2max)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N1), mean(N2), mean(log(BF))),
         N1, N2, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N1', 'N2', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 7, byrow = T)
    
    out.OCandASN$N1 = matrix(data = out.OCandASN$N1, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$N2 = matrix(data = out.OCandASN$N2, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN1', 'ASN2', 'avg.logBF')
  rownames(out.OCandASN$N1) = rownames(out.OCandASN$N2) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### two-sample t ####
SBFNAP_twoT = function(es = c(0, .2, .3, .5), n1min = 2, n2min = 2,
                       n1max = 5000, n2max = 5000,
                       tau.NAP = 0.3/sqrt(2),
                       RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                       batch1.size.increment, batch2.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
    # batch1.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
    # batch2.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.twosample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumSS1.not.reached.decision = cumSS2.not.reached.decision = 
      cumsum1.not.reached.decision = cumsum2.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N1 = rep(n1max, nReplicate)
    N2 = rep(n2max, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1min,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2min,
                   byrow = T)
        }
        
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1.increment,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2.increment,
                   byrow = T)
        }
      }
      
      # constant terms in BF
      m_n = (n1*n2)/(n1+n2)
      r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
      q_n = (r_n*(n1+n2-1))/(n1+n2-2)
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + rowSums(obs1.not.reached.decision)
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + rowSums(obs2.not.reached.decision)
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + rowSums(obs1.not.reached.decision^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + rowSums(obs2.not.reached.decision^2)
      
      # xbar and S until step n
      xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
      xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
      pooleds.not.reached.decision =
        sqrt((cumSS1.not.reached.decision - n1*(xbar1.not.reached.decision^2) +
                cumSS2.not.reached.decision - n2*(xbar2.not.reached.decision^2))/
               (n1+n2-2))
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(m_n)*(xbar2.not.reached.decision - xbar1.not.reached.decision))/
        pooleds.not.reached.decision
      G.not.reached.decision = 1 + (test.statistic.not.reached.decision^2)/(n1+n2-2)
      H.not.reached.decision = 1 + ((1-r_n)*(test.statistic.not.reached.decision^2))/(n1+n2-2)
      BF.not.reached.decision = ((m_n*(tau.NAP^2) + 1)^(-3/2))*
        ((G.not.reached.decision/H.not.reached.decision)^((n1+n2-1)/2))*
        (1 + (q_n*(test.statistic.not.reached.decision^2))/
           H.not.reached.decision)
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N1[not.reached.decision[reached.decision_n]] = n1
        N2[not.reached.decision[reached.decision_n]] = n2
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum1.not.reached.decision = cumsum1.not.reached.decision[-reached.decision_n]
        cumsum2.not.reached.decision = cumsum2.not.reached.decision[-reached.decision_n]
        cumSS1.not.reached.decision = cumSS1.not.reached.decision[-reached.decision_n]
        cumSS2.not.reached.decision = cumSS2.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n1==n1max)||(n2==n2max)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N1), mean(N2), mean(log(BF))),
         N1, N2, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N1', 'N2', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 7, byrow = T)
    
    out.OCandASN$N1 = matrix(data = out.OCandASN$N1, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$N2 = matrix(data = out.OCandASN$N2, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN1', 'ASN2', 'avg.logBF')
  rownames(out.OCandASN$N1) = rownames(out.OCandASN$N2) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### implement SBF+NAP for an observed data ####
#### one-sample z ####
implement.SBFNAP_oneZ = function(obs, sigma0 = 1, tau.NAP = 0.3/sqrt(2),
                                 RejectH1.threshold = exp(-3), 
                                 RejectH0.threshold = exp(3),
                                 batch.size, return.plot = T,
                                 until.decision.reached = T){
  
  if(missing(batch.size)) batch.size = rep(1, length(obs))
  
  obs.not.reached.decision = obs
  
  nAnalyses = min(which(length(obs)<=cumsum(batch.size)))
  
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS.not.reached.decision = cumsum.not.reached.decision = 0
  decision = 'I'
  N = length(obs)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[1:n])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[1:n]^2)
      
    }else{
      
      n.increment = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)]^2)
      
      n = n + n.increment
    }
    
    # constant terms in BF
    r_n = (n*(tau.NAP^2))/(n*(tau.NAP^2) + 1)
    
    # xbar and S until step n
    xbar.not.reached.decision = cumsum.not.reached.decision/n
    
    # BF at step n for those not reached decision
    test.statistic.not.reached.decision = 
      (sqrt(n)*xbar.not.reached.decision)/sigma0
    W.not.reached.decision = (r_n*(test.statistic.not.reached.decision^2))/2
    
    BF = c(BF,
           ((n*(tau.NAP^2) + 1)^(-3/2))*
             (1 + 2*W.not.reached.decision)*exp(W.not.reached.decision))
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N = n
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    if(decision=='A'){
      
      decision.plot = 'Accept the null'
      
    }else if(decision=='R'){
      
      decision.plot = 'Reject the null'
      
    }else if(decision=='I'){
      
      decision.plot = 'Inconclusive'
    }
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2, ylim = plot.range,
         main = paste('One-sample z-test'),
         sub = paste('Decision: ', decision.plot, 
                     ', N = ', N, sep = ''),
         xlab = 'Sequential order', ylab = 'Log(Bayes factor)')
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N' = N, 'BF' = BF, 'decision' = decision))
}

#### one-sample t ####
implement.SBFNAP_oneT = function(obs, tau.NAP = 0.3/sqrt(2),
                                 RejectH1.threshold = exp(-3), 
                                 RejectH0.threshold = exp(3),
                                 batch.size, return.plot = T,
                                 until.decision.reached = T){
  
  if(missing(batch.size)){
    
    batch.size = c(2, rep(1, length(obs)-2))
    
  }else{
    
    if(batch.size[1]<2) return("batch.size[1] (the size of the first batch) needs to be at least 2.")
  }
  
  obs.not.reached.decision = obs
  
  nAnalyses = min(which(length(obs)<=cumsum(batch.size)))
  
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS.not.reached.decision = cumsum.not.reached.decision = 0
  decision = 'I'
  N = length(obs)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[1:n])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[1:n]^2)
      
    }else{
      
      n.increment = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)]^2)
      
      n = n + n.increment
    }
    
    # constant terms in BF
    r_n = (n*(tau.NAP^2))/(n*(tau.NAP^2) + 1)
    q_n = (r_n*n)/(n-1)
    
    # xbar and S until step n
    xbar.not.reached.decision = cumsum.not.reached.decision/n
    s.not.reached.decision = sqrt((cumSS.not.reached.decision - 
                                     n*(xbar.not.reached.decision^2))/(n-1))
    
    # BF at step n for those not reached decision
    test.statistic.not.reached.decision = 
      (sqrt(n)*xbar.not.reached.decision)/s.not.reached.decision
    G.not.reached.decision = 1 + (test.statistic.not.reached.decision^2)/(n-1)
    H.not.reached.decision = 1 + ((1-r_n)*(test.statistic.not.reached.decision^2))/(n-1)
    
    BF = c(BF,
           ((n*(tau.NAP^2) + 1)^(-3/2))*
             ((G.not.reached.decision/H.not.reached.decision)^(n/2))*
             (1 + (q_n*(test.statistic.not.reached.decision^2))/
                H.not.reached.decision))
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N = n
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    if(decision=='A'){
      
      decision.plot = 'Accept the null'
      
    }else if(decision=='R'){
      
      decision.plot = 'Reject the null'
      
    }else if(decision=='I'){
      
      decision.plot = 'Inconclusive'
    }
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2, ylim = plot.range,
         main = paste('One-sample t-test'),
         sub = paste('Decision: ', decision.plot, 
                     ', N = ', N, sep = ''),
         xlab = 'Sequential order', ylab = 'Log(Bayes factor)')
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N' = N, 'BF' = BF, 'decision' = decision))
}

#### two-sample z ####
implement.SBFNAP_twoZ = function(obs1, obs2, sigma0 = 1, tau.NAP = 0.3/sqrt(2),
                                 RejectH1.threshold = exp(-3), 
                                 RejectH0.threshold = exp(3),
                                 batch1.size, batch2.size, return.plot = T,
                                 until.decision.reached = T){
  
  if(missing(batch1.size)) batch1.size = rep(1, length(obs1))
  if(missing(batch2.size)) batch2.size = rep(1, length(obs2))
  
  if(length(batch1.size)!=length(batch2.size)) return("Lengths of batch1.size and batch2.size needs to equal. They represent the number of steps in a sequential/group-sequential comparison.")
  
  obs1.not.reached.decision = obs1
  obs2.not.reached.decision = obs2
  
  nAnalyses = min(min(which(length(obs1)<=cumsum(batch1.size))),
                  min(which(length(obs2)<=cumsum(batch2.size))))
  
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS1.not.reached.decision = cumSS2.not.reached.decision = 
    cumsum1.not.reached.decision = cumsum2.not.reached.decision = 0
  decision = 'I'
  N1 = length(obs1)
  N2 = length(obs2)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n1 = batch1.size[seq.step]
      n2 = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[1:n1])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[1:n2])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[1:n1]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[1:n2]^2)
      
    }else{
      
      n1.increment = batch1.size[seq.step]
      n2.increment = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)]^2)
      
      n1 = n1 + n1.increment
      n2 = n2 + n2.increment
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    
    # xbar and S until step n
    xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
    xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
    
    # BF at step n for those not reached decision
    test.statistic.not.reached.decision = 
      (sqrt(m_n)*(xbar1.not.reached.decision - xbar2.not.reached.decision))/sigma0
    W.not.reached.decision = (r_n*(test.statistic.not.reached.decision^2))/2
    
    BF = c(BF,
           ((m_n*(tau.NAP^2) + 1)^(-3/2))*
             (1 + 2*W.not.reached.decision)*exp(W.not.reached.decision))
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N1 = n1
      N2 = n2
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    if(decision=='A'){
      
      decision.plot = 'Accept the null'
      
    }else if(decision=='R'){
      
      decision.plot = 'Reject the null'
      
    }else if(decision=='I'){
      
      decision.plot = 'Inconclusive'
    }
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2, ylim = plot.range,
         main = paste('Two-sample z-test'),
         sub = paste('Decision: ', decision.plot, 
                     ', N1 = ', N1, ', N2 = ', N2, sep = ''),
         xlab = 'Sequential order', ylab = 'Log(Bayes factor)')
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N1' = N1, 'N2' = N2, 'BF' = BF, 'decision' = decision))
}

#### two-sample t ####
implement.SBFNAP_twoT = function(obs1, obs2, tau.NAP = 0.3/sqrt(2),
                                 RejectH1.threshold = exp(-3), 
                                 RejectH0.threshold = exp(3),
                                 batch1.size, batch2.size, return.plot = T,
                                 until.decision.reached = T){
  
  if(missing(batch1.size)){
    
    batch1.size = c(2, rep(1, length(obs1)-2))
    
  }else{
    
    if(batch1.size[1]<2) return("batch1.size[1] (the size of the first batch from Group-1) needs to be at least 2.")
  }
  
  if(missing(batch2.size)){
    
    batch2.size = c(2, rep(1, length(obs2)-2))
    
  }else{
    
    if(batch2.size[1]<2) return("batch2.size[1] (the size of the first batch from Group-2) needs to be at least 2.")
  }
  
  if(length(batch1.size)!=length(batch2.size)) return("Lengths of batch1.size and batch2.size needs to equal. They represent the number of steps in a sequential/group-sequential comparison.")
  
  obs1.not.reached.decision = obs1
  obs2.not.reached.decision = obs2
  
  nAnalyses = min(min(which(length(obs1)<=cumsum(batch1.size))),
                  min(which(length(obs2)<=cumsum(batch2.size))))
  
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS1.not.reached.decision = cumSS2.not.reached.decision = 
    cumsum1.not.reached.decision = cumsum2.not.reached.decision = 0
  decision = 'I'
  N1 = length(obs1)
  N2 = length(obs2)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n1 = batch1.size[seq.step]
      n2 = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[1:n1])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[1:n2])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[1:n1]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[1:n2]^2)
      
    }else{
      
      n1.increment = batch1.size[seq.step]
      n2.increment = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)]^2)
      
      n1 = n1 + n1.increment
      n2 = n2 + n2.increment
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    r_n = (m_n*(tau.NAP^2))/(m_n*(tau.NAP^2) + 1)
    q_n = (r_n*(n1+n2-1))/(n1+n2-2)
    
    # xbar and S until step n
    xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
    xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
    
    # BF at step n for those not reached decision
    pooleds.not.reached.decision =
      sqrt((cumSS1.not.reached.decision - n1*(xbar1.not.reached.decision^2) +
              cumSS2.not.reached.decision - n2*(xbar2.not.reached.decision^2))/
             (n1+n2-2))
    test.statistic.not.reached.decision = 
      (sqrt(m_n)*(xbar1.not.reached.decision - xbar2.not.reached.decision))/
      pooleds.not.reached.decision
    G.not.reached.decision = 1 + (test.statistic.not.reached.decision^2)/(n1+n2-2)
    H.not.reached.decision = 1 + ((1-r_n)*(test.statistic.not.reached.decision^2))/(n1+n2-2)
    
    BF = c(BF,
           ((m_n*(tau.NAP^2) + 1)^(-3/2))*
             ((G.not.reached.decision/H.not.reached.decision)^((n1+n2-1)/2))*
             (1 + (q_n*(test.statistic.not.reached.decision^2))/
                H.not.reached.decision))
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N1 = n1
      N2 = n2
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2,
         ylim = plot.range)
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N1' = N1, 'N2' = N2, 'BF' = BF, 'decision' = decision))
}


#### SBF+composite alternative ####
#### OC and ASN ####
#### one-sample z ####
SBFHajnal_oneZ = function(es = c(0, .2, .3, .5), es1 = 0.3, nmin = 1, nmax = 5000,
                          sigma0 = 1,
                          RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                          batch.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
    # batch.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.onesample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumSS.not.reached.decision = 
      cumsum.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N = rep(nmax, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:nmin,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:nmin, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = nmin,
                   byrow = T)
        }
        
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:n.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:n.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n.increment,
                   byrow = T)
        }
      }
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + rowSums(obs.not.reached.decision)
      
      # xbar and S until step n
      xbar.not.reached.decision = cumsum.not.reached.decision/n
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(n)*xbar.not.reached.decision)/sigma0
      BF.not.reached.decision = 
        (dnorm(x = test.statistic.not.reached.decision,
               mean = sqrt(n)*abs(es1))/
           dnorm(x = test.statistic.not.reached.decision) +
           dnorm(x = test.statistic.not.reached.decision,
                 mean = sqrt(n)*abs(es1))/
           dnorm(x = test.statistic.not.reached.decision))/2
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N[not.reached.decision[reached.decision_n]] = n
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum.not.reached.decision = cumsum.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n==nmax)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N), mean(log(BF))),
         N, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 6, byrow = T)
    
    out.OCandASN$N = matrix(data = out.OCandASN$N, nrow = nReplicate,
                            ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN', 'avg.logBF')
  rownames(out.OCandASN$N) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### one-sample t ####
SBFHajnal_oneT = function(es = c(0, .2, .3, .5), es1 = 0.3, nmin = 2, nmax = 5000,
                          RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                          batch.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch.size.increment)){
    
    batch.size.increment = function(narg){20}
    # batch.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.onesample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumSS.not.reached.decision = 
      cumsum.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N = rep(nmax, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n = nmin
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:nmin,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:nmin, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = nmin,
                   byrow = T)
        }
        
        
      }else{
        
        n.increment = min(batch.size.increment(n), nmax-n)
        n = n + n.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs.not.reached.decision = 
            mapply(X = 1:n.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs.not.reached.decision = 
            matrix(mapply(X = 1:n.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n.increment,
                   byrow = T)
        }
      }
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + rowSums(obs.not.reached.decision)
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + rowSums(obs.not.reached.decision^2)
      
      # xbar and S until step n
      xbar.not.reached.decision = cumsum.not.reached.decision/n
      s.not.reached.decision = sqrt((cumSS.not.reached.decision - 
                                       n*(xbar.not.reached.decision^2))/(n-1))
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(n)*xbar.not.reached.decision)/s.not.reached.decision
      BF.not.reached.decision = 
        df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n-1,
           ncp = n*(es1^2))/
        df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n-1)
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N[not.reached.decision[reached.decision_n]] = n
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum.not.reached.decision = cumsum.not.reached.decision[-reached.decision_n]
        cumSS.not.reached.decision = cumSS.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n==nmax)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N), mean(log(BF))),
         N, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = as.data.frame(matrix(out.OCandASN$summary,
                                                nrow = 1, ncol = 6, byrow = T))
    
    out.OCandASN$N = matrix(data = out.OCandASN$N, nrow = nReplicate,
                            ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN', 'avg.logBF')
  rownames(out.OCandASN$N) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### two-sample z ####
SBFHajnal_twoZ = function(es = c(0, .2, .3, .5), es1 = 0.3, n1min = 1, n2min = 1,
                          n1max = 5000, n2max = 5000, sigma0 = 1,
                          RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                          batch1.size.increment, batch2.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
    # batch1.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
    # batch2.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.twosample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumsum1.not.reached.decision = cumsum2.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N1 = rep(n1max, nReplicate)
    N2 = rep(n2max, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1min,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2min,
                   byrow = T)
        }
        
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1.increment,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, sigma0)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2.increment,
                   byrow = T)
        }
      }
      
      # constant terms in BF
      m_n = (n1*n2)/(n1+n2)
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + rowSums(obs1.not.reached.decision)
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + rowSums(obs2.not.reached.decision)
      
      # xbar and S until step n
      xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
      xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(m_n)*(xbar2.not.reached.decision - xbar1.not.reached.decision))/sigma0
      BF.not.reached.decision = 
        (dnorm(x = test.statistic.not.reached.decision,
               mean = sqrt(m_n)*abs(es1))/
           dnorm(x = test.statistic.not.reached.decision) +
           dnorm(x = test.statistic.not.reached.decision,
                 mean = -sqrt(m_n)*abs(es1))/
           dnorm(x = test.statistic.not.reached.decision))/2
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N1[not.reached.decision[reached.decision_n]] = n1
        N2[not.reached.decision[reached.decision_n]] = n2
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum1.not.reached.decision = cumsum1.not.reached.decision[-reached.decision_n]
        cumsum2.not.reached.decision = cumsum2.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n1==n1max)||(n2==n2max)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N1), mean(N2), mean(log(BF))),
         N1, N2, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N1', 'N2', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 7, byrow = T)
    
    out.OCandASN$N1 = matrix(data = out.OCandASN$N1, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$N2 = matrix(data = out.OCandASN$N2, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN1', 'ASN2', 'avg.logBF')
  rownames(out.OCandASN$N1) = rownames(out.OCandASN$N2) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### two-sample t ####
SBFHajnal_twoT = function(es = c(0, .2, .3, .5), es1 = 0.3, n1min = 2, n2min = 2,
                          n1max = 5000, n2max = 5000,
                          RejectH1.threshold = exp(-3), RejectH0.threshold = exp(3),
                          batch1.size.increment, batch2.size.increment, nReplicate = 5e+4, nCore){
  
  if(missing(batch1.size.increment)){
    
    batch1.size.increment = function(narg){20}
    # batch1.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(batch2.size.increment)){
    
    batch2.size.increment = function(narg){20}
    # batch2.size.increment = function(narg){
    #   
    #   if(narg<100){
    #     return(1)
    #   }else if(narg<1000){
    #     return(5)
    #   }else if(narg<2500){
    #     return(10)
    #   }else if(narg<5000){
    #     return(20)
    #   }else{
    #     return(50)
    #   }
    # }
  }
  if(missing(nCore)) nCore = max(1, parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cores = nCore)
  out.OCandASN = foreach(es.gen = es, .combine = 'mycombine.seq.twosample', .multicombine = T) %dopar% {
    
    ## simulating data and calculating the BF
    
    # required storages
    cumSS1.not.reached.decision = cumSS2.not.reached.decision = 
      cumsum1.not.reached.decision = cumsum2.not.reached.decision = numeric(nReplicate)
    decision = rep('A', nReplicate)
    N1 = rep(n1max, nReplicate)
    N2 = rep(n2max, nReplicate)
    BF = rep(NA, nReplicate)
    not.reached.decision = 1:nReplicate
    nNot.reached.decision = nReplicate
    
    seq.step = 0
    terminate = F
    while(!terminate){
      
      # tracking sequential step
      seq.step = seq.step + 1
      
      # sample size used at this step
      if(seq.step==1){
        
        n1 = n1min
        n2 = n2min
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2min,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1min,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2min, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2min,
                   byrow = T)
        }
        
        
      }else{
        
        n1.increment = min(batch1.size.increment(n1), n1max-n1)
        n1 = n1 + n1.increment
        n2.increment = min(batch2.size.increment(n2), n2max-n2)
        n2 = n2 + n2.increment
        
        # simulating obs at this step
        set.seed(seq.step)
        if(nNot.reached.decision>1){
          
          obs1.not.reached.decision = 
            mapply(X = 1:n1.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                   })
          
          obs2.not.reached.decision = 
            mapply(X = 1:n2.increment,
                   FUN = function(X){
                     
                     rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                   })
          
        }else{
          
          obs1.not.reached.decision = 
            matrix(mapply(X = 1:n1.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, -es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n1.increment,
                   byrow = T)
          
          obs2.not.reached.decision = 
            matrix(mapply(X = 1:n2.increment, 
                          FUN = function(X){
                            
                            rnorm(nReplicate, es.gen/2, 1)[not.reached.decision]
                            
                          }), nrow = 1, ncol = n2.increment,
                   byrow = T)
        }
      }
      
      # constant terms in BF
      m_n = (n1*n2)/(n1+n2)
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + rowSums(obs1.not.reached.decision)
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + rowSums(obs2.not.reached.decision)
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + rowSums(obs1.not.reached.decision^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + rowSums(obs2.not.reached.decision^2)
      
      # xbar and S until step n
      xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
      xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
      pooleds.not.reached.decision =
        sqrt((cumSS1.not.reached.decision - n1*(xbar1.not.reached.decision^2) +
                cumSS2.not.reached.decision - n2*(xbar2.not.reached.decision^2))/
               (n1+n2-2))
      
      # BF at step n for those not reached decision
      test.statistic.not.reached.decision = 
        (sqrt(m_n)*(xbar2.not.reached.decision - xbar1.not.reached.decision))/
        pooleds.not.reached.decision
      BF.not.reached.decision = 
        df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n1+n2-2,
           ncp = m_n*((es1)^2))/
        df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n1+n2-2)
      
      # comparing with the thresholds
      AcceptedH0_n = which(BF.not.reached.decision<=RejectH1.threshold)
      RejectedH0_n = which(BF.not.reached.decision>=RejectH0.threshold)
      reached.decision_n = union(AcceptedH0_n, RejectedH0_n)
      
      # tracking those reaching/not reaching a decision at step n
      if(length(reached.decision_n)>0){
        
        # stopping BF, time and decision
        BF[not.reached.decision[reached.decision_n]] = BF.not.reached.decision[reached.decision_n]
        N1[not.reached.decision[reached.decision_n]] = n1
        N2[not.reached.decision[reached.decision_n]] = n2
        decision[not.reached.decision[RejectedH0_n]] = 'R'
        
        # tracking only those not reached decision yet
        cumsum1.not.reached.decision = cumsum1.not.reached.decision[-reached.decision_n]
        cumsum2.not.reached.decision = cumsum2.not.reached.decision[-reached.decision_n]
        cumSS1.not.reached.decision = cumSS1.not.reached.decision[-reached.decision_n]
        cumSS2.not.reached.decision = cumSS2.not.reached.decision[-reached.decision_n]
        not.reached.decision = not.reached.decision[-reached.decision_n]
        nNot.reached.decision = length(not.reached.decision)
      }
      
      terminate = (nNot.reached.decision==0)||(n1==n1max)||(n2==n2max)
    }
    
    # inconclusive
    if(nNot.reached.decision>0){
      
      decision[not.reached.decision] = 'I'
      
      if(length(reached.decision_n)>0){
        
        # BF at step n for those not reached decision
        BF.not.reached.decision = BF.not.reached.decision[-reached.decision_n]
      }
      
      # BF
      BF[not.reached.decision] = BF.not.reached.decision
    }
    
    decision.freq = table(decision)
    
    list(c(es.gen, 
           ifelse(is.na(decision.freq['A']), 0, as.numeric(decision.freq['A'])/nReplicate),
           ifelse(is.na(decision.freq['R']), 0, as.numeric(decision.freq['R'])/nReplicate),
           ifelse(is.na(decision.freq['I']), 0, as.numeric(decision.freq['I'])/nReplicate),
           mean(N1), mean(N2), mean(log(BF))),
         N1, N2, BF)
  }
  
  names(out.OCandASN) = c('summary', 'N1', 'N2', 'BF')
  
  if(length(es)==1){
    
    out.OCandASN$summary = matrix(out.OCandASN$summary,
                                  nrow = 1, ncol = 7, byrow = T)
    
    out.OCandASN$N1 = matrix(data = out.OCandASN$N1, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$N2 = matrix(data = out.OCandASN$N2, nrow = nReplicate,
                             ncol = 1, byrow = T)
    
    out.OCandASN$BF = matrix(data = out.OCandASN$BF, nrow = nReplicate,
                             ncol = 1, byrow = T)
  }
  
  out.OCandASN$summary = as.data.frame(out.OCandASN$summary)
  colnames(out.OCandASN$summary) = c('effect.size', 'acceptH0', 'rejectH0', 
                                     'inconclusive', 'ASN1', 'ASN2', 'avg.logBF')
  rownames(out.OCandASN$N1) = rownames(out.OCandASN$N2) = rownames(out.OCandASN$BF) = 
    as.character(es)
  
  return(out.OCandASN)
}

#### implementation for an observed data ####
#### one-sample z ####
implement.SBFHajnal_oneZ = function(obs, es1 = 0.3, sigma0 = 1,
                                    RejectH1.threshold = exp(-3), 
                                    RejectH0.threshold = exp(3),
                                    batch.size, return.plot = T,
                                    until.decision.reached = T){
  
  if(missing(batch.size)) batch.size = rep(1, length(obs))
  
  obs.not.reached.decision = obs
  
  nAnalyses = min(which(length(obs)<=cumsum(batch.size)))
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS.not.reached.decision = cumsum.not.reached.decision = 0
  decision = 'I'
  N = length(obs)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[1:n])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[1:n]^2)
      
    }else{
      
      n.increment = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)]^2)
      
      n = n + n.increment
    }
    
    # xbar and S until step n
    xbar.not.reached.decision = cumsum.not.reached.decision/n
    
    # BF at step n for those not reached decision
    test.statistic.not.reached.decision = 
      (sqrt(n)*xbar.not.reached.decision)/sigma0
    
    BF = c(BF,
           (dnorm(x = test.statistic.not.reached.decision,
                  mean = sqrt(n)*abs(es1))/
              dnorm(x = test.statistic.not.reached.decision) +
              dnorm(x = test.statistic.not.reached.decision,
                    mean = -sqrt(n)*abs(es1))/
              dnorm(x = test.statistic.not.reached.decision))/2)
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N = n
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    if(decision=='A'){
      
      decision.plot = 'Accept the null'
      
    }else if(decision=='R'){
      
      decision.plot = 'Reject the null'
      
    }else if(decision=='I'){
      
      decision.plot = 'Inconclusive'
    }
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2, ylim = plot.range,
         main = paste('One-sample z-test'),
         sub = paste('Decision: ', decision.plot, 
                     ', N = ', N, sep = ''),
         xlab = 'Sequential order', ylab = 'Log(Bayes factor)')
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N' = N, 'BF' = BF, 'decision' = decision))
}

#### one-sample t ####
implement.SBFHajnal_oneT = function(obs, es1 = 0.3,
                                    RejectH1.threshold = exp(-3), 
                                    RejectH0.threshold = exp(3),
                                    batch.size, return.plot = T,
                                    until.decision.reached = T){
  
  if(missing(batch.size)){
    
    batch.size = c(2, rep(1, length(obs)-2))
    
  }else{
    
    if(batch.size[1]<2) return("batch.size[1] (the size of the first batch) needs to be at least 2.")
  }
  
  obs.not.reached.decision = obs
  
  nAnalyses = min(which(length(obs)<=cumsum(batch.size)))
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS.not.reached.decision = cumsum.not.reached.decision = 0
  decision = 'I'
  N = length(obs)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[1:n])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[1:n]^2)
      
    }else{
      
      n.increment = batch.size[seq.step]
      
      # sum of observations until step n
      cumsum.not.reached.decision = 
        cumsum.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)])
      
      # sum of squares of observations until step n
      cumSS.not.reached.decision = 
        cumSS.not.reached.decision + sum(obs.not.reached.decision[(n+1):(n+n.increment)]^2)
      
      n = n + n.increment
    }
    
    # xbar and S until step n
    xbar.not.reached.decision = cumsum.not.reached.decision/n
    s.not.reached.decision = sqrt((cumSS.not.reached.decision - 
                                     n*(xbar.not.reached.decision^2))/(n-1))
    
    # BF at step n for those not reached decision
    test.statistic.not.reached.decision = 
      (sqrt(n)*xbar.not.reached.decision)/s.not.reached.decision
    
    BF = c(BF,
           df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n-1,
              ncp = n*((es1)^2))/
             df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n-1))
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N = n
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    if(decision=='A'){
      
      decision.plot = 'Accept the null'
      
    }else if(decision=='R'){
      
      decision.plot = 'Reject the null'
      
    }else if(decision=='I'){
      
      decision.plot = 'Inconclusive'
    }
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2, ylim = plot.range,
         main = paste('One-sample t-test'),
         sub = paste('Decision: ', decision.plot, 
                     ', N = ', N, sep = ''),
         xlab = 'Sequential order', ylab = 'Log(Bayes factor)')
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N' = N, 'BF' = BF, 'decision' = decision))
}

#### two-sample z ####
implement.SBFHajnal_twoZ = function(obs1, obs2, es1 = 0.3, sigma0 = 1,
                                    RejectH1.threshold = exp(-3), 
                                    RejectH0.threshold = exp(3),
                                    batch1.size, batch2.size, return.plot = T,
                                    until.decision.reached = T){
  
  if(missing(batch1.size)) batch1.size = rep(1, length(obs1))
  if(missing(batch2.size)) batch2.size = rep(1, length(obs2))
  
  if(length(batch1.size)!=length(batch2.size)) return("Lengths of batch1.size and batch2.size needs to equal. They represent the number of steps in a sequential/group-sequential comparison.")
  
  obs1.not.reached.decision = obs1
  obs2.not.reached.decision = obs2
  
  nAnalyses = min(min(which(length(obs1)<=cumsum(batch1.size))),
                  min(which(length(obs2)<=cumsum(batch2.size))))
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS1.not.reached.decision = cumSS2.not.reached.decision = 
    cumsum1.not.reached.decision = cumsum2.not.reached.decision = 0
  decision = 'I'
  N1 = length(obs1)
  N2 = length(obs2)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n1 = batch1.size[seq.step]
      n2 = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[1:n1])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[1:n2])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[1:n1]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[1:n2]^2)
      
    }else{
      
      n1.increment = batch1.size[seq.step]
      n2.increment = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)]^2)
      
      n1 = n1 + n1.increment
      n2 = n2 + n2.increment
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    
    # xbar and S until step n
    xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
    xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
    
    # BF at step n for those not reached decision
    test.statistic.not.reached.decision = 
      (sqrt(m_n)*(xbar1.not.reached.decision - xbar2.not.reached.decision))/sigma0
    
    BF = c(BF,
           (dnorm(x = test.statistic.not.reached.decision,
                  mean = sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic.not.reached.decision) +
              dnorm(x = test.statistic.not.reached.decision,
                    mean = -sqrt(m_n)*abs(es1))/
              dnorm(x = test.statistic.not.reached.decision))/2)
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N1 = n1
      N2 = n2
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    if(decision=='A'){
      
      decision.plot = 'Accept the null'
      
    }else if(decision=='R'){
      
      decision.plot = 'Reject the null'
      
    }else if(decision=='I'){
      
      decision.plot = 'Inconclusive'
    }
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2, ylim = plot.range,
         main = paste('Two-sample z-test'),
         sub = paste('Decision: ', decision.plot, 
                     ', N1 = ', N1, ', N2 = ', N2, sep = ''),
         xlab = 'Sequential order', ylab = 'Log(Bayes factor)')
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N1' = N1, 'N2' = N2, 'BF' = BF, 'decision' = decision))
}

#### two-sample t ####
implement.SBFHajnal_twoT = function(obs1, obs2, es1 = 0.3,
                                    RejectH1.threshold = exp(-3), 
                                    RejectH0.threshold = exp(3),
                                    batch1.size, batch2.size, return.plot = T,
                                    until.decision.reached = T){
  
  if(missing(batch1.size)){
    
    batch1.size = c(2, rep(1, length(obs1)-2))
    
  }else{
    
    if(batch1.size[1]<2) return("batch1.size[1] (the size of the first batch from Group-1) needs to be at least 2.")
  }
  
  if(missing(batch2.size)){
    
    batch2.size = c(2, rep(1, length(obs2)-2))
    
  }else{
    
    if(batch2.size[1]<2) return("batch2.size[1] (the size of the first batch from Group-2) needs to be at least 2.")
  }
  
  if(length(batch1.size)!=length(batch2.size)) return("Lengths of batch1.size and batch2.size needs to equal. They represent the number of steps in a sequential/group-sequential comparison.")
  
  obs1.not.reached.decision = obs1
  obs2.not.reached.decision = obs2
  
  nAnalyses = min(min(which(length(obs1)<=cumsum(batch1.size))),
                  min(which(length(obs2)<=cumsum(batch2.size))))
  
  ## random shuffling data and calculating the BF
  # required storages
  cumSS1.not.reached.decision = cumSS2.not.reached.decision = 
    cumsum1.not.reached.decision = cumsum2.not.reached.decision = 0
  decision = 'I'
  N1 = length(obs1)
  N2 = length(obs2)
  BF = NULL
  terminate = F
  
  seq.step = 0
  while(!terminate){
    
    # tracking sequential step
    seq.step = seq.step + 1
    
    # sample size used at this step
    if(seq.step==1){
      
      n1 = batch1.size[seq.step]
      n2 = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[1:n1])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[1:n2])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[1:n1]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[1:n2]^2)
      
    }else{
      
      n1.increment = batch1.size[seq.step]
      n2.increment = batch2.size[seq.step]
      
      # sum of observations until step n
      cumsum1.not.reached.decision = 
        cumsum1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)])
      cumsum2.not.reached.decision = 
        cumsum2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)])
      
      # sum of squares of observations until step n
      cumSS1.not.reached.decision = 
        cumSS1.not.reached.decision + sum(obs1.not.reached.decision[(n1+1):(n1+n1.increment)]^2)
      cumSS2.not.reached.decision = 
        cumSS2.not.reached.decision + sum(obs2.not.reached.decision[(n2+1):(n2+n2.increment)]^2)
      
      n1 = n1 + n1.increment
      n2 = n2 + n2.increment
    }
    
    # constant terms in BF
    m_n = (n1*n2)/(n1+n2)
    
    # xbar and S until step n
    xbar1.not.reached.decision = cumsum1.not.reached.decision/n1
    xbar2.not.reached.decision = cumsum2.not.reached.decision/n2
    
    # BF at step n for those not reached decision
    pooleds.not.reached.decision =
      sqrt((cumSS1.not.reached.decision - n1*(xbar1.not.reached.decision^2) +
              cumSS2.not.reached.decision - n2*(xbar2.not.reached.decision^2))/
             (n1+n2-2))
    test.statistic.not.reached.decision = 
      (sqrt(m_n)*(xbar1.not.reached.decision - xbar2.not.reached.decision))/
      pooleds.not.reached.decision
    
    BF = c(BF,
           df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n1+n2-2,
              ncp = m_n*((es1)^2))/
             df(x = test.statistic.not.reached.decision^2, df1 = 1, df2 = n1+n2-2))
    
    # comparing with the thresholds
    AcceptedH0_n = BF[seq.step]<=RejectH1.threshold
    RejectedH0_n = BF[seq.step]>=RejectH0.threshold
    reached.decision_n = AcceptedH0_n|RejectedH0_n
    
    # tracking those reaching/not reaching a decision at step n
    if(reached.decision_n){
      
      # stopping BF, time and decision
      N1 = n1
      N2 = n2
      if(AcceptedH0_n) decision = 'A'
      if(RejectedH0_n) decision = 'R'
    }
    
    if(until.decision.reached){
      
      terminate = reached.decision_n||(seq.step==nAnalyses)
      
    }else{terminate = seq.step==nAnalyses}
  }
  
  # sequential comparison plot
  if(return.plot){
    
    log.BF = log(BF)
    plot.range = range(range(log.BF), log(RejectH1.threshold),
                       log(RejectH0.threshold))
    plot(log.BF, type = 'l', lwd = 2,
         ylim = plot.range)
    points(log.BF, pch = 16)
    abline(h = log(RejectH1.threshold), lwd = 2, col = 2)
    abline(h = log(RejectH0.threshold), lwd = 2, col = 2)
  }
  
  return(list('N1' = N1, 'N2' = N2, 'BF' = BF, 'decision' = decision))
}

