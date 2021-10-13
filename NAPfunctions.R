

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
  
  abs.esdiff = abs(es1)
  
  if(!missing(obs)){
    
    nObs = length(obs)
    mean.obs = mean(obs)
    
    test.statistic = (sqrt(nObs)*mean.obs)/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = (sqrt(nObs)*(abs.esdiff))/sigma0)/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = (sqrt(nObs)*(-abs.esdiff))/sigma0)/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(nObs))&&(!missing(mean.obs))){
    
    test.statistic = (sqrt(nObs)*mean.obs)/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = (sqrt(nObs)*(abs.esdiff))/sigma0)/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = (sqrt(nObs)*(-abs.esdiff))/sigma0)/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(nObs))&&(!missing(test.statistic))){
    
    return((dnorm(x = test.statistic,
                  mean = (sqrt(nObs)*(abs.esdiff))/sigma0)/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = (sqrt(nObs)*(-abs.esdiff))/sigma0)/
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
  
  abs.esdiff = abs(es1)
  
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
                  mean = (sqrt(m_n)*(abs.esdiff))/sigma0)/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = (sqrt(m_n)*(-abs.esdiff))/sigma0)/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(mean.obs1))&&(!missing(mean.obs2))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    # BF at step n for those not reached decision
    test.statistic = (sqrt(m_n)*(mean.obs2 - mean.obs1))/sigma0
    
    return((dnorm(x = test.statistic,
                  mean = (sqrt(m_n)*(abs.esdiff))/sigma0)/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = (sqrt(m_n)*(-abs.esdiff))/sigma0)/
              dnorm(x = test.statistic))/2)
    
  }else if((!missing(n1Obs))&&(!missing(n2Obs))&&
           (!missing(test.statistic))){
    
    # constant terms in BF
    m_n = (n1Obs*n2Obs)/(n1Obs+n2Obs)
    
    return((dnorm(x = test.statistic,
                  mean = (sqrt(m_n)*(abs.esdiff))/sigma0)/
              dnorm(x = test.statistic) +
              dnorm(x = test.statistic,
                    mean = (sqrt(m_n)*(-abs.esdiff))/sigma0)/
              dnorm(x = test.statistic))/2)
  }
}

#### two-sample t ####
HajnalBF_twoT = function(obs1, obs2, n1Obs, n2Obs, 
                         mean.obs1, mean.obs2, sd.obs1, sd.obs2, test.statistic,
                         es1 = 0.3){
  
  abs.esdiff = abs(es1)
  
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

