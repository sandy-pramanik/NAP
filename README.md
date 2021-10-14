
# Non-local Alternative Priors (NAP) for Bayesian Hypothesis Tests in Psychology

# Article: Sandipan Pramanik and Valen E. Johnson (2021+)


## Description

Bayesian hypothesis testing procedures have gained increased acceptance in recent years.  A key advantage that Bayesian tests have over classical testing procedures is their potential to quantify information in support of true null hypotheses.  Ironically, default implementations of Bayesian tests prevent accumulation of strong evidence in favor of true null hypotheses because associated default alternative hypotheses assign high probability to data that are most consistent with a null effect. We propose the use of "non-local" alternative hypotheses to resolve this paradox. The resulting class of Bayesian hypothesis tests permit more rapid accumulation of evidence in favor of both true null hypotheses and alternative hypotheses that are compatible with standardized effect sizes of most interest in psychology. The functions available computes

* the Bayes factor in favor of the normal moment prior (the proposed NAP prior)
* the Hajnal's ratio (the ratio of the sampling distribution of the test-statistic) in favor of the alternative. The composite alternative places equal probabilities at two point masses symmetric about the point null.
             

## List of functions

The followings are the available functions in `NAPfunctions.R`:
  
* [**NAP**](#bf-nap)
  * [**NAPBF_oneZ**](#napbf_onez)
  * [**NAPBF_oneT**](#napbf_onet)
  * [**NAPBF_twoZ**](#napbf_twoz)
  * [**NAPBF_twoT**](#napbf_twot)
* [**Hajnal**](#bf-hajnal)
  * [**HajnalBF_oneZ**](#hajnalbf_onez)
  * [**HajnalBF_oneT**](#hajnalbf_onet)
  * [**HajnalBF_twoZ**](#hajnalbf_twoz)
  * [**HajnalBF_twoT**](#hajnalbf_twot)
             

## User's guide to the functions

### **NAPBF_oneZ**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In a ```math N(\mu,\sigma_0^2)``` population with known variance $\sigma_0^2$, consider the two-sided one-sample $z$-test for testing the point null hypothesis $H_0 : \mu = 0$ against $H_1 : \mu \neq 0$. Based on an observed data, this function calculates the Bayes factor in favor of $H_1$ when a *normal moment prior* is assumed on $\mu$ under the alternative [[Johnson & Rossell (2010)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2009.00730.x)].

* **Input:**
  
  * **obs:** Numeric vector. Observed vector of data.
  
  * **nObs:** Numeric or numeric vector. Sample size(s). Same as `length(obs)` when numeric.
  
  * **mean.obs:** Numeric or numeric vector. Sample mean(s). Same as `mean(obs)` when numeric.
    
  * **tau.NAP:** Positive numeric. Parameter in the moment prior. **Default:** $0.3/\sqrt{2}$. This places the prior modes of $\mu$ at $\pm0.3$.
  
  * **sigma0:** Positive numeric. Known standard deviation in the population. **Default:** 1.

* **Details:** 
  * Users can either specify `obs`, or `nObs` and `mean.obs` in the order of usage priorities.
  * If `obs` is provided, it returns the corresponding Bayes factor value.
  * If `nObs` and `mean.obs` are provided, the function is vectorized over both of the arguments.
    * If `nObs` is a numeric value and `mean.obs` is a numeric vector, it returns a numeric vector of the same length as `mean.obs`. These are Bayes factor values corresponding to different values in `mean.obs` for sample size `nObs`.
    * If `nObs` and `mean.obs` are both numeric vectors, both must be of the same length. In this case, it returns a numeric vector of the same length as `nObs` and `mean.obs`. These are Bayes factor values corresponding to each value in `nObs` and `mean.obs`.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **NAPBF_oneT**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In a $N(\mu,\sigma^2)$ population with unknown variance $\sigma^2$, consider the two-sided one-sample $t$-test for testing the point null hypothesis $H_0 : \mu = 0$ against $H_1 : \mu \neq 0$. Based on an observed data, this function calculates the Bayes factor in favor of $H_1$ when a *normal moment prior* is assumed on the standardized effect size $(\mu/\sigma)$ under the alternative [[Johnson & Rossell (2010)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2009.00730.x)]. Under both hypotheses, the Jeffrey's prior $\pi(\sigma^2) \propto 1/\sigma^2$ is assumed on $\sigma^2$.

* **Input:**
  
  * **obs:** Numeric vector. Observed vector of data.
  
  * **nObs:** Numeric or numeric vector. Sample size(s). Same as `length(obs)` when numeric.
  
  * **mean.obs:** Numeric or numeric vector. Sample mean(s). Same as `mean(obs)` when numeric.
  
  * **sd.obs:** Positive numeric or numeric vector. Sample standard deviation(s). Same as `sd(obs)` when numeric.
  
  * **test.statistic:** Numeric or numeric vector. Test-statistic value(s).
    
  * **tau.NAP:** Positive numeric. Parameter in the moment prior. **Default:** $0.3/\sqrt{2}$. This places the prior modes of $\mu$ at $\pm0.3 \sigma$.

* **Details:** 
  * A user can either specify `obs`; or `nObs`, `mean.obs` and `sd.obs`; or `nObs` and `test.statistic` in the order of usage priorities.
  * If `obs` is provided, it returns the corresponding Bayes factor value.
  * If `nObs`, `mean.obs` and `sd.obs` are provided, the function is vectorized over these arguments.
    * If `nObs` is a numeric value, and `mean.obs` and `sd.obs` are both numeric vectors, it returns a numeric vector of the same length as `mean.obs` and `sd.obs`. These are Bayes factor values corresponding to different values in `mean.obs` and `sd.obs` for sample size `nObs`.
    * If `nObs`, `mean.obs` and `sd.obs` are all numeric vectors, the function returns a numeric vector of the same length as `nObs`, `mean.obs` and `sd.obs`. These are Bayes factor values corresponding to the different values in `nObs`, `mean.obs and `sd.obs`.
  * If `nObs` and `test.statistic` are provided, the function is vectorized over both the arguments.
    * If `nObs` is a numeric value and `test.statistic` is a numeric vector, it returns a numeric vector of the same length as `test.statistic`. These are Bayes factor values corresponding to different values in `test.statistic` for sample size `nObs`.
    * If `nObs` and `test.statistic` are both numeric vectors, it returns a numeric vector of the same length as `nObs` and `test.statistic`. These are Bayes factor values corresponding to the different values in `nObs` and `test.statistic`.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **NAPBF_twoZ**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In case of two independent populations $N(\mu_1,\sigma_0^2)$ and $N(\mu_2,\sigma_0^2)$ with known common variance $\sigma_0^2$, consider the two-sample $z$-test for testing the point null hypothesis of difference in their means $H_0 : \mu_2 - \mu_1 = 0$ against $H_1 : \mu_2 - \mu_1 \neq 0$. Based on an observed data, this function calculates the Bayes factor in favor of $H_1$ when a *normal moment prior* is assumed on the difference $\mu_2 - \mu_1$ under the alternative [[Johnson & Rossell (2010)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2009.00730.x)].

* **Input:**
  
  * **obs1:** Numeric vector. Observed vector of data from Group-1.
  
  * **obs2:** Numeric vector. Observed vector of data from Group-2.
  
  * **n1Obs:** Numeric or numeric vector. Sample size(s) from Group-1. Same as `length(obs1)` when numeric.
  
  * **n2Obs:** Numeric or numeric vector. Sample size(s) from Group-2. Same as `length(obs2)` when numeric.
  
  * **mean.obs1:** Numeric or numeric vector. Sample mean(s) from Group-1. Same as `mean(obs1)` when numeric.
  
  * **mean.obs2:** Numeric or numeric vector. Sample mean(s) from Group-2. Same as `mean(obs2)` when numeric.
  
  * **test.statistic:** Numeric or numeric vector. Test-statistic value(s).
    
  * **tau.NAP:** Positive numeric. Parameter in the moment prior. **Default:** $0.3/\sqrt{2}$. This places the prior modes of $\mu_2 - \mu_1$ at $\pm0.3$.
  
  * **sigma0:** Positive numeric. Known common standard deviation of the populations. **Default:** 1.

* **Details:** 
  * A user can either specify `obs1` and `obs2`, or `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2`, or `n1Obs`, `n2Obs`, and `test.statistic` in the order of usage priorities.
  * If `obs1` and `obs2` are provided, it returns the corresponding Bayes factor value.
  * If `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `mean.obs1` and `mean.obs2` are numeric vectors, it returns a numeric vector of the same length as `mean.obs1` and `mean.obs2`. These are Bayes factor values corresponding to different values in `mean.obs1` and `mean.obs2` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2` are all numeric vectors, all must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
  *  If `n1Obs`, `n2Obs`, and `test.statistic` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `test.statistic` is a numeric vector, it returns a numeric vector of the same length as `test.statistic`. These are Bayes factor values corresponding to different values in `test.statistic` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, and `test.statistic` are all numeric vectors, they must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **NAPBF_twoT**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In case of two independent populations $N(\mu_1,\sigma^2)$ and $N(\mu_2,\sigma^2)$ with unknown common variance $\sigma^2$, consider the two-sample $t$-test for testing the point null hypothesis of difference in their means $H_0 : \mu_2 - \mu_1 = 0$ against $H_1 : \mu_2 - \mu_1 \neq 0$. Based on an observed data, this function calculates the Bayes factor in favor of $H_1$ when a *normal moment prior* is assumed on the standardized difference in effect sizes $(\mu_2 - \mu_1)/\sigma$ under the alternative [[Johnson & Rossell (2010)](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2009.00730.x)]. Under both hypotheses, the Jeffrey's prior $\pi(\sigma^2) \propto 1/\sigma^2$ is assumed on $\sigma^2$.

* **Input:**
  
  * **obs1:** Numeric vector. Observed vector of data from Group-1.
  
  * **obs2:** Numeric vector. Observed vector of data from Group-2.
  
  * **n1Obs:** Numeric or numeric vector. Sample size(s) from Group-1. Same as `length(obs1)` when numeric.
  
  * **n2Obs:** Numeric or numeric vector. Sample size(s) from Group-2. Same as `length(obs2)` when numeric.
  
  * **mean.obs1:** Numeric or numeric vector. Sample mean(s) from Group-1. Same as `mean(obs1)` when numeric.
  
  * **mean.obs2:** Numeric or numeric vector. Sample mean(s) from Group-2. Same as `mean(obs2)` when numeric.
  
  * **sd.obs1:** Numeric or numeric vector. Sample standard deviations(s) from Group-1. Same as `sd(obs1)` when numeric.
  
  * **sd.obs2:** Numeric or numeric vector. Sample standard deviations(s) from Group-2. Same as `sd(obs2)` when numeric.
  
  * **test.statistic:** Numeric or numeric vector. Test-statistic value(s).
    
  * **tau.NAP:** Positive numeric. Parameter in the moment prior. **Default:** $0.3/\sqrt{2}$. This places the prior modes of the difference $\mu_2 - \mu_1$ at $\pm0.3 \sigma$.
  
  * **sigma0:** Positive numeric. Known common standard deviation of the populations. **Default:** 1.

* **Details:** 
  * If `obs1` and `obs2` are provided, it returns the corresponding Bayes factor value.
  * If `n1Obs`, `n2Obs`, `mean.obs1`, `mean.obs2`, `sd.obs1` and `sd.obs2` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `mean.obs1`, `mean.obs2`, `sd.obs1` and `sd.obs2` are numeric vectors, it returns a numeric vector of the same length as `mean.obs1` and `mean.obs2`. These are Bayes factor values corresponding to different values in `mean.obs1`, `mean.obs2`, `sd.obs1` and `sd.obs2` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2` are all numeric vectors, all must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
  *  If `n1Obs`, `n2Obs`, and `test.statistic` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `test.statistic` is a numeric vector, it returns a numeric vector of the same length as `test.statistic`. These are Bayes factor values corresponding to different values in `test.statistic` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, and `test.statistic` are all numeric vectors, they must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **HajnalBF_oneZ**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In a $N(\mu,\sigma_0^2)$ population with known variance $\sigma_0^2$, consider the two-sided one-sample $z$-test for testing the point null hypothesis $H_0 : \mu = 0$ against $H_1 : \mu \neq 0$. Based on an observed data, this function calculates the Hajnal's ratio in favor of $H_1$ [[Hajnal (1961)](https://academic.oup.com/biomet/article/48/1-2/65/227215?login=true)]. The prior assumed on $\mu$ under the alternative places equal probability at two point masses $\pm\mu_1$ ($\mu_1>0$ prefixed).

* **Input:**
  
  * **obs:** Numeric vector. Observed vector of data.
  
  * **nObs:** Numeric or numeric vector. Sample size(s). Same as `length(obs)` when numeric.
  
  * **mean.obs:** Numeric or numeric vector. Sample mean(s). Same as `mean(obs)` when numeric.
  
  * **es1:** Numeric. Same as $\mu_1$ above. The prior on the population mean places equal probabilities at point masses $\pm$`abs(es1)`. **Default:** 0.3.
  
  * **sigma0:** Positive numeric. Known standard deviation in the population. **Default:** 1.

* **Details:** 
  * Users can either specify `obs`, or `nObs` and `mean.obs` in the order of usage priorities.
  * If `obs` is provided, it returns the corresponding Bayes factor value.
  * If `nObs` and `mean.obs` are provided, the function is vectorized over both of the arguments.
    * If `nObs` is a numeric value and `mean.obs` is a numeric vector, it returns a numeric vector of the same length as `mean.obs`. These are Bayes factor values corresponding to different values in `mean.obs` for sample size `nObs`.
    * If `nObs` and `mean.obs` are both numeric vectors, both must be of the same length. In this case, it returns a numeric vector of the same length as `nObs` and `mean.obs`. These are Bayes factor values corresponding to each value in `nObs` and `mean.obs`.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **HajnalBF_oneT**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In a $N(\mu,\sigma^2)$ population with unknown variance $\sigma^2$, consider the two-sided one-sample $t$-test for testing the point null hypothesis $H_0 : \mu = 0$ against $H_1 : \mu \neq 0$. Based on an observed data, this function calculates the Hajnal's ratio in favor of $H_1$ [[Hajnal (1961)](https://academic.oup.com/biomet/article/48/1-2/65/227215?login=true)]. The prior assumed on the standardized effect size $\mu/\sigma$ under the alternative places equal probabilities at two point masses $\pm\delta_1$ ($\delta_1>0$ prefixed). Under both hypotheses, the Jeffrey's prior $\pi(\sigma^2) \propto 1/\sigma^2$ is assumed on $\sigma^2$.

* **Input:**
  
  * **obs:** Numeric vector. Observed vector of data.
  
  * **nObs:** Numeric or numeric vector. Sample size(s). Same as `length(obs)` when numeric.
  
  * **mean.obs:** Numeric or numeric vector. Sample mean(s). Same as `mean(obs)` when numeric.
  
  * **sd.obs:** Positive numeric or numeric vector. Sample standard deviation(s). Same as `sd(obs)` when numeric.
  
  * **test.statistic:** Numeric or numeric vector. Test-statistic value(s).
  
  * **es1:** Numeric. Same as $\delta_1$ above. The prior on the population mean places equal probabilities at point masses $\pm$`es1`$\sigma$. **Default:** 0.3.

* **Details:** 
  * A user can either specify `obs`; or `nObs`, `mean.obs` and `sd.obs`; or `nObs` and `test.statistic` in the order of usage priorities.
  * If `obs` is provided, it returns the corresponding Bayes factor value.
  * If `nObs`, `mean.obs` and `sd.obs` are provided, the function is vectorized over these arguments.
    * If `nObs` is a numeric value, and `mean.obs` and `sd.obs` are both numeric vectors, it returns a numeric vector of the same length as `mean.obs` and `sd.obs`. These are Bayes factor values corresponding to different values in `mean.obs` and `sd.obs` for sample size `nObs`.
    * If `nObs`, `mean.obs` and `sd.obs` are all numeric vectors, the function returns a numeric vector of the same length as `nObs`, `mean.obs` and `sd.obs`. These are Bayes factor values corresponding to the different values in `nObs`, `mean.obs and `sd.obs`.
  * If `nObs` and `test.statistic` are provided, the function is vectorized over both the arguments.
    * If `nObs` is a numeric value and `test.statistic` is a numeric vector, it returns a numeric vector of the same length as `test.statistic`. These are Bayes factor values corresponding to different values in `test.statistic` for sample size `nObs`.
    * If `nObs` and `test.statistic` are both numeric vectors, it returns a numeric vector of the same length as `nObs` and `test.statistic`. These are Bayes factor values corresponding to the different values in `nObs` and `test.statistic`.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **HajnalBF_twoZ**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In case of two independent populations $N(\mu_1,\sigma_0^2)$ and $N(\mu_2,\sigma_0^2)$ with known common variance $\sigma_0^2$, consider the two-sample $z$-test for testing the point null hypothesis of difference in their means $H_0 : \mu_2 - \mu_1 = 0$ against $H_1 : \mu_2 - \mu_1 \neq 0$. Based on an observed data, this function calculates the Hajnal's ratio in favor of $H_1$ [[Hajnal (1961)](https://academic.oup.com/biomet/article/48/1-2/65/227215?login=true)]. In this case, the prior assumed on the difference $\mu_2-\mu_1$ under the alternative places equal probability at two point masses $\pm\delta_1$ ($\delta_1>0$ prefixed).

* **Input:**
  
  * **obs1:** Numeric vector. Observed vector of data from Group-1.
  
  * **obs2:** Numeric vector. Observed vector of data from Group-2.
  
  * **n1Obs:** Numeric or numeric vector. Sample size(s) from Group-1. Same as `length(obs1)` when numeric.
  
  * **n2Obs:** Numeric or numeric vector. Sample size(s) from Group-2. Same as `length(obs2)` when numeric.
  
  * **mean.obs1:** Numeric or numeric vector. Sample mean(s) from Group-1. Same as `mean(obs1)` when numeric.
  
  * **mean.obs2:** Numeric or numeric vector. Sample mean(s) from Group-2. Same as `mean(obs2)` when numeric.
  
  * **test.statistic:** Numeric or numeric vector. Test-statistic value(s).
  
  * **es1:** Numeric. $\delta_1$ as above. The prior on the difference $\mu_2-\mu_1$ places equal probabilities at point masses $\pm$`abs(es1)`. **Default:** 0.3.
  
  * **sigma0:** Positive numeric. Known common standard deviation of the populations. **Default:** 1.

* **Details:** 
  * A user can either specify `obs1` and `obs2`, or `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2`, or `n1Obs`, `n2Obs`, and `test.statistic` in the order of usage priorities.
  * If `obs1` and `obs2` are provided, it returns the corresponding Bayes factor value.
  * If `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `mean.obs1` and `mean.obs2` are numeric vectors, it returns a numeric vector of the same length as `mean.obs1` and `mean.obs2`. These are Bayes factor values corresponding to different values in `mean.obs1` and `mean.obs2` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2` are all numeric vectors, all must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
  *  If `n1Obs`, `n2Obs`, and `test.statistic` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `test.statistic` is a numeric vector, it returns a numeric vector of the same length as `test.statistic`. These are Bayes factor values corresponding to different values in `test.statistic` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, and `test.statistic` are all numeric vectors, they must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).

### **HajnalBF_twoT**

<font size="3"> [[List of all functions]](#list-of-functions) </font>

In case of two independent populations $N(\mu_1,\sigma^2)$ and $N(\mu_2,\sigma^2)$ with unknown common variance $\sigma^2$, consider the two-sample $t$-test for testing the point null hypothesis of difference in their means $H_0 : \mu_2 - \mu_1 = 0$ against $H_1 : \mu_2 - \mu_1 \neq 0$. Based on an observed data, this function calculates the Hajnal's ratio in favor of $H_1$ [[Hajnal (1961)](https://academic.oup.com/biomet/article/48/1-2/65/227215?login=true)]. In this case, the prior assumed on the standardized difference $(\mu_2-\mu_1)/\sigma$ under the alternative places equal probability at two point masses $\pm\delta_1$ ($\delta_1>0$ prefixed). Under both hypotheses, the Jeffrey's prior $\pi(\sigma^2) \propto 1/\sigma^2$ is assumed on $\sigma^2$.

* **Input:**
  
  * **obs1:** Numeric vector. Observed vector of data from Group-1.
  
  * **obs2:** Numeric vector. Observed vector of data from Group-2.
  
  * **n1Obs:** Numeric or numeric vector. Sample size(s) from Group-1. Same as `length(obs1)` when numeric.
  
  * **n2Obs:** Numeric or numeric vector. Sample size(s) from Group-2. Same as `length(obs2)` when numeric.
  
  * **mean.obs1:** Numeric or numeric vector. Sample mean(s) from Group-1. Same as `mean(obs1)` when numeric.
  
  * **mean.obs2:** Numeric or numeric vector. Sample mean(s) from Group-2. Same as `mean(obs2)` when numeric.
  
  * **sd.obs1:** Numeric or numeric vector. Sample standard deviations(s) from Group-1. Same as `sd(obs1)` when numeric.
  
  * **sd.obs2:** Numeric or numeric vector. Sample standard deviations(s) from Group-2. Same as `sd(obs2)` when numeric.
  
  * **test.statistic:** Numeric or numeric vector. Test-statistic value(s).
    
  * **es1:** Numeric. $\delta_1$ as above. The prior on the difference $\mu_2-\mu_1$ places equal probabilities at point masses $\pm$`abs(es1)`$\sigma$. **Default:** 0.3.
  
  * **sigma0:** Positive numeric. Known common standard deviation of the populations. **Default:** 1.

* **Details:** 
  * If `obs1` and `obs2` are provided, it returns the corresponding Bayes factor value.
  * If `n1Obs`, `n2Obs`, `mean.obs1`, `mean.obs2`, `sd.obs1` and `sd.obs2` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `mean.obs1`, `mean.obs2`, `sd.obs1` and `sd.obs2` are numeric vectors, it returns a numeric vector of the same length as `mean.obs1` and `mean.obs2`. These are Bayes factor values corresponding to different values in `mean.obs1`, `mean.obs2`, `sd.obs1` and `sd.obs2` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, `mean.obs1` and `mean.obs2` are all numeric vectors, all must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
  *  If `n1Obs`, `n2Obs`, and `test.statistic` are provided, the function is vectorized over each of these arguments.
    * If `n1Obs` and `n2Obs` are numeric values, and `test.statistic` is a numeric vector, it returns a numeric vector of the same length as `test.statistic`. These are Bayes factor values corresponding to different values in `test.statistic` for sample sizes `n1Obs` from Group-1 and `n2Obs` from Group-2.
    * If `n1Obs`, `n2Obs`, and `test.statistic` are all numeric vectors, they must have the same length. In this case, it returns a numeric vector of the same length. These are Bayes factor values corresponding to each value in these vectors.
    
* **Output:** Positive numeric or numeric vector. The Bayes factor value(s).


<font size="3"> [[Top]](#description) </font>
<font size="3"> [[List of all functions]](#list-of-functions) </font>

