
## Non-local Alternative Priors (NAP) for Bayesian Hypothesis Tests in Psychology

## Authors: Sandipan Pramanik and Valen E. Johnson (2021+)

## Description

Bayesian hypothesis testing procedures have gained increased acceptance in recent years.  A key advantage that Bayesian tests have over classical testing procedures is their potential to quantify information in support of true null hypotheses.  Ironically, default implementations of Bayesian tests prevent accumulation of strong evidence in favor of true null hypotheses because associated default alternative hypotheses assign high probability to data that are most consistent with a null effect. We propose the use of "non-local" alternative hypotheses to resolve this paradox. The resulting class of Bayesian hypothesis tests permit more rapid accumulation of evidence in favor of both true null hypotheses and alternative hypotheses that are compatible with standardized effect sizes of most interest in psychology.

The codes available here conducts Bayesian Hypothesis tests of a point null hypothesis against a two-sided alternative using Non-local Alternative Prior (NAP) for one- and two-sample z- and t-tests. Under the alternative, the NAP is assumed on the standardized effects size in one-sample tests and on their differences in two-sample tests. The package considers two types of NAP densities: (1) the normal moment prior, and (2) the composite alternative. In fixed design tests, the functions calculate the Bayes factors and the expected weight of evidence for varied effect size and sample size. The package also provides a sequential testing framework using the Sequential Bayes Factor (SBF) design. The functions calculate the operating characteristics (OC) and the average sample number (ASN), and also conducts sequential tests for a sequentially observed data.
