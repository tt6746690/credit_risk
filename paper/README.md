

#### Newer methods

+ constitutive relations
    + https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Constitutive_relations#Isotropic_materials

+ https://ac.els-cdn.com/S0377221710000111/1-s2.0-S0377221710000111-main.pdf?_tid=e891d830-3184-4093-a4e9-402ff58716ae&acdnat=1545954668_cf8cdc9e2f8368d525f000bc93592bc8
    + generalizes 2005 paper using t-copula model instead
    + an overview of what glasserman has done previously
        + 2005 Importance sampling for portfolio credit risk
            + initial paper on using IS to estimate credit risk
            + logarithmically efficient algorithm
        + 2007 Large deviations of multifactor portfolio credit risk.
            + generalizes to multifactor
        + 2008 fast simulation of multifactor portfolio credit risk 
            + fixes problem for previous papers

+ http://personal.cb.cityu.edu.hk/jeffhong/papers/HongHuLiuTOMACS2014.pdf
    + a good review to what glasserman has done over the years
        + 2005 paper: optimality for single-factor homogeneous models
        + 2008 paper: mixture of mean shifts for Z for multi-factor models

+ https://ac.els-cdn.com/S0167668715000955/1-s2.0-S0167668715000955-main.pdf?_tid=027b0be3-51a4-45da-94d9-17bcaf373ed6&acdnat=1545953738_d34fdb938168919b6e2d8482eeab94bb
    + normal copula model 2005_GL is best performing
    + experiment
        + compares GL and its algorithm (try to replicate result)
        + might need some time to do this ... since its hard coded and need to understand its use of symbols

#### Monte Carlo for Credit Risk


+ [2005_importance_sampling_for_portfolio_credit_risk](2005_importance_sampling_for_portfolio_credit_risk.pdf)
    + Glasserman&Li
    + two level IS for corpula credit model     
        + outer level: shift mean of normal distribution of systematic risk factors
            + works well when systematic risk highly dependent cross oligors
        + inner level: exponential twisting to default probabilities conditioned on the systematic factors
            + works well, and maybe sufficient when obligors are weakly dependent
    + think about using its metric, specifically, a plot of plain MC and IS where x-axis is the tail thresholds and the y-axis is the corresponding simulated probability. If we do this for many replications, we can see the effect of variance reduction with IS methods by plotting curves for the simulated probabilities and 95% CI. We can utilize this metric in evaluating the different models implemented in Adam Sturge's paper

+ [2008_fast_simulation_of_multifactor_portfolio_credit_risk](2008_fast_simulation_of_multifactor_portfolio_credit_risk.pdf)
    + revised GL method to correct the algorithm in 2005 paper

+ [2015_new_approaches_to_importance_sampling_for_portfolio_credit_risk_valuation](2015_new_approaches_to_importance_sampling_for_portfolio_credit_risk_valuation.pdf)
    + a better introduction to motivation & introduction

+ [2018_an_importance_sampling_scheme_for_multi_credit_state_portfolios](2018_an_importance_sampling_scheme_for_multi_credit_state_portfolios.pdf)
    + Adam Sturge Msc thesis



#### Importance Sampling 


+ [importance_sampling_a_review_duke](importance_sampling_a_review_duke.pdf)

+ [psu_astrostatistics_importance_sampling](psu_astrostatistics_importance_sampling.pdf)
    + a set of slides, intro to MC and importance sampling

+ [university_of_arizona_importance_sampling_tom_kennedy](university_of_arizona_importance_sampling_tom_kennedy.pdf)
    + chapter 6 of Tom Kennedy's book
    + good introduction on
        + importance sampling
        + strategy on how to pick proposal distribution `q`
        + self-normalizing IS
            + if `p` and `q` known up to an unknown constant factor
        + variance minimization
            + minimize `q` by restricting it to be in some family, `p(x, theta)`
                + exponential tiling: `q` is in exponential     family
            + consider IS as weighted MC
            + reference distribution

#### Optimization 

+ [2016_on_differentiating_parameterized_argmin_and_argmax_problems_with_applications_to_bi_level_optimization](2016_on_differentiating_parameterized_argmin_and_argmax_problems_with_applications_to_bi_level_optimization.pdf)
    + bi-level optimization


#### Sampling 

+ [2003_slice_sampling](2003_slice_sampling.pdf)
    + impls 
        + https://homepages.inf.ed.ac.uk/imurray2/teaching/09mlss/
        + http://www.mas.ncl.ac.uk/~ntwn/talks/slice.pdf
    + more efficient than metropolis & gibbs
    + but samples are correlated
    + 
