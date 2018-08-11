


#### Monte Carlo for Credit Risk


+ [2005_importance_sampling_for_portfolio_credit_risk](2005_importance_sampling_for_portfolio_credit_risk.pdf)
    + Glasserman&Li
    + two level IS for corpula credit model     
        + outer level: shift mean of normal distribution of systematic risk factors
            + works well when systematic risk highly dependent cross oligors
        + inner level: exponential twisting to default probabilities conditioned on the systematic factors
            + works well, and maybe sufficient when obligors are weakly dependent
    + think about using its metric, specifically, a plot of plain MC and IS where x-axis is the tail thresholds and the y-axis is the corresponding simulated probability. If we do this for many replications, we can see the effect of variance reduction with IS methods by plotting curves for the simulated probabilities and 95% CI. We can utilize this metric in evaluating the different models implemented in Adam Sturge's paper

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
