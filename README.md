# Master thesis

## Estimation and forecasting of multivariate high frequency volatility models: application to a leveraged risk-parity portfolio 

### Quick summary

The purpose of the thesis was to study a selection of bivariate high frequency volatility models and compare them using Hansens model confidence set (MCS) ([link](https://onlinelibrary.wiley.com/doi/pdf/10.3982/ECTA5771?casa_token=Zy0wE0wmFTQAAAAA:_OPaHsOkFadekMMWs7F6nmTxUmX7_ej5r6oBKDVzV8nWK849XUdTsb6gLeCaeTD9l-l3AC9x3yi6d1o)). A subset of the best models found via the MCS procedure, was then used in a risk-targeted risk-parity portfolio. The idea was to construct a simplified portfolio following the same strategy of [Man AHL TargetRisk](https://www.man.com/ahl-targetrisk). The general results for the thesis are described below.  

From 10 years of data on the ETFs, *SPY* and *TLT*, I had more than 250 GB of raw data which was cleaned following the procedure described in [Barndorff‐Nielsen et al. (2009)](https://academic.oup.com/ectj/article/12/3/C1/5061260?casa_token=JHbSnyQ9xa4AAAAA:oZ6WLWvC91FcyD9WKwB_JYrB4HEPpHQFj2sTSFDvBqmmoxowoHVD-ASOuo5nu_AnCXairSkzLb8K). This gave us enough usable data to construct sample schemes down to second frequencies. Therefore, we ended up with more than 350 volatility models across frequencies, realized measures and volatility models, to be tested using forecast comparison methods.  

*It needs to be noted that the portfolios was heavily affect by the stock-bond correlation, since I chose a time horizon with negative correlation, thus imposing better diversification in terms of how a bond position can hedge against an equity marrket sell-off.*  

---
### Results


The results from my thesis are described below.

**From the analysis of the realized measures:**

1. I found evidence that the literature-preferred realized covariance estimator on a 5-minute frequency, <img src="images/rcov5min.svg" />, was significantly out-performed by 4 out of 75 realized measures. 

2. In line with [Li et al. (2015)](https://www.sciencedirect.com/science/article/pii/S0304407615000329?casa_token=75e1b--a2ZoAAAAA:78E3kAVjp7QY3u8CKwHe-vEG7tprV26zJHp0W6VgIMiTNlWuKXElk98cCXTOPtKLSms4VJZ6bw), I concluded that <img src="images/rcov1min.svg" /> appeared more often in the superior set, and that noise-robust measures on second-frequencies tend to perform well for liquid equity-indices or proxies thereof.

3. An alteration of the MCS procedure with an economic loss function, favoured jump-robust estimators under a proxy of integrated covariation which contributes to the fact that measuring quadratic covariation induces more variability in each estimator, that might procure higher variances in the MVP framework.

4. Conclusively, we found evidence that many of the realized measures on minute-frequencies seems to provide much of the benefits of high frequency data without exposing the estimators to microstructure noise, and the empirical accuracy only slightly increases (on average) when considering noise-robust estimators on second frequencies. 

*In general, the analysis was done to investigate how the superior performance of the realized measures translated over to the volatility models.*

**From the analysis of the realized volatility models:**

1. I found evidence that continuous assymetrical models (accounting for leverage effect) could not 
capture the leverage effect for both models, due to a "inverse"-type leverage effect relationship 
for the assets. 

2. In general, the DRD-HAR(Q) type models performed the best and procurred the lowest statistical losses, while being very parsimonious. 

3. All of the high frequency volatility models (even the poorest ones) beat their open-to-close volatility model counterparts as noted by the forecast comparison analysis. 

4. I further observed that the superior performance of the realized measures was carried over to the volatility models. With alternating frequencies, it tells us that the choice of realized measure matters more than the frequency.

**From the portfolio analysis:**

1. Transitioning from volatility models on daily returns to any DRD-HAR type model under consideration, decreased the vol of vol substantially while retaining (approximately) the same annualized expected returns. This lead to a reduction in the expected shortfall under a time-horizon with no apparent crises, that would net the risk-averse investor a 2.8 million notional difference, if they invest 1 billion in the DRD-CHAR portfolio.

2. Contrary to [Bollerslev et al. (2018)](https://www.sciencedirect.com/science/article/pii/S0304407618301180?casa_token=gJ0dxNPfv2IAAAAA:LKASCIaIPZdJAR1CraKd3yu-8ljjr2fufoFVCrmlnOW9mJ9pNzZcTR2RfpI4QgvTwlre9KUAfA), we found evidence that the turnover rates of the ARQ type models increased substantially adverse of the HAR-type dynamics, which contributed to a reduction in performance (as seen by Sharpe ratios) under transaction costs. 

3. The most surprising result, was the universal performance of the CHAR model. Upon further investigation, we found that the DRD-CHAR model was in the top three across all cost-schemes in terms of Sharpe ratio, resulting 
from the reduced turnover rates.

4. In general, it lead me to believe, that good out-of-sample portfolio performance is not only dictated by models within the superior set of the MCS procedure, but in some cases, rejected parsimonious models might perform consistently well across realized measures with adequately chosen frequencies. 


## Code
The code for my master thesis was fractioned into small sniplets corresponding to different subsections of my thesis. You **cannot** download the repository and expect all of the code to be runable due to lacking the original cleaned stock data. The code has been provided for illustrational purposes. 


**```functions.R```** contains all of my self-written ad-hoc functions used throughout the entire thesis. 

**Summarystats.R**, **Summary stats for cleaning Podolskij method.R** and **jumpproportion.R** contains the code for the summary statistics and analysis of the assets, provided in the empirical analysis.

**DRD-HAR models.R** and **Volatility models.R** contains the code for the in-sample estimation and analysis of the volatility models. 

**mcs_realized_measures.R** and **out_of_sample_MCS_and_prelim_volmodels.R** contains the forecast comparison analysis for the realized measures and the volatility models. Morevover, the results was compared with the MCS function in [Kevin sheppards toolbox](https://www.kevinsheppard.com/code/matlab/mfe-toolbox/) and the code sniplet can be seen in **mainMCS.m**. 

**Out_of_sample_forecast_analysis.R** contains the code for the out-of-sample forecast analysis. 

**Out_of_sample_portfolio_analysis.R** contains the code for the out-of-sample portfolio analysis.

**riskparity.R** contains the code that produces the graphs provided in the theoretical section of the risk-parity portfolio. 

**Correlation check.R**, **IntradayACF.R**, **semicovplot.R**, **pre-averaging stability.R** and **realized semicovariance graph.R** contains the code that produces preliminary graphs before the empirical analysis analysis. 
