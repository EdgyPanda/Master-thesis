# Master thesis

## Estimation and forecasting of multivariate high frequency volatility models: application to a leveraged risk-parity portfolio 

The purpose of the thesis was to study a selection of bivariate high frequency volatility models and compare them using Hansens model confidence set (MCS) ([link](https://onlinelibrary.wiley.com/doi/pdf/10.3982/ECTA5771?casa_token=Zy0wE0wmFTQAAAAA:_OPaHsOkFadekMMWs7F6nmTxUmX7_ej5r6oBKDVzV8nWK849XUdTsb6gLeCaeTD9l-l3AC9x3yi6d1o)). From 10 years of data on the ETFs, SPY and TLT I had more than 250 GB of raw data which was cleaned following the procedure in [Barndorff‐Nielsen et al. (2009)](https://academic.oup.com/ectj/article/12/3/C1/5061260?casa_token=JHbSnyQ9xa4AAAAA:oZ6WLWvC91FcyD9WKwB_JYrB4HEPpHQFj2sTSFDvBqmmoxowoHVD-ASOuo5nu_AnCXairSkzLb8K). This gave us enough usable data to construct sample schemes down to second frequencies. Therefore, we ended up with more than 350 volatility models across frequencies, realized measures and volatility models, to be tested using forecast comparison methods.  

A selection of the best volatility models found within the superior set of the MCS procedure was then used in a risk-targeted risk-parity portfolio setup with surprising results. The results are described below:

**From the analysis of the realized measures:**

1. I found evidence that the realized covariance estimator on a 5-minute frequency, <img src="images/rcov5min.svg" />, was significantly out-performed by 4 out of 75 realized measures. 

2. In line with [Li et al. (2015)](https://www.sciencedirect.com/science/article/pii/S0304407615000329?casa_token=75e1b--a2ZoAAAAA:78E3kAVjp7QY3u8CKwHe-vEG7tprV26zJHp0W6VgIMiTNlWuKXElk98cCXTOPtKLSms4VJZ6bw), I concluded that <img src="images/rcov1min.svg" /> appeared more often in the superior set, and that noise-robust measures on second-frequencies tend to perform well for liquid equity-indices or proxies thereof.


It needs to be noted that the portfolios was heavily affect by the stock-bond correlation, since I chose a time horizon with negative correlation, thus imposing better diversification in terms of how a bond position can hedge against an equity marrket sell-off.  

The code for my master thesis. The code will be fractioned into small sniplets corresponding to different subsections of my thesis. This will provide a better overview and make things easier. The data cleaning code was only used once for time efficiency, and thus I stored the cleaned data as .Rdata of the assets that I'm working with.  