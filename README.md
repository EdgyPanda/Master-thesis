# Master thesis

## Estimation and forecasting of multivariate high frequency volatility models: application to a leveraged risk-parity portfolio 

The purpose of the thesis was to study a selection of bivariate high frequency volatility models and compare them using Hansens model confidence set (MCS) ([link](https://onlinelibrary.wiley.com/doi/pdf/10.3982/ECTA5771?casa_token=Zy0wE0wmFTQAAAAA:_OPaHsOkFadekMMWs7F6nmTxUmX7_ej5r6oBKDVzV8nWK849XUdTsb6gLeCaeTD9l-l3AC9x3yi6d1o)). From 10 years of data on the ETFs, SPY and TLT I had more than 250 GB of raw data which was cleaned following the procedure in [Barndorff‐Nielsen et al. (2009)](https://academic.oup.com/ectj/article/12/3/C1/5061260?casa_token=JHbSnyQ9xa4AAAAA:oZ6WLWvC91FcyD9WKwB_JYrB4HEPpHQFj2sTSFDvBqmmoxowoHVD-ASOuo5nu_AnCXairSkzLb8K). This gave us enough usable data to construct sample schemes down to second frequencies. Therefore, we ended up with more than 350 volatility models across frequencies, realized measures and volatility models, to be tested using forecast comparison methods.  

A selection of the best volatility models found within the superior set of the MCS procedure was then used in a risk-targeted risk-parity portfolio setup with surprising results. 

The code for my master thesis. The code will be fractioned into small sniplets corresponding to different subsections of my thesis. This will provide a better overview and make things easier. The data cleaning code was only used once for time efficiency, and thus I stored the cleaned data as .Rdata of the assets that I'm working with.  