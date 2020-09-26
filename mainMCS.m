%MCS procedure following Sheppards Toolbox. 

%losses computed in R and imported here.


losses = readmatrix("losses_transformed.csv");
losses = losses(:,2:end);

[includedR,pvalsR,excludedR,includedSQ,pvalsSQ,excludedSQ] = mcs(losses,0.05,1000,10,"Block")

includedR
includedSQ

losses(1:4,includedR)