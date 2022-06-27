%{
rng(18,'twister') % For reproducibility
lambda = 1.75;
n = 75;
x1 = poissrnd(lambda,n,1);
x1 = x1(x1 > 0);

histogram(x1,0:1:max(x1)+1)

pf_truncpoiss = @(x1,lambda) poisspdf(x1,lambda)./(1-poisscdf(0,lambda));

start = mean(x1)

[lambdaHat,lambdaCI] = mle(x1,'pdf',pf_truncpoiss,'Start',start, ...
    'LowerBound',0)

[lambdaHat2,lambdaCI2] = mle(x1,'Distribution','Poisson', ...
    'TruncationBounds',[0 Inf])
%}

fileID = fopen("integralData.txt", "r");
integralVals = fscanf(fileID, '%d', 100);
fclose(fileID);

fileID = fopen("integralBinData.txt");
xVals = fscanf(fileID, '%f', 100);
fclose(fileID);

data = [xVals, integralVals];
start = [0.05];

integralMethodFunction = @(xVals, lambda) (1 -lambda + xVals);

[phat1, phat2] = mle(data, 'pdf', integralMethodFunction, 'start', start);

disp(phat)
