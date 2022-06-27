fileID = fopen("integralData.txt", "r");
integralVals = fscanf(fileID, '%d', 100);
fclose(fileID);

fileID = fopen("integralBinData.txt");
xVals = fscanf(fileID, '%f', 100);
fclose(fileID);

data = plot(xVals, integralVals);

customLogPDF = @(data, N0, lambda) ... 
    N0(1-e^(-lambda*data));

phat = mle(data, 'nloglf', customLogPDF, 'start', 0.05);