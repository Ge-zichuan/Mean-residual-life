function [m,cumLambda] = mEstTSumMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,p,d,upperLimit,group)
% estSample = length(inputt(:,1));
% finer the grid of t
finerCoef = 100;
inputtFiner = [0:(upperLimit/finerCoef):upperLimit,inputt']';
inputtFiner = unique(inputtFiner);
inputgFiner = ones(length(inputtFiner),1)*inputg;
[cumLambdaFiner,lambdaEst] = lambdaCumulEstTSumMR(parasLower,inputtFiner,inputgFiner,xall,zall,deltaall,samplesize,bands,p,d,group);
inputtRever = inputtFiner(end)-inputtFiner(end:-1:1);
survCumSumRever = cumsum((exp(-cumLambdaFiner(end:-1:2))+exp(-cumLambdaFiner((end-1):-1:1)))/2.*diff(inputtRever));
intgral2 = [survCumSumRever(end:-1:1);0];
finerm = exp(cumLambdaFiner).*intgral2;
coarseInd = find(ismember(inputtFiner,inputt));
m = finerm(coarseInd);
cumLambda = cumLambdaFiner(coarseInd);
end