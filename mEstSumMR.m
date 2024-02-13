function [m,cumLambda] = mEstSumMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,p,d,upperLimit,group)
% estSample = length(inputt(:,1));
% finer the grid of t
finerCoef = 100;
inputtFiner = [0:(upperLimit/finerCoef):upperLimit,inputt']';
inputtFiner = unique(inputtFiner);
inputgFiner = ones(length(inputtFiner),1)*inputg;
[cumLambdaFiner,lambdaEst] = lambdaCumulEstSumMR(parasLower,inputtFiner,inputgFiner,xall,zall,deltaall,samplesize,bands,p,d,group);
inputtRever = inputtFiner(end)-inputtFiner(end:-1:1);
survCumSumRever = cumsum((exp(-cumLambdaFiner(end:-1:2))+exp(-cumLambdaFiner((end-1):-1:1)))/2.*diff(inputtRever));
intgral2 = [survCumSumRever(end:-1:1);0];
finerm = exp(cumLambdaFiner).*intgral2;
coarseInd = find(ismember(inputtFiner,inputt));
m = finerm(coarseInd);
cumLambda = cumLambdaFiner(coarseInd);
end

% %
% figure;
% plot(inputtFiner,inputtFiner*exp(inputg),inputtFiner,lambdaEst);
% figure;
% plot(inputtFiner,inputtFiner.^2*exp(inputg)/2,inputtFiner,cumLambdaFiner);
% figure;
% plot(inputtFiner,exp(inputtFiner.^2/2*exp(inputg)).*normcdf(-inputtFiner,0,sqrt(1/exp(inputg)))*sqrt(2*pi)*sqrt(1/exp(inputg)),inputtFiner,finerm);hold off;
