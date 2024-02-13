function [lambdaCumul,lambda] = lambdaCumulEstTSumMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,p,d,group)
lambda = lambdaEstMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,p,d,group);
lambdaCumul = 0;
lambdaCumul = [lambdaCumul; cumsum((lambda(1:(end-1))+lambda(2:end))/2.*diff(inputt))];
end