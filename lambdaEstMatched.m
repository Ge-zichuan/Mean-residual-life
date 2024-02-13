function lambdaVec = lambdaEstMatched(inputt,repeats,zall,aall,deltaall,samplesize,goodI)
% This is the estimation of lambda, rows are i, columns are t.
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
nt = length(inputt);
deltaallmat = repmat(deltaall,1,nt);
yMinusT = repmat(zall,1,nt)-repmat(inputt',samplesize,1);% {i,j} = Z_i - t_j
aallMat = aall*ones(1,nt);
repeatsMat = repeats*ones(1,nt);
% read bandwidth and bandwidth adjustment

% calculate nonparametric \lambda 
denoVec = sum(repeatsMat(goodI,:).*aallMat(goodI,:).*(yMinusT(goodI,:)>=0),1);

numVec = sum(repeatsMat(goodI,:).*aallMat(goodI,:).*(yMinusT(goodI,:)==0).*deltaallmat(goodI,:),1);
lambdaVec = numVec'./denoVec';
end