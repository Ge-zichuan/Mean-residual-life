function lambdaVec = lambdaTildeEst(inputt,inputg,zall,aall,deltaall,samplesize,goodI,bands)
% This is the estimation of lambda, rows are i, columns are t.
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
nt = length(inputt);
deltaallmat = repmat(deltaall,1,nt);
yMinusT = repmat(zall,1,nt)-repmat(inputt',samplesize,1);% {i,j} = Z_i - t_j
aallmat = aall*ones(1,nt);
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
inputgall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(inputg,[1 1 samplesize]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate nonparametric \lambda 
khbetax = kernel(inputgall/banddh/bandadjdhd,3,'norm')/banddh/bandadjdhd;
denoVec = khbetax(:,goodI)*(aallmat(goodI,:).*(yMinusT(goodI,:)>=0));
khbetax = kernel(inputgall/banddh/bandadjdhn,3,'norm')/banddh/bandadjdhn;

numVec = khbetax(:,goodI)*(aallmat(goodI,:).*(yMinusT(goodI,:)==0).*deltaallmat(goodI,:));
lambdaVec = numVec./denoVec;

end