function lambdaVec = lambdaEst(parasLower,inputt,inputg,xall,zall,aall,deltaall,samplesize,goodI,bands)
estSample = length(inputt);
% treatment matrix
aallmat = aall*ones(1,samplesize);
deltaallmat = repmat(deltaall,1,samplesize);
goodIndex1 = (goodI{1});
goodIndex0 = (goodI{2});
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
% calculate all Z_i-estZ
% zdiffall = inputt*ones(1,samplesize)-ones(estSample,1)*zall';% {i,j} = Z_i-Z_j
% kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
btransx = xall*[eye(d);parasLower]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate estimating equation, each term.

% calculate nonparametric \lambda and \lambda\prime
khbetax = kernel(btxall/banddh/bandadjdhd,3,'norm')'/banddh/bandadjdhd;
denoVec1 = khbetax*((zdiffall>=0).*aallmat);
denoVec1 = denoVec1(:,goodIndex1);
denoVec0 = khbetax*((zdiffall>=0).*(1-aallmat));
denoVec0 = denoVec0(:,goodIndex0);
khbetax = kernel(btxall/banddh/bandadjdhn,3,'norm')/banddh/bandadjdhn;
lambdaDen1 = sum(aallmat(goodIndex1,:)'.*kbz(:,goodIndex1).*deltaallmat(goodIndex1,:)'.*khbetax(:,goodIndex1)./denoVec1,2);
lambdaDen0 = sum((1-aallmat(goodIndex0,:))'.*kbz(:,goodIndex0).*deltaallmat(goodIndex0,:)'.*khbetax(:,goodIndex0)./denoVec0,2);
lambdaVec = aall.*lambdaDen1+(1-aall).*lambdaDen0;

end