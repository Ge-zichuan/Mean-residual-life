function lambdaVec = lambdaEstTMRSpec(parasLower,t,inputg,xall,zall,wall,deltaall,samplesize,bands,d,group)
% one t at several different \bb\trans\x
estSample = length(inputg(:,1));
inputt = t*ones(estSample,1);
% read bandwidth and bandwidth adjustment
if group == 1
    bandb = bands(1,2);
else
    bandb = bands(1,1);
end
bandadjb = bands(1,3);
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
% calculate all Z_i-estZ
% zdiffall = inputt*ones(1,samplesize)-ones(estSample,1)*zall';% {i,j} = Z_i-Z_j
% kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
zdiffest = zall*ones(1,estSample)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffest/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
btransx = [xall*[eye(d);parasLower],wall]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate nonparametric \lambda and \lambda\prime
khbetax = kernel(btxall/banddh/bandadjdhd,3,'norm')/(banddh*bandadjdhd)^(d+1);
denoVec = (zdiffall<=0)*khbetax;
khbetax = kernel(btxall/banddh/bandadjdhn,3,'norm')/(banddh*bandadjdhn)^(d+1);
lambdaVec = sum(kbz.*repmat(deltaall,1,estSample).*khbetax./denoVec,1)';

end