function lambdaVec = lambdaEstTMR(parasLower,inputt,inputg,xall,zall,wall,deltaall,samplesize,bands,d,group)
estSample = length(inputg(:,1));
% read bandwidth and bandwidth adjustment
if group == 1% transplant
    bandb = bands(1,2);
    bandadjb = bands(1,4);
    bandadjdhn = bands(2,2);
    bandadjdhd = bands(2,3);
    bandadjdhw = bands(2,4);
    banddh21 = bands(2,1);
    banddh22 = bands(2,2);
    bandadjdhn23 = bands(2,3);
    bandadjdhd24 = bands(2,4);
    bandadjdhn25 = bands(2,5);
    bandadjdhd26 = bands(2,6);

else% nontransplant
    bandb = bands(1,1);
    bandadjb = bands(1,3);
    banddh = bands(2,1);
    bandadjdhn = bands(2,2);
    bandadjdhd = bands(2,3);
    bandadjdhw = bands(2,4);
end
% calculate all Z_i-estZ
% zdiffall = inputt*ones(1,samplesize)-ones(estSample,1)*zall';% {i,j} = Z_i-Z_j
% kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
zdiffest = zall*ones(1,estSample)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffest/(bandb*bandadjb),1,'Epanechnikov')/(bandb*bandadjb);% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
btransx = [xall*[eye(d);parasLower],wall]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate nonparametric \lambda and \lambda\prime
bandTwo = permute(repmat(banddh21*[bandadjdhd bandadjdhw]',[1 estSample samplesize]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'norm')/(((banddh21*bandadjdhd)^d)*(banddh21*bandadjdhw));
denoVec = (zdiffall<=0)*khbetax;
bandTwo = permute(repmat(banddh21*[bandadjdhn bandadjdhw]',[1 estSample samplesize]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'Epanechnikov')/(((banddh21*bandadjdhn)^d)*(banddh21*bandadjdhw));
lambdaVec = sum(kbz.*repmat(deltaall,1,estSample).*khbetax./denoVec,1)';

end