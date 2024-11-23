function [stdVal,jaccobVal] = jaccobBetaMRW(parasLower,xall,wall,zdiffallN,zdiffallT,zdiffall,kbzN,kbzT...
    ,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands)
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
paran = d*(p-d);
jaccobValAll = zeros(samplesize,length(parasLower(:)));
jaccobVal = zeros(length(parasLower(:)),length(parasLower(:)));
numeVecT = zeros(samplesizeT,samplesizeT,d);
numeVecN = zeros(samplesizeN,samplesizeN,d);
lambdaDen = zeros(samplesize,1);
expectVec = zeros(samplesize,p-d);
% read bandwidth and bandwidth adjustment
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
bandnh = bands(3,1);
bandadjnhn = bands(3,2);
bandadjnhd = bands(3,3);
bandhe = bands(4,1);
bandadjed = bands(4,3);
% calculate \bb\trans\x
btransx = xall*[eye(d);reshape(parasLower,p-d,d)];
btransxT = [xall(groupIndex==1,:)*[eye(d);reshape(parasLower,p-d,d)],wall(groupIndex==1,:)]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    btxall = permute(repmat(btransx',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 samplesize]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    btxallT = permute(repmat(btransxT',[1 1 samplesizeT]),[3 2 1])-permute(repmat(btransxT,[1 1 samplesizeT]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    btxallN = btxall;% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 samplesize])-repmat(btransx',[samplesize 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    btxallT = repmat(btransxT,[1 samplesizeT])-repmat(btransxT',[samplesizeT 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    btxallN = btxall;% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
% calculate estimating equation, each term.

% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
khbetaxT = kernel(btxallT/(banddh*bandadjdhd),kerD,'Epanechnikov')/(banddh*bandadjdhd)^d;
khbetaxN = kernel(btxallN/(banddh*bandadjdhd),kerD,'Epanechnikov')/(banddh*bandadjdhd)^d;
denoVecT = khbetaxT*(zdiffallT>=0);
denoVecN = khbetaxN*(zdiffallN>=0);
khbetaxT = kernel(btxallT/(banddh*bandadjdhn),kerD,'Epanechnikov')/(banddh*bandadjdhn)^d;
khbetaxN = kernel(btxallN/(banddh*bandadjdhn),kerD,'Epanechnikov')/(banddh*bandadjdhn)^d;
lambdaDen(:,1) = sum(kbzN.*repmat(deltaallN,1,samplesizeN)'.*khbetaxN./denoVecN,2);
lambdaDen(groupIndex==1,1) = sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetaxT./denoVecT,2);
% calculate kernel values of K\prime(.)
khprimevecT = kernel2(btxallT/(bandnh*bandadjnhn),kerD,d+1,'Epanechnikov')/(bandnh*bandadjnhn)^(d+2);
khprimevecN = kernel2(btxallN/(bandnh*bandadjnhn),kerD,d,'Epanechnikov')/(bandnh*bandadjnhn)^(d+1);
for t = 1:d
    numeVecT(:,:,t) = squeeze(khprimevecT(:,:,t))*(zdiffallT>=0);
    numeVecN(:,:,t) = squeeze(khprimevecN(:,:,t))*(zdiffallN>=0);
end
lambdaNum = zeros(samplesize,d);
khbetaxT = kernel(btxallT/(bandnh*bandadjnhd),kerD,'Epanechnikov')/(bandnh*bandadjnhd)^d;
khbetaxN = kernel(btxallN/(bandnh*bandadjnhd),kerD,'Epanechnikov')/(bandnh*bandadjnhd)^d;
for t = 1:d
    lambdaNum(:,t) = sum(-kbzN.*repmat(deltaallN,1,samplesizeN)'.*squeeze(khprimevecN(:,:,t))./denoVecN,2)...
        +sum(kbzN.*repmat(deltaallN,1,samplesizeN)'.*khbetaxN.*squeeze(numeVecN(:,:,t))./(denoVecN.^2),2);
    lambdaNum(groupIndex==1,t) = sum(-kbzT.*repmat(deltaallT,1,samplesizeT)'.*squeeze(khprimevecT(:,:,t))./denoVecT,2)...
        +sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetaxT.*squeeze(numeVecT(:,:,t))./(denoVecT.^2),2);
end
lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
khbetax = kernel(btxall/(bandhe*bandadjed),3,'Epanechnikov')/(bandhe*bandadjed)^d;
expectDen = sum(khbetax.*(zdiffall<=0),2)./sum(khbetax,2);
khbetax = kernel(btxall/(bandhe*bandadjed),3,'Epanechnikov')/(bandhe*bandadjed)^d;
for t = 1:(p-d)
    expectVec(:,t) = xall(:,d+t)-sum(khbetax.*repmat(xall(:,d+t),1,samplesize)'.*(zdiffall<=0),2)./sum(khbetax,2)./expectDen;
end

for i = 1:d
    for j = 1:(p-d)
        jaccobValAll(:,(i-1)*(p-d)+j) = (lambdaVec(:,i).*expectVec(:,j));
    end
end
for i = 1:paran
    for j = 1:paran
        jaccobVal(i,j) = sum(deltaallN(groupIndex==0).*jaccobValAll(groupIndex==0,i).*jaccobValAll(groupIndex==0,j))+...
            sum(deltaallT.*jaccobValAll(groupIndex==1,i).*jaccobValAll(groupIndex==1,j));
    end
end

stdVal = sqrt(diag(inv(jaccobVal)));

end