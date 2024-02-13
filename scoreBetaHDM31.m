function scoreVal = scoreBetaHDM31(parasLower,xall,zdiffall,aall,aallmat,kbz1,kbz0,deltaall,samplesize,goodI,q,d,bands)
% This version of score function includes the treatment variable A.
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
numeVec1 = zeros(samplesize,samplesize,d);
numeVec0 = zeros(samplesize,samplesize,d);
scoreVal = length(parasLower(:));
expectVec1 = zeros(samplesize,q-d);
expectVec0 = zeros(samplesize,q-d);
deltaallmat = repmat(deltaall,1,samplesize);
goodIndex1 = (goodI{1});
goodIndex0 = (goodI{2});
goodIndexAll = goodI{3};
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
bandnh = bands(3,1);
bandadjnhn = bands(3,2);
bandadjnhd = bands(3,3);
bandhe = bands(4,1);
bandadjen = bands(4,2);
bandadjed = bands(4,3);
% calculate \bb\trans\x
btransx = xall*[eye(d);parasLower]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    btxall = permute(repmat(btransx',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 samplesize]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 samplesize])-repmat(btransx',[samplesize 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
% calculate estimating equation, each term.

% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
khbetax = kernel(btxall/(banddh*bandadjdhd),kerD,'norm')/(banddh*bandadjdhd)^d;
denoVec1 = khbetax*((zdiffall>=0).*aallmat);
denoVec1 = denoVec1(:,goodIndex1);
denoVec0 = khbetax*((zdiffall>=0).*(1-aallmat));
denoVec0 = denoVec0(:,goodIndex0);
khbetax = kernel(btxall/(banddh*bandadjdhn),kerD,'norm')/(banddh*bandadjdhn)^d;
lambdaDen1 = sum(aallmat(goodIndex1,:)'.*kbz1(:,goodIndex1).*deltaallmat(goodIndex1,:)'.*khbetax(:,goodIndex1)./denoVec1,2);
lambdaDen0 = sum((1-aallmat(goodIndex0,:))'.*kbz0(:,goodIndex0).*deltaallmat(goodIndex0,:)'.*khbetax(:,goodIndex0)./denoVec0,2);
lambdaDen = aall.*lambdaDen1+(1-aall).*lambdaDen0;
% calculate kernel values of K\prime(.)
khprimevec = kernel2(btxall/(bandnh*bandadjnhn),kerD,d,'quar')/(bandnh*bandadjnhn)^(d+1);
for t = 1:d
    numeVec1(:,:,t) = squeeze(khprimevec(:,:,t))*((zdiffall>=0).*aallmat);
    numeVec0(:,:,t) = squeeze(khprimevec(:,:,t))*((zdiffall>=0).*(1-aallmat));
end
lambdaNum1 = zeros(samplesize,d);
lambdaNum0 = zeros(samplesize,d);
khbetax = kernel(btxall/(bandnh*bandadjnhd),kerD,'norm')/(bandnh*bandadjnhd)^d;
for t = 1:d
    lambdaNum1(:,t) = sum(-aallmat(goodIndex1,:)'.*kbz1(:,goodIndex1).*deltaallmat(goodIndex1,:)'.*squeeze(khprimevec(:,goodIndex1,t))./denoVec1,2)...
        +sum(aallmat(goodIndex1,:)'.*kbz1(:,goodIndex1).*deltaallmat(goodIndex1,:)'.*khbetax(:,goodIndex1).*squeeze(numeVec1(:,goodIndex1,t))./(denoVec1.^2),2);
    lambdaNum0(:,t) = sum(-(1-aallmat(goodIndex0,:))'.*kbz0(:,goodIndex0).*deltaallmat(goodIndex0,:)'.*squeeze(khprimevec(:,goodIndex0,t))./denoVec0,2)...
        +sum((1-aallmat(goodIndex0,:))'.*kbz0(:,goodIndex0).*deltaallmat(goodIndex0,:)'.*khbetax(:,goodIndex0).*squeeze(numeVec0(:,goodIndex0,t))./(denoVec0.^2),2);
end
lambdaNum = repmat(aall,1,d).*lambdaNum1+(repmat(1-aall,1,d).*lambdaNum0);
lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
% lambdaVec = exp(btransx)./repmat(sum(exp(btransx),2),1,d);
% lambdaVec(isnan(lambdaVec)) = 0;
% calculate nonparametric expectations
khbetax = kernel(btxall/(bandhe*bandadjed),kerD,'norm')/(bandhe*bandadjed)^d;
expectDen1 = sum(aallmat'.*khbetax.*(zdiffall<=0),2);
expectDen0 = sum((1-aallmat)'.*khbetax.*(zdiffall<=0),2);
khbetax = kernel(btxall/(bandhe*bandadjen),kerD,'norm')/(bandhe*bandadjen)^d;
for t = 1:(q-d)
    expectVec1(:,t) = xall(:,d+t)-sum(aallmat'.*khbetax.*repmat(xall(:,d+t),1,samplesize)'.*(zdiffall<=0),2)...
        ./expectDen1;
    expectVec0(:,t) = xall(:,d+t)-sum((1-aallmat)'.*khbetax.*repmat(xall(:,d+t),1,samplesize)'.*(zdiffall<=0),2)...
        ./expectDen0;
end
expectVec = repmat(aall,1,q-d).*expectVec1+repmat(1-aall,1,q-d).*expectVec0;
for i = 1:d
    for j = 1:(q-d)
        scoreVal((i-1)*(q-d)+j) = sum(deltaall(goodIndexAll).*lambdaVec(goodIndexAll,i).*expectVec(goodIndexAll,j));
    end
end
scoreVal = scoreVal/samplesize;%;

end