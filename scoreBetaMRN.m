function scoreVal = scoreBetaMRN(parasLower,xall,zdiffallN,...
    kbzN,deltaallN,samplesizeN,groupIndex,p,d,bands)
% groupIndex: 1: transplant, 0: nontransplant
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
numeVecN = zeros(samplesizeN,samplesizeN,d);
scoreVal = zeros(length(parasLower(:)),1);
expectVecN = zeros(samplesizeN,p-d);
%% read bandwidth and bandwidth adjustment
banddh21 = bands(2,1);
bandadjdhn23 = bands(2,3);
bandadjdhd24 = bands(2,4);
bandnh31 = bands(3,1);
bandadjnhn32 = bands(3,2);
bandadjnhd33 = bands(3,3);% first 2 for Non-transplant group
bandhe41 = bands(4,1);
bandadjen42 = bands(4,2);
bandadjed43 = bands(4,3);
%% calculate the score from the non-transplant group
btransx = xall*[eye(d);parasLower]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    btxall = permute(repmat(btransx',[1 1 samplesizeN]),[3 2 1])-permute(repmat(btransx,[1 1 samplesizeN]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 samplesizeN])-repmat(btransx',[samplesizeN 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
khbetax = kernel(btxall/(banddh21*bandadjdhd24),kerD,'norm')/(banddh21*bandadjdhd24)^d;
denoVec = khbetax*(zdiffallN>=0);
khbetax = kernel(btxall/(banddh21*bandadjdhn23),kerD,'Epanechnikov')/(banddh21*bandadjdhn23)^d;
lambdaDen = sum(kbzN.*repmat(deltaallN,1,samplesizeN)'.*khbetax./denoVec,2);
% calculate kernel values of K\prime(.)
khprimevec = kernel2(btxall/(bandnh31*bandadjnhn32),kerD,d,'Epanechnikov')/(bandnh31*bandadjnhn32)^(d+1);
% a1 = ans(:,:,1);hist(a1(:));a2 = ans(:,:,2);s = waitforbuttonpress;hist(a2(:));
for t = 1:d
    numeVecN(:,:,t) = squeeze(khprimevec(:,:,t))*(zdiffallN>=0);
end

lambdaNum = zeros(samplesizeN,d);
khbetax = kernel(btxall/(bandnh31*bandadjnhd33),kerD,'Epanechnikov')/(bandnh31*bandadjnhd33)^d;
for t = 1:d
    lambdaNum(:,t) = sum(-kbzN.*repmat(deltaallN,1,samplesizeN)'.*squeeze(khprimevec(:,:,t))./denoVec,2)...
        +sum(kbzN.*repmat(deltaallN,1,samplesizeN)'.*khbetax.*squeeze(numeVecN(:,:,t))./(denoVec.^2),2);
end
lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
% lambdaVec(isnan(lambdaVec))=0;
% lambdaVec(isinf(lambdaVec))=0;
%     [temp1,~,~,temp2] = trueMRL(zall,btransx,'mimicN');
%     lambdaVec = temp2./temp1;
% calculate nonparametric expectations
khbetax = kernel(btxall/(bandhe41*bandadjed43),3,'norm')/(bandhe41*bandadjed43)^d;
expectDen = sum(khbetax.*(zdiffallN<=0),2)./sum(khbetax,2);
khbetax = kernel(btxall/(bandhe41*bandadjen42),3,'Epanechnikov')/(bandhe41*bandadjen42)^d;
for t = 1:(p-d)
    expectVecN(:,t) = xall(:,d+t)-sum(khbetax.*repmat(xall(:,d+t),1,samplesizeN)'.*(zdiffallN<=0),2)./sum(khbetax,2)./expectDen;
end
%% Combine the two groups as one score function
for i = 1:d
    for j = 1:(p-d)
        scoreVal((i-1)*(p-d)+j) = ...
            sum(deltaallN(groupIndex==0).*lambdaVec(groupIndex==0,i).*expectVecN(groupIndex==0,j));
    end
end
%     scoreVal = scoreVal/samplesize;

end
