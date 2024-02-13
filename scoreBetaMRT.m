function scoreVal = scoreBetaMRT(parasLower,xall,wall,zdiffallT,...
    kbzT,deltaallT,samplesizeT,groupIndex,p,d,bands)
% groupIndex: 1: transplant, 0: nontransplant
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
numeVecT = zeros(samplesizeT,samplesizeT,d);
scoreVal = zeros(length(parasLower(:)),1);
expectVecT = zeros(samplesizeT,p-d);
%% read bandwidth and bandwidth adjustment
banddh21 = bands(2,1);
banddh22 = bands(2,2);
bandadjdhn23 = bands(2,3);
bandadjdhd24 = bands(2,4);
bandadjdhn25 = bands(2,5);
bandadjdhd26 = bands(2,6);
bandnh31 = bands(3,1);
bandadjnhn34 = bands(3,4);
bandadjnhd35 = bands(3,5);
bandadjnhd36 = bands(3,6);
bandadjnhd37 = bands(3,7);% last 4 for Non-transplant group
bandhe41 = bands(4,1);
bandadjed44 = bands(4,4);
bandadjen45 = bands(4,5);
bandadjed46 = bands(4,6);
%% calculate the score from the transplant group
btransx = [xall(groupIndex==1,:)*[eye(d);parasLower] wall(groupIndex==1,end)]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
btxall = permute(repmat(btransx',[1 1 samplesizeT]),[3 2 1])-permute(repmat(btransx,[1 1 samplesizeT]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate estimating equation, each term.

% calculate nonparametric \lambda and \lambda\prime
kerDT = 3;
bandTwo = permute(repmat([banddh21*bandadjdhd24 banddh22*bandadjdhd26]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerDT,'norm')/(((banddh21*bandadjdhd24)^d)*(banddh22*bandadjdhd26));
denoVec = khbetax*(zdiffallT>=0);
bandTwo = permute(repmat([banddh21*bandadjdhn23 banddh22*bandadjdhn25]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((banddh21*bandadjdhn23)^d)*(banddh22*bandadjdhn25));
lambdaDen = sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetax./denoVec,2);
% calculate kernel values of K\prime(.)
bandTwo = permute(repmat([bandnh31*bandadjnhn34 banddh22*bandadjnhd36]',[1 samplesizeT samplesizeT]),[3 2 1]);
khprimevec = kernel2(btxall./bandTwo,kerDT,d+1,'Epanechnikov')/(((bandnh31*bandadjnhn34)^d)*(banddh22*bandadjnhd36));
% a1 = ans(:,:,1);a2 = ans(:,:,2);hist(a1(:));s = waitforbuttonpress;hist(a2(:));
for t = 1:d
    numeVecT(:,:,t) = squeeze(khprimevec(:,:,t))*(zdiffallT>=0);
end

lambdaNum = zeros(samplesizeT,d);
bandTwo = permute(repmat([bandnh31*bandadjnhd35 banddh22*bandadjnhd37]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((bandnh31*bandadjnhd35)^d)*(banddh22*bandadjnhd37));
for t = 1:d
    lambdaNum(:,t) = sum(-kbzT.*repmat(deltaallT,1,samplesizeT)'.*squeeze(khprimevec(:,:,t))./denoVec,2)...
        +sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetax.*squeeze(numeVecT(:,:,t))./(denoVec.^2),2);
end
lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
%     lambdaVec(isnan(lambdaVec))=0;
%     lambdaVec(isinf(lambdaVec))=0;

%     [temp1,~,~,temp2] = trueMRL(zall,btransx,'mimicT');
%     lambdaVec = temp2./temp1;

% calculate nonparametric expectations
bandTwo = permute(repmat([bandhe41*bandadjed44 banddh22*bandadjed46]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'norm')/(((bandhe41*bandadjed44)^d)*(banddh22*bandadjed46));
expectDen = sum(khbetax.*(zdiffallT<=0),2)./sum(khbetax,2);
bandTwo = permute(repmat([bandhe41*bandadjen45 banddh22*bandadjed46]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'Epanechnikov')/(((bandhe41*bandadjen45)^d)*(banddh22*bandadjed46));
for t = 1:(p-d)
    expectVecT(:,t) = xall(groupIndex==1,d+t)-sum(khbetax.*repmat(xall(groupIndex==1,d+t),1,samplesizeT)'.*(zdiffallT<=0),2)./sum(khbetax,2)./expectDen;
end
% lambdaVec(isnan(lambdaVec))=0;
% lambdaVec(isinf(lambdaVec))=0;
for i = 1:d
    for j = 1:(p-d)
        scoreVal((i-1)*(p-d)+j) = sum(deltaallT.*lambdaVec(:,i).*expectVecT(:,j));
    end
end
%     scoreVal = scoreVal/samplesize;

end