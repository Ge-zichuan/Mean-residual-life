function [stdVal,jaccobVal] = jaccobBetaMR(parasLower,xall,wall,zdiffallN,zdiffallT,kbzN,kbzT...
    ,deltaallN,deltaallT,deltaall,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands)
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
paran = d*(p-d);
jaccobValAll = zeros(samplesize,length(parasLower(:)));
jaccobVal = zeros(length(parasLower(:)),length(parasLower(:)));
numeVecT = zeros(samplesizeT,samplesizeT,d);
numeVecN = zeros(samplesizeN,samplesizeN,d);
expectVecT = zeros(samplesizeT,p-d);
expectVecN = zeros(samplesizeN,p-d);
% read bandwidth and bandwidth adjustment
% read bandwidth and bandwidth adjustment
%% read bandwidth and bandwidth adjustment
banddh21 = bands(2,1);
banddh22 = bands(2,2);
bandadjdhn23 = bands(2,3);
bandadjdhd24 = bands(2,4);
bandadjdhn25 = bands(2,5);
bandadjdhd26 = bands(2,6);
bandnh31 = bands(3,1);
bandadjnhn32 = bands(3,2);
bandadjnhd33 = bands(3,3);% first 2 for Non-transplant group
bandadjnhn34 = bands(3,4);
bandadjnhd35 = bands(3,5);
bandadjnhd36 = bands(3,6);
bandadjnhd37 = bands(3,7);% last 4 for Non-transplant group
bandhe41 = bands(4,1);
bandadjen42 = bands(4,2);
bandadjed43 = bands(4,3);
bandadjed44 = bands(4,4);
bandadjen45 = bands(4,5);
bandadjed46 = bands(4,6);

%% calculate the score from the transplant group
btransx = [xall(groupIndex==1,:)*[eye(d);reshape(parasLower,p-d,d)] wall(groupIndex==1,end)]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
btxall = permute(repmat(btransx',[1 1 samplesizeT]),[3 2 1])-permute(repmat(btransx,[1 1 samplesizeT]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate estimating equation, each term.

% calculate nonparametric \lambda and \lambda\prime
kerDT = 3;
bandTwo = permute(repmat([banddh21*bandadjdhd24*ones(1,d) banddh22*bandadjdhd26]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((banddh21*bandadjdhd24)^d)*(banddh22*bandadjdhd26));
denoVec = khbetax*(zdiffallT>=0);
bandTwo = permute(repmat([banddh21*bandadjdhn23*ones(1,d) banddh22*bandadjdhn25]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((banddh21*bandadjdhn23)^d)*(banddh22*bandadjdhn25));
lambdaDen = sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetax./denoVec,2);
% calculate kernel values of K\prime(.)
bandTwo = permute(repmat([bandnh31*bandadjnhn34*ones(1,d) banddh22*bandadjnhd36]',[1 samplesizeT samplesizeT]),[3 2 1]);
khprimevec = kernel2(btxall./bandTwo,kerDT,d+1,'Epanechnikov')/(((bandnh31*bandadjnhn34)^d)*(banddh22*bandadjnhd36));
% a1 = ans(:,:,1);a2 = ans(:,:,2);hist(a1(:));s = waitforbuttonpress;hist(a2(:));
for t = 1:d
    numeVecT(:,:,t) = squeeze(khprimevec(:,:,t))*(zdiffallT>=0);
end

lambdaNum = zeros(samplesizeT,d);
bandTwo = permute(repmat([bandnh31*bandadjnhd35*ones(1,d) banddh22*bandadjnhd37]',[1 samplesizeT samplesizeT]),[3 2 1]);
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
bandTwo = permute(repmat([bandhe41*bandadjed44*ones(1,d) banddh22*bandadjed46]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'Epanechnikov')/(((bandhe41*bandadjed44)^d)*(banddh22*bandadjed46));
expectDen = sum(khbetax.*(zdiffallT<=0),2)./sum(khbetax,2);
bandTwo = permute(repmat([bandhe41*bandadjen45*ones(1,d) banddh22*bandadjed46]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'Epanechnikov')/(((bandhe41*bandadjen45)^d)*(banddh22*bandadjed46));
for t = 1:(p-d)
    expectVecT(:,t) = xall(groupIndex==1,d+t)-sum(khbetax.*repmat(xall(groupIndex==1,d+t),1,samplesizeT)'.*(zdiffallT<=0),2)./sum(khbetax,2)./expectDen;
end
% lambdaVec(isnan(lambdaVec))=0;
% lambdaVec(isinf(lambdaVec))=0;
for i = 1:d
    for j = 1:(p-d)
        jaccobValAll(groupIndex==1,(i-1)*(p-d)+j) = ...
            lambdaVec(:,i).*expectVecT(:,t);
    end
end
%% calculate the score from the non-transplant group
btransx = xall*[eye(d);reshape(parasLower,p-d,d)]; % paras1 include upper d-by-d block.
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
khbetax = kernel(btxall/(banddh21*bandadjdhd24),kerD,'Epanechnikov')/(banddh21*bandadjdhd24)^d;
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
khbetax = kernel(btxall/(bandhe41*bandadjed43),3,'Epanechnikov')/(bandhe41*bandadjed43)^d;
expectDen = sum(khbetax.*(zdiffallN<=0),2)./sum(khbetax,2);
khbetax = kernel(btxall/(bandhe41*bandadjen42),3,'Epanechnikov')/(bandhe41*bandadjen42)^d;
for t = 1:(p-d)
    expectVecN(:,t) = xall(:,d+t)-sum(khbetax.*repmat(xall(:,d+t),1,samplesizeN)'.*(zdiffallN<=0),2)./sum(khbetax,2)./expectDen;
end
%% Combine the two groups as one score function
for i = 1:d
    for j = 1:(p-d)
        jaccobValAll(groupIndex==0,(i-1)*(p-d)+j) = ...
            lambdaVec(groupIndex==0,i).*expectVecN(groupIndex==0,j);
    end
end
for i = 1:paran
    for j = 1:paran
        jaccobVal(i,j) = sum(deltaall.*jaccobValAll(:,i).*jaccobValAll(:,j));
    end
end

stdVal = sqrt(diag(inv(jaccobVal)));

end