function scoreVal = scoreBetaMRH(parasLower,xall,zall,zdiffallN,zdiffallT,zdiffall,kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands)
% groupIndex: 1: transplant, 0: nontransplant
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
numeVecT = zeros(samplesizeT,samplesizeT,d);
numeVecN = zeros(samplesizeN,samplesizeN,d);
lambdaDen = zeros(samplesize,1);
scoreVal = length(parasLower(:));
expectVec = zeros(samplesize,p-d);
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
btransxT = btransx(groupIndex==1,:); % paras1 include upper d-by-d block.
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
% khbetaxT = kernel(btxallT/banddh/bandadjdhd,kerD,'norm')/banddh/bandadjdhd;
% khbetaxN = kernel(btxallN/banddh/bandadjdhd,kerD,'norm')/banddh/bandadjdhd;
% denoVecT = khbetaxT*(zdiffallT>=0);
% denoVecN = khbetaxN*(zdiffallN>=0);
% khbetaxT = kernel(btxallT/banddh/bandadjdhn,kerD,'norm')/banddh/bandadjdhn;
% khbetaxN = kernel(btxallN/banddh/bandadjdhn,kerD,'norm')/banddh/bandadjdhn;
% lambdaDen(:,1) = sum(kbzN.*repmat(deltaallN,1,samplesizeN)'.*khbetaxN./denoVecN,2);
% lambdaDen(groupIndex==1,1) = sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetaxT./denoVecT,2);
% % calculate kernel values of K\prime(.)
% khprimevecT = kernel2(btxallT/bandnh/bandadjnhn,kerD,d,'quar')/bandnh/bandadjnhn;
% khprimevecN = kernel2(btxallN/bandnh/bandadjnhn,kerD,d,'quar')/bandnh/bandadjnhn;
% for t = 1:d
%     numeVecT(:,:,t) = squeeze(khprimevecT(:,:,t))*(zdiffallT>=0);
%     numeVecN(:,:,t) = squeeze(khprimevecN(:,:,t))*(zdiffallN>=0);
% end
% lambdaNum = zeros(samplesize,d);
% khbetaxT = kernel(btxallT/bandnh/bandadjnhd,kerD,'norm')/bandnh/bandadjnhd;
% khbetaxN = kernel(btxallN/bandnh/bandadjnhd,kerD,'norm')/bandnh/bandadjnhd;
% for t = 1:d
%     lambdaNum(:,t) = sum(-kbzN.*repmat(deltaallN,1,samplesizeN)'.*squeeze(khprimevecN(:,:,t))./denoVecN,2)...
%         +sum(kbzN.*repmat(deltaallN,1,samplesizeN)'.*khbetaxN.*squeeze(numeVecN(:,:,t))./(denoVecN.^2),2);
%     lambdaNum(groupIndex==1,t) = sum(-kbzT.*repmat(deltaallT,1,samplesizeT)'.*squeeze(khprimevecT(:,:,t))./denoVecT,2)...
%         +sum(kbzT.*repmat(deltaallT,1,samplesizeT)'.*khbetaxT.*squeeze(numeVecT(:,:,t))./(denoVecT.^2),2);
% end
% % lambdaVec(groupIndex==0,i) = ones(sum(groupIndex==0),d);
% % lambdaVec(groupIndex==1,i) = ones(sum(groupIndex==1),d);
lambdaVec = (-2*sum(exp(btransx)/2,2)./(sum(exp(btransx)/2,2).^2+zall.^2))*ones(1,2).*exp(btransx)/2;
tempstar = log(zall(groupIndex==1))-5+0.1*sum((1-sqrt(2)*btransxT).^2,2);
lambdaVec(groupIndex==1,:) = 2*0.1*sqrt(2)*(1-sqrt(2)*btransxT).*(((tempstar.*normpdf(tempstar,0,1)./zall(groupIndex==1).*(1-normcdf(tempstar))-...
    normpdf(tempstar).^2./zall(groupIndex==1))./(normcdf(-tempstar,0,1).^2))*ones(1,2));
% calculate nonparametric expectations
% khbetax = kernel(btxallT/bandhe/bandadjed,3,'norm')/bandhe/bandadjed;
% expectDenT = sum(khbetax.*(zdiffallT<=0),2)./sum(khbetax,2);
% khbetaxT = kernel(btxallT/bandhe/bandadjen,3,'norm')/bandhe/bandadjen;
% khbetax = kernel(btxallN/bandhe/bandadjed,3,'norm')/bandhe/bandadjed;
% expectDenN = sum(khbetax.*(zdiffallN<=0),2)./sum(khbetax,2);
% khbetaxN = kernel(btxallN/bandhe/bandadjen,3,'norm')/bandhe/bandadjen;
khbetax = kernel(btxall/bandhe/bandadjed,3,'norm')/bandhe/bandadjed;
expectDen = sum(khbetax.*(zdiffall<=0),2)./sum(khbetax,2);
khbetax = kernel(btxall/bandhe/bandadjen,3,'norm')/bandhe/bandadjen;
for t = 1:(p-d)
    expectVec(:,t) = xall(:,d+t)-sum(khbetax.*repmat(xall(:,d+t),1,samplesize)'.*(zdiffall<=0),2)./sum(khbetax,2)./expectDen;
%     expectVec(:,t) = xallN(:,d+t)-sum(khbetaxN.*repmat(xallN(:,d+t),1,samplesizeN)'.*(zdiffallN<=0),2)./sum(khbetaxN,2)./expectDenN;
%     expectVec(groupIndex==1,t) = xallT(:,d+t)-sum(khbetaxT.*repmat(xallT(:,d+t),1,samplesizeT)'.*(zdiffallT<=0),2)./sum(khbetaxT,2)./expectDenT;
end

for i = 1:d
    for j = 1:(p-d)
        scoreVal((i-1)*(p-d)+j) = sum(deltaallT.*lambdaVec(groupIndex==1,i).*expectVec(groupIndex==1,j))+...
            sum(deltaallN(groupIndex==0).*lambdaVec(groupIndex==0,i).*expectVec(groupIndex==0,j));
    end
end
scoreVal = scoreVal/samplesizeN;

end