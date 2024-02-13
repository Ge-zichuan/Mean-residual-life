function scoreVal = scoreBetaMRWTrue(parasLower,xall,wall,zdiffallN,zdiffallT,...
    zallN,zallT,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands)
% groupIndex: 1: transplant, 0: nontransplant
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
scoreVal = zeros(length(parasLower(:)),1);
expectVecT = zeros(samplesizeT,p-d);
expectVecN = zeros(samplesizeN,p-d);
%% read bandwidth and bandwidth adjustment
banddh22 = bands(2,2);
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
% calculate kernel values of K\prime(.)

% [temp1,~,~,temp2] = trueMRL(zallT,btransx,'LinLi22W');
% lambdaVec = temp2./repmat(temp1,[1 d]);
lambdaVec = exp(btransx(:,1:(end-1)))./repmat(sum(exp(btransx(:,1:(end-1))),2),[1 d]);

% calculate nonparametric expectations
bandTwo = permute(repmat([bandhe41*bandadjed44*ones(1,d) banddh22*bandadjed46]',[1 samplesizeT samplesizeT]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,3,'norm')/(((bandhe41*bandadjed44)^d)*(banddh22*bandadjed46));
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
        scoreVal((i-1)*(p-d)+j) = sum(deltaallT.*lambdaVec(:,i).*expectVecT(:,j));
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

% [temp1,~,~,temp2] = trueMRL(zallN,btransx,'LinLi21');
% lambdaVec = temp2./repmat(temp1,[1 d]);
lambdaVec = exp(btransx)./repmat(sum(exp(btransx),2),[1 d]);

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
        scoreVal((i-1)*(p-d)+j) = scoreVal((i-1)*(p-d)+j)+...
            sum(deltaallN(groupIndex==0).*lambdaVec(groupIndex==0,i).*expectVecN(groupIndex==0,j));
    end
end
%     scoreVal = scoreVal/samplesize;

end