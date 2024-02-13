function scoreVal = scoreBetaHazardOnly(parasLower,xall,zdiffall,zall,...
    kbz,deltaall,samplesize,p,d,bands,group)
% groupIndex: 1: transplant, 0: nontransplant
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
numeVec = zeros(samplesize,samplesize,d+group);
lambdaDen = zeros(samplesize,1);
scoreVal = length(parasLower(:));
expectVec = zeros(samplesize,p-d);
% read bandwidth and bandwidth adjustment
banddh21 = bands(2,1);
banddh25 = bands(2,5);
bandadjdhn22 = bands(2,2);
bandadjdhd23 = bands(2,3);
bandadjdhd24 = bands(2,4);
bandnh31 = bands(3,1);
bandadjnhn32 = bands(3,2);
bandadjnhd33 = bands(3,3);
bandadjnhd34 = bands(3,4);
bandhe41 = bands(4,1);
bandadjen42 = bands(4,2);
bandadjed43 = bands(4,3);
bandadjed44 = bands(4,4);
% calculate \bb\trans\x
if group == 1
    btransx = [xall(:,1:(end-1))*[eye(d);parasLower] xall(:,end)]; % paras1 include upper d-by-d block.
    % calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
    btxall = permute(repmat(btransx',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 samplesize]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    % calculate estimating equation, each term.

    % calculate nonparametric \lambda and \lambda\prime
    kerDT = 3;
    bandTwo = permute(repmat([banddh21*bandadjdhd23 banddh25*bandadjdhd24]',[1 samplesize samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((banddh21*bandadjdhd23)^d)*(banddh25*bandadjdhd24));
    denoVec = khbetax*(zdiffall>=0);
    bandTwo = permute(repmat([banddh21*bandadjdhn22 banddh25*bandadjdhd24]',[1 samplesize samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((banddh21*bandadjdhn22)^d)*(banddh25*bandadjdhd24));
    lambdaDen = sum(kbz.*repmat(deltaall,1,samplesize)'.*khbetax./denoVec,2);
    % calculate kernel values of K\prime(.)
    bandTwo = permute(repmat([bandnh31*bandadjnhn32 banddh25*bandadjnhd34]',[1 samplesize samplesize]),[3 2 1]);
    khprimevec = kernel2(btxall./bandTwo,kerDT,d+1,'Epanechnikov')/(((bandnh31*bandadjnhn32)^d)*(banddh25*bandadjnhd34));
    % a1 = ans(:,:,1);a2 = ans(:,:,2);hist(a1(:));s = waitforbuttonpress;hist(a2(:));
    for t = 1:d
        numeVec(:,:,t) = squeeze(khprimevec(:,:,t))*(zdiffall>=0);
    end

    lambdaNum = zeros(samplesize,d);
    bandTwo = permute(repmat([bandnh31*bandadjnhd33 banddh25*bandadjnhd34]',[1 samplesize samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((bandnh31*bandadjnhd33)^d)*(banddh25*bandadjnhd34));
    for t = 1:d
        lambdaNum(:,t) = sum(-kbz.*repmat(deltaall,1,samplesize)'.*squeeze(khprimevec(:,:,t))./denoVec,2)...
            +sum(kbz.*repmat(deltaall,1,samplesize)'.*khbetax.*squeeze(numeVec(:,:,t))./(denoVec.^2),2);
    end
    lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
%     lambdaVec(isnan(lambdaVec))=0;
%     lambdaVec(isinf(lambdaVec))=0;

%     [temp1,~,~,temp2] = trueMRL(zall,btransx,'mimicT');
%     lambdaVec = temp2./temp1;

    % calculate nonparametric expectations
    bandTwo = permute(repmat([bandhe41*bandadjed43 banddh25*bandadjed44]',[1 samplesize samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,3,'Epanechnikov')/(((bandhe41*bandadjed43)^d)*(banddh25*bandadjed44));
    expectDen = sum(khbetax.*(zdiffall<=0),2)./sum(khbetax,2);
    bandTwo = permute(repmat([bandhe41*bandadjen42 banddh25*bandadjed44]',[1 samplesize samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,3,'Epanechnikov')/(((bandhe41*bandadjen42)^d)*(banddh25*bandadjed44));
else
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
    khbetax = kernel(btxall/(banddh21*bandadjdhd23),kerD,'Epanechnikov')/(banddh21*bandadjdhd23)^d;
    denoVec = khbetax*(zdiffall>=0);
    khbetax = kernel(btxall/(banddh21*bandadjdhn22),kerD,'Epanechnikov')/(banddh21*bandadjdhn22)^d;
    lambdaDen = sum(kbz.*repmat(deltaall,1,samplesize)'.*khbetax./denoVec,2);
    % calculate kernel values of K\prime(.)
    khprimevec = kernel2(btxall/(bandnh31*bandadjnhn32),kerD,d,'Epanechnikov')/(bandnh31*bandadjnhn32)^(d+1);
    % a1 = ans(:,:,1);hist(a1(:));a2 = ans(:,:,2);s = waitforbuttonpress;hist(a2(:));
    for t = 1:d
        numeVec(:,:,t) = squeeze(khprimevec(:,:,t))*(zdiffall>=0);
    end

    lambdaNum = zeros(samplesize,d);
    khbetax = kernel(btxall/(bandnh31*bandadjnhd33),kerD,'Epanechnikov')/(bandnh31*bandadjnhd33)^d;
    for t = 1:d
        lambdaNum(:,t) = sum(-kbz.*repmat(deltaall,1,samplesize)'.*squeeze(khprimevec(:,:,t))./denoVec,2)...
            +sum(kbz.*repmat(deltaall,1,samplesize)'.*khbetax.*squeeze(numeVec(:,:,t))./(denoVec.^2),2);
    end
    lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
%     [temp1,~,~,temp2] = trueMRL(zall,btransx,'mimicN');
%     lambdaVec = temp2./temp1;
    % calculate nonparametric expectations
    khbetax = kernel(btxall/(bandhe41*bandadjed43),3,'Epanechnikov')/(bandhe41*bandadjed43)^d;
    expectDen = sum(khbetax.*(zdiffall<=0),2);%./sum(khbetax,2);
    khbetax = kernel(btxall/(bandhe41*bandadjen42),3,'Epanechnikov')/(bandhe41*bandadjen42)^d;
end
for t = 1:(p-d)
    expectVec(:,t) = xall(:,d+t)-sum(khbetax.*repmat(xall(:,d+t),1,samplesize)'.*(zdiffall<=0),2)./expectDen;
end
% lambdaVec(isnan(lambdaVec))=0;
% lambdaVec(isinf(lambdaVec))=0;
for i = 1:d
    for j = 1:(p-d)
        scoreVal((i-1)*(p-d)+j) = sum(deltaall.*lambdaVec(:,i).*expectVec(:,j));
    end
end
% scoreVal = scoreVal/samplesize;

end