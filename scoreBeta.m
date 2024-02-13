function scoreVal = scoreBeta(parasLower,xall,zdiffall,kbz,deltaall,samplesize,q,d,bands)
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
numeVec = zeros(samplesize,samplesize,d);
scoreVal = length(parasLower(:));
expectVec = zeros(samplesize,q-d);
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
khbetax = kernel(btxall/banddh/bandadjdhd,kerD,'norm')/banddh/bandadjdhd;
denoVec = khbetax*(zdiffall>=0);
khbetax = kernel(btxall/banddh/bandadjdhn,kerD,'norm')/banddh/bandadjdhn;
lambdaDen = sum(kbz.*repmat(deltaall,1,samplesize)'.*khbetax./denoVec,2);
% calculate kernel values of K\prime(.)
khprimevec = kernel2(btxall/bandnh/bandadjnhn,kerD,d,'quar')/bandnh/bandadjnhn;
for t = 1:d
    numeVec(:,:,t) = squeeze(khprimevec(:,:,t))*(zdiffall>=0);
end
lambdaNum = zeros(samplesize,d);
khbetax = kernel(btxall/bandnh/bandadjnhd,kerD,'norm')/bandnh/bandadjnhd;
for t = 1:d
    lambdaNum(:,t) = sum(-kbz.*repmat(deltaall,1,samplesize)'.*squeeze(khprimevec(:,:,t))./denoVec,2)...
        +sum(kbz.*repmat(deltaall,1,samplesize)'.*khbetax.*squeeze(numeVec(:,:,t))./(denoVec.^2),2);
end
lambdaVec = lambdaNum./repmat(lambdaDen,1,d);
% lambdaVec = exp(btransx)./repmat(sum(exp(btransx),2),1,d);
% lambdaVec(isnan(lambdaVec)) = 0;
% calculate nonparametric expectations
khbetax = kernel(btxall/bandhe/bandadjed,3,'norm')/bandhe/bandadjed;
expectDen = sum(khbetax.*(zdiffallT<=0),2)./sum(khbetax,2);
khbetax = kernel(btxall/bandhe/bandadjen,3,'norm')/bandhe/bandadjen;
for t = 1:(p-d)
    expectVec(:,t) = xall(:,d+t)-sum(khbetax.*repmat(xall(:,d+t),1,samplesize)'.*(zdiffall<=0),2)./sum(khbetax,2)./expectDen;
end

for i = 1:d
    for j = 1:(q-d)
        scoreVal((i-1)*(q-d)+j) = sum(deltaall.*lambdaVec(:,i).*expectVec(:,j));
    end
end
scoreVal = scoreVal/samplesize;

end