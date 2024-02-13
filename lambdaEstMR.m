function lambdaVec = lambdaEstMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,p,d,group)
estSample = length(inputg(:,1));
% read bandwidth and bandwidth adjustment
if group == 1
    bandb = bands(1,2);
    bandadjb = bands(1,4);
    banddh21 = bands(2,1);
    banddh22 = bands(2,2);
    bandadjdhn23 = bands(2,3);
    bandadjdhd24 = bands(2,4);
    bandadjdhn25 = bands(2,5);
    bandadjdhd26 = bands(2,6);
    zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
    zdiffest = zall*ones(1,estSample)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
    kbz = kernel(zdiffest/(bandb*bandadjb),1,'Epanechnikov')/(bandb*bandadjb);% {i,j} = Z_i-Z_j
    btransx = [xall(:,1:(end-1))*[eye(d);reshape(parasLower,p-d,d)] xall(:,end)]; % paras1 include upper d-by-d block.
    % calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
    btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    kerDT = 3;
    bandTwo = permute(repmat([banddh21*bandadjdhd24*ones(1,d) banddh22*bandadjdhd26]',[1 estSample samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,kerDT,'norm')/(((banddh21*bandadjdhd24)^d)*(banddh22*bandadjdhd26));
    denoVec = (zdiffall<=0)*khbetax;
    bandTwo = permute(repmat([banddh21*bandadjdhn23*ones(1,d) banddh22*bandadjdhn25]',[1 estSample samplesize]),[3 2 1]);
    khbetax = kernel(btxall./bandTwo,kerDT,'Epanechnikov')/(((banddh21*bandadjdhn23)^d)*(banddh22*bandadjdhn25));
% a1 = ans(:,:,1);hist(a1(:));a2 = ans(:,:,2);s = waitforbuttonpress;hist(a2(:));
    lambdaVec = sum(kbz.*repmat(deltaall,1,estSample).*khbetax./denoVec,1)';
    
else
    bandb = bands(1,1);
    bandadjb = bands(1,3);
    banddh21 = bands(2,1);
    bandadjdhn23 = bands(2,3);
    bandadjdhd24 = bands(2,4);
    zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
    zdiffest = zall*ones(1,estSample)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
    kbz = kernel(zdiffest/(bandb*bandadjb),1,'Epanechnikov')/(bandb*bandadjb);% {i,j} = Z_i-Z_j
    btransx = xall*[eye(d);reshape(parasLower,p-d,d)]; % paras1 include upper d-by-d block.
    % calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
    if (d>1)
        btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    else
        btxall = repmat(btransx,[1 estSample])-repmat(inputg',[samplesize 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
    end
    % calculate nonparametric \lambda and \lambda\prime
    if (d>1)
        kerD = 3;
    else
        kerD = 1;
    end
    khbetax = kernel(btxall/(banddh21*bandadjdhd24),kerD,'norm')/(banddh21*bandadjdhd24)^d;
    denoVec = (zdiffall<=0)*khbetax;
    khbetax = kernel(btxall/(banddh21*bandadjdhn23),kerD,'Epanechnikov')/(banddh21*bandadjdhn23)^d;
    lambdaVec = sum(kbz.*repmat(deltaall,1,estSample).*khbetax./denoVec,1)';
end

end