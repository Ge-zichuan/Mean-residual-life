function lambdaVec = lambdaD1EstMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,d,group)
estSample = length(inputg(:,1));
% read bandwidth and bandwidth adjustment
numeVec = zeros(samplesize,estSample,d);
if group == 1
    bandb = bands(1,2);
    bandadjb = bands(1,4);
else
    bandb = bands(1,1);
    bandadjb = bands(1,3);
end
bandnh = bands(3,1);
bandadjnhn = bands(3,2);
bandadjnhd = bands(3,3);
% calculate all Z_i-estZ
% zdiffall = inputt*ones(1,samplesize)-ones(estSample,1)*zall';% {i,j} = Z_i-Z_j
% kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
zdiffest = zall*ones(1,estSample)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffest/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
btransx = xall*[eye(d);parasLower]; % paras1 include upper d-by-d block.
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
khbetax = kernel(btxall/(bandnh*bandadjnhd),kerD,'norm')/(bandnh*bandadjnhd);
denoVec = (zdiffall<=0)*khbetax;
% calculate kernel values of K\prime(.)
khprimevec = kernel2(btxall/(bandnh*bandadjnhn),kerD,d,'quar')/(bandnh*bandadjnhn)^(d+1);
for t = 1:d
    numeVec(:,:,t) = (zdiffall>=0)*squeeze(khprimevec(:,:,t));
end

lambdaVec = zeros(estSample,d);
% khbetax = kernel(btxall/bandnh/bandadjnhd,kerD,'norm')/bandnh/bandadjnhd;
for t = 1:d
    lambdaVec(:,t) = sum(-kbz.*repmat(deltaall,1,estSample).*squeeze(khprimevec(:,:,t))./denoVec,1)'...
        +sum(kbz.*repmat(deltaall,1,estSample).*khbetax.*squeeze(numeVec(:,:,t))./(denoVec.^2),1)';
end

end