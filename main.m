rng('shuffle');
%%% basic information of numerical study
simulations = 1000; % number of simulations
samplesize = 2000;
p = 6;
d = 2;
% model setting:p, d, samplesize, covariates, coefficient \beta, transplant
% model, waitlist model, nontransplant model, censoring, censoring rate
options = {num2str(p);num2str(d);num2str(samplesize);'norm';'ma62';'LinLi21';'unif';'LinLi22W';'gamm';'0.0'};

%% Assigning stores
xallsave = zeros(simulations,samplesize,p);
zallsave = zeros(simulations,samplesize);
tallsave = zeros(simulations,samplesize);
cenallsave = zeros(simulations,samplesize);
wallsave = zeros(simulations,samplesize);
deltaallsave = zeros(simulations,samplesize);
groupIndexsave = zeros(simulations,samplesize);
censorInfo = zeros(simulations,1);
vicResult = zeros(simulations,1);
betaEstResult = zeros(simulations,d*(p-d));
betaStdResult = zeros(simulations,d*(p-d));
rhoResult = zeros(simulations,1);
flags = zeros(simulations,2);
% bandwidth tuning for kernel estimation
bands = zeros(4,7);
bands(1,3:4) =    [10,10];
bands(2,3:6) =    [50,50,50,20];
bands(3,2:7) = [50,50,50,50,50,20];
bands(4,2:6) = [8,5,5,5,10];

% the grid for predict mean residual life
inputt = 0:0.1:1;
inputbbtx = repmat((-2:(2+2)/10:2)',[1 d]);
combOfVarsN = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsN = combvec(combOfVarsN',inputbbtx(:,i)')';
    end
end
inputt = 0:0.5:5;
inputw = 0:0.1:1;
inputbbtx = repmat((-2:(2+2)/10:2)',[1 d]);
combOfVarsT = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsT = combvec(combOfVarsT',inputbbtx(:,i)')';
    end
end
combOfVarsT = combvec(combOfVarsT',inputw)';
mResultNsave1 = zeros(simulations,length(combOfVarsN(:,1)));
mResultTsave1 = zeros(simulations,length(combOfVarsT(:,1)));

%% iterations
tt = [1 1];
tic
while tt(1) <= simulations
    %% data generation
    %%% generate survival data
    xall                           = gendata_x(samplesize,p,options{4});% genearte covariates
    betaLower                      = gendata_beta(p,d,options{5},[]);% generate parameter
    [tall,lambdaNtrue]             = gendata_t(xall,samplesize,options{6},[eye(d);betaLower]);% generate nontransplant survival time
    wall                           = gendata_w(samplesize,options{7},1);% generate transplant time
    transIndex                     = find(wall<tall);% find the patients whose transplant happens earlier than survival (transplanted)
    [tallt,lambdaTtrue]            = gendata_t([xall,wall],samplesize,options{8},[eye(d);betaLower]);% generate transplanted survival
    tall(transIndex)               = wall(transIndex)+tallt(transIndex);% transplanted group has new survival time, w+t_T
    groupIndex                     = zeros(samplesize,1);
    groupIndex(transIndex)         = 1;
    zall                           = ones(samplesize,1);
    deltaall                       = ones(samplesize,1);
    [cenallT,deltas,~,~]    = gendata_cen(tallt(groupIndex==1),xall(groupIndex==1,:),options{9},str2num(options{10}),sum(groupIndex==1));% generate censoring for all
    zall(groupIndex==1)            = wall(groupIndex==1)+min(tallt(groupIndex==1),cenallT);
    deltaall(groupIndex==1)        = deltas;
    [cenallN,deltas,~,~]    = gendata_cen(tall(groupIndex==0),xall(groupIndex==0,:),options{9},str2num(options{10}),sum(groupIndex==0));% generate censoring for all
    zall(groupIndex==0)            = min(tall(groupIndex==0),cenallN);
    deltaall(groupIndex==0)        = deltas;
    cenall = zeros(samplesize,1);cenall(groupIndex==1) = cenallT;cenall(groupIndex==0) = cenallN;

    %% separate data as transplant and nontransplant groups
    xallN = xall;xallT = xall(groupIndex==1,:);
    zallN = zall;zallN(groupIndex==1) = wall(groupIndex==1);zallT = zall(groupIndex==1)-wall(groupIndex==1);
    deltaallN = deltaall;deltaallN(groupIndex==1) = 0;deltaallT = deltaall(groupIndex==1);
    samplesizeN = samplesize;samplesizeT = sum(groupIndex);

    %% Step: Semiparametric Survival (from Ge, Ma, and Lu 2021).
    % bandwidth selection
    ntemp = samplesize^(-1/(10/3));
    bandez1 = mean(mad(sum(xall(:,1:(end-1)),2)))*ntemp;
    bandew = mean(mad(wall(groupIndex==1)))*ntemp;
    ntemp = samplesize^(-1/12);
    bandbN = mad(zallN)*ntemp;
    bandbT = mad(zallT)*ntemp;
    bands(1,1) = bandbN; % bandwidth in estimating \Z
    bands(1,2) = bandbT; % bandwidth in estimating \Z
    bands(2,1) = bandez1; % bandwidth in estimating \lambda, denominator
    bands(2,2) = bandew; % bandwidth in estimating E(YX), numerator, denominator
    bands(3,1) = bandez1; % bandwidth in estimating \lambda, numerator
    bands(4,1) = bandez1; % bandwidth in estimating E(YX), numerator, denominator
    % calculate all z_i-z_j
    zdiffall = zall-zall';
    zdiffallT = zallT-zallT';
    zdiffallN = zallN-zallN';
    bandadjbN = bands(1,3);
    bandadjbT = bands(1,4);
    kbzT = kernel(zdiffallT/(bandbT*bandadjbT),1,'Epanechnikov')/(bandbT*bandadjbT);
    kbzN = kernel(zdiffallN/(bandbN*bandadjbN),1,'Epanechnikov')/(bandbN*bandadjbN);
    options1 = optimoptions('fsolve','Display','off',...
        'Algorithm','trust-region-dogleg');
    beta0 = (rand((p-d)*d,1)-0.5); % initial value
    [betaEst,~,eflag,output] = fsolve(@(x)scoreBetaMRW(x,xall,wall,zdiffallN,zdiffallT,...
        kbzN,kbzT,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands),beta0,options1);%
    vic = calvicMRW(reshape(betaEst,p-d,d),xall,wall,zdiffallN,zdiffallT,...
        kbzN,kbzT,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands);
    [stdVal,jaccobMat] = jaccobBetaMRW(reshape(betaEst,p-d,d),xall,wall,zdiffallN,zdiffallT,zdiffall,kbzN,kbzT...
        ,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands);
    %% Estimate the mean residual life function
    % estimate the non-transplant group
    lambdaEstN = lambdaEstMR(betaEst,combOfVarsN(:,1),combOfVarsN(:,2:end),xallN,zallN,deltaallN,samplesizeN,bands,p,d,0);
    mResultN = 1:length(combOfVarsN(:,1));
    uniqueGN = unique(combOfVarsN(:,2:end),'rows');
    for i = 1:length(uniqueGN(:,1))
        temp = find(ismember(combOfVarsN(:,2:end), uniqueGN(i,:), 'rows'));
        idx = temp(1:length(unique(combOfVarsN(:,1))));
        upperLimit = min(quantile(zallN,0.99),max(combOfVarsN(:,1)));
        [m,cl] = mEstSumMR(betaEst,combOfVarsN(idx,1),uniqueGN(i,:),xallN,zallN,deltaallN,samplesizeN,bands,p,d,upperLimit,0);
        mResultN(idx) = m;
    end
    % estimate the transplant group
    lambdaEstT = lambdaEstMR(betaEst,combOfVarsT(:,1),combOfVarsT(:,2:end),[xallT wall(groupIndex==1)],zallT,deltaallT,samplesizeT,bands,p,d,1);
    mResultT = 1:length(combOfVarsT(:,1));
    uniqueGT = unique(combOfVarsT(:,2:end),'rows');
    for i = 1:length(uniqueGT(:,1))
        temp = find(ismember(combOfVarsT(:,2:end), uniqueGT(i,:), 'rows'));
        idx = temp(1:length(unique(combOfVarsT(:,1))));
        upperLimit = min(quantile(zallT,0.99),max(combOfVarsT(:,1)));
        [m,cl] = mEstTSumMR(betaEst,combOfVarsT(idx,1),uniqueGT(i,:),[xallT,wall(groupIndex==1)],zallT,deltaallT,samplesizeT,bands,p,d,upperLimit,1);
        mResultT(idx) = m;
    end

    %% Wrap up everything
    betaEstResult(tt(1),:) = betaEst(:);
    xallsave(tt(1),:,:) = xall;
    zallsave(tt(1),:) = zall;
    tallsave(tt(1),:) = tall;
    cenallsave(tt(1),:) = cenall;
    wallsave(tt(1),:) = wall;
    deltaallsave(tt(1),:) = deltaall;
    groupIndexsave(tt(1),:) = groupIndex;
    flags(tt(1),1) = eflag;
    betaStdResult(tt(1),:) = stdVal;
    mResultNsave1(tt(1),:) = mResultN;
    mResultTsave1(tt(1),:) = mResultT;
    vicResult(tt(1)) = vic;

    tt = tt+1;
end
%% output the result
dlmwrite(strcat('bands.txt'),bands,'delimiter','\t');
dlmwrite(strcat('options.txt'),char(options),'delimiter','');
dlmwrite(strcat('trueBeta.txt'),betaLower(:),'delimiter','\t');
save(strcat('xall.mat'),'xallsave');
dlmwrite(strcat('zall.txt'),zallsave,'delimiter','\t');
dlmwrite(strcat('tall.txt'),tallsave,'delimiter','\t');
dlmwrite(strcat('cenall.txt'),cenallsave,'delimiter','\t');
dlmwrite(strcat('deltaall.txt'),deltaallsave,'delimiter','\t');
dlmwrite(strcat('wall.txt'),wallsave,'delimiter','\t');
dlmwrite(strcat('groupIndex.txt'),groupIndexsave,'delimiter','\t');
dlmwrite(strcat('censorInfo.txt'),censorInfo,'delimiter','\t');
dlmwrite(strcat('flag.txt'),flags,'delimiter','\t');
dlmwrite(strcat('resultBeta.txt'),betaEstResult,'delimiter','\t');
dlmwrite(strcat('resultBetaStd.txt'),betaStdResult,'delimiter','\t');
dlmwrite(strcat('estHazardGridNon.txt'),combOfVarsN,'delimiter','\t');
dlmwrite(strcat('estMNon.txt'),mResultNsave1,'delimiter','\t');
dlmwrite(strcat('estHazardGridTrans.txt'),combOfVarsT,'delimiter','\t');
dlmwrite(strcat('estMTrans.txt'),mResultTsave1,'delimiter','\t');
dlmwrite(strcat('vic.txt'),vicResult,'delimiter','\t');

