rng('shuffle');
%%% basic information of numerical study
filenum = 422;
simulations = 1050/75;
samplesize = 2000;
p = 6;
d = 2;
% options = {num2str(p);num2str(d);num2str(samplesize);'norm';'ma91';'lnl2';'unif';'xia3WW';'gamm';'0.0'};
options = {num2str(p);num2str(d);num2str(samplesize);'norm';'ma62';'LinLi21';'unif';'LinLi22W';'gamm';'0.0'};

%% Assigning stores
jump = 0;
if exist(strcat('xall',num2str(filenum),'_75','.mat'))
    paraIndStart = textread('parallel_file.txt','%d')
    dlmwrite('parallel_file.txt',paraIndStart+1,'delimiter','\t');
    load(strcat('xall',num2str(filenum),'_',num2str(paraIndStart),'.mat'));
    xallsave = reshape(xallsave,simulations,samplesize,p);
    zallsave = dlmread(strcat('zall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    tallsave = dlmread(strcat('tall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    cenallsave = dlmread(strcat('cenall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    deltaallsave = dlmread(strcat('deltaall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    wallsave = dlmread(strcat('wall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    groupIndexsave = dlmread(strcat('groupIndex',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    if exist(strcat('resultBeta',num2str(filenum),'_',num2str(paraIndStart),'.txt'))
        jump = 1;
    end
else
    xallsave = zeros(simulations,samplesize,p);
    zallsave = zeros(simulations,samplesize);
    tallsave = zeros(simulations,samplesize);
    cenallsave = zeros(simulations,samplesize);
    wallsave = zeros(simulations,samplesize);
    deltaallsave = zeros(simulations,samplesize);
    groupIndexsave = zeros(simulations,samplesize);
    censorInfo = zeros(simulations,1);
end

if jump == 1
    display('existed');
else
    bands = zeros(4,7);
    betaEstResult = zeros(simulations,d*(p-d));
    betaStdResult = zeros(simulations,d*(p-d));
    rhoResult = zeros(simulations,1);
    flags = zeros(simulations,2);
    
                bands(1,3:4) =    [10,10]; % tuning the bandwidth in estimating \Z
                bands(2,3:6) =    [50,50,50,20]; % tuning the bandwidth in estimating \lambda, denominator
                bands(3,2:7) = [50,50,50,50,50,20]; % tuning the bandwidth in estimating \lambda, numerator
                bands(4,2:6) = [8,5,5,5,10]; % tuning the bandwidth in estimating E(YX), numerator, denominator
    % bands(1,:) 5 seems to be the good choice. positive related to est.
    % std. The larger, the better? stable at 10
    % bands(2,:) larger -> larger biases. Seems 8 (or 7) is the best. Skew
    % or not doesn't change too much. Change all variables, positively, a
    % little. larger -> smaller est std.
    % bands(3,:) around 15 is the best. larger or smaller leads to large
    % biase and strange shape, std estimation. Large -> large std. Seems 10
    % is better now. It changes estimated std a lot. The larger, the
    % smaller emp std?
    % bands(4,:) smaller, biases seem to get smaller. 2 seems to be a good
    % choice, larger value will increase the sample sd a little bit.
    % Negative to est. std. Positive to empirical std.
    
    %% iterations
    tt = [1 1];
    tic
    while tt(1) <= simulations
        %% data generation
        if exist(strcat('xall',num2str(filenum),'_75','.mat'))
            xall = squeeze(xallsave(tt(1),:,:));
            tall = squeeze(tallsave(tt(1),:))';
            wall = squeeze(wallsave(tt(1),:))';
            zall = squeeze(zallsave(tt(1),:))';
            cenall = squeeze(cenallsave(tt(1),:))';
            deltaall = squeeze(deltaallsave(tt(1),:))';
            groupIndex = squeeze(groupIndexsave(tt(1),:))';
            betaLower  = gendata_beta(p,d,options{5},[]);% generate parameter
        else
            %%% generate survival data
            xall                           = gendata_x(samplesize,p,options{4});% genearte covariates
%             [xall,c,s]                     = gendata_normalize(xall,'mimic');
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
            [cenallT,deltas,kval,randx]    = gendata_cen(tallt(groupIndex==1),xall(groupIndex==1,:),options{9},str2num(options{10}),sum(groupIndex==1));% generate censoring for all
            zall(groupIndex==1)            = wall(groupIndex==1)+min(tallt(groupIndex==1),cenallT);
            deltaall(groupIndex==1)        = deltas;
            [cenallN,deltas,kval,randx]    = gendata_cen(tall(groupIndex==0),xall(groupIndex==0,:),options{9},str2num(options{10}),sum(groupIndex==0));% generate censoring for all
            zall(groupIndex==0)            = min(tall(groupIndex==0),cenallN);
            deltaall(groupIndex==0)        = deltas;
            cenall = zeros(samplesize,1);cenall(groupIndex==1) = cenallT;cenall(groupIndex==0) = cenallN;
        end
        
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
        bands(1,1) = bandbN; % tuning the bandwidth in estimating \Z
        bands(1,2) = bandbT; % tuning the bandwidth in estimating \Z
        bands(2,1) = bandez1; % tuning the bandwidth in estimating \lambda, denominator
        bands(2,2) = bandew; % tuning the bandwidth in estimating E(YX), numerator, denominator
        bands(3,1) = bandez1; % tuning the bandwidth in estimating \lambda, numerator
        bands(4,1) = bandez1; % tuning the bandwidth in estimating E(YX), numerator, denominator
        % calculate all z_i-z_j
        if tt(1)==1
            display(bands);
        end
        zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
        zdiffallT = zallT*ones(1,samplesizeT)-ones(samplesizeT,1)*zallT';% {i,j} = Z_i-Z_j
        zdiffallN = zallN*ones(1,samplesizeN)-ones(samplesizeN,1)*zallN';% {i,j} = Z_i-Z_j
        bandadjbN = bands(1,3);
        bandadjbT = bands(1,4);
        kbzT = kernel(zdiffallT/(bandbT*bandadjbT),1,'Epanechnikov')/(bandbT*bandadjbT);% {i,j} = Z_i-Z_j
        % a1 = ans(:,:,1);hist(a1(:));a2 = ans(:,:,2);s = waitforbuttonpress;hist(a2(:));
        kbzN = kernel(zdiffallN/(bandbN*bandadjbN),1,'Epanechnikov')/(bandbN*bandadjbN);% {i,j} = Z_i-Z_j
        try
            options1 = optimoptions('fsolve','Display','off',...
                'Algorithm','trust-region-dogleg');%,'TolFun',1e-5,'MaxIter', 100,'MaxFunEvals',1000
            beta0 = betaLower(:)+(rand((p-d)*d,1)-0.5)*0.2;
            eflag = 98;
            [betaEst,~,eflag,output] = fsolve(@(x)scoreBetaMRW(x,xall,wall,zdiffallN,zdiffallT,...
                kbzN,kbzT,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands),beta0,options1);%
            algorithmtest = 1;
        catch
            algorithmtest = 0;
            eflag = 99;
        end
        
        
        %% Wrap up everything
        if (algorithmtest==0 || eflag<1)
            tt(2) = tt(2)+1;
        else
            betaEstResult(tt(1),:) = betaEst(:);
            xallsave(tt(1),:,:) = xall;
            zallsave(tt(1),:) = zall;
            tallsave(tt(1),:) = tall;
            cenallsave(tt(1),:) = cenall;
            wallsave(tt(1),:) = wall;
            deltaallsave(tt(1),:) = deltaall;
            groupIndexsave(tt(1),:) = groupIndex;
            
            %             [stdVal,jaccobMat] = jaccobBetaMR(betaEst,xall,zdiffallN,zdiffallT,zdiffall,...
            %                 kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands);
            %             betaStdResult(tt(1),:) = stdVal;
            flags(tt(1),1) = eflag;
            tt = tt+1;
            % summarize bandwidth adjustment
            bandwidthCheck = [bands(1,1)*bands(1,2) bands(2,1)*[bands(2,2:3) bands(3,2:3) bands(4,2:3)]];
        end
        if (mod(tt(1),1)==0 || tt(1)<5)
            disp([tt algorithmtest eflag mean(deltaall)]);
        end
    end
    toc
    %% Store the results
    dlmwrite(strcat('bands',num2str(filenum),'.txt'),bands,'delimiter','\t');
    dlmwrite(strcat('options',num2str(filenum),'.txt'),char(options),'delimiter','');
    dlmwrite(strcat('trueBeta',num2str(filenum),'.txt'),betaLower(:),'delimiter','\t');
    
    if exist(strcat('xall',num2str(filenum),'_75','.mat'))
    else
        paraIndStart = textread('parallel_file.txt','%d');
        dlmwrite('parallel_file.txt',paraIndStart+1,'delimiter','\t');
        save(strcat('xall',num2str(filenum),'_',num2str(paraIndStart),'.mat'),'xallsave');
        dlmwrite(strcat('zall',num2str(filenum),'_',num2str(paraIndStart),'.txt'),zallsave,'delimiter','\t');
        dlmwrite(strcat('tall',num2str(filenum),'_',num2str(paraIndStart),'.txt'),tallsave,'delimiter','\t');
        dlmwrite(strcat('cenall',num2str(filenum),'_',num2str(paraIndStart),'.txt'),cenallsave,'delimiter','\t');
        dlmwrite(strcat('deltaall',num2str(filenum),'_',num2str(paraIndStart),'.txt'),deltaallsave,'delimiter','\t');
        dlmwrite(strcat('wall',num2str(filenum),'_',num2str(paraIndStart),'.txt'),wallsave,'delimiter','\t');
        dlmwrite(strcat('groupIndex',num2str(filenum),'_',num2str(paraIndStart),'.txt'),groupIndexsave,'delimiter','\t');
        dlmwrite(strcat('censorInfo',num2str(filenum),'_',num2str(paraIndStart),'.txt'),censorInfo,'delimiter','\t');
    end
    
    dlmwrite(strcat('flag',num2str(filenum),'_',num2str(paraIndStart),'.txt'),flags,'delimiter','\t');
    dlmwrite(strcat('resultBeta',num2str(filenum),'_',num2str(paraIndStart),'.txt'),betaEstResult,'delimiter','\t');
    dlmwrite(strcat('resultBetaStd',num2str(filenum),'_',num2str(paraIndStart),'.txt'),betaStdResult,'delimiter','\t');
end

