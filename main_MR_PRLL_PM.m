rng('shuffle');
%%% basic information of numerical study
filenum = 220;
simulations = 1050/75;
samplesize = 1000;
p = 9;
d = 1;
options = {num2str(p);num2str(d);num2str(samplesize);'norm';'ma91';'lnl2';'unif';'xia3W';'gamm';'0.0'};
inputt = 0:2:50;
inputbbtx = repmat((-2:(2+2)/10:2)',[1 d]);
combOfVarsN = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsN = combvec(combOfVarsN',inputbbtx(:,i)')';
    end
end
inputt = 0:5:100;
inputw = 0:10:200;
inputbbtx = repmat((-2:(2+2)/10:2)',[1 d]);
combOfVarsT = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsT = combvec(combOfVarsT',inputbbtx(:,i)')';
    end
end
combOfVarsT = combvec(combOfVarsT',inputw)';

%% Assigning stores
jump = 0;
    betaTrue = dlmread(strcat('trueBeta',num2str(filenum),'.txt'));
if exist(strcat('xall',num2str(filenum),'_75','.mat'))
    paraIndStart = textread('parallel_file4.txt','%d')
    dlmwrite('parallel_file4.txt',paraIndStart+1,'delimiter','\t');
    load(strcat('xall',num2str(filenum),'_',num2str(paraIndStart),'.mat'));
    xallsave = reshape(xallsave,simulations,samplesize,p);
    zallsave = dlmread(strcat('zall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    tallsave = dlmread(strcat('tall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    cenallsave = dlmread(strcat('cenall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    deltaallsave = dlmread(strcat('deltaall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    wallsave = dlmread(strcat('wall',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    groupIndexsave = dlmread(strcat('groupIndex',num2str(filenum),'_',num2str(paraIndStart),'.txt'));
    if exist(strcat('estBetaPM',num2str(filenum),'_',num2str(paraIndStart),'.txt'))
        jump = 1;
    end
end
betaResultPM = zeros(simulations,p-d+1);%dlmread(strcat('estBetaPM',num2str(filenum),'.txt'));%
mResultPMNsave1 = zeros(simulations,length(combOfVarsN(:,1)));
mResultPMTsave1 = zeros(simulations,length(combOfVarsT(:,1)));

if jump == 1
    display('existed');
else
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
        else
            xall                           = gendata_x(samplesize,p,options{4});% genearte covariates
            % [xall,c,s]                     = gendata_normalize(xall,'mimic');
            betaLower                      = gendata_beta(p,d,options{5},[]);% generate parameter
            [tall,lambdaNtrue]             = gendata_t(xall,samplesize,options{6},[eye(d);betaLower]);% generate nontransplant survival time
            wall                           = gendata_w(samplesize,options{7},200);% generate transplant time
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
        
        %% Solve the problem by Chen's methods
        options1 = optimoptions('fsolve','Display','off','Algorithm','trust-region');
        beta0 = [1;betaTrue+(rand((p-d),d)-0.5)*0.5;1];
        [betaEstPM,~,eflag,output] = fsolve(@(x)scoreBetaPM1(x,[xallN zeros(samplesizeN,1)],zallN,deltaallN,samplesizeN,[xallT wall(groupIndex==1)],zallT,deltaallT,...
            samplesizeT,groupIndex==1,p+1),beta0,options1);%
        betaResultPM(tt(1),:) = betaEstPM(2:(end))/betaEstPM(1);
        mResultPMN = estMPM1(betaEstPM(2:(end-1))/betaEstPM(1),unique(combOfVarsN(:,1)),unique(combOfVarsN(:,2:end)),xallN,zallN,deltaallN,samplesizeN);
        combOfVarsTTemp = unique(combOfVarsT(:,2:end),'rows');
        mResultPMT = estMPM1T(betaEstPM(2:end)/betaEstPM(1),unique(combOfVarsT(:,1)),(combOfVarsTTemp(:,1)+combOfVarsTTemp(:,2)/5-10),[xallT wall(groupIndex==1)],zallT,deltaallT,samplesizeT);
        mResultPMTsave1(tt(1),:) = mResultPMT;
        mResultPMNsave1(tt(1),:) = mResultPMN;
        tt = tt+1
    end
    toc
    %% Store the results
    paraIndStart = textread('parallel_file4.txt','%d')
    dlmwrite('parallel_file4.txt',paraIndStart+1,'delimiter','\t');
    dlmwrite(strcat('estBetaPM',num2str(filenum),'_',num2str(paraIndStart),'.txt'),betaResultPM,'delimiter','\t');
    dlmwrite(strcat('estMPMNon',num2str(filenum),'_',num2str(paraIndStart),'.txt'),mResultPMNsave1,'delimiter','\t');
    dlmwrite(strcat('estMPMTrans',num2str(filenum),'_',num2str(paraIndStart),'.txt'),mResultPMTsave1,'delimiter','\t');
end

