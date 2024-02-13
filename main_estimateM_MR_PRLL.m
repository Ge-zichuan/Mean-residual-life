paraIndStart = textread('parallel_file.txt','%d');
dlmwrite('parallel_file.txt',paraIndStart+1,'delimiter','\t');
%%% basic information of numerical study
%%% basic information of numerical study
filenum = 424;
simulations = 1050/75;
[filefieldname,filefieldval] = optionRead(filenum);
p = str2num(filefieldval{1});
d = str2num(filefieldval{2});

betaEstResult = dlmread(strcat('resultBeta',num2str(filenum),'.txt'));
betaTrue = dlmread(strcat('trueBeta',num2str(filenum),'.txt'));
bands = dlmread(strcat('bands',num2str(filenum),'.txt'));

load(strcat('xall',num2str(filenum),'.mat'));
zallsave = dlmread(strcat('zall',num2str(filenum),'.txt'));
tallsave = dlmread(strcat('tall',num2str(filenum),'.txt'));
cenallsave = dlmread(strcat('cenall',num2str(filenum),'.txt'));
deltaallsave = dlmread(strcat('deltaall',num2str(filenum),'.txt'));
wallsave = dlmread(strcat('wall',num2str(filenum),'.txt'));
groupIndexsave = dlmread(strcat('groupIndex',num2str(filenum),'.txt'));

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

%% Assigning stores
betaResultPM = zeros(simulations,p-d+1);%dlmread(strcat('estBetaPM',num2str(filenum),'.txt'));%
mResultPMNsave1 = zeros(simulations,length(combOfVarsN(:,1)));
mResultPMTsave1 = zeros(simulations,length(combOfVarsT(:,1)));
mResultNsave1 = zeros(simulations,length(combOfVarsN(:,1)));
mResultTsave1 = zeros(simulations,length(combOfVarsT(:,1)));
lambdaResultTsave1 = zeros(simulations,length(combOfVarsT(:,1)));
lambdaResultNsave1 = zeros(simulations,length(combOfVarsN(:,1)));
tausave = zeros(simulations,2);
betaStdResult = zeros(simulations,d*(p-d));
% %% iterations
% for ibandPos = 14%[9 10 13 14 18 22]
%
bands(1,3:4) = [10,10]; % tuning the bandwidth in estimating \Z
bands(2,3:6) = [20,50,30,80]; % tuning the bandwidth in estimating \lambda, denominator
bands(3,2:7) = [50,80,50,50,100,50]; % tuning the bandwidth in estimating \lambda, numerator
bands(4,2:6) = [1,15,20,30,30]; % tuning the bandwidth in estimating E(YX), numerator, denominator
%
% % ibandPos = 18;
% switch ibandPos
%     case 10
%         tempibands = 0.5:0.5:10;
%     case 14
%         tempibands = 100:50:500;
% end
% tempibands = 30:5:80;
% for iband = tempibands
%     bands(ibandPos) = iband;

tt = [1 1];
tic
while (tt(1) <= simulations)
    transIndex = find(groupIndexsave(tt(1),:)==1);
    xall = squeeze(xallsave(tt(1)+(paraIndStart-1)*simulations,:,:));
    zall = zallsave(tt(1)+(paraIndStart-1)*simulations,:)';
    tall = tallsave(tt(1)+(paraIndStart-1)*simulations,:)';
    cenall = cenallsave(tt(1)+(paraIndStart-1)*simulations,:)';
    deltaall = deltaallsave(tt(1)+(paraIndStart-1)*simulations,:)';
    wall = wallsave(tt(1)+(paraIndStart-1)*simulations,:)';
    groupIndex = groupIndexsave(tt(1)+(paraIndStart-1)*simulations,:)';
    samplesize = length(tall);
    %% separate data as transplant and nontransplant groups
    xallN = xall;xallT = xall(groupIndex==1,:);
    zallN = zall;zallN(groupIndex==1) = wall(groupIndex==1);zallT = zall(groupIndex==1)-wall(groupIndex==1);
    deltaallN = deltaall;deltaallN(groupIndex==1) = 0;deltaallT = deltaall(groupIndex==1);
    samplesizeN = samplesize;samplesizeT = sum(groupIndex);
    betaEst = betaEstResult(tt(1)+(paraIndStart-1)*simulations,:)';%betaLower;%
    
    %     btransx = hist(xall*[eye(d);betaEstResult(1,:)'])
    %% Solve the problem by Chen's methods
%     options1 = optimoptions('fsolve','Display','off','Algorithm','trust-region');
%     beta0 = [1;betaTrue+(rand((p-d),d)-0.5)*0.0;1]*0.0;
%     [betaEstPM,~,eflag,output] = fsolve(@(x)scoreBetaPM1(x,[xallN zeros(samplesizeN,1)],zallN,deltaallN,samplesizeN,[xallT wall(groupIndex==1)],zallT,deltaallT,...
%         samplesizeT,groupIndex==1,p+1),beta0,options1);%
%     betaResultPM(tt(1),:) = betaEstPM(2:(end))/betaEstPM(1);
%     mResultPMN = estMPM1(betaEstPM(2:(end-1))/betaEstPM(1),unique(combOfVarsN(:,1)),unique(combOfVarsN(:,2:end)),xallN,zallN,deltaallN,samplesizeN);
%     combOfVarsTTemp = unique(combOfVarsT(:,2:end),'rows');
%     mResultPMT = estMPM1T(betaEstPM(2:end)/betaEstPM(1),unique(combOfVarsT(:,1)),(combOfVarsTTemp(:,1)+combOfVarsTTemp(:,2)/5-10),[xallT wall(groupIndex==1)],zallT,deltaallT,samplesizeT);
%     mResultPMTsave1(tt(1),:) = mResultPMT;
%     mResultPMNsave1(tt(1),:) = mResultPMN;
    
    %% Estimate the hazard and mean residual function
    
    bands(1,3:4) = [1,10]; % good N.  tuning the bandwidth in estimating \Z
    bands(2,3:6) = [200,200,10,20]; % tuning the bandwidth in estimating \lambda, denominator
    lambdaEstN = lambdaEstMR(betaEst,combOfVarsN(:,1),combOfVarsN(:,2:end),xallN,zallN,deltaallN,samplesizeN,bands,p,d,0);
    mResultN = 1:length(combOfVarsN(:,1));
    cumulLambdaN = mResultN;
    uniqueGN = unique(combOfVarsN(:,2:end),'rows');
    for i = 1:length(uniqueGN(:,1))
        temp = find(ismember(combOfVarsN(:,2:end), uniqueGN(i,:), 'rows'));
        idx = temp(1:length(unique(combOfVarsN(:,1))));
        upperLimit = min(quantile(zallN,0.99),max(combOfVarsN(:,1)));
        tausave(tt(1),1) = upperLimit;
        [m,cl] = mEstSumMR(betaEst,combOfVarsN(idx,1),uniqueGN(i,:),xallN,zallN,deltaallN,samplesizeN,bands,p,d,upperLimit,0);
        cumulLambdaN(idx) = cl;
        mResultN(idx) = m;
    end
    mResultNsave1(tt(1),:) = mResultN;
    lambdaResultNsave1(tt(1),:) = lambdaEstN;
    
    bands(1,3:4) = [10,10]; % good N.  tuning the bandwidth in estimating \Z
    bands(2,3:6) = [10,20,100,100]; % tuning the bandwidth in estimating \lambda, denominator
    lambdaEstT = lambdaEstMR(betaEst,combOfVarsT(:,1),combOfVarsT(:,2:end),[xallT wall(groupIndex==1)],zallT,deltaallT,samplesizeT,bands,p,d,1);
    mResultT = 1:length(combOfVarsT(:,1));
    cumulLambdaT = mResultT;
    uniqueGT = unique(combOfVarsT(:,2:end),'rows');
    for i = 1:length(uniqueGT(:,1))
        temp = find(ismember(combOfVarsT(:,2:end), uniqueGT(i,:), 'rows'));
        idx = temp(1:length(unique(combOfVarsT(:,1))));
        upperLimit = min(quantile(zallT,0.99),max(combOfVarsT(:,1)));
        tausave(tt(1),2) = upperLimit;
        [m,cl] = mEstTSumMR(betaEst,combOfVarsT(idx,1),uniqueGT(i,:),[xallT,wall(groupIndex==1)],zallT,deltaallT,samplesizeT,bands,p,d,upperLimit,1);
        cumulLambdaT(idx) = cl;
        mResultT(idx) = m;
    end
    mResultTsave1(tt(1),:) = mResultT;
    lambdaResultTsave1(tt(1),:) = lambdaEstT;
    
    %% Estimate the standard deviation of \beta
% %     % calculate all z_i-z_j
% %     for jbandPosition = 7:27%[7 8 10 11 12 14] %
% %         if bands(jbandPosition) ~= 0
% %             bandVar = [];
% %             for jband = 1:1:160
% %                 bands(jbandPosition) = jband;
%                 bands(1,3:4) =    [10,10]; % tuning the bandwidth in estimating \Z
%                 bands(2,3:6) =    [50,50,50,50]; %
%                 bands(3,2:7) = [18,50,50,50,50,50]; % (3,2) small;(3,3) two ends;
%                 bands(4,2:6) = [50,50,50,50,50]; %
%                 zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
%                 zdiffallT = zallT*ones(1,samplesizeT)-ones(samplesizeT,1)*zallT';% {i,j} = Z_i-Z_j
%                 zdiffallN = zallN*ones(1,samplesizeN)-ones(samplesizeN,1)*zallN';% {i,j} = Z_i-Z_j
%                 bandbN = bands(1,1);
%                 bandbT = bands(1,2);
%                 bandadjbN = bands(1,3);
%                 bandadjbT = bands(1,4);
%                 kbzT = kernel(zdiffallT/(bandbT*bandadjbT),1,'Epanechnikov')/(bandbT*bandadjbT);% {i,j} = Z_i-Z_j
%                 kbzN = kernel(zdiffallN/(bandbN*bandadjbN),1,'Epanechnikov')/(bandbN*bandadjbN);% {i,j} = Z_i-Z_j
%                 [stdVal,jaccobVal] = jaccobBetaMR(betaEst,xall,wall,zdiffallN,zdiffallT,kbzN,kbzT...
%                     ,deltaallN,deltaallT,deltaall,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands);
%                   bandsVar = bands;
% %                 bandVar = [bandVar;[jband stdVal']];
% %             end
% %             bands(1,3:4) =    [10,10]; % tuning the bandwidth in estimating \Z
% %             bands(2,3:6) =    [10,30,60,70]; %
% %             bands(3,2:7) = [10,100,80,50,40,60]; % (3,2) small;(3,3) two ends;
% %             bands(4,2:6) = [5,5,60,50,60]; %
% %             t = figure('visible','on');
% %             plot(bandVar(:,1),bandVar(:,2:end));
% %             hold on
% %             indy = sub2ind(size(bandVar), round(linspace(1,length(bandVar(:,1)),length(bandVar(1,2:end)))), 2:length(bandVar(1,:)));
% %             text_y = bandVar(indy);
% %             text_x = bandVar(round(linspace(1,length(bandVar(:,1)),length(bandVar(1,2:end)))), 1);
% %             labels = cellstr(strcat('beta', int2str((1:8).')));
% %             tttt = text(text_x, text_y, labels, 'Color', 'k');
% %             saveas(t,[pwd,strcat('\plots\std',num2str(filenum),'_',regexprep(num2str(jbandPosition),'\s+','_'),'.png')],'png')
% %         end
% %     end
%     betaStdResult(tt(1),:) = stdVal;
%     
    %% Wrap up everything
    %     tau(tt(1),:) = [max(inputt),max(inputt)];
    disp(tt);
    tt = tt+1;
end
toc

%% Store the results
% % % dlmwrite(strcat('estLambdaNon',num2str(filenum),'_',num2str(paraIndStart),'.txt'),lambdaResultNsave1,'delimiter','\t');
% % % dlmwrite(strcat('estLambdaTrans',num2str(filenum),'_',num2str(paraIndStart),'.txt'),lambdaResultTsave1,'delimiter','\t');
% % % dlmwrite(strcat('tau',num2str(filenum),'_',num2str(paraIndStart),'.txt'),tausave,'delimiter','\t');
% % dlmwrite(strcat('estBetaPM',num2str(filenum),'_',num2str(paraIndStart),'.txt'),betaResultPM,'delimiter','\t');
% % dlmwrite(strcat('estMPMNon',num2str(filenum),'_',num2str(paraIndStart),'.txt'),mResultPMNsave1,'delimiter','\t');
% % dlmwrite(strcat('estMPMTrans',num2str(filenum),'_',num2str(paraIndStart),'.txt'),mResultPMTsave1,'delimiter','\t');
dlmwrite(strcat('estHazardGridNon',num2str(filenum),'.txt'),combOfVarsN,'delimiter','\t');
dlmwrite(strcat('estMNon',num2str(filenum),'_',num2str(paraIndStart),'.txt'),mResultNsave1,'delimiter','\t');
dlmwrite(strcat('estHazardGridTrans',num2str(filenum),'.txt'),combOfVarsT,'delimiter','\t');
dlmwrite(strcat('estMTrans',num2str(filenum),'_',num2str(paraIndStart),'.txt'),mResultTsave1,'delimiter','\t');
% dlmwrite(strcat('resultBetaStd',num2str(filenum),'_',num2str(paraIndStart),'.txt'),betaStdResult,'delimiter','\t');
% dlmwrite(strcat('bandsVar',num2str(filenum),'.txt'),bandsVar,'delimiter','\t');

