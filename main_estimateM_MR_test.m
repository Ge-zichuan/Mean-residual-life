%%% basic information of numerical study
filenum = 122;
simulations = 1050;
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
subplot(1,3,1);
hist(tallsave(groupIndexsave==1),50);xlabel('T');
subplot(1,3,2);
hist(tallsave(groupIndexsave==0),50);xlabel('N');
subplot(1,3,3);
hist(wallsave(wallsave~=0),50);xlabel('W');
inputt = 0:0.1:2;
inputbbtx = repmat((-2:(2+2)/10:2)',[1 d]);
combOfVarsN = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsN = combvec(combOfVarsN',inputbbtx(:,i)')';
    end
end
inputt = 0:0.05:0.5;
inputw = 0:0.1:2;
inputbbtx = repmat((-2:(2+2)/10:2)',[1 d]);
combOfVarsT = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsT = combvec(combOfVarsT',inputbbtx(:,i)')';
    end
end
combOfVarsT = combvec(combOfVarsT',inputw)';

%% Assigning stores
mResultTsaveMany = zeros(10,simulations,length(combOfVarsT(:,1)));
mResultNsaveMany = zeros(10,simulations,length(combOfVarsN(:,1)));
lambdaResultTsaveMany = zeros(10,simulations,length(combOfVarsT(:,1)));
lambdaResultNsaveMany = zeros(10,simulations,length(combOfVarsN(:,1)));
tausave = zeros(simulations,2);
betaStdResult = zeros(simulations,d*(p-d));
bands(1,3:4) =   [10,10]; % tuning the bandwidth in estimating \Z
bands(2,3:6) =   [10,100,30,80]; % tuning the bandwidth in estimating \lambda, denominator
bands(3,2:7) = [1,80,50,50,100,50]; % tuning the bandwidth in estimating \lambda, numerator
bands(4,2:6) = [1,15,20,30,30]; % tuning the bandwidth in estimating E(YX), numerator, denominator
% %% iterations
for ibandPos = [9 10 13 14 18 22]
    tempibands = bands(ibandPos)/20:(bands(ibandPos)*5-bands(ibandPos)/20)/10:bands(ibandPos)*5;
    for iband = 1:length(tempibands)
        bands(ibandPos) = tempibands(iband);
        
        tt = [1 1];
        tic
        while (tt(1) <= simulations)
            transIndex = find(groupIndexsave(tt(1),:)==1);
            xall = squeeze(xallsave(tt(1),:,:));
            p = length(xall(1,:));
            zall = zallsave(tt(1),:)';
            tall = tallsave(tt(1),:)';
            cenall = cenallsave(tt(1),:)';
            deltaall = deltaallsave(tt(1),:)';
            wall = wallsave(tt(1),:)';
            groupIndex = groupIndexsave(tt(1),:)';
            samplesize = length(tall);
            %% separate data as transplant and nontransplant groups
            xallN = xall;xallT = xall(groupIndex==1,:);
            zallN = zall;zallN(groupIndex==1) = wall(groupIndex==1);zallT = zall(groupIndex==1)-wall(groupIndex==1);
            deltaallN = deltaall;deltaallN(groupIndex==1) = 0;deltaallT = deltaall(groupIndex==1);
            samplesizeN = samplesize;samplesizeT = sum(groupIndex);
            betaEst = betaEstResult(tt(1),:)';%betaLower;%
            
            %         btransx = hist(xall*[eye(d);betaEstResult(1,:)'])
            %% Solve the problem by Chen's methods
%             options1 = optimoptions('fsolve','Display','off');
%             beta0 = [1;betaTrue+(rand((p-d),d)-0.5)*0.2;1];
%             [betaEstPM,~,eflag,output] = fsolve(@(x)scoreBetaPM1(x,[xallN zeros(samplesizeN,1)],zallN,deltaallN,samplesizeN,[xallT wall(groupIndex==1)],zallT,deltaallT,...
%                 samplesizeT,transIndex,p+1),beta0,options1);%
%             mResultPMN = estMPM1(betaEstPM(2:(end-1))/betaEstPM(1),unique(combOfVarsN(:,1)),unique(combOfVarsN(:,2:end)),xallN,zallN,deltaallN,samplesizeN);
%             combOfVarsTTemp = unique(combOfVarsT(:,2:end),'rows');
%             mResultPMT = estMPM1T(betaEstPM(2:end)/betaEstPM(1),unique(combOfVarsT(:,1)),(combOfVarsTTemp(:,1)+combOfVarsTTemp(:,2)/5-10),[xallT wall(groupIndex==1)],zallT,deltaallT,samplesizeT);

            %% Estimate the hazard and mean residual function
            bandsN = bands;
            bandsN(1,3:4) = [5,10]; % tuning the bandwidth in estimating \Z
            bandsN(2,3:6) = [10,5,100,70]; % tuning the bandwidth in estimating \lambda, denominator
            lambdaEstN = lambdaEstMR(betaEst,combOfVarsN(:,1),combOfVarsN(:,2:end),xallN,zallN,deltaallN,samplesizeN,bandsN,d,0);
            mResultN = 1:length(combOfVarsN(:,1));
            cumulLambdaN = mResultN;
            uniqueGN = unique(combOfVarsN(:,2:end),'rows');
            for i = 1:length(uniqueGN(:,1))
                temp = find(ismember(combOfVarsN(:,2:end), uniqueGN(i,:), 'rows'));
                idx = temp(1:length(unique(combOfVarsN(:,1))));
                upperLimit = min(quantile(zallN,0.99),max(combOfVarsN(:,1)));
                tausave(tt(1),1) = upperLimit;
                [m,cl] = mEstSumMR(betaEst,combOfVarsN(idx,1),uniqueGN(i,:),xallN,zallN,deltaallN,samplesizeN,bandsN,d,upperLimit,0);
                cumulLambdaN(idx) = cl;
                mResultN(idx) = m;
            end
            
            bandsT = bands;
            bandsT(1,3:4) = [20,50]; % tuning the bandwidth in estimating \Z
            bandsT(2,3:6) = [10,2,100,50]; % tuning the bandwidth in estimating \lambda, denominator
            lambdaEstT = lambdaEstMR(betaEst,combOfVarsT(:,1),combOfVarsT(:,2:end),[xallT wall(groupIndex==1)],zallT,deltaallT,samplesizeT,bandsT,d,1);
            mResultT = 1:length(combOfVarsT(:,1));
            cumulLambdaT = mResultT;
            uniqueGT = unique(combOfVarsT(:,2:end),'rows');
            for i = 1:length(uniqueGT(:,1))
                temp = find(ismember(combOfVarsT(:,2:end), uniqueGT(i,:), 'rows'));
                idx = temp(1:length(unique(combOfVarsT(:,1))));
                upperLimit = min(quantile(zallT,0.99),max(combOfVarsT(:,1)));
                tausave(tt(1),2) = upperLimit;
                [m,cl] = mEstTSumMR(betaEst,combOfVarsT(idx,1),uniqueGT(i,:),[xallT,wall(groupIndex==1)],zallT,deltaallT,samplesizeT,bandsT,d,upperLimit,1);
                cumulLambdaT(idx) = cl;
                mResultT(idx) = m;
            end
            %% Estimate the standard deviation of \beta
%             bands(1,3:4) = [20,10]; % tuning the bandwidth in estimating \Z
%             bands(2,3:6) = [20,50,30,80]; % tuning the bandwidth in estimating \lambda, denominator
%             bands(3,2:7) = [1,80,50,50,100,50]; % tuning the bandwidth in estimating \lambda, numerator
%             bands(4,2:6) = [1,15,20,30,30]; % tuning the bandwidth in estimating E(YX), numerator, denominator
% %             for jbandPosition = [7:8 10:12 14 17 21 24]
% %                 if bands(jbandPosition) ~= 0
% %                     bandVar = [];
% %                     for jband = 1:1:160
% %                         bands(jbandPosition) = jband;
%                         zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
%                         zdiffallT = zallT*ones(1,samplesizeT)-ones(samplesizeT,1)*zallT';% {i,j} = Z_i-Z_j
%                         zdiffallN = zallN*ones(1,samplesizeN)-ones(samplesizeN,1)*zallN';% {i,j} = Z_i-Z_j
%                         bandbN = bands(1,1);
%                         bandbT = bands(1,2);
%                         kbzN = kernel(zdiffallN/(bandbN*bandadjbN),1,'Epanechnikov')/(bandbN*bandadjbN);% {i,j} = Z_i-Z_j
%                         kbzT = kernel(zdiffallT/(bandbT*bandadjbT),1,'Epanechnikov')/(bandbT*bandadjbT);% {i,j} = Z_i-Z_j
%                         [stdVal,jaccobVal] = jaccobBetaMR(betaEst,xall,wall,zdiffallN,zdiffallT,kbzN,kbzT...
%                             ,deltaallN,deltaallT,deltaall,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands);
% %                         bandVar = [bandVar;[jband stdVal']];
% %                     end
% %                     bands(1,3:4) = [20,10]; % tuning the bandwidth in estimating \Z
% %                     bands(2,3:6) = [20,50,30,80]; % tuning the bandwidth in estimating \lambda, denominator
% %                     bands(3,2:7) = [1,80,50,50,100,50]; % tuning the bandwidth in estimating \lambda, numerator
% %                     bands(4,2:6) = [1,15,20,30,30]; % tuning the bandwidth in estimating E(YX), numerator, denominator
% %                     t = figure('visible','off','WindowState', 'maximized','Position', get(0, 'Screensize'));
% %                     plot(bandVar(:,1),bandVar(:,2:end));
% %                     saveas(t,[pwd,strcat('\plots\std',num2str(filenum),'_',regexprep(num2str(jbandPosition),'\s+','_'),'.png')],'png')
% %                 end
% %             end
%             betaStdResult(tt(1),:) = stdVal;
            
            %% Wrap up everything
            lambdaResultNsaveMany(iband,tt(1),:) = lambdaEstN;
            lambdaResultTsaveMany(iband,tt(1),:) = lambdaEstT;
            mResultNsaveMany(iband,tt(1),:) = mResultN;
            mResultTsaveMany(iband,tt(1),:) = mResultT;
            tau(tt(1),:) = [max(inputt),max(inputt)];
            if (mod(tt(1),500)==0 || tt(1)<5)
                disp(tt);
            end
            tt = tt+1;
        end
    end
    disp(ibandPos);
    [trueHazardN,trueMRLN,~] = trueMRL(combOfVarsN(:,1),combOfVarsN(:,2:end),filefieldval{6});
    [trueHazardT,trueMRLT,~] = trueMRL(combOfVarsT(:,1),combOfVarsT(:,2:end),filefieldval{8});
    %     plotTestTrend = plotResult(combOfVarsN,[],lambdaResultNsaveTest,[],trueHazardN,[],...
    %         [filenum,ibandPos,tempibands],'on','biasChange','mrlTest',bands);
    plotT = plotResult(combOfVarsN,[],squeeze(mResultNsaveMany(iband,:,:)),[],trueMRLN,[],[filenum,ibandPos,iband],'on','N','mrl1',bandsN);
    plotT = plotResult(combOfVarsT,[],squeeze(mResultTsaveMany(iband,:,:)),[],trueMRLT,[],[filenum,ibandPos,iband],'on','T','mrl1',bandsT);
%     plotTestTrend = plotResult(combOfVarsN,[],lambdaResultNsaveMany,[],trueHazardN,[],...
%         [filenum,ibandPos,combOfVarsN(2,1),tempibands],'on','biasChange','mrlTestCloser',bandsN);
%     plotTestTrend = plotResult(combOfVarsT,[],lambdaResultTsaveMany,[],trueHazardT,[],...
%         [filenum,ibandPos,combOfVarsT(2,1),tempibands],'on','biasChange','mrlTestCloser',bandsT);
    bands(1,3:4) =    [10,10]; % tuning the bandwidth in estimating \Z
    bands(2,3:6) =    [10,100,30,80]; % tuning the bandwidth in estimating \lambda, denominator
    bands(3,2:7) = [50,80,50,50,100,50]; % tuning the bandwidth in estimating \lambda, numerator
    bands(4,2:6) = [1,15,20,30,30]; % tuning the bandwidth in estimating E(YX), numerator, denominator
    %% plot hazard and m
    % hist(xall*[eye(d);betaEst]); % paras1 include upper d-by-d block.
    %     plotT = plotResult(combOfVarsN,[],mResultNsave1,[],trueMRLN,[],[filenum,ibandPos,iband],'off','N','mrl1',bands);
    %     plotT = plotResult(combOfVarsT,[],mResultTsave1,[],trueMRLT,[],[filenum,ibandPos,iband],'off','T','mrl1',bands);
    %     plotImprove = plotResult(combOfVarsT,combOfVarsN,mResultTsave1,mResultNsave1,trueMRLT,trueMRLN,...
    %         [filenum,0,0],'on','improve','mrlimprove',bands);
end
toc
%% Store the results
% dlmwrite(strcat('estHazardGridNon',num2str(filenum),'.txt'),combOfVarsN,'delimiter','\t');
% dlmwrite(strcat('estHazardGridTrans',num2str(filenum),'.txt'),combOfVarsT,'delimiter','\t');
% dlmwrite(strcat('estMNon',num2str(filenum),'.txt'),mResultNsave1,'delimiter','\t');
% dlmwrite(strcat('estMTrans',num2str(filenum),'.txt'),mResultTsave1,'delimiter','\t');
% dlmwrite(strcat('estLambdaNon',num2str(filenum),'.txt'),lambdaResultNsave1,'delimiter','\t');
% dlmwrite(strcat('estLambdaTrans',num2str(filenum),'.txt'),lambdaResultTsave1,'delimiter','\t');
% dlmwrite(strcat('tau',num2str(filenum),'.txt'),tausave,'delimiter','\t');

