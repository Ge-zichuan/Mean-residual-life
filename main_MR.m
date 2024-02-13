%%% basic information of numerical study
filenum = 2029;
simulations = 1050;
samplesize = 1000;
p = 9;
d = 1;
options = {num2str(p);num2str(d);num2str(samplesize);'norm';'ma91';'lnl2';'unif';'xia3';'gamm';'0.2'};

%% created the estimated grid of hazard function
inputt = 0:0.2:4;
inputbbtx = repmat((-3:(3+3)/10:3)',d);
combOfVarsN = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsN = combvec(combOfVarsN',inputbbtx(:,i)')';
    end
end
inputt = 0:0.2:4;
inputbbtx = repmat((-3:(3+3)/10:3)',d);
combOfVarsT = combvec(inputt,inputbbtx(:,1)')';
if (d>1)
    for i = 2:(d)
        combOfVarsT = combvec(combOfVarsT',inputbbtx(:,i)')';
    end 
end

%% Assigning stores
betaEstResult = zeros(simulations,d*(p-d));
betaStdResult = zeros(simulations,d*(p-d));
rhoResult = zeros(simulations,1);
flags = zeros(simulations,2);
lambdaResultT = zeros(simulations,length(combOfVarsT(:,1)));
lambdaResultN = zeros(simulations,length(combOfVarsN(:,1)));
lambdaD1ResultT = zeros(simulations,length(combOfVarsT(:,1)),d);
lambdaD1ResultN = zeros(simulations,length(combOfVarsN(:,1)),d);
xallsave = zeros(simulations,samplesize,p);
zallsave = zeros(simulations,samplesize);
tallsave = zeros(simulations,samplesize);
cenallsave = zeros(simulations,samplesize);
wallsave = zeros(simulations,samplesize);
deltaallsave = zeros(simulations,samplesize);
groupIndexsave = zeros(simulations,samplesize);
censorInfo = zeros(simulations,1);

estLambdaAll = zeros(samplesize,simulations);
bands(1,:) = [1,1,0.05,0.5]; % tuning the bandwidth in estimating \Z
bands(2,:) = [1,5,5,0]; % tuning the bandwidth in estimating \lambda, denominator
bands(3,:) = [1,10,10,0]; % tuning the bandwidth in estimating \lambda, numerator
bands(4,:) = [1,1,1,0]; % tuning the bandwidth in estimating E(YX), numerator, denominator
display(bands);
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
while (tt(1) <= simulations)
    %% data generation
    %%% generate survival data
    xall                         = gendata_x(samplesize,p,options{4});% genearte covariates
    betaLower                    = gendata_beta(p,d,options{5});% generate parameter
    [tall,lambdaNtrue,l]         = gendata_t(xall,samplesize,options{6},[eye(d);betaLower]);% generate nontransplant survival time
    wall                         = gendata_w(samplesize,options{7},quantile(tall,0.9));% generate transplant time
    transIndex                   = find(wall<tall);% find the patients whose transplant happens earlier than survival (transplanted)
    [tallt,lambdaTtrue,~]        = gendata_t(xall,samplesize,options{8},[eye(d);betaLower]);% generate transplanted survival
    tall(transIndex)             = wall(transIndex)+tallt(transIndex);% transplanted group has new survival time, w+t_T
    [cenall,deltaall,kval,randx] = gendata_cen(tall,xall,options{9},str2num(options{10}),samplesize);% generate censoring for all
    zall                         = min(tall,cenall);
    transIndex                   = setdiff(transIndex,find(zall<=wall));% exclude the case of which the censoring 
                                                                                %happened before transplantation, i.e. it should be considered as the nontransplant group
    nontransIndex                = setdiff(1:samplesize,transIndex)';
    wall(nontransIndex)          = 0;
    groupIndex                   = zeros(samplesize,1);
    groupIndex(transIndex)       = 1;sum(groupIndex);
    
    %% separate data as transplant and nontransplant groups
    xallN = xall;xallT = xall(transIndex,:);
    zallN = zall;zallN(transIndex) = wall(transIndex);zallT = zall(transIndex)-wall(transIndex);
    deltaallN = deltaall;deltaallN(transIndex) = 0;deltaallT = deltaall(transIndex);
    samplesizeN = samplesize;samplesizeT = length(transIndex);

    %% Step: Semiparametric Survival (from Ge, Ma, and Lu 2021).
    % bandwidth selection
    ntemp = samplesize^(-1/(10/3));
    bandez1 = mean(mad(xall))*ntemp;
    ntemp = samplesize^(-1/12);
    bandbN = mad(zallN)*ntemp;
    bandbT = mad(zallT)*ntemp;
    bands(1,1) = bandbN; % tuning the bandwidth in estimating \Z
    bands(1,2) = bandbT; % tuning the bandwidth in estimating \Z
    bands(2,1) = bandez1; % tuning the bandwidth in estimating \lambda, denominator
    bands(3,1) = bandez1; % tuning the bandwidth in estimating \lambda, numerator
    bands(4,1) = bandez1; % tuning the bandwidth in estimating E(YX), numerator, denominator
    % calculate all z_i-z_j
    bandadjbN = bands(1,3);
    bandadjbT = bands(1,4);
    zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
    zdiffallT = zallT*ones(1,samplesizeT)-ones(samplesizeT,1)*zallT';% {i,j} = Z_i-Z_j
    zdiffallN = zallN*ones(1,samplesizeN)-ones(samplesizeN,1)*zallN';% {i,j} = Z_i-Z_j
    kbzT = kernel(zdiffallT/bandbT/bandadjbT,1,'norm')/bandbT/bandadjbT;% {i,j} = Z_i-Z_j
    kbzN = kernel(zdiffallN/bandbN/bandadjbN,1,'norm')/bandbN/bandadjbN;% {i,j} = Z_i-Z_j
    kbzTD = kernel(zdiffallT/bandbT/bandadjbT,1,'norm')/bandbT/bandadjbT;% {i,j} = Z_i-Z_j
    kbzND = kernel(zdiffallN/bandbN/bandadjbN,1,'norm')/bandbN/bandadjbN;% {i,j} = Z_i-Z_j

    try
%         options1 = optimoptions('fsolve','Display','off',...
%             'Algorithm','trust-region-dogleg','MaxIter',100,'MaxFunEvals',10000,...
%             'FunctionTolerance',1e-5,'OptimalityTolerance',1e-5);
        options1 = optimoptions('fsolve','Display','off',...
            'Algorithm','trust-region-dogleg','MaxIter', 1000,'MaxFunEvals',10000,...
            'TolFun',1e-5);
        beta0 = betaLower+(rand((p-d),d)-0.5)*0.02;
        betaEst = beta0;
        eflag = 98;
        [betaEst,~,eflag,output] = fsolve(@(x)scoreBetaMR(x,xall,zdiffallN,zdiffallT,zdiffall,...
            kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands),beta0,options1);%
%         [betaEst,~,eflag,output] = fsolve(@(x)scoreBetaMRH(x,xall,zall,zdiffallN,zdiffallT,zdiffall,...
%             kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands),beta0,options1);%
%         [betaEst,~,eflag,~] = fsolve(@(x)scoreBeta(x,estF,zdiffall,kbz,deltaall,samplesize,p,d,bands),beta0,options1);%
        algorithmtest = 1;
    catch
        algorithmtest = 0;
        eflag = 99;
    end
    %% Estimate the hazard and mean residual function
    lambdaEstN = lambdaEstMR(betaEst,combOfVarsN(:,1),combOfVarsN(:,2:end),xallN,zallN,deltaallN,samplesizeN,bands,d,0);
    lambdaEstT = lambdaEstMR(betaEst,combOfVarsT(:,1),combOfVarsT(:,2:end),xallT,zallT,deltaallT,samplesizeT,bands,d,1);
    lambdaD1EstN = lambdaD1EstMR(betaEst,combOfVarsN(:,1),combOfVarsN(:,2:end),xallT,zallT,deltaallT,samplesizeT,bands,d,0);
    lambdaD1EstT = lambdaD1EstMR(betaEst,combOfVarsT(:,1),combOfVarsT(:,2:end),xallN,zallN,deltaallN,samplesizeN,bands,d,1);
    
    %% Wrap up everything
    if (algorithmtest==0)
        tt(2) = tt(2)+1;
    else
        betaEstResult(tt(1),:) = betaEst(:);
        lambdaResultT(tt(1),:) = lambdaEstT;
        lambdaResultN(tt(1),:) = lambdaEstN;
        lambdaD1ResultT(tt(1),:,:) = lambdaD1EstT;
        lambdaD1ResultN(tt(1),:,:) = lambdaD1EstN;
        xallsave(tt(1),:,:) = xall;
        zallsave(tt(1),:) = zall;
        tallsave(tt(1),:) = tall;
        cenallsave(tt(1),:) = cenall;
        wallsave(tt(1),:) = wall;
        deltaallsave(tt(1),:) = deltaall;
        groupIndexsave(tt(1),:) = groupIndex;
        censorInfo(tt(1),1) = kval;

        [stdVal,jaccobMat] = jaccobBetaMR(betaEst,xall,zdiffallN,zdiffallT,zdiffall,...
                kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands);
        betaStdResult(tt(1),:) = stdVal;
        flags(tt(1),1) = eflag;
        tt = tt+1;
        % summarize bandwidth adjustment
        bandwidthCheck = [bands(1,1)*bands(1,3) bands(1,2)*bands(1,4) bands(2,1)*[bands(2,2:3) bands(3,2:3) bands(4,2:3)]];
    end
    if (mod(tt(1),50)==0 || tt(1)<50)
        disp([tt algorithmtest eflag mean(deltaall)]);
    end
end
toc;
    display(bands);
%% plot hazard
% hist(xall*[eye(d);betaEst]);hist(zallN);hist(zallT) % paras1 include upper d-by-d block.

lambdaTrueN = lambdaTrueMR(combOfVarsN(:,1),combOfVarsN(:,2:end),options{6});
lambdaTrueT = lambdaTrueMR(combOfVarsT(:,1),combOfVarsT(:,2:end),options{8});
bbx = zeros(1,d);
% plotsubindex = combOfVars(:,2)==bbx(1);% & combOfVars(:,3)==bbx(2);
% if d>1
%     for i = 1:(d-1)
%         plotsubindex = plotsubindex & combOfVars(:,2+i)==bbx(1+i);
%     end
% end
plotsubindex = find(ismember(combOfVarsN(:,2:end), bbx,'rows'));
t = figure;
subplot(1,2,1)
scatter(combOfVarsN(plotsubindex,1),lambdaTrueN(plotsubindex),'filled');hold on;
scatter(combOfVarsN(plotsubindex,1),nanmedian(lambdaResultN(:,plotsubindex),1),'filled');
scatter(combOfVarsN(plotsubindex,1),quantile(lambdaResultN(:,plotsubindex),0.975,1),'+');
scatter(combOfVarsN(plotsubindex,1),quantile(lambdaResultN(:,plotsubindex),0.025,1),'o');
title(num2str(bbx));hold off;
plotsubindex = find(ismember(combOfVarsT(:,2:end), bbx,'rows'));
subplot(1,2,2)
scatter(combOfVarsT(plotsubindex,1),lambdaTrueT(plotsubindex),'filled');hold on;
scatter(combOfVarsT(plotsubindex,1),nanmedian(lambdaResultT(:,plotsubindex),1),'filled');
scatter(combOfVarsT(plotsubindex,1),quantile(lambdaResultT(:,plotsubindex),0.975,1),'+');
scatter(combOfVarsT(plotsubindex,1),quantile(lambdaResultT(:,plotsubindex),0.025,1),'o');
title(num2str([samplesize,bands(1,3:4),bands(2,2:3)]));hold off;
saveas(t,[pwd,strcat('/plots/hazard',num2str(filenum),regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bbx]),'\s+','_'),'.jpg')])
% % plot3(combOfVars(combOfVars(:,3)==0 & combOfVars(:,2)==0,1),combOfVars(combOfVars(:,3)==0,2),lambdaTrueN(combOfVars(:,3)==0));
% x = unique(combOfVars(:,1));
% y = unique(combOfVars(:,2));
% [X,Y] = meshgrid(x,y);
% trueLambdamesh = X.*exp(Y);
% t = figure;
% subplot(1,2,1);
% contour(X,Y,trueLambdamesh,'ShowText','on');
% subplot(1,2,2);
% contour(X,Y,reshape(median(lambdaResultN,1),length(x),length(y))','ShowText','on');hold off;
%% plot hazard derivatives 
lambdaD1TrueN = lambdaD1TrueMR(combOfVarsN(:,1),combOfVarsN(:,2:end),d,options{6});
lambdaD1TrueT = lambdaD1TrueMR(combOfVarsT(:,1),combOfVarsT(:,2:end),d,options{8});
bbx = zeros(1,d);
plotsubindex = find(ismember(combOfVarsN(:,2:end), bbx,'rows'));
t = figure;
for i = 1:d
    subplot(2,d,i)
    scatter(combOfVarsN(plotsubindex,1),lambdaD1TrueN(plotsubindex,i),'filled');hold on;
    scatter(combOfVarsN(plotsubindex,1),median(lambdaD1ResultN(:,plotsubindex,i),1),'*');
    scatter(combOfVarsN(plotsubindex,1),quantile(lambdaD1ResultN(:,plotsubindex,i),0.975,1),'+');
    scatter(combOfVarsN(plotsubindex,1),quantile(lambdaD1ResultN(:,plotsubindex,i),0.025,1),'o');
    title(num2str(bbx));hold off;
end
plotsubindex = find(ismember(combOfVarsT(:,2:end), bbx,'rows'));
for i = 1:d
    subplot(2,d,i+d)
    scatter(combOfVarsT(plotsubindex,1),lambdaD1TrueT(plotsubindex),'filled');hold on;
    scatter(combOfVarsT(plotsubindex,1),median(lambdaD1ResultT(:,plotsubindex,i),1),'*');
    scatter(combOfVarsT(plotsubindex,1),quantile(lambdaD1ResultT(:,plotsubindex,i),0.975,1),'+');
    scatter(combOfVarsT(plotsubindex,1),quantile(lambdaD1ResultT(:,plotsubindex,i),0.025,1),'o');
    title(num2str(bbx));
end
title(num2str([samplesize,bands(1,3:4),bands(2,2:3)]));hold off;

saveas(t,[pwd,strcat('/plots/hazardD1',num2str(filenum),regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bbx]),'\s+','_'),'.jpg')])
% plot3(combOfVars(combOfVars(:,3)==0 & combOfVars(:,2)==0,1),combOfVars(combOfVars(:,3)==0,2),lambdaTrueN(combOfVars(:,3)==0));

%% plot parameter
removedNum = 1;
[~,I] = sort(sum((betaEstResult(1:(simulations-50),:)-repmat(betaLower(:)',(simulations-50),1)).^2,2));
orderGood = [I(1:(simulations-50)-removedNum);((simulations-50+1):((simulations-50)+removedNum))'];
coverage = sum(((ones((simulations-50),1)*betaLower(:)')<=(betaEstResult(orderGood,:)+1.96*betaStdResult(orderGood,:))) & ...
    ((ones((simulations-50),1)*betaLower(:)')>=(betaEstResult(orderGood,:)-1.96*betaStdResult(orderGood,:))),1);
% t = tiledlayout(2,3,'TileSpacing','Compact');
% nexttile;
% histogram(betaEstResult(:,1))
% nexttile;
% qqplot(betaEstResult(:,1));
% nexttile;
% histogram(betaStdResult(:,1))
% nexttile;
% histogram(betaEstResult(:,2))
% title(num2str(bands));
% nexttile;
% qqplot(betaEstResult(:,2));
% nexttile;
% histogram(betaStdResult(:,2))
% title(num2str([min(rhoResult) mean(rhoResult) std(rhoResult) max(rhoResult)]));
% title(t,num2str([betaLower(:)';mean(betaEstResult);std(betaEstResult);mean(betaStdResult);coverage],'%5.4f'))
t = figure;
for i = 1:(d*(p-d))
    subplot(d,p-d,i);
    qqplot(betaEstResult(orderGood,i));xlabel('');ylabel('');
    title([num2str(mean(betaEstResult(orderGood,i))) '  ' num2str(betaLower(i))]);
end
saveas(t,[pwd,strcat('/plots/figQQ',num2str(filenum),regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bands(3,2:3),bands(4,2:3)]),'\s+','_'),'.jpg')])

t = figure;
subplot(2,4,1)
hist(betaEstResult(orderGood,1));
subplot(2,4,2)
qqplot(betaEstResult(orderGood,1));xlabel('');ylabel('');
title(num2str([betaLower(:)';nanmean(betaEstResult(orderGood,:));nanstd(betaEstResult(orderGood,:));nanmedian(betaStdResult(orderGood,:));coverage/10]))
subplot(2,4,3)
hist(betaStdResult(orderGood,1))
subplot(2,4,4)
plot(mean(estLambdaAll,2))
subplot(2,4,5)
hist(betaEstResult(orderGood,3))
title(num2str(bands));
subplot(2,4,6)
qqplot(betaEstResult(orderGood,3));xlabel('');ylabel('');
title([samplesize p d]);
subplot(2,4,7)
hist(betaStdResult(orderGood,3))
title(num2str([min(rhoResult) mean(rhoResult) std(rhoResult) max(rhoResult)]));
saveas(t,[pwd,strcat('/plots/figTemp',num2str(filenum),regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bands(3,2:3),bands(4,2:3)]),'\s+','_'),'.jpg')])
beep

%% Store the results
dlmwrite(strcat('flag',num2str(filenum),'.txt'),flags,'delimiter','\t');
dlmwrite(strcat('bands',num2str(filenum),'.txt'),bands,'delimiter','\t');
dlmwrite(strcat('options',num2str(filenum),'.txt'),char(options),'delimiter','');

save(strcat('xall',num2str(filenum),'.mat'),'xallsave');
dlmwrite(strcat('zall',num2str(filenum),'.txt'),zallsave,'delimiter','\t');
dlmwrite(strcat('tall',num2str(filenum),'.txt'),tallsave,'delimiter','\t');
dlmwrite(strcat('cenall',num2str(filenum),'.txt'),cenallsave,'delimiter','\t');
dlmwrite(strcat('deltaall',num2str(filenum),'.txt'),deltaallsave,'delimiter','\t');
dlmwrite(strcat('wall',num2str(filenum),'.txt'),wallsave,'delimiter','\t');
dlmwrite(strcat('groupIndex',num2str(filenum),'.txt'),groupIndexsave,'delimiter','\t');
dlmwrite(strcat('censorInfo',num2str(filenum),'.txt'),censorInfo,'delimiter','\t');

dlmwrite(strcat('trueBeta',num2str(filenum),'.txt'),betaLower(:),'delimiter','\t');
dlmwrite(strcat('resultBeta',num2str(filenum),'.txt'),betaEstResult,'delimiter','\t');
dlmwrite(strcat('resultBetaStd',num2str(filenum),'.txt'),betaStdResult,'delimiter','\t');
dlmwrite(strcat('estHazardGrid',num2str(filenum),'.txt'),combOfVarsN,'delimiter','\t');

