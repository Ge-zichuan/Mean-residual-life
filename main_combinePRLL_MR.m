%%% basic information of numerical study
filenum = 420;
numOfPRLL = 75;
simulations = 1050;
[filefieldname,filefieldval] = optionRead(filenum);
p           = str2num(filefieldval{1});
d           = str2num(filefieldval{2});
samplesize  = str2num(filefieldval{3});
eachsim = simulations/numOfPRLL;

if ~exist(strcat('xall',num2str(filenum),'.mat'))
    betaEstResult = zeros(simulations,d*(p-d));
    betaStdResult = zeros(simulations,d*(p-d));
    rhoResult = zeros(simulations,1);
    flags = zeros(simulations,2);
    xallsave1 = zeros(simulations,samplesize,p);
    zallsave = zeros(simulations,samplesize);
    tallsave = zeros(simulations,samplesize);
    cenallsave = zeros(simulations,samplesize);
    wallsave = zeros(simulations,samplesize);
    deltaallsave = zeros(simulations,samplesize);
    groupIndexsave = zeros(simulations,samplesize);
    censorInfo = zeros(simulations,1);
    display(simulations/numOfPRLL);
    test = [];
    for i = 1:numOfPRLL
        betaEstResult1 = dlmread(strcat('resultBeta',num2str(filenum),'_',num2str(i),'.txt'));
        test = [test;[i length(betaEstResult1(:,1))]];
    end
    display(test);
    for i = 1:numOfPRLL

        display(i);
        betaEstResult1 = dlmread(strcat('resultBeta',num2str(filenum),'_',num2str(i),'.txt'));
        zall = dlmread(strcat('zall',num2str(filenum),'_',num2str(i),'.txt'));
        tall = dlmread(strcat('tall',num2str(filenum),'_',num2str(i),'.txt'));
        cenall = dlmread(strcat('cenall',num2str(filenum),'_',num2str(i),'.txt'));
        deltaall = dlmread(strcat('deltaall',num2str(filenum),'_',num2str(i),'.txt'));
        wall = dlmread(strcat('wall',num2str(filenum),'_',num2str(i),'.txt'));
        groupIndex = dlmread(strcat('groupIndex',num2str(filenum),'_',num2str(i),'.txt'));
        censorInfo1 = dlmread(strcat('censorInfo',num2str(filenum),'_',num2str(i),'.txt'));
        betaEstStdResult1 = dlmread(strcat('resultBetaStd',num2str(filenum),'_',num2str(i),'.txt'));
        flag = dlmread(strcat('flag',num2str(filenum),'_',num2str(i),'.txt'));
        load(strcat('xall',num2str(filenum),'_',num2str(i),'.mat'));
        %         if (length(betaEstResult1(:,1))<eachsim*i)
        %             lackNum = eachsim*i-length(betaEstResult1(:,1));
        %             betaEstResult1 = [betaEstResult1;betaEstResult1(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             zall = [zall;zall(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             tall = [tall;tall(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             cenall = [cenall;cenall(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             deltaall = [deltaall;deltaall(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             wall = [wall;wall(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             groupIndex = [groupIndex;groupIndex(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             censorInfo1 = [censorInfo1;censorInfo1(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             betaEstStdResult1 = [betaEstStdResult1;betaEstStdResult1(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             flag = [flag;flag(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:)];
        %             load(strcat('xall',num2str(filenum),'_',num2str(i),'.mat'));
        %             temp = [xallsave;xallsave(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:,:)];
        %             [size(xallsave) size(xallsave(((i-1)*eachsim+1):((i-1)*eachsim+lackNum),:,:))]
        %             xallsave1(((i-1)*eachsim+1):(i*eachsim),:,:) = temp(((i-1)*eachsim+1):(i*eachsim),:,:);
        %         else
        %             load(strcat('xall',num2str(filenum),'_',num2str(i),'.mat'));
        %             xallsave1(((i-1)*eachsim+1):(i*eachsim),:,:) = xallsave(((i-1)*eachsim+1):(i*eachsim),:,:);
        %         end
        xallsave1(((i-1)*eachsim+1):(i*eachsim),:,:) = xallsave;
        betaEstResult(((i-1)*eachsim+1):(i*eachsim),:) = betaEstResult1;
        betaStdResult(((i-1)*eachsim+1):(i*eachsim),:) = betaEstStdResult1;
        zallsave(((i-1)*eachsim+1):(i*eachsim),:) = zall;
        tallsave(((i-1)*eachsim+1):(i*eachsim),:) = tall;
        cenallsave(((i-1)*eachsim+1):(i*eachsim),:) = cenall;
        deltaallsave(((i-1)*eachsim+1):(i*eachsim),:) = deltaall;
        wallsave(((i-1)*eachsim+1):(i*eachsim),:) = wall;
        groupIndexsave(((i-1)*eachsim+1):(i*eachsim),:) = groupIndex;
        censorInfo(((i-1)*eachsim+1):(i*eachsim),:) = censorInfo1;
        flags(((i-1)*eachsim+1):(i*eachsim),:) = flag;
    end
    mean(groupIndexsave,2)
    xallsave = xallsave1;
    save(strcat('xall',num2str(filenum),'.mat'),'xallsave');
    dlmwrite(strcat('zall',num2str(filenum),'.txt'),zallsave,'delimiter','\t');
    dlmwrite(strcat('tall',num2str(filenum),'.txt'),tallsave,'delimiter','\t');
    dlmwrite(strcat('cenall',num2str(filenum),'.txt'),cenallsave,'delimiter','\t');
    dlmwrite(strcat('deltaall',num2str(filenum),'.txt'),deltaallsave,'delimiter','\t');
    dlmwrite(strcat('wall',num2str(filenum),'.txt'),wallsave,'delimiter','\t');
    dlmwrite(strcat('groupIndex',num2str(filenum),'.txt'),groupIndexsave,'delimiter','\t');
    dlmwrite(strcat('censorInfo',num2str(filenum),'.txt'),censorInfo,'delimiter','\t');
    dlmwrite(strcat('flags',num2str(filenum),'.txt'),flags,'delimiter','\t');
    
    dlmwrite(strcat('resultBeta',num2str(filenum),'.txt'),betaEstResult,'delimiter','\t');
    dlmwrite(strcat('resultBetaStd',num2str(filenum),'.txt'),betaStdResult,'delimiter','\t');
    
    betaTrue = dlmread(strcat('trueBeta',num2str(filenum),'.txt'));
    dlmwrite(strcat('trueBeta',num2str(filenum),'.txt'),betaTrue,'delimiter','\t');
    bands = dlmread(strcat('bands',num2str(filenum),'.txt'));
    dlmwrite(strcat('bands',num2str(filenum),'.txt'),bands,'delimiter','\t');
    options = filefieldval;
    dlmwrite(strcat('options',num2str(filenum),'.txt'),char(options),'delimiter','');
    
    
else
    display('existed');
    options = filefieldval;
    betaTrue = dlmread(strcat('trueBeta',num2str(filenum),'.txt'));
    betaEstResult = dlmread(strcat('resultBeta',num2str(filenum),'.txt'));
    betaStdResult = dlmread(strcat('resultBetaStd',num2str(filenum),'.txt'));
    bands = dlmread(strcat('bands',num2str(filenum),'.txt'));
    load(strcat('xall',num2str(filenum),'.mat'));
    zallsave = dlmread(strcat('zall',num2str(filenum),'.txt'));
    tallsave = dlmread(strcat('tall',num2str(filenum),'.txt'));
    cenallsave = dlmread(strcat('cenall',num2str(filenum),'.txt'));
    deltaallsave = dlmread(strcat('deltaall',num2str(filenum),'.txt'));
    wallsave = dlmread(strcat('wall',num2str(filenum),'.txt'));
    groupIndexsave = dlmread(strcat('groupIndex',num2str(filenum),'.txt'));
    flags = dlmread(strcat('flags',num2str(filenum),'.txt'));
end

[betaTrue';median(betaEstResult,1);mean(betaStdResult,1)]
removedNum = 1;
[~,I] = sort(sum((betaEstResult(1:1000,:)-repmat(betaTrue',1000,1)).^2,2));
orderGood = [I(1:1000-removedNum);(1001:(1000+removedNum))'];
coverage = sum(((ones(1000,1)*betaTrue')<=(betaEstResult(orderGood,:)+1.96*betaStdResult(orderGood,:))) & ...
    ((ones(1000,1)*betaTrue')>=(betaEstResult(orderGood,:)-1.96*betaStdResult(orderGood,:))),1);
t = figure;
for i = 1:(d*(p-d))
    subplot(d,p-d,i);
    qqplot(betaEstResult(orderGood,i));xlabel('');ylabel('');
    title([num2str(mean(betaEstResult(orderGood,i))) '  ' num2str(betaTrue(i))]);
end
saveas(t,[pwd,strcat('\plots\figQQ',num2str(filenum),'_',regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bands(3,2:3),bands(4,2:3)]),'\s+','_'),'.jpg')])

t = figure;
subplot(1,2,1)
hist(zallsave(find(groupIndexsave==0)));xlabel('group=0');ylabel('');title('');
subplot(1,2,2)
hist(zallsave(find(groupIndexsave==1)));xlabel('group=1');ylabel('');title('');
saveas(t,[pwd,strcat('/plots/figZHist',num2str(filenum),'_',regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bands(3,2:3),bands(4,2:3)]),'\s+','_'),'.jpg')])

t = figure;
subplot(2,4,1)
hist(betaEstResult(orderGood,1));xlabel('');ylabel('');title('');
subplot(2,4,2)
qqplot(betaEstResult(orderGood,1));xlabel('');ylabel('');
title(num2str([betaTrue';nanmedian(betaEstResult(orderGood,:));nanstd(betaEstResult(orderGood,:));nanmedian(betaStdResult(orderGood,:));coverage/10]))
subplot(2,4,3)
hist(betaStdResult(orderGood,1));xlabel('');ylabel('');title('');
subplot(2,4,4)
%    plot(mean(estLambdaAll,2))
subplot(2,4,5)
hist(betaEstResult(orderGood,3))
title(num2str(bands));
subplot(2,4,6)
qqplot(betaEstResult(orderGood,3));xlabel('');ylabel('');title('');
subplot(2,4,7)
hist(betaStdResult(orderGood,3));
title([samplesize p d mean(flags)]);
saveas(t,[pwd,strcat('/plots/figTemp',num2str(filenum),'_',regexprep(num2str([samplesize,bands(1,3:4),bands(2,2:3),bands(3,2:3),bands(4,2:3)]),'\s+','_'),'.jpg')])


%% Assigning stores
if exist(strcat('estMNon',num2str(filenum),'_1.txt'))
    combOfVarsN = dlmread(strcat('estHazardGridNon',num2str(filenum),'.txt'));
    combOfVarsT = dlmread(strcat('estHazardGridTrans',num2str(filenum),'.txt'));
    mResultTsave = zeros(simulations,length(combOfVarsT(:,1)));
    mResultNsave = zeros(simulations,length(combOfVarsN(:,1)));
    tausave = zeros(simulations,2);
    %% iterations
    for i = 1:numOfPRLL
        mResultN = dlmread(strcat('estMNon',num2str(filenum),'_',num2str(i),'.txt'));
        mResultT = dlmread(strcat('estMTrans',num2str(filenum),'_',num2str(i),'.txt'));
        mResultTsave(((i-1)*eachsim+1):(i*eachsim),:) = mResultT(((i-1)*eachsim+1):(i*eachsim),:);
        mResultNsave(((i-1)*eachsim+1):(i*eachsim),:) = mResultN(((i-1)*eachsim+1):(i*eachsim),:);
%         tau = dlmread(strcat('tau',num2str(filenum),'_',num2str(i),'.txt'));
%         tausave(((i-1)*eachsim+1):(i*eachsim),:) = tau(((i-1)*eachsim+1):(i*eachsim),:);
    end
    %% Store the results
    dlmwrite(strcat('estMNon',num2str(filenum),'.txt'),mResultNsave,'delimiter','\t');
    dlmwrite(strcat('estMTrans',num2str(filenum),'.txt'),mResultTsave,'delimiter','\t');
%     dlmwrite(strcat('tau',num2str(filenum),'.txt'),tausave,'delimiter','\t');
    display('m finished');
end


if exist(strcat('estBetaPM',num2str(filenum),'_1.txt'))
    disp('PM');
    combOfVarsN = dlmread(strcat('estHazardGridNon',num2str(filenum),'.txt'));
    combOfVarsT = dlmread(strcat('estHazardGridTrans',num2str(filenum),'.txt'));
    betaResultPM = zeros(simulations,p-d+1);
    mResultPMNsave1 = zeros(simulations,length(combOfVarsN(:,1)));
    mResultPMTsave1 = zeros(simulations,length(combOfVarsT(:,1)));
    %% iterations
    for i = 1:numOfPRLL
        betaResultPM1 = dlmread(strcat('estBetaPM',num2str(filenum),'_',num2str(i),'.txt'));
        mResultN = dlmread(strcat('estMPMNon',num2str(filenum),'_',num2str(i),'.txt'));
        mResultT = dlmread(strcat('estMPMTrans',num2str(filenum),'_',num2str(i),'.txt'));
        betaResultPM(((i-1)*eachsim+1):(i*eachsim),:) = betaResultPM1;
        mResultPMNsave1(((i-1)*eachsim+1):(i*eachsim),:) = mResultN;
        mResultPMTsave1(((i-1)*eachsim+1):(i*eachsim),:) = mResultT;
    end
    %% Store the results
    dlmwrite(strcat('estBetaPM',num2str(filenum),'.txt'),betaResultPM,'delimiter','\t');
    dlmwrite(strcat('estMPMNon',num2str(filenum),'.txt'),mResultPMNsave1,'delimiter','\t');
    dlmwrite(strcat('estMPMTrans',num2str(filenum),'.txt'),mResultPMTsave1,'delimiter','\t');
    display('PM finished');
end