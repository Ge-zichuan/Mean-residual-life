filenum1 = 101; % file number to store different numeric result
simulations = 1000; % number of repeat times
p = 10;d = 1; % dimension of parameter \gamma (q-by-d, top d-by-d is identity matrix) (model logit(S)~\gamma\Z)
samplen = 500; % sample size
censors = (1-0.2)*samplen;
tt = 1;
bands = zeros(3); % bandwidth adjustment in nonparametric regression
resultBeta = zeros(simulations,(p-d)*d);
resultBetaAdd = zeros(simulations,p);
resultBetaProp = zeros(simulations,p);
resultBetaProp2 = zeros(simulations,p);
estmAdd = zeros(simulations,samplen);
estmProp = zeros(simulations,samplen);
estmProp2 = zeros(simulations,samplen);
% resultBetaSd = zeros(simulations,(p-d)*d);
betaTrue = [1,-0.75,0.4,-0.2,0.15,0,-0.15,0.2,-0.4,0.75]'; % study 1
xallAll = zeros(simulations,samplen,p);
zallAll = zeros(simulations,samplen);
tallAll = zeros(simulations,samplen);
callAll = zeros(simulations,samplen);
deltaallAll = zeros(simulations,samplen);
censoringDist = zeros(simulations,samplen);
censoringK = zeros(simulations,1);
xtype = 'norm';ttype = 'lnl2';ctype = 'gamm';
while (tt<=simulations)
  %set initial parameters
  % generate data
  xall = gendata_x(samplen,p,xtype);
  tall = gendata_t(xall,samplen,ttype,betaTrue);
  [call,deltaall,kval,randx] = gendata_cen(tall,xall,ctype,censors,samplen);
  zall = min(tall,call);
%   sum(deltaall)/samplen
  %%%%%%% calculating bandwidth %%%%%%%%%%%%%
  ntemp = samplen^(-1/4);
  bandez1 = sqrt(sum(std(xall))/p)*ntemp;
  bandb = sqrt(sum(std(zall)))*ntemp;
  bands(1,:) = [bandez1,1,1]; % tuning the bandwidth in estimating \beta
  bands(2,:) = [bandb,1,1]; % tuning the bandwidth in estimating \beta
  bands(3,:) = [bandez1,1,1]; % tuning the bandwidth in estimating \beta
  %%% solving estimating equation
  %%% beta
  options1 = optimoptions('fsolve','Display','off',...
    'Algorithm','trust-region-dogleg','MaxIterations',1000,'MaxFunctionEvaluations',50000,'OptimalityTolerance',1e-6);
  paras0 = betaTrue(((d+1):end),:)+(rand((p-d),d,1)-0.5)*0.1;
  paras0additon = betaTrue;
  [betaEst,~,eflag,~] = fsolve(@(x)scoreBeta(x,xall,zall,deltaall,samplen,p,d,bands),paras0,options1);
  [betaEstAdd,~,~,~] = fsolve(@(x)scoreBetaAdd(x,xall,zall,deltaall,samplen,p),paras0additon,options1);
  [betaEstProp,~,~,~] = fsolve(@(x)scoreBetaProp(x,xall,zall,deltaall,samplen,p),paras0additon,options1);
  [betaEstProp2,~,~,~] = fsolve(@(x)scoreBetaProp2(x,xall,zall,deltaall,samplen,p),paras0additon,options1);
  io = 1;
  while (eflag > 2 && io<1)
    paras0 = betaTrue+(rand((p-d),d)-0.5)*0.01;
    [betaEst,~,eflag,~] = fsolve(@(x)scoreBeta(x,xall,zall,deltaall,samplen,p,d,bands),paras0(:),options1);
    io = io+1;
  end
  %%% estimate mean residual life function
  mAdd = estMhatAdd(paras0additon,xall,zall,deltaall,samplen);
  mProp = estMhatProp(paras0additon,xall,zall,deltaall,samplen);
  mProp2 = estMhatProp2(paras0additon,xall,zall,deltaall,samplen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% save result
  resultBeta(tt,:) = betaEst(:);
  resultBetaAdd(tt,:) = betaEstAdd(:);
  resultBetaProp(tt,:) = betaEstProp(:);
  resultBetaProp2(tt,:) = betaEstProp2(:);
  estmAdd(tt,:) = mAdd;
  estmProp(tt,:) = mProp;
  estmProp2(tt,:) = mProp2;
  xallAll(tt,:,p) = xall;
  zallAll(tt,:) = zall;
  tallAll(tt,:) = tall;
  callAll(tt,:) = call;
  deltaallAll(tt,:) = deltaall;
  censoringDist(tt,:) = randx;
  censoringK(tt,1) = kval;
  tt = tt+1
end
save(strcat('betaResult',num2str(filenum1),'new.mat'),'resultBeta');
save(strcat('betaAddResult',num2str(filenum1),'new.mat'),'resultBetaAdd');
save(strcat('betaPropResult',num2str(filenum1),'new.mat'),'resultBetaProp');
save(strcat('betaProp2Result',num2str(filenum1),'new.mat'),'resultBetaProp2');
save(strcat('mAddResult',num2str(filenum1),'new.mat'),'estmAdd');
save(strcat('mPropResult',num2str(filenum1),'new.mat'),'estmProp');
save(strcat('mProp2Result',num2str(filenum1),'new.mat'),'estmProp2');
save(strcat('X',num2str(filenum1),'new.mat'),'xallAll');
save(strcat('Z',num2str(filenum1),'new.mat'),'zallAll');
save(strcat('T',num2str(filenum1),'new.mat'),'tallAll');
save(strcat('C',num2str(filenum1),'new.mat'),'callAll');
save(strcat('Delta',num2str(filenum1),'new.mat'),'deltaallAll');
save(strcat('censoringDist',num2str(filenum1),'new.mat'),'censoringDist');
save(strcat('censoringK',num2str(filenum1),'new.mat'),'censoringK');
save(strcat('bandsAdj',num2str(filenum1),'new.mat'),'bands');

