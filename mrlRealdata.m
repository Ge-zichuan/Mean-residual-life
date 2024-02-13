filenum1 = 401; % file number to store different numeric result
tonum = 10; % number of repeat times
p = 217;d = 2; % dimension of parameter \gamma (q-by-d, top d-by-d is identity matrix) (model logit(S)~\gamma\Z)
samplen = 400; % sample size
tt = 1;
bands = zeros(3); % bandwidth adjustment in nonparametric regression
resultBeta = zeros(tonum,p+1);
resultMX = zeros(tonum,d*samplen2);
resultMY = zeros(tonum,samplen2);
resultBetaSd = zeros(tonum,p+1);
etaXEst = zeros(tonum,100,2);
etaYEst = zeros(tonum,100,100);
etaest = zeros(samplen2,d+1);
etaprime = ones(samplen1,d);
betaTrue = rand((p-d),d,1); % initial guess of \gamma or just random number, excluding upper d-by-d block.
% you can use [1,2,3,4] to generate sequence, [1,2;3,4;5,6;7,8] is 4 by 2 matrix. ";" starts a new row.
while (tt<=tonum)
  %set initial parameters
  % generate data
  fileID = fopen('xApp.txt','r');
  xall = fscanf(fileID,'%f',[samplen1,p]); % read covariates x of validation data set
  fclose(fileID);
  fileID = fopen('zApp.txt','r');
  zall = fscanf(fileID,'%f',[samplen1,1]); % read covariates z of validation data set
  fclose(fileID);
  fileID = fopen('deltaApp.txt','r');
  deltaall = fscanf(fileID,'%f',[samplen1,1]); % read index s of validation data set
  fclose(fileID);
  %%%%%%% calculating bandwidth %%%%%%%%%%%%%
  ntemp = samplen1^(-1/3);
  bandez1 = sqrt(sum(std(xall))/p)*ntemp;
  bandb = sqrt(sum(std(zall)))*ntemp;
  bands(1,:) = [bandez1,1,1]; % tuning the bandwidth in estimating \beta
  bands(2,:) = [bandb,1,1]; % tuning the bandwidth in estimating \beta
  bands(3,:) = [bandez1,1,1]; % tuning the bandwidth in estimating \beta
  %%% solving estimating equation
  %%% gamma-->eta-->gamma(optional)-->beta
  %%% gamma
  options1 = optimoptions('fsolve','Display','off',...
    'Algorithm','trust-region-dogleg','MaxIterations',1000,'MaxFunctionEvaluations',50000,'OptimalityTolerance',1e-6);

  paras0 = betaTrue+(rand((p-d),d,1)-0.5)*0.1;
  [betaEst,~,eflag,~] = fsolve(@(x)scoreBeta(x,xall,zall,deltaall,samplen,p,d,bands),paras0(:),options1);
  io = 1;
  while (eflag > 2 && io<1)
    paras0 = betaTrue+(rand((p-d),d)-0.5)*0.01;
    [betaEst,~,eflag,~] = fsolve(@(x)scoreBeta(x,xall,zall,deltaall,samplen,p,d,bands),paras0(:),options1);
    io = io+1;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% save result
  resultBeta(tt,:) = betaEst(:);
%   resultGamma(tt,:) = gammaEst;
  tt = tt+1
end
save(strcat('betaResult',num2str(filenum1),'.mat'),'resultBeta');
