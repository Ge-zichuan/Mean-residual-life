filenum1 = 101; % file number to store different numeric result
simulations = 1000; % number of repeat times
samplen = 500; % sample size
censors = (1-0.95)*samplen;
bands = zeros(3); % bandwidth adjustment in nonparametric regression
% resultBetaSd = zeros(simulations,(p-d)*d);
% betaTrue = [1,-0.75,0.4,-0.2,0.15,0,-0.15,0.2,-0.4,0.75]'; % study 1
% betaTrue = [1,1.3,-1.3,1,-0.5,0.5,-0.5]'; % study 2
betaTrue = [1,0,2.75,-0.75,-1,2;0,1,-3.125,-1.125,1,-2]'; % study 3
dims = size(betaTrue);
p = dims(1);d = dims(2); % dimension of parameter \gamma (q-by-d, top d-by-d is identity matrix) (model logit(S)~\gamma\Z)
resultBeta = zeros(simulations,(p-d)*d);
resultBetaAdd = zeros(simulations,p);
resultBetaProp = zeros(simulations,p);
resultBetaProp2 = zeros(simulations,p);
estmAdd = zeros(simulations,samplen);
estmProp = zeros(simulations,samplen);
estmProp2 = zeros(simulations,samplen);
xallAll = zeros(simulations,samplen,p);
zallAll = zeros(simulations,samplen);
tallAll = zeros(simulations,samplen);
callAll = zeros(simulations,samplen);
deltaallAll = zeros(simulations,samplen);
censoringDist = zeros(simulations,samplen);
censoringK = zeros(simulations,1);
xtype = 'norm';ttype = 'lnl2';ctype = 'gamm';
tt = 1;
while (tt<=simulations)
  %set initial parameters
  % generate data
  xall = gendata_x(samplen,p,xtype);
  tall = gendata_t(xall,samplen,ttype,betaTrue);
  [call,deltaall,kval,randx] = gendata_cen(tall,xall,ctype,censors,samplen);
  zall = min(tall,call);
  censoringK(tt) = kval;
  tt = tt+1
end
[min(censoringK) max(censoringK) mean(censoringK) std(censoringK)]