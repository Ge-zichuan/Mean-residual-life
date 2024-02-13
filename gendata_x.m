function xout = gendata_x(samplesize,xdimp,xtype)
% generate data X, n: sample size, xdim: dimension of each observation
% xtype---
%       norm: standard normal ranform variable
%       unif: standard uniform ranform variable
switch xtype
    case 'unif'
        xout = rand(samplesize,xdimp);
    case 'norm'
        xout = normrnd(0,1,samplesize,xdimp);
    case 'pm1a'
        xout(:,1) = binornd(1,0.5,samplesize);
        xout(:,2) = rand(samplesize,1);
    case 'pm2a'
        xout = binornd(1,0.5,samplesize,ximp);
    case 'ma62'
        matSigma = 0.5.^abs(repmat(0:(2-1),2,1)-repmat(0:(2-1),2,1)');
        dataF12 = mvnrnd(zeros(2,1),matSigma,samplesize);
        errorxi1 = normrnd(0,1,samplesize,1);
        dataF3 = abs(sum(dataF12,2))+dataF12(:,1).*errorxi1;
        errorxi2 = normrnd(0,1,samplesize,1);
        dataF4 = sum(dataF12,2).^2+abs(dataF12(:,2)).*errorxi2;
        dataF5 = binornd(1,logsig(dataF12(:,2)),samplesize,1);
        dataF6 = binornd(1,normcdf(dataF12(:,2)),samplesize,1);
%         dataF7 = binornd(1,normcdf(dataF3),samplesize,1);
%         dataF8 = normrnd(0,1,samplesize,1);
%         dataF9 = normrnd(0,1,samplesize,1);
        xout = [dataF12, dataF3, dataF4, dataF5, dataF6];
    case 'mimic'
        temp = rand(samplesize,1);
        X1 = 1*(temp<=0.6)+2*(temp>0.6);
        temp = rand(samplesize,1);
        X2 = 1*(temp<=0.5)+2*(temp>0.5).*(temp<=0.8)...
            +3*(temp>0.8).*(temp<=.95)+4*(temp>0.95);
        X3 = rand(samplesize,1)*100;
        temp = rand(samplesize,1);
        X4 = 1*(temp<=0.65)+2*(temp>0.65);
        X5 = exprnd(18,samplesize,1);
        X6 = randi([1,11],samplesize,1);
        X7 = normrnd(0,20,samplesize,1);
        temp = rand(samplesize,1);
        X8 = 1*(temp<=0.9)+2*(temp>0.9).*(temp<=0.95)+3*(temp>0.95);
        temp = rand(samplesize,1);
        X9 = 0*(temp<=0.3)+1*(temp>0.3);
       
        xout = [X1,X2,X3,X4,X5,X6,X7,X8,X9];
end
end