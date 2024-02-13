function wout = gendata_w(n,wtype,optionw)
% generate data X, n: sample size, xdim: dimension of each observation
% xtype---
%       norm: standard normal ranform variable
%       unif: standard uniform ranform variable
switch wtype
    case 'unif'
        wout = rand(n,1)*optionw;
    case 'unif10'
        wout = rand(n,1)*10;
    case 'mimic'
        wout = rand(n,1)*optionw;
end
end