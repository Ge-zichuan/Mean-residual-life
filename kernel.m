function kernelVal1 = kernel(xVal,whichd,typeKert)
% whichd: which index in X represent the dimension of X, not the smaples we want to
% compute. If dimension of X is 1, whichd = 1;
switch (typeKert)
    case 'quar' % this is quaratic kernel
        xVals1 = 0.9375*(1-xVal.^2).^2.*(abs(xVal)<=1);
    case 'norm' % this use pdf of normal distribution as kernel
        xVals1 = 1/sqrt(2*pi)*exp(-xVal.^2/2);
    case 'Epanechnikov' % this is Epanechnikov kernel
        xVals1 = 0.75*(1-xVal.^2).*(abs(xVal)<=1);
end
if (whichd>1)
    kernelVal1 = prod(xVals1,whichd);
else
    kernelVal1 = xVals1;
end
if (length(xVals1(:,1))==1)
    kernelVal1 = xVals1';
end
end