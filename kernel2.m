function kernelVal2 = kernel2(xVal,whichd,d,typeKert)
% First derivative of kernel function.
% whichd: which index in X represent the dimension of X, not the smaples we want to
% compute. If dimension of X is 1, whichd = 1;
switch (typeKert)
    case 'quar' % this is quaratic kernel
        xVals1 = 0.9375*(1-xVal.^2).^2.*(abs(xVal)<=1);
        xVals2 = 4*0.9375*xVal.*(xVal.^2-1).*(abs(xVal)<=1);
    case 'norm' % this use pdf of normal distribution as kernel
        xVals1 = 1/sqrt(2*pi)*exp(-xVal.^2/2);
        xVals2 = 1/sqrt(2*pi)*exp(-xVal.^2/2).*(-xVal);
    case 'Epanechnikov' % this is Epanechnikov kernel
        xVals1 = 0.75*(1-xVal.^2).*(abs(xVal)<=1);
        xVals2 = -1.5*xVal.*(abs(xVal)<=1);
end
if (whichd > 1)
    xVals3 = repmat(prod(xVals1,whichd),[ones(1,length(size(xVal))-1),d])./xVals1;
else
    xVals3 = ones(size(xVals1));
end
if d == 1
    kernelVal2 = xVals2;
else
    kernelVal2 = xVals3.*xVals2;
    kernelVal2(isnan(kernelVal2)) = 0;
end

end