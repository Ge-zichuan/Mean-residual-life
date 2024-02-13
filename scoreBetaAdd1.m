function scoreVal = scoreBetaAdd1(paras1,xallN,zallN,deltaallN,samplesizeN,xallT,zallT,deltaallT,samplesizeT,transIndex,p)
% method from Chen 2007
% paras: parameter \beta in equation (5)
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1
parasall = [eye(1);reshape(paras1,p-1,1)];
scoreVal = 1:p;
integral1N = ones(samplesizeN,p);
integral1T = ones(samplesizeT,p);

[zallordN,sortIndex] = sort(zallN);
xallordN = xallN(sortIndex,:);
deltaordN = deltaallN(sortIndex);
[zallordT,sortIndex] = sort(zallT);
xallordT = xallT(sortIndex,:);
deltaordT = deltaallT(sortIndex);
% calculate \bb\trans\x
btransxN = xallordN*parasall; % paras1 include upper d-by-d block.
btransxT = xallordT*parasall; % paras1 include upper d-by-d block.
% calculate all z_i-z_j
mdiffallN = (zallordN*ones(1,samplesizeN))<=(ones(samplesizeN,1)*zallordN');% T_i<=T_j
denominatorN = sum(mdiffallN,2);
mdiffallT = (zallordT*ones(1,samplesizeT))<=(ones(samplesizeT,1)*zallordT');% T_i<=T_j
denominatorT = sum(mdiffallT,2);
% calculate estimating equation, each term.
xbarN = (mdiffallN*xallordN)./repmat(denominatorN,1,p);
xbarT = (mdiffallT*xallordT)./repmat(denominatorT,1,p);
SNAinsideT = exp(-(deltaordT'*mdiffallT)'./denominatorT);
QinsideT = mdiffallT*exp(-btransxT(:,1))./denominatorT;
SNAinsideN = exp(-(deltaordN'*mdiffallN)'./denominatorN);
QinsideN = mdiffallN*exp(-btransxN(:,1))./denominatorN;

for i = 1:p
    innerTemp = repmat(xallordN(:,i)-xbarN(:,i),1,samplesizeN).*mdiffallN;
    integral1N(:,i) = sum((innerTemp(1:(end-1),:)+innerTemp(2:end,:))/2.*...
        repmat(diff(zallordN),[1 samplesizeN]),1);
    innerTemp = repmat(xallordT(:,i)-xbarT(:,i),1,samplesizeT).*mdiffallT;
    integral1T(:,i) = sum((innerTemp(1:(end-1),:)+innerTemp(2:end,:))/2.*...
        repmat(diff(zallordT),[1 samplesizeT]),1);
end
% calculate \wh m(t) first at each z_i
innerTemp = SNAinsideN.*QinsideN;
integral2N = sum((innerTemp(1:(end-1))+innerTemp(2:end)).*diff(zallordN)/2)-...
    [0 ; cumsum((innerTemp(1:(end-1))+innerTemp(2:end))/2.*diff(zallordN))];
mhatN = integral2N./SNAinsideN;
innerTemp = SNAinsideT.*QinsideT;
integral2T = sum((innerTemp(1:(end-1))+innerTemp(2:end)).*diff(zallordT)/2)-...
    [0 ; cumsum((innerTemp(1:(end-1))+innerTemp(2:end))/2.*diff(zallordT))];
mhatT = integral2T./SNAinsideT;
for i = 1:p
    scoreValall = (xallordN(:,i)-xbarN(:,i)).*mhatN.*deltaordN...
        -exp(-btransxN(:,1)).*integral1N(:,i);
    scoreValall(transIndex) = (xallordT(:,i)-xbarT(:,i)).*mhatT.*deltaordT...
        -exp(-btransxT(:,1)).*integral1T(:,i);
    scoreVal(i) = mean(scoreValall);
end
end