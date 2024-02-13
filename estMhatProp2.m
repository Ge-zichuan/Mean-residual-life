function mhat = estMhatProp2(paras1,xall,zall,deltaall,samplen)
% mean residual life function estimation by proportional hazard model 2
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1
  
  mhat = zeros(samplen,1);
  [zallord,sortIndex] = sort(zall);
  xallord = xall(sortIndex,:);
  deltaord = deltaall(sortIndex);
  xbinside = zeros(samplen,1);
  denominator1 = zeros(samplen,1);
  numerator = zeros(samplen,1);
  
  btransx = xallord*paras1; % paras1 include upper d-by-d block.
  zdiffall = zallord*ones(1,samplen)-ones(samplen,1)*zallord';
  mdiffall = zdiffall<=0;
  finside = 1-(1-deltaord)./sum(mdiffall,2);
  fhat = cumprod(finside);
  fhat(end) = (fhat(end)==0)*fhat(end-1)+fhat(end)*(fhat(end)~=0);

  for i = 1:samplen
    xbinside(i) = sum(exp(-2.0*btransx(:,1)).*deltaord...
      .*zdiffall(:,i).*mdiffall(i,:)'./fhat)/samplen;
  end
  for i = 1:samplen
    denominator1(i) = sum(exp(-btransx(:,1)).*mdiffall(i,:)'.*...
      deltaord./fhat)/samplen;
    numerator(i) = sum(exp(-btransx(:,1)).*deltaord.*zallord...
        ./fhat)/samplen;
    if (denominator1(i)==0)
      denominator1(i) = denominator1(i-1);
    end
  end
  
  mhat = xbinside.*numerator./denominator1/samplen;
end