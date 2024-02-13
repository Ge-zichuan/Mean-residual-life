function mhat = estMhatPropNew(paras1,xall,zall,deltaall,samplen)
% mean residual life function estimation by proportional model 1
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1
  
  mhat = zeros(samplen,1);
  [zallord,sortIndex] = sort(zall);
  xallord = xall(sortIndex,:);
  deltaord = deltaall(sortIndex);
  xbinside = zeros(samplen,1);
  ainside = zeros(samplen,1);
  
  btransx = xallord*paras1; % paras1 include upper d-by-d block.
  zdiffall = zallord*ones(1,samplen)-ones(samplen,1)*zallord';
  mdiffall = zdiffall<=0;
  zdiff = [zallord(1);zallord(2:end)-zallord(1:(end-1))];
  denompart = sum(mdiffall,2);
  
  for t = 1:samplen
    ainside(t) = exp(-sum(mdiffall(:,t).*deltaord./denompart));
    xbinside(t) = sum(exp(-btransx(:,1)).*mdiffall(t,:)')/sum(mdiffall(t,:));
  end

  for i = 1:(samplen-1)
    mhat(i) = 1/ainside(i);
    mhat(i) = mhat(i)*sum(zdiff((i+1):samplen)/2.*(ainside((i+1):samplen).*...
      xbinside((i+1):samplen)+ainside(i:(samplen-1)).*xbinside(i:(samplen-1))));
  end
end