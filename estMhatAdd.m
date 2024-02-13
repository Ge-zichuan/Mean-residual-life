function mhat = estMhatAdd(paras1,xall,zall,deltaall,samplen)
% mean residual life function estimation by add model
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1
  [zallord,sortIndex] = sort(zall);
  xallord = xall(sortIndex,:);
  deltaord = deltaall(sortIndex);
  firstinside = zeros(samplen,1);
  denominator1 = zeros(samplen,1);
  secinside = zeros(samplen,1);
  
  btransx = xallord*paras1; % paras1 include upper d-by-d block.
  zdiffall = zallord*ones(1,samplen)-ones(samplen,1)*zallord';
  mdiffall = zdiffall<=0;
  zdiff = [zallord(1);zallord(2:end)-zallord(1:(end-1))];
  denominator1 = sum(mdiffall,2);
  ainside = exp(-cumsum(deltaord./denominator1));

  for i = 1:samplen
    firstinside(i) = sum(zdiff((i+1):samplen).*ainside((i):(samplen-1)));
    secinside(i) = sum(btransx(i:samplen,1).*deltaord(i:samplen).*ainside(i:samplen)...
      ./denominator1(i:samplen));
  end
  mhat = (firstinside-secinside)./ainside;

end