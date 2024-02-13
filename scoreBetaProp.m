function scoreVal = scoreBetaProp(paras1,xall,zall,deltaall,samplen,p)
% paras: parameter \beta in equation (5)
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1

  % khbetax: K_h(\bb\trans\x)
  % denoVev: \sumi I()K_h()
	scoreVal = 1:p;
  xbar = zeros(samplen,p);
  ainside = ones(samplen,1);
  xbinside = ones(samplen,1);
  secinside = ones(samplen,p);
  mhat = 1:samplen;
  
  [zallord,sortIndex] = sort(zall);
  xallord = xall(sortIndex,:);
  deltaord = deltaall(sortIndex);
  % calculate \bb\trans\x
  parasAll = paras1;
  btransx = xallord*parasAll; % paras1 include upper d-by-d block.
  % calculate \wh m(t) first at each z_i
  mdiffall = (zallord*ones(1,samplen))<=(ones(samplen,1)*zallord');

  % calculate all z_i-z_j
  zdiffall = [zallord(1);zallord(2:end)-zallord(1:(end-1))];
  denominator1 = sum(mdiffall,2);
  % calculate estimating equation, each term.
  for i = 1:samplen
    for j = 1:p
      xbar(i,j) = sum(xallord(:,j).*mdiffall(i,:)')/sum(mdiffall(i,:));
    end
    ainside(i,1) = exp(-sum(deltaord(1:i)./denominator1(1:i)));
    xbinside(i) = sum(exp(-btransx(:,1)).*mdiffall(i,:)')/sum(mdiffall(i,:));
  end
  for i = 1:samplen
    for j = 1:p
      secinside(i,j) = sum(zdiffall(1:i).*xbar(1:i,j));
    end
  end
  for i = 1:(samplen-1)
    mhat(i) = sum(zdiffall((i+1):samplen).*ainside((i):(samplen-1)).*xbinside((i+1):samplen))/ainside(i);
  end
  mhat(samplen) = mhat(samplen-1);
  for i = 1:p
    scoreVal(i) = sum((xallord(:,i)-xbar(:,i)).*mhat'.*deltaord...
      -xallord(:,i).*exp(-btransx(:,1)).*zallord+secinside(:,i).*exp(-btransx(:,1)))/samplen;
  end
end