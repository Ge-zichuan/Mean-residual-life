function scoreVal = scoreBetaAdd(paras1,xall,zall,deltaall,samplen,p)
% paras: parameter \beta in equation (5)
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1

  % khbetax: K_h(\bb\trans\x)
  % denoVev: \sumi I()K_h()
	scoreVal = 1:p;
  xbar = zeros(samplen,p);
  ainside = ones(samplen,1);
  firstinside = ones(samplen,1);
  secinside = ones(samplen,1);
  [zallord,sortIndex] = sort(zall);
  xallord = xall(sortIndex,:);
  deltaord = deltaall(sortIndex);
  % calculate \bb\trans\x
  btransx = xallord*paras1; % paras1 include upper d-by-d block.
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
  end
  for i = 1:samplen
    firstinside(i,1) = sum(zdiffall((i+1):samplen).*ainside((i):(samplen-1),1));
    secinside(i,1) = sum(btransx(i:samplen,1).*deltaord(i:samplen).*ainside(i:samplen,1)...
      ./denominator1(i:samplen));
  end
  mhat = (firstinside-secinside)./ainside;
  for i = 1:p
    scoreVal(i) = sum((xallord(:,i)-xbar(:,i)).*(mhat+btransx(:,1)).*deltaord);
  end
  scoreVal = scoreVal/samplen;
end