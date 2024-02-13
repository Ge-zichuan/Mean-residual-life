function scoreVal = scoreBetaProp2(paras1,xall,zall,deltaall,samplen,p)
% paras: parameter \beta in equation (5)
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1

  % khbetax: K_h(\bb\trans\x)
  % denoVev: \sumi I()K_h()
	scoreVal = 1:p;
  xbar = zeros(samplen,p);
  intex = ones(samplen*10,1);
  xbinside = ones(samplen,1);
  firstinside = ones(samplen,p);
  
  [zallord,sortIndex] = sort(zall);
  xallord = xall(sortIndex,:);
  deltaord = deltaall(sortIndex);
  % calculate \bb\trans\x
  parasAll = paras1;
  btransx = xallord*parasAll; 
  % calculate \wh m(t) first at each z_i
  mdiffall = (zallord*ones(1,samplen))<=(ones(samplen,1)*zallord');

  % calculate all z_i-z_j
  finside = 1-(1-deltaord)./sum(mdiffall,2);
  fhat = cumprod(finside);
  if (fhat(samplen)==0) 
    fhat(samplen) = fhat(samplen-1);
  end

  for i = 1:samplen
    xbinside(i) = sum(exp(-btransx(:,1)).*deltaord.*zallord./fhat)/samplen;
  end

  for i = 1:samplen
    for j = 1:(samplen*10)
      intex(j) = zallord(samplen)/(samplen*10)/2+(j-1)*zallord(samplen)/(samplen*10);
    end
    integration1 = inteProp2(paras1,samplen*10,samplen,intex,zallord,xallord,deltaord,fhat,zallord(i));
    firstinside(i) = sum(integration1)*zallord(samplen)/(samplen*10);
  end
  for j = 1:p
      xbar(:,j) = xallord(:,j).*exp(-2*btransx(:,1)).*firstinside(:,j)./xbinside;
  end
  for i = 1:p
    scoreVal(i) = sum(deltaord.*(xallord(:,i)-xbar(:,i))./fhat);
  end
  scoreVal = scoreVal/(samplen);
end