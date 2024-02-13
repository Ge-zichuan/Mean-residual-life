function outputy = inteProp2(paras1,n,samplen,inputx,zallord,xallord,deltaord,fhat,tend)
  numerator = zeros(n,1);
  denominator1 = zeros(n,1);
  diffplus = zeros(samplen,1);
  indplus = zeros(samplen,1);
  tendarray = zeros(n,1);
  btransx = xallord*paras1; 
  for i = 1:n
    for j = 1:samplen
      if (zallord(j)>=inputx(i))
        diffplus(j) = zallord(j)-inputx(i);
        indplus(j) = 1;
      end
    end
    numerator(i) = sum(exp(-btransx(:,1)).*deltaord./fhat.*indplus);
    denominator1(i) = sum(exp(-2*btransx(:,1)).*deltaord./fhat.*diffplus);
    if (tend>=inputx(i))
      tendarray(i) = tend-inputx(i);
    end
  end
  outputy = tendarray.*numerator./denominator1;
end
