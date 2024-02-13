function fw = fwfun(x,a,b,c)
% fw = x.^(a)/b.*exp(-x/c);
fw = (1-exp(-x/a)/b)*2.5;
end