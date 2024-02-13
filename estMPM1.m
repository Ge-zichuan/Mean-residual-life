function mhatfinal = estMPM1(paras1,inputt,inputg,xall,zall,deltaall,samplesize)
% method from Chen and Chen 2005
% paras: parameter \beta in equation (5)
% xinput, deltaall, zinput: matrix X, \Delat vector, Z vector in equation (5)
% Y is indicator function which calculated in this function
% samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1

% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
inputn = length(inputt);
parasall = [1;paras1];

[zallord,sortIndex] = sort(zall);
xallord = xall(sortIndex,:);
deltaord = deltaall(sortIndex);
% calculate \bb\trans\x
btransx = xallord*parasall; % paras1 include upper d-by-d block.
% calculate all z_i-z_j
mdiffall = (inputt*ones(1,samplesize))<=(ones(inputn,1)*zallord');% T_i<=T_j
denominator = sum(mdiffall,2);
% calculate estimating equation, each term.
SNAinside = exp(-(mdiffall*deltaord)./denominator);
Qinside = mdiffall*exp(-btransx(:,1))./denominator;
SNAinside(isnan(SNAinside)) = 1;
Qinside(isnan(Qinside)) = 0;
% calculate \wh m(t) first at each z_i
innerTemp = SNAinside.*Qinside;
integral2 = sum((innerTemp(1:(end-1))+innerTemp(2:end)).*diff(inputt)/2)-...
    [0 ; cumsum((innerTemp(1:(end-1))+innerTemp(2:end))/2.*diff(inputt))];
mhat = integral2./SNAinside;
% mhat(isnan(mhat)) = 0;
matTemp = zeros(inputn*length(inputg),length(inputg));
for i = 1:length(inputg)
    matTemp(((i-1)*inputn+1):(i*inputn),i) = 1;
end
mhatfinal = repmat(mhat,length(inputg),1).*(matTemp*exp(inputg));

end