function lambdaVec = lambdaD1TrueMR(inputt,inputg,d,ttype)
lambdaVec = zeros(length(inputt),d);
switch ttype
    case 'cox1' % exp(1)*\bb\X\trans
    case 'cox2' % exp(1)*\bb\X\trans
    case 'sm1h' % hazard from paper SM1
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
    case 'coxt' % cox model:t*exp(\bb\X\trans), study 3
        lambdaVec = exp(inputg).*repmat(inputt,1,d);
    case 'reciprocal' % 10*1/t*exp(\bb\X\trans), study 3 cont.
        templong = sum(exp(inputg),2);
        for i = 1:d
            lambdaVec(:,i) = templong.*10./(inputt+1);
        end
    case 'xia1' % Model 1 from Xia's paper
    case 'xia2' % Model 2 from Xia's paper
    case 'xia3' % Model 3 from Xia's paper
        tempstar = log(inputt)-5+0.1*sum((1-sqrt(2)*inputg).^2,2);
        for i = 1:d
            lambdaVec(:,i) = -sqrt(2)*2*0.1*(1-sqrt(2)*inputg(:,i)).*(-tempstar.*normpdf(tempstar,0,1).*normcdf(-tempstar,0,1)...
                ./inputt+normpdf(tempstar,0,1).^2./inputt)./(normcdf(-tempstar,0,1).^2);
        end
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 2
    case 'add1' % additive model from Chen 2007
    case 'prp1' % proportional mean residual life model from Chen 2005a
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputg)*10,2);
        lambdaVec = repmat(-4*templong.*inputt./(templong.^2+inputt.^2).^2,1,d).*(exp(inputg)*10);
    case 'LinLi21' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        lambdaVec = exp(inputg).*repmat(inputt.^(2/5),1,d);
    case 'LinLi22' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        lambdaVec = exp(inputg).*repmat(inputt.^(7/5),1,d);
end

end