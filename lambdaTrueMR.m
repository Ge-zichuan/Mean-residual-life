function [lambdaVec,lambdaCumVec,m] = lambdaTrueMR(inputt,inputg,ttype)
switch ttype
    case 'cox1' % exp(1)*\bb\X\trans
    case 'cox2' % exp(1)*\bb\X\trans
    case 'sm1h' % hazard from paper SM1
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
    case 'coxt' % cox model:t*exp(\bb\X\trans), study 3
        templong = sum(exp(inputg),2);
        lambdaVec = templong.*inputt;
        lambdaCumVec = sum(exp(inputg),2).*inputt.^2/2;
        m = exp(inputt.^2/2.*sum(exp(inputg),2)).*normcdf(-inputt,0,sqrt(1./sum(exp(inputg),2))).*...
            sqrt(2*pi).*sqrt(1/sum(exp(inputg),2));
    case 'reciprocal' % 10*1/t*exp(\bb\X\trans), study 3 cont.
        const = 10;
        templong = sum(exp(inputg),2);
        lambdaVec = (templong.*const+1)./(inputt+1);
        lambdaCumVec = (templong.*const+1)./log(inputt+1);
        m = (inputt+1)./(templong.*const);
    case 'xia1' % Model 1 from Xia's paper
    case 'xia2' % Model 2 from Xia's paper
    case 'xia3' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*inputg).^2,2);
        tempstar = log(inputt)-5+0.1*templong;
        lambdaVec = normpdf(tempstar,0,1)./inputt./normcdf(-tempstar,0,1);
        lambdaCumVec = -log(normpdf(-tempstar,0,1));
        tempm = templong;
        uniqueT = unique(inputt);
        for i = 1:length(uniqueT)
            idx = find(ismember(inputt,uniqueT(i)));
            tempm(idx) = integral(@(x)funXia3(x,templong(idx)),uniqueT(i),Inf,'ArrayValued',true);
        end
        m = tempm./normcdf(-tempstar,0,1);
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 2
        templong = sum(exp(inputg),2);
        lambdaVec = 2*inputt./(templong.^2+inputt.^2);
        lambdaCumVec = log(1+inputt.^2./templong.^2);
        m = (1+inputt.^2./templong.^2)./templong.*(pi/2-atan(inputt./templong));
    case 'add1' % additive model from Chen 2007
    case 'prp1' % proportional mean residual life model from Chen 2005a
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputg)*10,2);
        lambdaVec = 2*inputt./(templong.^2+inputt.^2);
        lambdaCumVec = log(1+inputt.^2./templong.^2);
        m = (1+inputt.^2./templong.^2).*templong.*(pi/2-atan(inputt./templong));
    case 'LinLi21' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputg),2);
        lambdaVec = inputt.^(2/5).*templong;
        lambdaCumVec = 5/7*inputt.^(7/5).*templong;
        uniquet = unique(inputt);
        m1 = inputt;
        for i = 1:length(uniquet)
            m1(inputt == uniquet(i)) = integral(@(x)funLinLi21(x,templong(inputt==uniquet(i))),uniquet(i),inf,'ArrayValued',true);
        end
        m = exp(5/7*inputt.^(7/5).*templong).*m1;
    case 'LinLi22' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputg),2);
        lambdaVec = inputt.^(7/5).*templong;
        lambdaCumVec = 5/12*inputt.^(12/5).*templong;
        uniquet = unique(inputt);
        m1 = inputt;
        for i = 1:length(uniquet)
            m1(inputt == uniquet(i)) = integral(@(x)funLinLi22(x,templong(inputt==uniquet(i))),uniquet(i),inf,'ArrayValued',true);
        end
        m = exp(5/12*inputt.^(12/5).*templong).*m1;
end

end