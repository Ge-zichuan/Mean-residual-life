function [hazard,MRL,hazardPrime,hazard2bt] = trueMRL(inputt,inputcov,modeltype)
% samplesize = parameter(1);
% d = parameter(2);
hazardPrime = 0;
hazard2bt = 0;
switch modeltype
    case 'cox1' % exp(1)*\bb\X\trans
    case 'cox2' % exp(1)*\bb\X\trans
    case 'sm1h' % hazard from paper SM1
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
    case 'coxt' % cox model:t*exp(\bb\X\trans), study 3
        templong = sum(exp(inputcov),2);
        hazard = templong.*inputt;
%         lambdaCumVec = sum(exp(covinput),2).*tinput.^2/2;
        MRL = exp(inputt.^2/2.*sum(exp(inputcov),2)).*normcdf(-inputt,0,sqrt(1./sum(exp(inputcov),2))).*...
            sqrt(2*pi).*sqrt(1./sum(exp(inputcov),2));
    case 'reciprocal' % 10*1/t*exp(\bb\X\trans), study 3 cont.
        const = 10;
        templong = sum(exp(inputcov),2);
        hazard = (templong.*const+1)./(inputt+1);
%         lambdaCumVec = (templong.*const+1)./log(tinput+1);
        MRL = (inputt+1)./(templong.*const);
    case 'reciprocalW' % 10*1/t*exp(\bb\X\trans), study 3 cont.
        const = 10;
        templong = sum(exp(inputcov(:,1:(end-1))+inputcov(:,end)),2);
        hazard = (templong.*const+1)./(inputt+1);
        MRL = (inputt+1)./(templong.*const);
    case 'xia1' % Model 1 from Xia's paper
    case 'xia2' % Model 2 from Xia's paper
    case 'xia3' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*inputcov).^2,2);
        tempstar = log(inputt)-5+0.1*templong;
        hazard = normpdf(tempstar,0,1)./inputt./normcdf(-tempstar,0,1);
%         lambdaCumVec = -log(normpdf(-tempstar,0,1));
        tempm = templong;
        uniqueT = unique(inputt);
        for i = 1:length(uniqueT)
            idx = find(ismember(inputt,uniqueT(i)));
            tempm(idx) = integral(@(x)funXia3(x,templong(idx)),uniqueT(i),Inf,'ArrayValued',true);
        end
        MRL = tempm./normcdf(-tempstar,0,1);
    case 'xia3W' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*inputcov(:,1:(end-1))).^2,2);
        tempstar = log(inputt)-5+0.1*templong-inputcov(:,end)/500;
        hazard = normpdf(tempstar,0,1)./inputt./normcdf(-tempstar,0,1);
%         lambdaCumVec = -log(normpdf(-tempstar,0,1));
        tempm = templong;
        uniqueT = unique(inputt);
        for i = 1:length(uniqueT)
            idx = find(ismember(inputt,uniqueT(i)));
            tempm(idx) = integral(@(x)funXia3(x,templong(idx)-inputcov(:,end)/50),uniqueT(i),Inf,'ArrayValued',true);
        end
        MRL = tempm./normcdf(-tempstar,0,1);
    case 'xia3WW' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*inputcov(:,1:(end-1))).^2,2);
        tempstar = log(inputt)-(3+inputcov(:,end)/100)+0.1*templong;
        hazard = normpdf(tempstar,0,1)./inputt./normcdf(-tempstar,0,1);
%         lambdaCumVec = -log(normpdf(-tempstar,0,1));
        tempm = templong;
        uniqueT = unique(inputt);
        for i = 1:length(uniqueT)
            idx = find(ismember(inputt,uniqueT(i)));
            tempm(idx) = integral(@(x)funXia3WW(x,templong(idx),inputcov(:,end)/100),uniqueT(i),Inf,'ArrayValued',true);
        end
        MRL = tempm./normcdf(-tempstar,0,1);
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 2
        templong = sum(exp(inputcov),2);
        hazard = 2*inputt./(templong.^2+inputt.^2);
%         lambdaCumVec = log(1+tinput.^2./templong.^2);
        MRL = (1+inputt.^2./templong.^2)./templong.*(pi/2-atan(inputt./templong));
    case 'add1' % additive model from Chen 2007
    case 'prp1' % proportional mean residual life model from Chen 2005a
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputcov)*10,2);
        hazard = 2*inputt./(templong.^2+inputt.^2);
%         lambdaCumVec = log(1+tinput.^2./templong.^2);
        MRL = (1+inputt.^2./templong.^2).*templong.*(pi/2-atan(inputt./templong));
    case 'LinLi21' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputcov),2);
        hazard = inputt.^(2/5).*templong;
%         lambdaCumVec = 5/7*tinput.^(7/5).*templong;
        uniquet = unique(inputt);
        m1 = inputt;
        for i = 1:length(uniquet)
            m1(inputt == uniquet(i)) = integral(@(x)funLinLi21(x,templong(inputt==uniquet(i))),uniquet(i),inf,'ArrayValued',true);
        end
        MRL = exp(5/7*inputt.^(7/5).*templong).*m1;
        hazardPrime = 2/5*inputt.^(-3/5).*templong;
        hazard2bt = repmat(inputt.^(2/5),[1 length(inputcov(1,:))]).*exp(inputcov);
    case 'LinLi22' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputcov),2);
        hazard = inputt.^(7/5).*templong;
%         lambdaCumVec = 5/12*tinput.^(12/5).*templong;
        uniquet = unique(inputt);
        m1 = inputt;
        for i = 1:length(uniquet)
            m1(inputt == uniquet(i)) = integral(@(x)funLinLi22(x,templong(inputt==uniquet(i))),uniquet(i),inf,'ArrayValued',true);
        end
        MRL = exp(5/12*inputt.^(12/5).*templong).*m1;
    case 'LinLi22W' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(inputcov(:,1:(end-1))),2);
        hazard = inputt.^(7/5).*templong.*inputcov(:,end);
        uniquet = unique(inputt);
        m1 = inputt;
        for i = 1:length(uniquet)
            m1(inputt == uniquet(i)) = integral(@(x)funLinLi22(x,templong(inputt==uniquet(i))).*inputcov(inputt==uniquet(i),end),uniquet(i),inf,'ArrayValued',true);
        end
        MRL = exp(5/12*inputt.^(12/5).*templong.*inputcov(:,end)).*m1;
        hazardPrime = 7/5*inputt.^(2/5).*templong.*inputcov(:,end);
        hazard2bt = repmat(inputt.^(7/5).*inputcov(:,end),[1 length(inputcov(1,:))-1]).*exp(inputcov(:,1:(end-1)));
    case 'mimicN' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(atan(inputcov)+pi/2,2);
        a = 200;c = 1/200;
        hazard = exp(c*inputt+templong)/a-c;
        hazardPrime = c*exp(c*inputt+templong)/a;
        hazard2bt = exp(c*inputt+templong)./(a.*(1+inputcov.^2));
        MRL = a*exp(-c*inputt-templong);
    case 'mimicT' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(atan(inputcov(:,1:(end-1))+10-inputcov(:,end)/5)+pi/2,2);
        a = 300;c = 1/300;
        hazard = exp(c*inputt+templong)./(a)-c;
        hazardPrime = c*exp(c*inputt+templong)./(a);
        hazard2bt = exp(c*inputt+templong)./((a).*(1+(inputcov(:,1:(end-1))-10+inputcov(:,end)/5).^2));
        MRL = a*exp(-c*inputt-templong);
%     case 'mimicT' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
%         templong = sum(atan(inputcov)+pi/2,2);
%         a = 300;c = 1/300;
%         hazard = exp(c*inputt+templong)./(a)-c;
%         hazardPrime = c*exp(c*inputt+templong)./(a);
%         hazard2bt = exp(c*inputt+templong)./((a).*(1+inputcov.^2));
%         MRL = a*exp(-c*inputt-templong);
end

end