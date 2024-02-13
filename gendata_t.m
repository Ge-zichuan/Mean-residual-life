function [tout,lambdaout,templong] = gendata_t(xinput,n,ttype,parainput)
% generate data T, n: sample size, xdim: dimension of each observation
% ttype---
%     cox1: ~exp(1)*\bb\X\trans;
%     sm1h: hazard from paper SM1;
%     texp: truncated exp(1)*\bb\X\trans at 100;
%     coxt: cox model with t*exp(\bb\X\trans);
%     xia1: Model 1 from Xia's paper;
%     xia2: Model 2 from Xia's paper;
%     xia3: Model 3 from Xia's paper;
switch ttype
    case 'cox1' % exp(1)*\bb\X\trans
        tout = -log(rand(n,1)).*sum(exp(xinput*parainput),2);
    case 'cox2' % exp(1)*\bb\X\trans
        tout = -log(rand(n,1)).*sum(exp(xinput*parainput),2);
    case 'sm1h' % hazard from paper SM1
        r1 = 0.5;
        templong = sum(exp(xinput*parainput),2);
        tout = templong.*exp(log(((1+r1)*(rand(n,1)).^(-templong*r1)-1)/r1)./templong);
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
        templong = sum(exp(xinput*parainput),2);
        c = 1/(1-exp(-100*templong));
        tout = (-log(rand(n,1)./c)./templong);
        while (sum(tout>100) > 0)
            tout(tout>100) = -log(1-rand(sum(tout>100),1)./c)./templong(tout>100);
        end
    case 'coxt' % cox model:t*exp(\bb\X\trans)
        templong = sum(exp(xinput*parainput),2);
        tout = sqrt(-2*log(rand(n,1))./templong);
        lambdaout = templong.*tout;
    case 'reciprocal' % 3*1/(t+1)*exp(\bb\X\trans)
        const = 10;
        templong = sum(exp(xinput*parainput),2);
        tout = rand(n,1).^(-1./(const*templong+1))-1;
        lambdaout = (const*templong+1)./(tout+1);
    case 'reciprocalW' % 3*1/(t+1)*exp(\bb\X\trans)
        const = 10;
        templong = sum(exp(xinput(:,1:(end-1))*parainput+xinput(:,end)),2);
        tout = rand(n,1).^(-1./(const*templong+1))-1;
        lambdaout = (const*templong+1)./(tout+1);
    case 'xia1' % Model 1 from Xia's paper
        templong = sum(exp(xinput*parainput),2);
        tout = normcdf(5*(-log(rand(n,1)).*(templong+1)-2));
        lambdaout = 1./(5*templong.*normpdf(norminv(tout,0,1),0,1));
    case 'xia2' % Model 2 from Xia's paper
        templong = sum((xinput*parainput/2),2);
        tout = exp(-log(rand(n,1))+templong);
    case 'xia3' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*xinput*parainput).^2,2);
        tout = exp(5-0.1*templong+normrnd(0,1,n,1));
        lambdaout = normpdf(log(tout)-5+0.1*templong,0,1)./tout./(1-normcdf(log(tout)-5+0.1*templong,0,1));
    case 'xia3W' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*xinput(:,1:(end-1))*parainput).^2,2);
        tout = exp(5-0.1*templong+xinput(:,end)/500+normrnd(0,1,n,1));
        lambdaout = normpdf(log(tout)-5+0.1*templong-xinput(:,end)/500,0,1)./tout./(1-normcdf(log(tout)-5+0.1*templong-xinput(:,end)/500,0,1));
    case 'xia3WW' % Model 3 from Xia's paper
        templong = sum((1-sqrt(2)*xinput(:,1:(end-1))*parainput).^2,2);
        tout = exp((3+xinput(:,end)/100)-0.1*templong+normrnd(0,1,n,1));
        lambdaout = normpdf(log(tout)-(3+xinput(:,end)/100)+0.1*templong,0,1)./tout./(1-normcdf(log(tout)-(3+xinput(:,end)/100)+0.1*templong,0,1));
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 2
        templong = sum(exp(xinput*parainput),2);
        parabeta = 25;
        tout = templong.*(1./rand(n,1)-1.0).^(1/parabeta);
        %         while (sum(tout>20) > 0)
        %             tout(tout>20) = templong(tout>20).*(1./rand(sum(tout>20),1)-1).^(1/8);
        %         end
        lambdaout = parabeta*tout.^(parabeta-1)./templong.^parabeta./(1+(tout./templong).^parabeta);
    case 'add1' % additive model from Chen 2007
        templong = sum((xinput*parainput),2);
        tout = (1+templong)./sqrt(rand(n,1))-1-templong;
    case 'prp1' % proportional mean residual life model from Chen 2005a
        templong = exp(sum(xinput,2));
        tout = (rand(n,1)).^(-templong./(1+templong))-1;
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x),2), study 1
        templong = sum(exp(xinput*parainput)*10,2);
        tout = (templong).*sqrt((1./rand(n,1)-1));
        lambdaout = 2*tout./(templong.^2+tout.^2);
    case 'LinLi1' % simulation 1 from Lin and Li 2016
        templong = sum(exp(xinput*parainput),2);
        tout = (templong).*sqrt((1./rand(n,1)-1));
        lambdaout = 2*tout./(templong.^2+tout.^2);
    case 'LinLi21' % simulation 2 from Lin and Li 2016
        templong = sum(exp(xinput*parainput),2);
        tout = (-7/5*log(rand(n,1))./templong).^(5/7);
        lambdaout = tout.^(2/5).*templong;
    case 'LinLi22' % simulation 2 from Lin and Li 2016
        templong = sum(exp(xinput*parainput),2);
        tout = (-12/5*log(rand(n,1))./templong).^(5/12);
        lambdaout = tout.^(7/5).*templong;
    case 'LinLi22W' % simulation 2 from Lin and Li 2016
        templong = sum(exp(xinput(:,1:(end-1))*parainput),2);
        tout = (-12/5*log(rand(n,1))./(templong.*xinput(:,end))).^(5/12);
        lambdaout = tout.^(7/5).*templong.*xinput(:,end);
    case 'mimicN' % mimic real data, N group
        templong = sum(atan(xinput*parainput)+pi/2,2);
        a = 200;c = 1/200;
        options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt');
        check = 0;
        while check == 0
            u = rand(n,1);
            tout = fsolve(@(x)tFromHazard(x,'N',u,a,c,templong,999),ones(n,1),options);
            if (sum(tout<=0)==0)
                check = 1;
            end
        end
        lambdaout = exp(c*tout+templong)/a-c;
    case 'mimicT' % mimic real data, N group
        templong = sum(atan(xinput(:,1:(end-1))*parainput+10-xinput(:,end)/5)+pi/2,2);
        a = 300;c = 1/300;
        options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt');
        check = 0;
        while check == 0
            u = rand(n,1);
            tout = fsolve(@(x)tFromHazard(x,'T',u,a,c,templong,1),ones(n,1),options);
            if (sum(tout<=0)==0)
                check = 1;
            end
        end
        lambdaout = exp(c*tout+templong)/a-c;
%     case 'mimicT' % mimic real data, N group
%         templong = sum(atan(xinput(:,1:(end-1))*parainput)+pi/2,2);
%         a = 300;c = 1/300;
%         fw = fwfun(xinput(:,end),10,1,999);
%         options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt');
%         check = 0;
%         while check == 0
%             u = rand(n,1);
%             tout = fsolve(@(x)tFromHazard(x,'T',u,a,c,templong,fw),ones(n,1),options);
%             if (sum(tout<=0)==0)
%                 check = 1;
%             end
%         end
%         lambdaout = exp(c*tout./fw+templong)/a./fw-c./fw;
end
end