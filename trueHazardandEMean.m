function [Hazard,CumulHazard,EMean] = trueHazardandEMean(tinput,covinput,parameter,modeltype)
samplesize = parameter(1);
d = parameter(2);
switch modeltype
    case 'cox1' % exp(1)*\bb\X\trans
        Hazard = -log(rand(n,1)).*sum(exp(xinput*parainput),2);
    case 'cox2' % exp(1)*\bb\X\trans
        Hazard = -log(rand(n,1)).*sum(exp(xinput*parainput),2);
    case 'sm1h' % hazard from paper SM1
        r1 = 0.5;
        templong = sum(exp(xinput*parainput),2);
        Hazard = templong.*exp(log(((1+r1)*(rand(n,1)).^(-templong*r1)-1)/r1)./templong);
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
        templong = sum(exp(xinput*parainput),2);
        c = 1/(1-exp(-100*templong));
        Hazard = (-log(1-rand(n,1)./c)./templong);
        while (sum(Hazard>100) > 0)
            Hazard(Hazard>100) = -log(1-rand(sum(Hazard>100),1)./c)./templong(Hazard>100);
        end
    case 'coxt' % cox model:t*exp(\bb\X\trans, study 3
        Hazard = repmat(tinput,[1 ones(1,d)*samplesize]).*repmat(exp(covinput'),[samplesize 1]);
        CumulHazard = repmat(tinput.^2,[1 ones(1,d)*samplesize]).*repmat(exp(covinput'),[samplesize 1]);
        EMean = sqrt(pi./exp(covinput)).*(0.5-normcdf(tauinput,0,sqrt(1./2./exp(covinput))));
    case 'xia1' % Model 1 from Xia's paper
        templong = sum(exp(6*xinput*parainput),2);
        Hazard = normcdf(5*(-log(rand(n,1)).*(templong+1)-2));
    case 'xia2' % Model 2 from Xia's paper
        templong = sum((xinput*parainput/2),2);
        Hazard = exp(-log(rand(n,1))+templong);
    case 'xia3' % Model 3 from Xia's paper
        Hazard = exp(5-10*templong+normrnd(0,1,n,1));
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 2
        templong = sum(exp(xinput*parainput),2);
        Hazard = templong.*(1./rand(n,1)-1.0).^(1/8);
        while (sum(Hazard>20) > 0)
            Hazard(Hazard>20) = templong(Hazard>20).*(1./rand(sum(Hazard>20),1)-1).^(1/8);
        end
    case 'add1' % additive model from Chen 2007
        templong = sum((xinput*parainput),2);
        Hazard = (1+templong)./sqrt(rand(n,1))-1-templong;
    case 'prp1' % proportional mean residual life model from Chen 2005a
        templong = exp(sum(xinput,2));
        Hazard = (rand(n,1)).^(-templong./(1+templong))-1;
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(xinput*parainput),2);
        Hazard = sqrt((1./templong+1).*(1./rand(n,1)-1));
    case 'mimicN' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(xinput*parainput),2);
        Hazard = sqrt((1./templong+1).*(1./rand(n,1)-1));
    case 'mimicT' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(xinput*parainput),2);
        Hazard = sqrt((1./templong+1).*(1./rand(n,1)-1));
end

end