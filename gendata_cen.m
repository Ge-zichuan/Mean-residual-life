function [cout,deltaout,kval,randx] = gendata_cen(tinput,xinput,ctype,censorRate,n)
kval = 20.0;
adjrate = 0.01;
% ttype----
%     xia1: Sensoring 1 from Xia's paper, CDF of N(0,1)+r.v.;
%     xia2: Sensoring 2 from Xia's paper, \|x\|^2+uniform;
%     unif: uniform
censorRate = round((1-censorRate)*n);
switch ctype
    case 'xia1'
        ctemp = normcdf(2*xinput(:,2)+2*xinput(:,3));
        randx = rand(n,1);
        deltaout = tinput<=ctemp;
        countt2 = 0;
        while (~(sum(deltaout)==(censorRate)) && countt2<=100000)
            if (sum(deltaout)<censorRate)
                kval = kval+adjrate;
            elseif (sum(deltaout)>censorRate)
                kval = kval-adjrate;
            end
            cout = ctemp+kval*randx;
            deltaout = tinput<=cout;
            countt2 = countt2 + 1;
        end
    case  'xia2'
        if (censorRate > n-1)
            deltaout = ones(n,1);
            cout = ones(n,1)*maxval(tinput);
        else
            temp = median(tinput);
            randx = rand(n,1);
            ctemp = temp*ones(n,1);
            deltaout = tinput<=ctemp;
            cout = ctemp;
            countt2 = 0;
            while (~(sum(deltaout)==(censorRate)) && countt2<=100000)
                if (sum(deltaout)<censorRate)
                    kval = kval+adjrate;
                elseif (sum(deltaout)>censorRate)
                    kval = kval-adjrate;
                end
                cout = ctemp+kval*randx;
                deltaout = tinput<=cout;
                countt2 = countt2 + 1;
            end
        end
    case 'gamm'
        cout = gamrnd(25,kval,n,1);
        randx = zeros(n,1);
        deltaout = tinput<=cout;
        countt2 = 0;
        while (~(sum(deltaout)==(censorRate)) && countt2<=100000)
            if (sum(deltaout)<censorRate)
                kval = kval+adjrate;
            elseif (sum(deltaout)>censorRate)
                kval = kval-adjrate;
            end
            cout = gamrnd(15,kval,n,1);
            deltaout = tinput<=cout;
            countt2 = countt2 + 1;
        end
    case 'expl'
        cout = exprnd(kval,n,1);
        deltaout = tinput<=cout;
        countt2 = 0;
        while (~(sum(deltaout)==(censorRate)) && countt2<=100000)
            if (sum(deltaout)<censorRate)
                kval = kval+adjrate;
            elseif (sum(deltaout)>censorRate)
                kval = kval-adjrate;
            end
            cout = exprnd(kval,n,1);
            deltaout = tinput<=cout;
            countt2 = countt2 + 1;
        end
    case 'unif'
        kval = max(tinput);
        randx = rand(n,1);
        cout = randx*kval;
        deltaout = tinput<=cout;
        countt2 = 0;
        while (~(sum(deltaout)==(censorRate)) && countt2<=100000)
            if (sum(deltaout)<censorRate)
                kval = kval+adjrate;
            elseif (sum(deltaout)>censorRate)
                kval = kval-adjrate;
            end
            cout = randx*kval;
            deltaout = tinput<=cout;
            countt2 = countt2 + 1;
        end
end
end