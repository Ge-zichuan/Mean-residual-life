function y = tFromHazard(x,scenario,u,a,c,bx,fw)
switch scenario
    case 'N'
        %y = -log(u)-(1/a-a/c)*x+a*log(1+exp(-x/c+bx))-c/a*exp(x/c-bx)+c/a*exp(-bx)-a*log(1+exp(bx));
        y = log(u)-c*x-exp(bx)/(a*c)+exp(c*x+bx)./(a*c);
%     case 'T'
% %         y = -log(u)-(1./a-a./c).*x+a.*log(1+exp(-x/c+bx))-c./a.*exp(x/c-bx)+c./a.*exp(-bx)-a.*log(1+exp(bx));
%         y = log(u)-c.*x./fw-exp(bx)./(a.*c)+exp(c.*x./fw+bx)./(a.*c);
    case 'T'
        y = log(u)-c.*x-exp(bx)./(a.*c)+exp(c.*x+bx)./(a.*c);
end
end