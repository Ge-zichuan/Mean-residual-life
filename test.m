samplesize = 2000;
x = 0:500;
w = 300;
fw = ((2*w/10-5)./((w/10).^2-w/10+5)+1);
plot(w,fw);
bx = -3;
yN = exp(0.04*x.^2./(2*exp(bx))).*normcdf(-0.2*x/sqrt(exp(bx)),0,1)*sqrt(2*pi*exp(bx)/0.04);
yT = fw*exp(0.01*x.^2./(2*exp(bx))).*normcdf(-0.1*x/sqrt(exp(bx)),0,1)*sqrt(2*pi*exp(bx)/0.01);
% hold on;
% plot(x,yN,'color','r');
% plot(x,yT,'color','b');
% hold off
w = 0:300;
fw = fwfun(w,10,1,100);
plot(w,fw);
wAll = [0.1,1,10,50,100,200];
for i = 1:6
    w = wAll(i);
    fw = fwfun(w,10,1,100);%w*100/(1+w*100);%1/(w/10+1);%
    a = 200;c = 1/200;
    yN = a*exp(-x*c-(atan(bx)+pi/2));
    u = rand(samplesize,1);
    t = fsolve(@(x)tFromHazard(x,'N',u,a,c,atan(bx)+pi/2,999),zeros(samplesize,1));
    hist(t);
    s = waitforbuttonpress;
    u = rand(samplesize,1);
%     y = -log(u)-(1/a-a/c)*x+a*log(1+exp(-x/c+bx))-c/a*exp(x/c-bx)+c/a*exp(-bx)-a*log(1+exp(bx))
    a = 300;c = 1/300;
    t = fsolve(@(x)tFromHazard(x,'T',u,a,c,atan(bx+10-w/5)+pi/2,fw),zeros(samplesize,1));
    hist(t);
    s = waitforbuttonpress;
    yT = a*exp(-x*c-(atan(bx+10-w/5)+pi/2));
    plot(x,yN);hold on;plot(x,yT);hold off;
    s = waitforbuttonpress;
    ydiff = yT(1:(500-w))-yN((w+1):500);
    plot(x(1:(500-w)),ydiff);title(num2str(w));
    s = waitforbuttonpress;
end


bx = 5;a = 200;c = 1/200;
bx = -100:100;
plot(bx,atan(bx/10));
t = 0:0.1:100;
ft = -c*t-exp(bx)/(a*c)+exp(c*t+bx)./(a*c);
plot(t,ft);title(min(ft));