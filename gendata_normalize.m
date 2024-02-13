function [xout,c,s] = gendata_normalize(x,normaltype)
% normalize data X (columnwise)
switch normaltype
    case 'mimic'
        xout = x;
        xout(:,3) = x(:,3)/max(x(:,3));
        xout(:,5) = x(:,5)/max(x(:,5));
        xout(:,6) = x(:,6)/max(x(:,6));
        xout(:,7) = x(:,7)/max(x(:,7));
        c = zeros(1,length(x(1,:)));
        s = ones(1,length(x(1,:)));
        s([3 5 6 7]) = max(x(:,[3 5 6 7]));
end
end