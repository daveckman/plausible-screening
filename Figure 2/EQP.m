
function [d,lambda]= EQP(A,b,Am,bm,xh)
p = size(Am,2);
if ~isempty(bm)
    xht = [A Am;
        Am' 10^(-12)*eye(p)]\[-b;bm];
    d  =xht(1:end-p);
    lambda = xht(end-p+1:end);
else
    xht = A\[-A*xh+b];
    d  =xht(1:end);
    lambda = [];
end
end