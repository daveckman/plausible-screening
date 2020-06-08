function [Lval,lambda] = LP_sub(o,A,b,x0,lambda0,CUTOFF)
A = -A';
b = -b;

n= size(A,1);
m= size(A,2);
x= x0;
lambda = lambda0;
s = (A'*x-b);
s(s <= 10^(-4)) = 10^(-4);
lambda(lambda <= 10^(-4)) = 10^(-4);
lambda(s>10) = 10^(-4);
%lambda(s<10^(-1)) = mean(lambda > 10^(-3));

UB_sub = inf;
UB = inf;
LB = -inf;
Lval = nan;

G = zeros(length(o));
g = o;

dv = lambda./s;
mu = sum(dv)/m;
% [~,REE]= sort(-lambda./s);
% if mu^(1/4)*m <= n
%     QEE = n;
% elseif mu^(1/4)*m < m
%     QEE = ceil(mu^(1/4)*m);
% else
%     QEE = m;
% end
%QEE = m;
IVS =1:m;


Ar = A(:,IVS);

rd = G*x+g-Ar*lambda(IVS);
resid = -A'*x+b;
rp = s+resid;
rsl = s.*lambda;


max_nu = nan;
min_nu = nan;
LB = -inf;

for k = 1:100
    rd = G*x+g-Ar*lambda(IVS);
    rp = s+resid;
    rsl = s.*lambda;
    
    rbar = Ar*(lambda(IVS)-dv(IVS).*rp(IVS));
    D = sparse(1:length(IVS),1:length(IVS),sqrt(dv(IVS)));
    F = sparse(Ar*D);
    Gbar2 = G+(F*F');
    
    
    
    dx_aff = Gbar2\(-rd-rbar);
    ds_aff = -rp+A'*dx_aff;
    dl_aff = (1./s).*(-rsl-lambda.*ds_aff);
    
    Il = find(dl_aff<0);
    Is = find(ds_aff<0);
    if ~isempty(Il)
        al_aff = 1*min(min(-lambda(Il)./dl_aff(Il)),1);
    else
        al_aff =1;
    end
    if ~isempty(Is)
        as_aff = 1*min(min(-s(Is)./ds_aff(Is)),1);
    else
        as_aff = 1;
    end
    alpha_aff = min(al_aff,as_aff);
    
    mu_aff =sum( (s+alpha_aff*ds_aff).*(lambda+alpha_aff*dl_aff))/m;
    sigma = (mu_aff/mu)^3;
    
    rsln = s.*lambda+dl_aff.*ds_aff-sigma*mu;
    rbarn = Ar*((1./s(IVS)).*(rsln(IVS)-lambda(IVS).*rp(IVS)));
    
    dx = Gbar2\(-rd-rbarn);
    Atdx = A'*dx;
    ds = -rp+Atdx;
    
    dl = (1./s).*(-rsln-lambda.*ds);
    
    
    Il = find(dl<0);
    Is = find(ds<0);
    if ~isempty(Il)
        al = 0.98*min(min(-lambda(Il)./dl(Il)),1);
    else
        al = 0.98;
    end
    if ~isempty(Is)
        as = 0.98*min(min(-s(Is)./ds(Is)),1);
    else
        as = 0.98;
    end
    alpha = min(al,as);
    
    x = x+alpha*dx;
    lambda = lambda+alpha*dl;
    s = s+alpha*ds;
    dv = lambda./s;
    mu = sum(s.*lambda)/m;
    
    IVS =1:m;
    Ar = A(:,IVS);
    
    resid = resid-alpha*Atdx;
    qterm = x'*G*x;
    rd = G*x+g-Ar*lambda(IVS);
    rp = s+resid;
    rsl = s.*lambda;
    
    if max(abs(rd)) < 10^(-2)
        LB = max( g'*x+lambda'*resid,LB);
    end
    
    if (max(resid) < 10^(-2))
        UB = g'*x;
    end
    
    if LB > CUTOFF
        Lval = LB;
        break;
    end
    if UB < CUTOFF
        Lval = UB;
        break;
    end
    if abs(UB-LB) < 10^(-3)
        Lval = mean([UB,LB]);
        break;
    end
end
if LB < CUTOFF
    Lval = inf;
end

end

