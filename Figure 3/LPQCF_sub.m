function [Lval,x,nu,lambda] = LPQCF_sub(o,A,b,G0,g0,c,x0,nu0,lambda0)
A = -A';
b = -b;

n= size(A,1);
m= size(A,2);
x= x0;
nu = nu0;
lambda = lambda0;
s = (A'*x-b);
s(s <= 10^(-1)) = 10^(-1);
lambda(lambda <= 10^(-1)) = 10^(-1);

LB_sub = inf;
UB_sub = inf;
UB = inf;
LB = -inf;
Lval = nan;

G = nu0*G0;
g = nu0*g0+o;

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

qterm = x'*G*x;

max_nu = nan;
min_nu = nan;
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
    
    
    if max(abs(rd))< 10^(-3)
        LB = 1/2*qterm+g'*x+lambda'*resid+nu*c;
    end
    
    if (max(resid) < 10^(-3))
        UB_sub = 1/2*qterm+g'*x+nu*c;
    end
    
    
    if UB_sub-LB < 10^(-3)
        if  (max(resid) < 10^(-3)) && (qterm/2/nu+g0'*x + c < 10^(-1))
            UB = o'*x;
        end
        
        if UB-LB < 10^(-3)
            break;
        end
        
        if UB < -10^9
            LB = -inf;
            break;
        end
        
        if x'*G0*x/2+g0'*x + c < 0
            max_nu = nu;
            if isnan(min_nu)
                nu = nu/2;
            else
                nu = exp((log(max_nu)+log(min_nu))/2);
            end
        else
            min_nu = nu;
            if isnan(max_nu)
                nu = nu*2;
            else
                nu = exp((log(max_nu)+log(min_nu))/2);
            end
        end
        
        
        
        G = nu*G0;
        g = nu*g0+o;
        qterm = x'*G*x;
        rd = G*x+g-A*lambda;
        UB_sub = inf;
    end
    
end
if k >200
[LB UB]
ada
end
if ~isnan(LB)
    Lval = LB;
else
    Lval = -inf;
end
