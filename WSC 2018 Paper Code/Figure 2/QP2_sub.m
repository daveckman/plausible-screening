function [Lval,xv] = QP2_sub(G,g,A,b,varargin)
A = -A';
b = -b;

n= size(A,1);
m= size(A,2);

if nargin < 5
    x = zeros(n,1);
    lambda = ones(m,1);
    s = ones(m,1);
else
    if ~isempty(varargin{1})
        x= varargin{1}(1:n);
        lambda = varargin{1}((n+1):(n+m));
        s = varargin{1}((n+m+1):end);
        
        lambda(lambda < 10^(-2)) = 10^(-2);
        s(s < 10^(-2)) = 10^(-2);
    else
        x = zeros(n,1);
        lambda = ones(m,1);
        s = ones(m,1);
    end
end
if nargin > 5
    Lim = varargin{2};
else
    Lim =0;
end

rd = G*x+g-A*lambda;
resid = -A'*x+b;
rp = s+resid;

LB = -inf;
UB = inf;

for k = 1:1000
    mu = sum(s.*lambda)/m;
    
    if mu < 10^(-8)
        xv =[x;lambda;s];
        Lval = qterm*1/2 + x'*g;
        break;
    end
    
    
    dv= lambda./s;
    D = sparse(1:m,1:m,sqrt(dv));
    rbar = A*(lambda-dv.*rp);
    
    F = sparse(A*D);
    Gbar2 = (G+(F*F'));
    dx_aff = Gbar2\(-rd-rbar);
    
    qvo = A'*dx_aff;
    
    dl_aff = -lambda+dv.*rp-dv.*qvo;
    ds_aff = -rp+qvo;
    Il = find(dl_aff<0);
    Is = find(ds_aff<0);
    if ~isempty(Il)
        al_aff = min(min(-lambda(Il)./dl_aff(Il)),1);
    else
        al_aff = 1;
    end
    if ~isempty(Is)
        as_aff = min(min(-s(Is)./ds_aff(Is)),1);
    else
        as_aff =1;
    end
    alpha_aff = min(al_aff,as_aff);
    
    mu_aff =sum( (s+alpha_aff*ds_aff).*(lambda+alpha_aff*dl_aff))/m;
    sigma = (mu_aff/mu)^3;
    
    rsln = s.*lambda+dl_aff.*ds_aff-sigma*mu;
    rbarn = A*((1./s).*rsln-dv.*rp);
    dx = Gbar2\(-rd-rbarn);
    qv = A'*dx;
    dl = (1./s).*(-rsln)+dv.*rp-dv.*qv;
    ds = -rp+qv;
    
    Il = find(dl<0);
    Is = find(ds<0);
    if ~isempty(Il)
        al = 0.99*min(min(-lambda(Il)./dl(Il)),1);
    else
        al = 0.99;
    end
    if ~isempty(Is)
        as = 0.99*min(min(-s(Is)./ds(Is)),1);
    else
        as = 0.99;
    end
    alpha = min(al,as);
    
    x = x+alpha*dx;
    qterm = x'*G*x;
    lambda = lambda+alpha*dl;
    s = s+alpha*ds;
    resid = resid-alpha*qv;
    rd = G*x+g-A*lambda;
    rp = s+resid;
    
    if max(abs(rd)) < 10^(-8) || max(abs(rp))< 10^(-8)
        xv =[x;lambda;s];
        Lval = qterm*1/2 + x'*g;
        break;
    end
    
    
    if nargin > 5
        if max(s.*lambda) < 10^(-6)
            LB = 1/2*qterm+g'*x+lambda'*resid;
            if (max(resid) < 10^(-6))
                UB = 1/2*qterm+g'*x;
            else
                UB = inf;
            end
            
            if LB-Lim > 0
                Lval = LB;
                xv =[x;lambda;s];
                break
            elseif UB-Lim < 0
                Lval = UB;
                xv =[x;lambda;s];
                break
            else
                xv =[x;lambda;s];
                Lval = qterm*1/2 + x'*g;
            end
        else
            xv =[x;lambda;s];
            Lval = qterm*1/2 + x'*g;
        end
    else
        qterm = x'*G*x;
        
        xv =[x;lambda;s];
        Lval = qterm*1/2 + x'*g;
    end
end
xv =[x;lambda;s];
qterm = x'*G*x;
Lval = qterm*1/2 + x'*g;
end
