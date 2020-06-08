function [Lval,xv] = QP_sub(G,g,A,b,x0,Lim)
A = -A';
b = -b;

n= size(A,1);
m= size(A,2);
x= x0;
lambda = ones(size(b));
s = (A'*x-b);
s(s <= 10^(-3)) = 10^(-3);
lambda(s >= 10^(-3)) = 10^(-3);


rd = G*x+g-A*lambda;
resid = -A'*x+b;
rp = s+resid;
rsl = s.*lambda;
qterm = x'*G*x;

for k = 1:20
    
    mu = sum(s.*lambda)/m;
    
    dv= lambda./s;
    D = sparse(1:m,1:m,sqrt(dv));
    rbar = A*(lambda-dv.*rp);
    
    F = sparse(A*D);
    Gbar2 = ((G+(F*F')));
    
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
    sigma = (mu_aff/mu)^4;
    
    rsln = s.*lambda+dl_aff.*ds_aff-sigma*mu;
    rbarn = A*((1./s).*(rsln-lambda.*rp));
    
    dx = Gbar2\(-rd-rbarn);
    ds = -rp+A'*dx;
    
    dl = (1./s).*(-rsln-lambda.*ds);
    
    
    Il = find(dl<0);
    Is = find(ds<0);
    if ~isempty(Il)
        al = 0.95*min(min(-lambda(Il)./dl(Il)),1);
    else
        al = 0.95;
    end
    if ~isempty(Is)
        as = 0.95*min(min(-s(Is)./ds(Is)),1);
    else
        as = 0.95;
    end
    alpha = min(al,as);
    
    x = x+alpha*dx;
    qterm = x'*G*x;
    lambda = lambda+alpha*dl;
    s = s+alpha*ds;
    resid = -A'*x+b;
    rd = G*x+g-A*lambda;
    
    rp = s+resid;
    rsl = s.*lambda;
    
    if max(s.*lambda) < 10^(-3)
        LB = 1/2*qterm+g'*x+lambda'*resid;
        if (max(resid) < 10^(-3))
            UB = 1/2*qterm+g'*x;
        else
            UB = inf;
        end
        [LB-Lim UB-Lim]
    else
        xv =[x;lambda;s];
        Lval = qterm*1/2 + x'*g;
    end
end
