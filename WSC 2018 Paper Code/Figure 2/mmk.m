function profit = mmk(c_val)
n = length(c_val);
profit  = zeros(1,n);

lambda = 1;
mu = 1000; %muc(q);
lcv = 1;
for q = 1:n
%     if (q/n) >= 0.1*lcv
%         fprintf('\n %d %% complete',round(100*q/n));
%         lcv = lcv+1;
%     end
    c = c_val(q);
    rho = mu*lambda/c;
    
    SS = 1:5000;
    w = zeros(size(SS));
    w(SS<c) = SS(SS<c)*log(c)+SS(SS<c)*log(rho)-cumsum(log(min(SS):(c-1)));
    w(SS>=c) = c*log(c)+SS(SS>=c)*log(rho)-sum(log(min(SS):c));
    w = w - max(w)+8;
    piv = exp(w);
    piv = piv/sum(piv);
    
    Uv = cumsum(piv);
    num = min(SS(rand()<Uv));
    U = exprnd(mu,[c,1]);   
    
    
    if num <= c
        U(1:(c-num)) = 0;
    else
        for k = 1:(num-c)
            [~,j] = min(U);
            U(j) = U(j)+exprnd(mu);
        end
    end
    
    T = 0;
    W = zeros(10000,1);
    i = 1;
    while T < 10000
        wv = exprnd(1);
        [~,j] = min(U);
        Up = U-wv;
        S = exprnd(mu);
        W(i) = U(j)+S;
        i = i + 1;
        Up(j) = Up(j)+S;
        U = Up.*(Up>=0);
        T = T + wv;
    end
    profit(q) =  0.01*c+mean(sqrt(W(1:(i-1))));
end