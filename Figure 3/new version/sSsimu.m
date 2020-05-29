function y = sSsimu(x)

n = size(x,1);
y  = zeros(1,n);

for q = 1:n
    s = x(q,1);
    S = x(q,2);
    I = S;
    cost = 0;
    for k = 1:100
        J = I;
        D = poissrnd(25);
        if J <= s
            cost = cost + 32 + 3*(S-J);
            J = S;
        end
        if J >= D
            cost = cost + (J-D);
        else
            cost = cost + 5*(D-J);
        end
        
        I = J-D;
    end
    y(q) = cost/100;
end
        
        
        