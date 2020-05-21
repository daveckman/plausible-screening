function y = sSsimuCRN(x,reps)

if length(reps) == 1
    reps = reps*ones(size(x,1),1);
end

n = size(x,1);
y  = zeros(n,max(reps));

Ds = poissrnd(25,[100,max(reps)]);
for q = 1:n
    for g = 1:reps(q)
        s = x(q,1);
        S = x(q,2);
        I = S;
        cost = 0;
        for k = 1:100
            
            J = I;
            D = Ds(k,g);
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
        y(q,g) = cost/100;
    end
end


