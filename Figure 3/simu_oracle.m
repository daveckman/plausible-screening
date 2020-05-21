function y = simu_oracle(x)

n = size(x,1);
y  = zeros(1,n);

pvec =@(x) (poisspdf(x,25)/(sum(poisspdf(0:50,25)))).*(x<=50).*(x>=0);
for q = 1:n
    s = x(q,1);
    S = x(q,2);
    if s < S
        
        PI = [(s-50):S];
        TranMat = zeros(length(PI),length(PI));
        for k = 1:length(PI)
            if PI(k) <= s
                TranMat(k,:) = pvec(S-PI);
            else
                TranMat(k,:) = pvec(PI(k)-PI);
            end
        end
        e = zeros(size(PI))';
        e(PI==S) = 1;
        cost_vec = zeros(size(PI));
        
        for k = 1:length(PI)
            cost_vec(k) = (PI(k)<=s).*(32+3*(S-PI(k)))+sum( TranMat(k,:).*((PI<0).*abs(PI)*5+(PI>0).*abs(PI)));
        end
        
        avg_cost = 0;
        for k = 1:100
            avg_cost = avg_cost+ 1/100*(cost_vec*(TranMat')^(k-1)*e);
        end
        
        y(q) = avg_cost;
    else
        y(q) = nan;
    end
    
    
end

