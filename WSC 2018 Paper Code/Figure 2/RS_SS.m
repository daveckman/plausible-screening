function L = RS_SS(Xtry,X,Y,alpha)
Xu = unique(X);
for k =1:length(Xu)
    Ym = Y(X == Xu(k));
    n(k) =length(Ym);
    sigma2hat(k) = var(Ym);
    muhat(k) = mean(Ym);
end
spos = sigma2hat(sigma2hat>(10^(-10)*max(sigma2hat)));
sigma2hat(sigma2hat<=quantile(spos,0.1)) = quantile(spos,0.1);

K = length(Xu);

L = ones(size(Xtry));
for k = 1:length(Xtry)
    x0 = Xtry(k);
    if ismember(x0,Xu)
        i = find(Xu==x0);
        if n(i) >= 1.5
            ti = tinv((1-alpha)^(1/(K-1)),n(i)-1);
            for l = 1:K
                if n(l) >= 1.5
                    tj = tinv((1-alpha)^(1/(K-1)),n(l)-1);
                    W = (ti^2*sigma2hat(i)/n(i)+tj^2*sigma2hat(l)/n(l))^(1/2);
                    if muhat(i) > muhat(l) + W
                        L(k) = 0;
                    end
                end
            end
        else
            L(k) = 1;
        end
    end

end


% end
% 
% function p = ttest_me(m1,v1,n1,m2,v2,n2)
% 
% t_stat = (m1-m2)/sqrt(v1/n1+v2/n2);
% 
% dof = (v1/n1+v2/n2)^2/(v1^2/n1^2/(n1-1)+v2^2/n2^2/(n2-1));
% 
% p = 2*tcdf(-abs(t_stat),dof);
% 
% if n1 == 0 || v1 < 10^(-14) || n2 == 0 || v2 < 10^(-14)
%     p = 1;
% end
% end