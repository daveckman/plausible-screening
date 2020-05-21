alpha = 0.05;

z0 = norminv(1-alpha/2);
for k = 1:100
    z(k) = norminv(1-alpha/k/2);
end
    