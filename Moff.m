function [ yoff, ytotal ] = Moff ( M, N, m, Possort )
% compute off diagonal values and total values
yfro2 = sum(M.^2,2);
ydiag2 = zeros(N,1);
for n = 1 : N
    ydiag2(n) = M(n, Possort(n))^2;
end

ytotal = sum(yfro2);
yoff = ytotal - sum(ydiag2);

end