% test JointOrtho.m
N = 50;
m = 10;
% eigenvalues
MaxPos = randi(m,N,1);
Omega = zeros(N,m);
for n = 1 : N
    Omega(n,MaxPos(n)) = 1;
end
U = randn(m,m);
[U,~] = qr(U);
M = Omega*U';

% joint ortho
eps = 1e-14;
[ V ] = JointOrtho( eps, M );
Sigma = M*V;
[~, MaxPos2] = max(abs(Sigma),[],2);
for i = 1 : N     
    Sigma(i,MaxPos2(i)) = 0;
end
error = sum(Sigma.^2,2);

plot(error);