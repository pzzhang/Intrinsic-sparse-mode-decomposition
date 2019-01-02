% test JointDiag.m
K = 20;
n = 10;
% eigenvalues
D = randn(n,K);
U = randn(n,n);
[U,~] = qr(U);

% matrices
M = zeros(n,n,K);
for k = 1 : K
    M(:,:,k) = U*diag(D(:,k))*U';
end

% joint diag
eps = 1e-14;
[ V, Sigma ] = JointDiag( eps, M );
error = zeros(1,K);
for k = 1 : K     
    error(k) = norm(V*Sigma(:,:,k)*V' - M(:,:,k), 'fro');
end

plot(error);