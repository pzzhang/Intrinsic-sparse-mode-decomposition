Cov = kron(Cov,ones(5));
tic
[U,Sigma] = eig(Cov);
toc

tic
[U,Sigma] = eigs(Cov, K, 'lm', opts);
toc