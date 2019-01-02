% Cov = kron(Cov,ones(5));
tic
[U,Sigma] = eig(Cov);
toc

% set eigs options
opts.issym = 1;
opts.isreal = 1;
tic
[U,Sigma] = eigs(Cov, K, 'lm', opts);
toc

tic
[U,Sigma] = eigs(Cov, Nx*Ny-2, 'lm', opts);
toc