% Given Cov, now do sparse EIG
% we have Nx = Ny = 32; 
patchL = 2;
nfx = patchL; nfy = patchL;
npx = Nx/nfx; npy = Ny/nfy;

% set eigs options
opts.issym = 1;
opts.isreal = 1;
% local dimension upper bound estimate and tolerance
dloc = 10;
eigtol = 1e-4;

%% locally do eigenvalue decomposition
tic
Covloc = cell(npx, npy);
Hloc = cell(npx, npy);
Indloc = cell(npx, npy);
IndTotal = reshape(1:Nx*Ny, Nx, Ny);
Kloc = zeros(npx,npy);
for i = 1 : npx
    indx = (i-1)*nfx+1 : i*nfx;
    for j = 1 : npy
        indy = (j-1)*nfy+1 : j*nfy;
        Indloc{i,j} = IndTotal(indx,indy);
        Indloc{i,j} = reshape(Indloc{i,j}, 1, []);
        Covloc{i,j} = Cov(Indloc{i,j},Indloc{i,j});
%         [Vloc,Dloc] = eigs(Covloc{i,j}, nfx*nfy-2, 'lm', opts);
        [Vloc,Dloc] = eigs(Covloc{i,j}, nfx*nfy, 'lm', opts);
        indpos = (diag(Dloc) > eigtol);
        Kloc(i,j) = sum(indpos);
        [Vloc,Dloc] = eigs(Covloc{i,j}, Kloc(i,j), 'lm', opts);
        Hloc{i,j} = Vloc*sqrt(Dloc);
    end
end
clear Vloc Dloc;
toc

% 2-d patch index
M = npx*npy;
Pindex = zeros(M,2);
Pindex(:,1) = repmat((1:npx)',npy,1);
Pindex(:,2) = kron((1:npy)', ones(npx,1));
Kpatch = reshape(Kloc,M,1);
Dindex = ones(M+1,1);
for m = 2 : M+1
    Dindex(m) = Dindex(m-1) + Kpatch(m-1);
end
Dtotal = Dindex(M+1)-1;

%% assemble Lambda
tic
Lambda = eye(Dtotal);
for m = 1 : M
    indm = Indloc{Pindex(m,1),Pindex(m,2)};
    Hm = Hloc{Pindex(m,1),Pindex(m,2)};
    for n = m+1 : M
        indn = Indloc{Pindex(n,1),Pindex(n,2)};
        Hn = Hloc{Pindex(n,1),Pindex(n,2)};
        Lambda(Dindex(m):Dindex(m+1)-1,Dindex(n):Dindex(n+1)-1)...
            = Hm\(Cov(indm,indn)/Hn');
        Lambda(Dindex(n):Dindex(n+1)-1,Dindex(m):Dindex(m+1)-1)...
            = Lambda(Dindex(m):Dindex(m+1)-1,Dindex(n):Dindex(n+1)-1)';
    end
end
toc

% eigenvalue decomposition of Lambda
tic
% [Vlambda, Dlambda] = eigs(Lambda,Dtotal-2,'lm',opts);
[Vlambda, Dlambda] = eigs(Lambda,Dtotal,'lm',opts);
toc
filename1 = ['data/dl',num2str(patchL),'/EigenvalueOfLambda.eps'];
h1 = figure(1);
plot(diag(Dlambda),'r*');
title('Eigenvalues of Correlation Matrix \Lambda');
pause;
print(gcf,'-depsc2',filename1);
close(h1);
% determine K
indpos = (diag(Dlambda) > 0.5);
K = sum(indpos);
clear Vlambda;

%% simultaneous diagonalization
tic
epsilon = 1e-14;
Dloc = cell(1,M);
Gloc = cell(1,M);
for m = 1 : M
    indm = Dindex(m):Dindex(m+1)-1;
    Sigma = zeros(Kpatch(m),Kpatch(m),M);
    for n = 1 : M
        indn = Dindex(n):Dindex(n+1)-1;
        Sigma(:,:,n) = Lambda(indm,indn)*Lambda(indm,indn)';
    end
    Dloc{1,m} = JointDiag( epsilon, Sigma );
    Gloc{1,m} = Hloc{Pindex(m,1),Pindex(m,2)}*Dloc{1,m};
end
toc

%% recover Omega, LL^T
tic
Omega = eye(Dtotal);
for m = 1 : M
    indm = Dindex(m):Dindex(m+1)-1;
    for n = m+1 : M
        indn = Dindex(n):Dindex(n+1)-1;
        Omega(indm,indn) = Dloc{1,m}'*Lambda(indm,indn)*Dloc{1,n};
        Omega(indn,indm) = Omega(indm,indn)';
    end
end
toc
% Deal with the case when it is not identifiable
Dloc2 = cell(1,M);
Dloc2{1,1} = eye(Kpatch(1));
for m = 2 : M
    indm = Dindex(m):Dindex(m+1)-1;
    % get the right rotation
    Dloc2{1,m} = JointOrtho(epsilon, Omega(1:Dindex(m)-1,indm));
    % change Omega accordingly
    Omega(indm,:) = Dloc2{1,m}'*Omega(indm,:);
    Omega(:,indm) = Omega(:,indm)*Dloc2{1,m};
    % change Gm accordingly
    Gloc{1,m} = Gloc{1,m}*Dloc2{1,m};
end

%% recover L
tic
Omega = IntAdjust(Omega);
[Vomega,Domega] = eigs(Omega,K,'lm',opts);
L = Vomega*sqrt(Domega);
U = JointOrtho(epsilon,L);
L = IntAdjust(L*U);
toc

%% recover g
tic
gISMD = zeros(Nx*Ny,K);
for m = 1 : M
    indm = Indloc{Pindex(m,1),Pindex(m,2)};
    Grec = Gloc{1,m}*L(Dindex(m):Dindex(m+1)-1,:);
    gISMD(indm,:) = Grec;
end
toc
%% plot the result
% for k = 1 : K
%     filename2 = ['data/dl',num2str(patchL),'/media_recover_',num2str(k),'.eps'];
%     h2 = figure(2);
%     set(gca,'fontsize',20);
%     contour(reshape(gISMD(:,k),Nx,Ny),5);
%     xlim([0 32]);ylim([0 32]);
%     title(['Sparse recovery']);
%     set(gca,'xtick',[0:nfx:Nx])
%     set(gca,'ytick',[0:nfy:Ny])
%     grid on;
%     print(gcf,'-depsc2',filename2);
%     pause;
%     close(h2);
% end