% Given Cov, now do sparse EIG
% we have Nx = Ny = 50; 
npx = 10; npy = 10;
nfx = 5; nfy = 5;

% set eigs options
opts.issym = 1;
opts.isreal = 1;
% local dimension upper bound estimate and tolerance
dloc = 10;
eigtol = 1e-4;

%% locally do eigenvalue decomposition
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
        [Vloc,Dloc] = eigs(Covloc{i,j}, dloc, 'lm', opts);
        indpos = (diag(Dloc) > eigtol);
        Kloc(i,j) = sum(indpos);
        Hloc{i,j} = Vloc(:,indpos)*sqrt(Dloc(indpos,indpos));
    end
end

% total dimension and assemble H
Dtoty = sum(Kloc);
Dtotal = sum(Dtoty);
H = zeros(Nx*Ny, Dtotal);
for j = 1 : npy
    nystart = sum(Dtoty(1:j-1));
    for i = 1 : npx
        nxstart = nystart+sum(Kloc(1:i-1,j))+1;
        nxend = nystart+sum(Kloc(1:i,j));
        H(Indloc{i,j}, nxstart:nxend) = Hloc{i,j};
    end
end

% compute Gamma
condH = cond(H);
Gamma = H\(Cov/H');
errH = norm(Cov-H*Gamma*H','fro');
[Vglo, Dglo] = eigs(Gamma,nModes, 'lm', opts);
Veff = Vglo*sqrt(Dglo);

% collect Gamma's eigenspace
dglo = diag(Dglo);
dglo = floor(dglo + 0.1);
uniqdglo = unique(dglo);
uniqdglo = sort(uniqdglo, 'ascend');
Nd = length(uniqdglo);
numD = zeros(1,Nd);
eigenspaceG = cell(1,Nd);
Covdglo = cell(1,Nd);
eigenspaceG2 = cell(1,Nd);
for i = 1 : Nd
    indtemp = (dglo == uniqdglo(i));
    numD(1,i) = sum(indtemp);
    eigenspaceG{1,i} = Veff(:,indtemp);
    Covdglo{1,i} = eigenspaceG{1,i}*eigenspaceG{1,i}';
%     [Vdglo,Ddglo] = eigs(Covdglo{1,i}, numD(1,i), 'lm', opts);
    [Vdglo,Ddglo] = eig(Covdglo{1,i});
    [ddglo, indpos] = sort(diag(Ddglo), 'descend');
    Vdglo = Vdglo(:,indpos);
    eigenspaceG2{1,i} = Vdglo(:,1:numD(1,i))*diag(sqrt(ddglo(1:numD(1,i))));
end

%% solve L
