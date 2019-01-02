% Given Cov, now do sparse EIG
% Originally, we have Nx = Ny = 32;
load('data/mediaData.mat');

% amplify ratio
AR = 3;
Nx = AR*Nx; Ny = AR*Ny;

% amplify the physical modes
charModes0 = charModes;
charModes = zeros(Nx,Ny,nModes);
charModes2 = zeros(Nx*Ny, nModes);
for i = 1 : nModes
    charModes(:,:,i) = kron(charModes0(:,:,i),ones(AR));
    charModes2(:,i) = reshape(charModes(:,:,i), Nx*Ny, 1);
end
clear charModes0;

% form amplified COV
% number of samples
Ns = 10;
% scale
c = 1;
sigma = 24;
weights = c*(1+sigma*rand(nModes, Ns));
fieldSamples = charModes2*weights;
fieldCha = charModes2*ones(nModes,1);

% variance field
Cov = charModes2*charModes2';
Cov = c^2*sigma^2*Cov/12;

hx = 1/Nx; hy = 1/Ny;
xx = hx/2 : hx : 1-hx/2;
yy = hy/2 : hy : 1-hy/2;