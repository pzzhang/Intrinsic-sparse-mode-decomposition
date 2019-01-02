% media covariance matrix generator
%% random modes generate
Nx = 50;
Ny = 50;
nModes = 50;
% characteristics generating
charModes = zeros(Nx,Ny,nModes);
charModes2 = zeros(Nx*Ny, nModes);
indModes = randi([1 3],1,nModes);
xCenter = randperm(Nx, nModes);
yCenter = randperm(Ny, nModes);
for i = 1 : nModes
    if indModes(i) == 1
        % inclusions
        xmin = max(1, xCenter(i)); xmax = min(Nx, xCenter(i)+1);
        ymin = max(1, yCenter(i)); ymax = min(Ny, yCenter(i)+1);
        charModes(xmin:xmax,ymin:ymax,i) = ones(xmax-xmin+1,ymax-ymin+1);
    else if indModes(i) == 2
            % x direction channels
            Ltemp = 2*randi([3 8]);
            xmin = max(1, xCenter(i)-Ltemp/2); xmax = min(Nx, xCenter(i)+Ltemp/2);
            ymin = max(1,yCenter(i)); ymax = min(Ny, yCenter(i));
            charModes(xmin:xmax,ymin:ymax,i) = ones(xmax-xmin+1,ymax-ymin+1);
        else
            % y direction channels
            Ltemp = 2*randi([3 8]);
            xmin = max(1, xCenter(i)); xmax = min(Nx, xCenter(i));
            ymin = max(1, yCenter(i)-Ltemp/2); ymax = min(Ny, yCenter(i)+Ltemp/2);
            charModes(xmin:xmax,ymin:ymax,i) = ones(xmax-xmin+1,ymax-ymin+1);
        end
    end
    charModes2(:,i) = reshape(charModes(:,:,i), Nx*Ny, 1);
end
% number of samples
Ns = 10000;
% scale
c = 1;
sigma = 24;
weights = c*(1+sigma*rand(nModes, Ns));
fieldSamples = charModes2*weights;
fieldCha = charModes2*ones(nModes,1);

% variance field
Cov = charModes2*charModes2';
Cov = c^2*sigma^2*Cov/12;

% plot
% figure(1)
% surf(reshape(fieldCha,Nx,Ny))
% title('Characteristics')
% 
% figure(2)
% surf(reshape(fieldSamples(:,1),Nx,Ny))
% title('fieldSamples')
% 
% figure(3)
% mesh(Cov)
% title('Covariance')