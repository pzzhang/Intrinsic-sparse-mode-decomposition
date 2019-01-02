% % media covariance matrix generator
% %% random modes generate
% Nx = 32;
% Ny = 32;
% nChannels = 16;
% nInclusions = 18;
% nOthers = 1;
% nModes = nChannels + nInclusions + nOthers;
% % characteristics generating
% charModes = zeros(Nx,Ny,nModes);
% charModes2 = zeros(Nx*Ny, nModes);
% % Channels
% indChannels = randi([3 4],1,nChannels);
% xCenter = randperm(Nx,nChannels);
% yCenter = randperm(Ny,nChannels);
% for i = 1 : nChannels
%     switch indChannels(i)
%         case 3
%             % x direction channels
%             Ltemp = 2*randi([3 14]);
%             xmin = max(1, xCenter(i)-Ltemp/2); xmax = min(Nx, xCenter(i)+Ltemp/2);
%             ymin = max(1,yCenter(i)); ymax = min(Ny, yCenter(i));
%             charModes(xmin:xmax,ymin:ymax,i) = ones(xmax-xmin+1,ymax-ymin+1);
%         case 4
%             % y direction channels
%             Ltemp = 2*randi([3 14]);
%             xmin = max(1, xCenter(i)); xmax = min(Nx, xCenter(i));
%             ymin = max(1, yCenter(i)-Ltemp/2); ymax = min(Ny, yCenter(i)+Ltemp/2);
%             charModes(xmin:xmax,ymin:ymax,i) = ones(xmax-xmin+1,ymax-ymin+1);
%     end
%     charModes2(:,i) = reshape(charModes(:,:,i), Nx*Ny, 1);
% end
% % Inclusions
% indInclusions = randi([1 2],1,nInclusions);
% xCenter = randperm(Nx,nInclusions);
% yCenter = randperm(Ny,nInclusions);
% for i = 1 : nInclusions
%     switch indInclusions(i)
%         case 1
%            % inclusions
%             xmin = max(1, xCenter(i)); xmax = min(Nx, xCenter(i)+1);
%             ymin = max(1, yCenter(i)); ymax = min(Ny, yCenter(i)+1);
%             charModes(xmin:xmax,ymin:ymax,nChannels+i) = ones(xmax-xmin+1,ymax-ymin+1);
%         case 2
%             % cross
%             xmin = max(1, xCenter(i)-1); xmax = min(Nx, xCenter(i)+1);
%             ymin = max(1, yCenter(i)-1); ymax = min(Ny, yCenter(i)+1);
%             charModes(xmin:xmax,yCenter(i),nChannels+i) = ones(xmax-xmin+1,1);
%             charModes(xCenter(i),ymin:ymax,nChannels+i) = ones(1,ymax-ymin+1);
%     end
%     charModes2(:,nChannels+i) = reshape(charModes(:,:,nChannels+i), Nx*Ny, 1);
% end
% % add another smiling face
% load('data/smileFace32.mat');
% charModes(:,:,nModes) = smileFace32;
% charModes2(:,nModes) = reshape(charModes(:,:,nModes), Nx*Ny, 1);
% % number of samples
% Ns = 10;
% % scale
% c = 1;
% sigma = 24;
% weights = c*(1+sigma*rand(nModes, Ns));
% fieldSamples = charModes2*weights;
% fieldCha = charModes2*ones(nModes,1);
% 
% % variance field
% Cov = charModes2*charModes2';
% Cov = c^2*sigma^2*Cov/12;

% load media data
% load('data/mediaData.mat');
% hx = 1/Nx; hy = 1/Ny;
% xx = hx/2 : hx : 1-hx/2;
% yy = hy/2 : hy : 1-hy/2;
% plot
h1 = figure(1);
set(gca,'fontsize',20);
surf(xx, yy, reshape(fieldCha,Nx,Ny))
xlim([0 1]);ylim([0 1]);
grid on;
title('Characteristics')
pause;
print(gcf,'-depsc2','data/Characteristics.eps');
% 
h2 = figure(2);
set(gca,'fontsize',20);
surf(xx, yy, reshape(fieldSamples(:,2),Nx,Ny))
xlim([0 1]);ylim([0 1]);
grid on;
title('fieldSample')
pause;
print(gcf,'-depsc2','data/fieldSample.eps');
% 
% h3 = figure(3);
% set(gca,'fontsize',20);
% contour(Cov)
% xlim([0 96^2]);ylim([0 96^2]);
% grid on;
% title('Covariance')
% pause;
% print(gcf,'-depsc2','data/Covariance.eps');
% 
h4 = figure(4);
set(gca,'fontsize',20);
contour(xx, yy, reshape(fieldSamples(:,1),Nx,Ny),5)
xlim([0 1]);ylim([0 1]);
grid on;
title('fieldSample')
pause;
print(gcf,'-depsc2','data/fieldSampleContour.eps');

close(h1);
close(h2);
% close(h3);
close(h4)

% eigen-decomposition
% set eigs options
opts.issym = 1;
opts.isreal = 1;
% compute the first nmodes eigenfunctions
[VCov,DCov] = eigs(Cov, nModes, 'lm', opts);

for k = 1:nModes
    filename1 = ['data/Realmodes/media_real_',num2str(k),'.eps'];
    filename2 = ['data/KLmodes/media_svd_',num2str(k),'.eps'];
    
    h1 = figure(1);
    set(gca,'fontsize',20);
    contour(xx, yy, reshape(charModes2(:,k),Nx,Ny),5);
    xlim([0 1]);ylim([0 1]);
    title(['Real mode']);
    grid on;
    print(gcf,'-depsc2',filename1);
    
    h2 = figure(2);
    set(gca,'fontsize',20);
    contour(xx, yy, reshape(VCov(:,k),Nx,Ny),5);
    xlim([0 1]);ylim([0 1]);
    title(['KL recovery']);
    grid on;
    print(gcf,'-depsc2',filename2);
    
    pause;
    close(h1);
    close(h2);
end