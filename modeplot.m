% plot modes
for i = 1 : nModes
    mesh(reshape(V(:,i),Nx,Ny));
    pause;
end