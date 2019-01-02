% plot error
NoiseLevel = [];
ErrorMat = [];

for n = 16 : -2 : 2
    filename4 = ['ErrorModes_e',num2str(n),'.mat'];
    load(filename4);
    NoiseLevel = [NoiseLevel; noiselevel];
    ErrorMat = [ErrorMat; ErrorModes];
end

% plot 
filename5 = ['ErrorPlotLInfinity.eps'];
h1 = figure(1);
loglog(NoiseLevel,max(ErrorMat,[],2),'*-');
xlabel('NoiseLevel')
ylabel('RecoveredError')
title('Stability Results in L^{\infty} norm');
pause;
print(gcf,'-depsc2',filename5);
close(h1);

filename6 = ['ErrorPlotL2.eps'];
h2 = figure(2);
loglog(NoiseLevel,sqrt(sum(ErrorMat.^2,2)),'*-');
xlabel('NoiseLevel')
ylabel('RecoveredError')
title('Stability Results in L^{2} norm');
pause;
print(gcf,'-depsc2',filename6);
close(h2);

filename7 = ['ErrorPlot.eps'];
h3 = figure(3);
loglog(NoiseLevel,max(ErrorMat,[],2),'*-');
hold on
loglog(NoiseLevel,sqrt(sum(ErrorMat.^2,2)),'*-');
xlabel('NoiseLevel')
ylabel('RecoveredError')
legend('L^{\infty} error', 'L^{2} error')
title('Stability Results');
pause;
print(gcf,'-depsc2',filename7);
close(h3);