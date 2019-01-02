% HH and CPU time
HH= [1/2, 1/3, 1/4, 1/6, 1/8, 1/12, 1/16, 1/24];
CPUtime = [7.2961, 1.9438, 1.0621, 0.7846,  0.8988, 1.9623, 3.2092, 5.7554];
CPUtime2 = [1.0942, 0.8856, 0.6204, 0.6652,  0.8304, 1.9059, 3.1999, 5.7468];
EigenTime = 95.908625*ones(size(HH));
EigenTime2 = 2*ones(size(HH));

% plot
filename1 = ['data/CPUtime.eps'];
h1 = figure(1);
loglog(HH, CPUtime,'*-');
hold on
loglog(HH, EigenTime, 'r--','LineWidth',3);
xlim([1/32 3/4]); ylim([0.5 200]);
legend('ISMD', 'Eigen Decomposition');
xlabel('subdomain size: H');
ylabel('CPU time');
set(gca,'xtick',fliplr(HH));
title('Comparison of CPU time');
pause;
print(gcf,'-depsc2',filename1);
close(h1);

filename1 = ['data/CPUtime2.eps'];
h1 = figure(1);
loglog(HH, CPUtime2,'*-');
hold on
loglog(HH, EigenTime2, 'k--','LineWidth',3);
xlim([1/32 3/4]); ylim([0.5 8]);
legend('ISMD', 'Eigen Decomposition');
xlabel('subdomain size: H');
ylabel('CPU time');
set(gca,'xtick',fliplr(HH));
title('Comparison of CPU time with known rank');
pause;
print(gcf,'-depsc2',filename1);
close(h1);
