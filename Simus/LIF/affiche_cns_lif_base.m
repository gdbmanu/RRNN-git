clf
axe_tau = 1:41
%[tau,axe_tau]=
hist(nonzeros(net.tau{1}{1})*net.delta_t,axe_tau/2);%);
colormap('default')
axis([0 20 0 30000])
%bar(axe_tau,tau);
%hold on;
%plot(axe_tau+net.tau_min(1,1)*net.delta_t,pdf('Poisson',axe_tau/net.delta_t,net.tau_moy(1,1)),'-o')
%plot(axe_tau*net.delta_t,pdf('Poisson',axe_tau-net.tau_min(1,1),net.tau_moy(1,1)),'-o')

xlabel('Delay (ms)')
ylabel('# axons')