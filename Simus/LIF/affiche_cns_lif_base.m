clf

subplot(3,1,1)

[w_plus,axe_w_plus]=hist(nonzeros(net.J_ref{1}{1}),0.025:0.05:1.975);
bar(axe_w_plus,w_plus);
hold on;
[w_moins,axe_w_moins]=hist(nonzeros(net.J_ref{1}{2}),-1.975:0.05:0.025);
bar(axe_w_moins,w_moins);
title('INDIVIDUAL SYNAPTIC WEIGHTS')

subplot(3,1,2)

[W_plus,axe_W_plus]=hist(sum(net.J_ref{1}{1}'),0.5:1:15.5);
bar(axe_W_plus,W_plus/net.N(1));
hold on;
[W_moins,axe_W_moins]=hist(sum(net.J_ref{1}{2}'),-15.5:1:5);
bar(axe_W_moins,W_moins/net.N(1));
[W_sum,axe_W_sum]=hist(sum([net.J_ref{1}{1},net.J_ref{1}{2}]'),-15.5:1:15.5);
bar(axe_W_sum,W_sum/net.N(1));
axis([-15 15 0 0.5]);
axe=-15:0.1:15;
plot(axe,pdf('norm',axe,net.J_barre(1,2),net.sigma_J_eff(1,2)));
plot(axe,pdf('norm',axe,net.J_barre(1,1),net.sigma_J_eff(1,1)));
plot(axe,pdf('norm',axe,net.J_barre(1,1)+net.J_barre(1,2),sqrt(net.sigma_J_eff(1,1)^2+net.sigma_J_eff(1,2)^2)));
title('WEIGHTS SUM');


subplot(3,1,3)
axe_tau = 1:41
%[tau,axe_tau]=
hist(nonzeros(net.tau{1}{1})*net.delta_t,axe_tau/2);%);
colormap('default')
axis([0 20 0 net.K(1,1)^2])
%bar(axe_tau,tau);
%hold on;
%plot(axe_tau+net.tau_min(1,1)*net.delta_t,pdf('Poisson',axe_tau/net.delta_t,net.tau_moy(1,1)),'-o')
%plot(axe_tau*net.delta_t,pdf('Poisson',axe_tau-net.tau_min(1,1),net.tau_moy(1,1)),'-o')

xlabel('Delay (ms)')
ylabel('# axons')