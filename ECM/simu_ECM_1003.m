
clear all;
FLAG_CORR=1;


% net.J_barre = 0;%-1;
% net.sigma_J = 3;
% net.nb_pop = 1;
% net.I = {0};


I_ref = 1; %0;%3.5; % 
net.I= {I_ref,I_ref};
net.nb_pop = 2;
J_ref= 10; %3; %3; % 20; %20; %15; %10; %
g = 0;
sigma_J_ref = 1; %2.2;  % sigma_J_total 
d = J_ref * sqrt(1 + (1 + g/J_ref)^2) / sigma_J_ref;

tt = 1; %net.tau_m_ref(1) / net.tau_m_ref(2);
net.J_barre =[J_ref        -(J_ref - g) ;
             J_ref * tt   -(J_ref - g) * tt];  
      
net.sigma_J = abs(net.J_barre/d);

%net.sigma_J=FLAG_CORR*abs (net.J_barre/d);
%net.nb_pop = 2;

%net.theta_barre=1;
%net.sigma_theta=0; 

net.delta_t=0.5;
net.theta=1;
net.tau_m=10;
net.nbp_m=net.tau_m/net.delta_t;
net.tau_r=1;
net.tau_d=3;
net.nbp_d=net.tau_d/net.delta_t;

nbp=500;

load tab_stat_091228;
net.data=data;

prop_H = 3 % Horizon

net=init_dyn_ECM_s(net,prop_H);

%net.ECM.m(2,:)=0*net.ECM.m(2,:);
%net.ECM.q(2,:)=0*net.ECM.q(2,:);
%net.ECM.C(2,:)=0*net.ECM.C(2,:);

net=iter_dyn_ECM_win_s(net,nbp,FLAG_CORR,prop_H);
%net=iter_dyn_ECM_s_1003(net,nbp,FLAG_CORR);

figure(1);clf;
subplot(2,1,1);
plot((1:nbp)*net.delta_t,net.ECM.DYN_M*1000);
title('MEAN ACTIVITY (Hz)');
xlabel('Time (ms)');
%hold on;
%plot((1:nbp)*net.delta_t,net.ECM.DYN_Q*net.tau_m/net.delta_t,'r');
subplot(2,1,2);
%plot((1:nbp)*net.delta_t,net.ECM.DYN_MU_H);
hold on;
plot((1:nbp)*net.delta_t,sqrt(net.ECM.DYN_NU_H),'r');
title('SYNAPTC POTENTIAL STANDARD DEVIATION');
xlabel('Time (ms)');

