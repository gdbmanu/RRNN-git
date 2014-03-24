% Initialisation champ moyen P populations
% On suppose connus stat_theta_m stat_theta_std stat_J_m stat_J_std

function net=init_dyn_ECM(net,prop_H);

nb_pop=length(net.J_barre);

net.ECM.DYN_MU=[];
net.ECM.DYN_NU=[];
net.ECM.DYN_DELTA=[];
net.ECM.DYN_ALPHA=[];

net.ECM.DYN_MU_H=[];
net.ECM.DYN_NU_H=[];
net.ECM.DYN_DELTA_H=[];
net.ECM.DYN_DELTA_H_K=[];
net.ECM.DYN_ALPHA_H=[];

net.ECM.DYN_M=[];
net.ECM.DYN_Q=[];
net.ECM.DYN_C=[];

tau_ref = net.tau_min(1,1) + net.tau_moy(1,1);
H = prop_H * net.nbp_m;

f0 =  0.1; %0.01; %

tau_max = max(H+1, 2 * tau_ref - net.tau_min(1,1));

net.ECM.m=f0*ones(nb_pop,tau_max);
net.ECM.q=f0^2*ones(nb_pop,tau_max);
net.ECM.C=f0^2*ones(nb_pop,tau_max);

net.ECM.mu=0*ones(nb_pop,tau_max);
net.ECM.nu=0.03^2*ones(nb_pop,tau_max);
net.ECM.Delta=0.03^2*ones(nb_pop,tau_max);
net.ECM.alpha=0*ones(nb_pop,tau_max);

net.ECM.mu_h=zeros(nb_pop,tau_max);
net.ECM.nu_h=f0^2 * ones(nb_pop,tau_max);
net.ECM.Delta_h=f0^2 * ones(nb_pop,tau_max);
net.ECM.Delta_h_K=f0^2 * ones(nb_pop,tau_max);
net.ECM.alpha_h=zeros(nb_pop,tau_max);

net.ECM.sigma=0.1*ones(nb_pop,tau_max);
