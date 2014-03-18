path(path, genpath('~/RRNN-git'));
clear all;

% Parametres ECM
FLAG_CORR=1;
prop_H = 1 % Horizon

% Reseau : general
net.fichier_I={''};  
net.script_out='iter_env_cns_lif';
net.script_init='init_env_cns_lif';

net.nb_pop = 2;
net.delta_t=0.5;
net.theta=1;
net.tau_m=10;
net.nbp_m=net.tau_m/net.delta_t;
net.tau_r=1;

% Poids
J_ref= 10;
net.J_barre = [J_ref -J_ref;
               J_ref -J_ref];

% Delais (uniformes)
sigma_tau =  19
net.tau_min= (20 - sigma_tau) * ones(net.nb_pop);
net.tau_moy= sigma_tau * ones(2);

nbp=500;

load tab_stat_091228;
net.data=data;

mem_M_1 = zeros(31,31);
mem_Q_1 = zeros(31,31);
mem_C_1 = zeros(31,31);
mem_M_2 = zeros(31,31);
mem_Q_2 = zeros(31,31);
mem_C_2 = zeros(31,31);

for i = 1:31
    I_ref = 0.1 * (i -1);
    net.I{1} = I_ref;%0;%
    net.I{2} = I_ref;%0;%
    for j = 1:31
        sigma_J_ref = 0.1 * (j - 1);
        d = J_ref * sqrt(2) / sigma_J_ref;
        net.sigma_J = abs(net.J_barre/d);        
        net=init_dyn_ECM_s(net,prop_H);
        net=iter_dyn_ECM_win_s(net,nbp,FLAG_CORR,prop_H);
        mem_M_1(i,j) = net.ECM.DYN_M(1,nbp);
        mem_Q_1(i,j) = net.ECM.DYN_Q(1,nbp);
        mem_C_1(i,j) = net.ECM.DYN_C(1,nbp);
        mem_M_2(i,j) = net.ECM.DYN_M(2,nbp);
        mem_Q_2(i,j) = net.ECM.DYN_Q(2,nbp);
        mem_C_2(i,j) = net.ECM.DYN_C(2,nbp);
        save simu_ECM_1403_J10_tau19 mem_M_1 mem_Q_1 mem_C_1 mem_M_2 mem_Q_2 mem_C_2
    end
end

