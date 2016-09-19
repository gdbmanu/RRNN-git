path(path, genpath('~/RRNN-git'));

%mem = {}
nb_app = 10000

for SEED = 1:1

    net = init_param_cns_lif(0.5); net = init_systeme_lif(net,SEED); % 8

    net = init_dyn_lif(net);

    net.ENV_PERIOD = 30;
    net.ENV_RENEWAL_RATE = 11/12 % %1; %2/3;

    net=iter_dyn_lif(net,400,0);

    net=iter_dyn_lif(net,nb_app,1);

    net=iter_dyn_lif(net,1000,0);
    
    %mem{SEED} = net
    
    %save tmp mem
end;

%save simu_cns_lif_1403_seed60_tau10_period60ms_RENEW066 net
