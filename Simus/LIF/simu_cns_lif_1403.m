path(path, genpath('/donnees/edauce/RRNN-git'));

SEED = 60

net = init_param_cns_lif(0.5); net = init_systeme_lif(net,SEED); % 8

net = init_dyn_lif(net);

net.ENV_PERIOD = 60;
net.ENV_RENEWAL_RATE = 2/3;

net=iter_dyn_lif(net,2000,0);

net=iter_dyn_lif(net,6000,1);

net=iter_dyn_lif(net,2000,0);

save simu_cns_lif_1403_seed60_tau10_period60ms_RENEW066 net
