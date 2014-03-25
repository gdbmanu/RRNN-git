path(path,genpath('/donnees/RRNN-git'))

SEED = 245

net = init_param_cns_lif(0.5); net = init_systeme_lif(net,SEED); % 8

net = init_dyn_lif(net);

net.ENV_PERIOD = 60;
net.ENV_RENEWAL_RATE = 0;

net=iter_dyn_lif(net,2000,0);

net=iter_dyn_lif(net,6000,1);

net=iter_dyn_lif(net,2000,0);

save simu_cns_lif_1403_seed245_tau10_period60ms_RENEW00 net
