
net = init_param_cns_lif(0.5); net = init_systeme_lif(net,60); % 8

net = init_dyn_lif(net);

net.FLAG_RENEWAL_RATE = 0;

net=iter_dyn_lif(net,2000,0);

net=iter_dyn_lif(net,6000,1);

net=iter_dyn_lif(net,2000,0);

save simu_cns_lif_1403_seed60_period60ms_RENEW00 net
