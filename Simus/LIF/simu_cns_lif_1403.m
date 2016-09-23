%path(path, genpath('~/RRNN-git'));
path(path, genpath('/donnees/edauce/RRNN-git'));

mem = {}
for SEED = 2 % 8%60 %

    net = init_param_cns_lif(0.5); 
    
    % taille 
    %net.N = [100;
    %         40];

    % type de neurone
    net.FLAG_REV = 1; 
    net.FLAG_SFA = 1;
    
    % type de delai
    net.tau_d = 10;   % Délai moyen (ms)
    net.nbp_d = floor(net.tau_d/net.delta_t);  % D�lai effectif (en pas de temps)
    sigma_tau =  8;  % ecart-type (en pas de temps!!)
    net.tau_min= (net.nbp_d - sigma_tau) * ones(net.nb_pop);             % net.tau_min : delai de transmission minimal                
    net.tau_moy= sigma_tau * ones(2);   % parametre distri Poisson      
    net.FLAG_UNIF = 0;                  % delais uniformes / Poisson    

    % type d'input
    net.ENV_PERIOD = 30;
    net.ENV_RENEWAL_RATE = 1; %2/3; %2/3; %1; %2/3;
    net.ENV_FLAG_NON_FLAT = 0;
    
    % Param apprentissage
    net.alpha = 1;
    net.REWARD = 1;
    % Optim si pas de plasticité
    net.FLAG_DASHBOARD=1;

    net = init_systeme_lif(net,SEED); % 8
    net = init_dyn_lif(net);
    net=iter_dyn_lif(net,600,0);
    

    %net=iter_dyn_lif(net,10000,1);

    %net=iter_dyn_lif(net,1000,0);

    %mem{SEED} = net;
    
    %save tmp mem
    
end
%save simu_cns_lif_1403_seed60_tau10_period60ms_RENEW066 net
