
clear all;

% !! d = 2.8, J_ref = 5;

load 090629.mat
net.tau_m_ref = [10;
                 10];
nb_app = 1;


% MODIFS
% alpha
%net.alpha = 0.3 * 200 / 1000; % !!??
net.alpha = 10 * 200 / 1000; % !!??
% modif num_res (initialement 3)
net.num_reseau = 3;
num_res = net.num_reseau;    
% modif d (initialement 2.8)
d = 1.5;
net.sigma_J_eff = abs(net.connex/d);
net.sigma_J = sqrt(abs(net.sigma_J_eff.^2-(1-net.densite).*(net.J_barre.^2)./net.K));
% norm_alpha
%net.norm_alpha = -[1 -1; 1 -1];
net.norm_alpha = [1 -1; 1 -1];
net.FLAG_SFA = 1;

randn('seed',2);
for i=1:20
    I1{i}=randn(net.N(1),1) * 4;
    I2{i}=randn(net.N(2),1) * 4;
    %I1{i}=(rand(net.N(1),1) < 0.2) * 4;
    %I2{i}=(rand(net.N(2),1) < 0.2) * 4;
end;

net=init_systeme_lif(net,num_res);
                
for cpt_app = [1] %[2 2 2] %[4 4 4 4 4] %2 %:2:6%nb_pat

    net=init_dyn_lif(net);  
    
%%% INPUT   
    net.I{1} = zeros(net.N(1),1);
    net=iter_dyn_lif(net,t_max/delta_t,0);

    if nb_app>0
        for i=1:nb_app 
            net.I{1} = I1{cpt_app};
            net=iter_dyn_lif(net,0.5 * t_max/delta_t,1);                            
            net.I{1} = zeros(net.N(1),1);
            net=iter_dyn_lif(net,0.5 * t_max/delta_t,0);                            
        end;
    end;

    net.I{1} = I1{cpt_app} ;
    net=iter_dyn_lif(net,0.5 * t_max/delta_t,0);
    net.I{1} = zeros(net.N(1),1);
    net=iter_dyn_lif(net,0.5 * t_max/delta_t,0);
    
    %pause;
    affiche_nips_lif(net,t_max,nb_pat,nb_app);
    affiche_cns_lif;
    %pause;
    
    for num_fig=1:5
        figure(num_fig);
        %eval(['print -dpng 071112-num-',num2str(num_res),'-nb_pat-',num2str(cpt_app),'-fig-',num2str(num_fig)]);
    end;
    
    
end;
