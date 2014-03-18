
clear all;

nb_pat = 5;%
t_max = 300*nb_pat;%2000;   % ms
nb_app = 1;%
num_res = 5;    

cpt=1;
mem_net={};

delta_t=0.5;
net=init_param_cns_lif(delta_t);
net.app='stdp';
%m=net.theta_barre(1)/(1-exp(-net.tau_r/net.tau_m));
net.alpha = 0.1; %0.1;%0.003;%0.5;

net.FLAG_DASHBOARD=0;
net.FLAG_STD = 0;
net.FLAG_SFA = 0;
net.FLAG_REV=1;
net.FLAG_SCALING=0;
net.FLAG_STDP_EPS=0;
net.FLAG_DELTA=1;
net.FLAG_BAYES=0;
net.FLAG_THETA=0;

randn('seed',2);
for i=1:20
    I1{i}=randn(net.N(1),1) * 4;
    I2{i}=randn(net.N(2),1) * 4;
    %I1{i}=(rand(net.N(1),1) < 0.2) * 4;
    %I2{i}=(rand(net.N(2),1) < 0.2) * 4;
end;

net=init_systeme_lif(net,num_res);

%for p=1:net.nb_pop
%    for q=1:net.nb_pop
%        if net.connex(p,q)~=0
%            net.J{p}{q}=sign(net.connex(p,q))*abs(net.J{p}{q});
%            net.J_ref{p}{q}=sign(net.connex(p,q))*abs(net.J_ref{p}{q});
%        end;
%    end;
%end;
                

for cpt_app = [2 2 2] %[4 4 4 4 4] %2 %:2:6%nb_pat

    net=init_dyn_lif(net);  
    
%%% INPUT
    
    FLAG_P1=1;
    FLAG_P2=1;
    
    for num_pat=1:nb_pat
        net.I{1}=I1{num_pat}*FLAG_P1;net.I{2}=I2{num_pat}*FLAG_P2;%net.I{3}=net.delta_t/net.tau_mb{3};
        net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
    end;


    if nb_app>0
        for i=1:nb_app
            for num_pat=ones(1,6)*cpt_app%[4 4 4 4 4 4]
               net.I{1}=I1{num_pat}*FLAG_P1;net.I{2}=I2{num_pat}*FLAG_P2;%net.I{3}=net.delta_t/net.tau_mb{3};
               %net=iter_dyn_lif(net,0.05/6*t_max/delta_t,0);
               %net=iter_dyn_lif(net,0.95/6*t_max/delta_t,1);                            
               net=iter_dyn_lif(net,1/6*t_max/delta_t,1);                            
            end;
        end;

    end;

    for num_pat=1:nb_pat
        net.I{1}=I1{num_pat}*FLAG_P1;net.I{2}=I2{num_pat}*FLAG_P2;%net.I{3}=net.delta_t/net.tau_mb{3};
        net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
    end;
    
    %pause;
    affiche_nips_lif(net,t_max,nb_pat,nb_app);
    affiche_cns_lif;
    %pause;
    
    for num_fig=1:5
        figure(num_fig);
        %eval(['print -dpng 071112-num-',num2str(num_res),'-nb_pat-',num2str(cpt_app),'-fig-',num2str(num_fig)]);
    end;
    
    
end;
