
clear all;

nb_pat=3%7;%
t_max=1000*nb_pat;%2000;   % ms
nb_app=1;%10%


cpt=1;
mem_net={};

delta_t=0.5;
net=init_param_rrnn_lif(delta_t);
net.app='stdp';
%m=net.theta_barre(1)/(1-exp(-net.tau_r/net.tau_m));
net.alpha=1;%0.5;
net.FLAG_DASHBOARD=1;
net.FLAG_STD=0;
net.FLAG_REV=1;
net.FLAG_SCALING=1;
net.FLAG_STDP_EPS=0;
net.FLAG_DELTA=1;
net.FLAG_BAYES=0;
net.FLAG_THETA=0;


randn('seed',1);
for i=1:10
    buf_I{i}=randn(net.N(1),1)*2*net.delta_t/10;%*0.15;%*2;
    
    %buf_I{i}=-ones(net.N(1),1)*0.5;
    %buf_I{i}=zeros(net.N(1),1);
    
    %buf_I{i}=2*((rand(net.N(1),1)<0.5)-0.5);
    %buf_I{i}=rand(net.N(1),1)<0.5;
    
    %nb_in=20;
    %b_inf=(i-1)*nb_in+1;b_sup=i*nb_in;
    %buf_I{i}(b_inf:b_sup)=ones(nb_in,1);
    
    %motif = rand(net.N(1),1)<0.1;    
    %buf_I{i}=buf_I{i}+motif*2.5;
end;
%buf=buf_I{2};buf_I{2}=buf_I{nb_pat};buf_I{nb_pat}=buf;

for i=1:nb_pat
    I{i}=buf_I{i};
end;

for i=1:9
    %buf_I{i}=randn(net.N(1),1);%*2;
    buf_I{i}=-ones(net.N(1),1);
    %buf_I{i}=(rand(net.N(1),1)<0.5);
    nb_in=20;
    b_inf=(i-1)*nb_in+11;b_sup=i*nb_in+10;
    buf_I{i}(b_inf:b_sup)=ones(nb_in,1);
end;
for i=1:nb_pat-1
    I_test{i}=buf_I{i};
end;

num_res=2;    
net=init_systeme_lif(net,num_res);

N2=net.N(1)/2;
net.J{1}{1}(:,1:N2)=abs(net.J{1}{1}(:,1:N2));
net.J{1}{1}(:,(N2+1):net.N(1))=-abs(net.J{1}{1}(:,(N2+1):net.N(1)));
net.J_ref{1}{1}=net.J{1}{1};
net.masque_app=[ones(N2),zeros(N2);zeros(N2),zeros(N2)];
                

%net.I{3}=0.5*(1-net.gamma{3});

for cpt_app=1:1 %nb_pat

    net=init_dyn_lif(net);  
    
%%% INPUT
    
    for i=1:nb_pat
        net.I{1}=I{i};
        net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
    end;


    if nb_app>0
        %net=iter_dyn_lif(net,t_max*nb_app/delta_t,1);
        for i=1:nb_app
            %net.REWARD=0;
            net.I{1}=I{3};
            net=iter_dyn_lif(net,0.02*t_max/delta_t,0);
            net=iter_dyn_lif(net,0.18*t_max/delta_t,1);
            %net.REWARD=1;
            %net=iter_dyn_lif(net,1,1);net.REWARD=0;%
            %net=iter_dyn_lif(net,0.1*t_max/delta_t,0);
            %net.REWARD=1;
            %net=iter_dyn_lif(net,1,1);net.REWARD=0;%
            %net=iter_dyn_lif(net,0.1*t_max/delta_t,0);
            %net.REWARD=1;
            %net=iter_dyn_lif(net,1,1);net.REWARD=0;%
            %net=iter_dyn_lif(net,0.1*t_max/delta_t,0);
            %net.REWARD=1;
            %net=iter_dyn_lif(net,1,1);net.REWARD=0;%
            %net=iter_dyn_lif(net,0.1*t_max/delta_t,0);
            %net.REWARD=1;
            %net=iter_dyn_lif(net,1,1);net.REWARD=0;%
                    
            if 1
            %net=iter_dyn_lif(net,0.4*t_max/delta_t,0);           
             net.I{1}=I{3};net.I{2}=0;%I{3};
             net=iter_dyn_lif(net,0.02*t_max/delta_t,0);
             net=iter_dyn_lif(net,0.18*t_max/delta_t,1);
%             net.REWARD=1;
%             net=iter_dyn_lif(net,1,1);net=iter_dyn_lif(net,1,0);
%             net.REWARD=0;
             net.I{1}=I{3};net.I{2}=0;%I{3};
             net=iter_dyn_lif(net,0.02*t_max/delta_t,0);
             net=iter_dyn_lif(net,0.18*t_max/delta_t,1);
%             net.REWARD=1;
%             net=iter_dyn_lif(net,1,1);net=iter_dyn_lif(net,1,0);
%             net.REWARD=0;
             net.I{1}=I{3};net.I{2}=0;%I{3};
             net=iter_dyn_lif(net,0.02*t_max/delta_t,0);
             net=iter_dyn_lif(net,0.18*t_max/delta_t,1);
%             net.REWARD=1;
%             net=iter_dyn_lif(net,1,1);net=iter_dyn_lif(net,1,0);
%             net.REWARD=0;
             net.I{1}=I{3};net.I{2}=0;%I{3};
             net=iter_dyn_lif(net,0.02*t_max/delta_t,0);
             net=iter_dyn_lif(net,0.18*t_max/delta_t,1);
%             net.REWARD=1;
%             net=iter_dyn_lif(net,1,1);net=iter_dyn_lif(net,1,0);
%             net.REWARD=0; 
            end;
        end;
%        net=iter_dyn_lif(net,0.05*t_max*nb_app/delta_t,0);         
%         net.I{1}=I{5};
%         net=iter_dyn_lif(net,0.1*t_max*nb_app/delta_t,0);
%         net=iter_dyn_lif(net,0.15*t_max*nb_app/delta_t,1);
%         net=iter_dyn_lif(net,0*t_max*nb_app/delta_t,0);
%         net.I{1}=I{3};
%         net=iter_dyn_lif(net,0.1*t_max*nb_app/delta_t,0);
%         net=iter_dyn_lif(net,0.15*t_max*nb_app/delta_t,1);
%         net=iter_dyn_lif(net,0*t_max*nb_app/delta_t,0);
%         net.I{1}=I{1};
%         net=iter_dyn_lif(net,0.1*t_max*nb_app/delta_t,0);
%         net=iter_dyn_lif(net,0.15*t_max*nb_app/delta_t,1);
%         net=iter_dyn_lif(net,0*t_max*nb_app/delta_t,0);
        %save 090223_simu_nips;
    end;

    for i=1:nb_pat
        net.I{1}=I{i};
        net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
    end;
        %save 090223_simu_nips;

%    net.I{1}=0;
%    net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
%    net.I{1}=I{3};
%    net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
%    net.I{1}=0;
%    net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
%    net.I{1}=I{8};
%    net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
%    net.I{1}=0;
%    net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);

%    for i=1:nb_pat-1
%        net.I{1}=I_test{i};
%        net=iter_dyn_lif(net,t_max/(nb_pat*delta_t),0);
%    end;
    

    %save 081126;
    %pause;
    affiche_nips_lif(net,t_max,nb_pat,nb_app);
    %pause;
    
    for num_fig=1:5
        figure(num_fig);
        %eval(['print -dpng 071112-num-',num2str(num_res),'-nb_pat-',num2str(cpt_app),'-fig-',num2str(num_fig)]);
    end;
    %pause;

    %eval(['save 071011-',num2str(cpt_app),' net']);
    %net=clear_dyn_lif(net);
    
    
end;
%     
% DYN_Y=lissage(net,net.dyn_lif{1}(:,:));
% 
% figure(cpt)
% clf;
% 
% subplot(4,2,[1;3;5]);
% imagesc(-net.DYN_I{1});colormap(gray);
% 
% subplot(4,2,[2;4;6]);
% imagesc(-DYN_Y);colormap(gray);
% subplot(4,4,[13:15]);
% nb_iter=nb_pat/2+nb_app*nb_pat+nb_pat/2;
% %plot(100+delta_t:delta_t:t_max*nb_iter-100,mean(DYN_Y(:,100/delta_t+1:t_max*nb_iter/delta_t-100/delta_t))*net.tau_r/net.tau_m);
% subplot(4,4,16);
% %plot2(mean(DYN_Y(:,100/delta_t:t_max*nb_iter/delta_t))*net.tau_r/net.tau_m);
% mem_net{cpt}=net;
% cpt=cpt+1;
% 
% %net=iter_dyn_lif(net,1000,0);
% 
% % net.sigma_theta=1;%0.3/0.4;
% % net=init_dyn_ECM(net);
% % net=iter_dyn_ECM_s(net,t_max/delta_t*nb_pat,0);
% % for i=1:2
% %     figure(i);
% %     subplot(4,4,[13:15]);
% %     hold on;
% %     plot(delta_t:delta_t:t_max*nb_pat,net.ECM.DYN_M*net.tau_r/net.tau_m,'r');
% % end;
% 
% 
% trace=mean(DYN_Y);
% 
% nb_iter=nb_pat*2+nb_app;
% 
% mem_corr=[];%zeros(nb_iter*50,500);
% axe_t=[];
% cpt=0;
% 
% t_max=(nb_pat-1)*2000+1
% for t0=1:2000:t_max
%     for delta_t=20:20:980
%         t=t0+delta_t;
%         cpt=cpt+1;
%         for tau=1:500
%             c=corrcoef(trace(t:t+500),trace(t+tau:t+tau+500));
%             mem_corr(cpt,tau)=c(1,2);
%         end;
%         disp(t+500);
%         axe_t=[axe_t,t+500];
%     end;
% end;
% 
% if nb_app>0
%     t_min=t_max+1500
%     t_max=t_min+2000*nb_app
%     for t=t_min:40:(t_max-40)
%         cpt=cpt+1;
%         for tau=1:500
%             c=corrcoef(trace(t:t+500),trace(t+tau:t+tau+500));
%             mem_corr(cpt,tau)=c(1,2);
%         end;
%         disp(t+500);
%         axe_t=[axe_t,t+500];
%     end;
%     
%   t_min=t_max+500
%   t_max=t_min+(nb_pat-1)*2000
%   for t0=t_min:2000:t_max
%     for delta_t=0:20:980
%         t=t0+delta_t;
%         cpt=cpt+1;
%         for tau=1:500
%             c=corrcoef(trace(t:t+500),trace(t+tau:t+tau+500));
%             mem_corr(cpt,tau)=c(1,2);
%         end;
%         disp(t+500);
%         axe_t=[axe_t,t+500];
%     end;
%   end;
% end;
%    
% figure(2);    
% clf    
% 
% subplot(3,1,1:2);
% axe_tau = 0.5:0.5:250;
% h = imagesc(axe_t*net.delta_t,axe_tau,mem_corr');
% axis('xy');
% set(gca,'xtick',0);
% title('A U T O C O R R E L O G R A M');
% ylabel('\tau (ms)','FontSize', 14);
% hold on;
% h = plot([500 500],[0 250],':r');
% set(h,'linewidth',5);
% 
% subplot(3,1,3);
% axe=501:5500;
% plot(axe*net.delta_t,trace(axe)/5);
% axis([250 2750 0.25 0.45]);
% hold on;
% h = plot([500 500],[0.1 0.4],':r');
% set(h,'linewidth',5);
% title('Mean activity','FontSize', 14);
% xlabel('Simulation time (ms)','FontSize', 14);
% 
% 
% %DYN_Y=DYN_Y(:,1001:25000);
% %trace=mean(DYN_Y);
