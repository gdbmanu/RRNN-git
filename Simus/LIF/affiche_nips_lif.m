

function affiche_nips_lif(net,t_max,nb_pat,nb_app)

Delta_t=t_max/net.delta_t;
if nb_app>0
    t_start=Delta_t+1;
    t_stop=t_start+Delta_t*nb_app;
    t_final=(2+nb_app)*Delta_t;
else
    t_final=Delta_t;
    t_start=t_final;
    t_stop=t_final;
end;


figure(1);
clf;

DYN_Y=lissage_lif(net,net.DYN_S{1}(:,:));
trace=mean(DYN_Y);

subplot(3,1,1:2);
imagesc(net.delta_t:net.delta_t:t_final*net.delta_t,1:net.N(1),-DYN_Y);axis('xy');colormap(gray);
ylabel('Neuron #','FontSize', 14);
hold on;
h=plot([t_start*net.delta_t t_start*net.delta_t],[1 net.N(1)],'r:');
set(h,'linewidth',5);
h=plot([t_stop*net.delta_t t_stop*net.delta_t],[1 net.N(1)],'r:');
set(h,'linewidth',5);
Title('Individual Activities','FontSize', 14);

subplot(3,1,3);
axe=1:t_final;
plot(axe*net.delta_t,trace(axe)/5);
axis([net.delta_t t_final*net.delta_t 0 0.65]);
hold on;
h=plot([t_start*net.delta_t t_start*net.delta_t],[0 0.65],'r:');
set(h,'linewidth',5);
h=plot([t_stop*net.delta_t t_stop*net.delta_t],[0 0.65],'r:');
set(h,'linewidth',5);
title('Mean activity','FontSize', 14);
xlabel('Simulation time (ms)','FontSize', 14);



figure(2)
clf;

subplot(2,2,1)
[h,axe] = hist(nonzeros(net.J_ref{1}{1}),-0.9:0.03:0.9);
bar(axe,h);
axis([-0.9 0.9 0 12000])
title('Initial weights distribution','FontSize', 14)

subplot(2,2,2)
[h,axe] = hist(nonzeros(net.J{1}{1}-net.J_ref{1}{1}),-0.9:0.03:0.9);
bar(axe,h,'r');
axis([-0.9 0.9 0 12000])
title('\Delta W_{ij} distribution','FontSize', 14)

subplot(2,2,3)
[h,axe] = hist(nonzeros(net.J{1}{1}),-0.9:0.03:0.9);
bar(axe,h);
axis([-0.9 0.9 0 12000])
title('Final weights distribution','FontSize', 14)

subplot(2,2,4);
sigma=std(nonzeros(net.J{1}{1}-net.J_ref{1}{1}));
sigma_ini=std(nonzeros(net.J_ref{1}{1}));
sigma_final=std(nonzeros(net.J{1}{1}));
cr=sigma/sigma_ini*100;
cr=round(cr*10)/10;
sigma_ini_aff=round(sigma_ini*100)/100;
sigma_final_aff=round(sigma_final*100)/100;
sigma_aff=round(sigma*100)/100;
hold off;
plot([t_start*net.delta_t t_stop*net.delta_t],[sigma_ini,sigma_final]);
hold on;
plot([t_start*net.delta_t t_stop*net.delta_t],[sigma_ini,sigma_final],'d');
plot([t_start*net.delta_t t_stop*net.delta_t],[0,sigma],'r');
plot([t_start*net.delta_t t_stop*net.delta_t],[0,sigma],'rd');
text(t_start*net.delta_t-150,0.05,['\sigma (\Delta W_{ij}) = ',num2str(sigma_aff)]);
text(t_start*net.delta_t-150,sigma_ini+0.05,['\sigma (W_{ij}^0) = ',num2str(sigma_ini_aff)]);
text(t_stop*net.delta_t-500,sigma_final+0.05,['\sigma (W_{ij}^f) = ',num2str(sigma_final_aff)]);
title(['Change rate = ',num2str(cr),' %']);
axis([t_start*net.delta_t-250 t_stop*net.delta_t+250 -0.05 0.4]);
plot([t_start*net.delta_t-250 t_stop*net.delta_t+250],[0 0],':');
plot([t_start*net.delta_t-250 t_stop*net.delta_t+250],[sigma_ini sigma_ini],':');
plot([t_start*net.delta_t t_start*net.delta_t],[0 0.4],':');
plot([t_stop*net.delta_t t_stop*net.delta_t],[0 0.4],':');
xlabel('Simulation time (ms)');


figure(3)
clf;

subplot(2,2,1)
[h,axe] = hist(net.DYN_H{1}{1}(:,t_start),-30:4:30);
bar(axe,h);colormap('default');
title('Initial Input Repartition','FontSize', 14);
axis([-30 30 0 60]);

subplot(2,2,2)
[h,axe] = hist(net.DYN_H{1}{1}(:,t_stop),-30:4:30);
bar(axe,h);colormap('default');
title('Final Input Repartition','FontSize', 14);
axis([-30 30 0 60]);

subplot(2,2,3)
%plot(net.delta_t:net.delta_t:t_final*net.delta_t,net.DYN_H{1}{1}(101:200,1:t_final)','color',[0.9 0.9 0.9]);
hold on;
plot(net.delta_t:net.delta_t:t_final*net.delta_t,net.DYN_H{1}{1}(51:100,1:t_final)','color',[0.8 0.8 0.8]);
plot(net.delta_t:net.delta_t:t_final*net.delta_t,net.DYN_H{1}{1}(21:50,1:t_final)','color',[0.7 0.7 0.7]);
plot(net.delta_t:net.delta_t:t_final*net.delta_t,net.DYN_H{1}{1}(6:20,1:t_final)','color',[0.5 0.5 0.5]);
plot(net.delta_t:net.delta_t:t_final*net.delta_t,net.DYN_H{1}{1}(1:5,1:t_final)','k');
axis([0 t_final*net.delta_t -30 30]);
xlabel('Simulation time (ms)');
hold on;
h=plot([t_start*net.delta_t t_start*net.delta_t],[-30 30],'r:');
set(h,'linewidth',5);
h=plot([t_stop*net.delta_t t_stop*net.delta_t],[-30 30],'r:');
set(h,'linewidth',5);
Title('Individual Inputs','FontSize', 14);

subplot(2,2,4);
plot(net.delta_t:net.delta_t:t_final*net.delta_t,mean(net.DYN_H{1}{1}(:,1:t_final)));
hold on;
plot(net.delta_t:net.delta_t:t_final*net.delta_t,std(net.DYN_H{1}{1}(:,1:t_final)),'r');
axis([0 t_final*net.delta_t -5 15]);
plot(net.delta_t:net.delta_t:t_final*net.delta_t,mean(DYN_Y(:,1:t_final)),'k');
xlabel('Simulation time (ms)');
h=plot([t_start*net.delta_t t_start*net.delta_t],[-5 15],'r:');
set(h,'linewidth',5);
h=plot([t_stop*net.delta_t t_stop*net.delta_t],[-5 15],'r:');
set(h,'linewidth',5);
plot([0 t_final*net.delta_t],[0 0],':');
text(500,-3,'Input mean');
text(500,10,'Input dispersion');
text(500,3,'Mean activity');

figure(4);    
clf    

mem_corr=[];cpt_corr=0;
axe_t=[];

t_max=Delta_t+1;
%for t0=1:delta_t:t_max
t0=1;
    for delta_t=40:40:Delta_t-1040
        t=t0+delta_t;
        cpt_corr=cpt_corr+1;
        for tau=1:500
            c=corrcoef(trace(t:t+500),trace(t+tau:t+tau+500));
            mem_corr(cpt_corr,tau)=c(1,2);
        end;
        disp(t+500);
        axe_t=[axe_t,t+500];
    end;
%end;

if nb_app>0
    t_min=t_max+500
    t_max=t_min+Delta_t*nb_app-1000
    for t=t_min:40:(t_max-40)
        cpt_corr=cpt_corr+1;
        for tau=1:500
            c=corrcoef(trace(t:t+500),trace(t+tau:t+tau+500));
            mem_corr(cpt_corr,tau)=c(1,2);
        end;
        disp(t+500);
        axe_t=[axe_t,t+500];
    end;
    
  %t_min=t_max+500
  %t_max=t_min+(nb_pat-1)*2000
  %for t0=t_min:2000:t_max
  t0=t_max+500;
    for delta_t=40:40:Delta_t-1040
        t=t0+delta_t;
        cpt_corr=cpt_corr+1;
        for tau=1:500
            c=corrcoef(trace(t:t+500),trace(t+tau:t+tau+500));
            mem_corr(cpt_corr,tau)=c(1,2);
        end;
        disp(t+500);
        axe_t=[axe_t,t+500];
    end;
  %end;
end;


subplot(3,1,1:2);
axe_tau = 0.5:0.5:250;
h = imagesc(axe_t*net.delta_t,axe_tau,mem_corr');
axis('xy');
set(gca,'xtick',0);
title('A U T O C O R R E L O G R A M');
ylabel('\tau (ms)','FontSize', 14);
hold on;
h=plot([t_start*net.delta_t t_start*net.delta_t],[0 250],'r:');
set(h,'linewidth',5);
h=plot([t_stop*net.delta_t t_stop*net.delta_t],[0 250],'r:');
set(h,'linewidth',5);


subplot(3,1,3);
axe=1:t_final;
plot(axe*net.delta_t,trace(axe)/5);
axis([0.5 t_final*net.delta_t 0 0.65]);
hold on;
h=plot([t_start*net.delta_t t_start*net.delta_t],[0 0.65],'r:');
set(h,'linewidth',5);
h=plot([t_stop*net.delta_t t_stop*net.delta_t],[0 0.65],'r:');
set(h,'linewidth',5);
title('Mean activity','FontSize', 14);
xlabel('Simulation time (ms)','FontSize', 14);

if net.FLAG_DASHBOARD
    figure(5);
    clf;
    warning off;
    plot(net.dash.dim{1}(1,:),net.dash.dim{1}(2,:));
    hold on;
    plot(net.dash.dim{1}(1,:),net.dash.dim{1}(2,:),'d');
    plot(net.dash.dim{1}(1,:),max(net.dash.periode{1}),'k');
    plot(net.dash.dim{1}(1,:),max(net.dash.periode{1}),'kd');
    warning on;
    axis([0 t_final*net.delta_t 0 125]);
    for num_pat=1:nb_pat
        plot([1000*num_pat 1000*num_pat],[0 125],':');
    end;
    if nb_app>0
        h=plot([t_start*net.delta_t t_start*net.delta_t],[0 125],'r:');
        set(h,'linewidth',5);
        h=plot([t_stop*net.delta_t t_stop*net.delta_t],[0 125],'r:');
        set(h,'linewidth',5);
        for num_pat=1:nb_pat
            plot([t_stop*net.delta_t+1000*num_pat t_stop*net.delta_t+1000*num_pat],[0 125],':');
        end;
    end;
    xlabel('Simulation time (ms)','FontSize', 14);
    title('Network activity : Embedded Dimension (blue) and Internal Period (black)');
end;
