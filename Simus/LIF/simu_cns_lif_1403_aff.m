%load simu_cns_lif_1403_seed60_tau10_period60ms_RENEW066_noSTDP.mat

DYN_Y = lissage_lif(net,net.DYN_S{1}); 

win = 120;

cpt_0 = 0;
mem_0 = zeros(1,win);
mem_Y_0 = zeros(net.N(1), win);
mem_I_0 = zeros(net.N(1), win);

cpt_1 = 0;
mem_1 = zeros(1,win);
mem_Y_1 = zeros(net.N(1), win);
mem_I_1 = zeros(net.N(1), win);


for t = 5011 : 60: (net.t_abs - win) % 31 : 60: (net.t_abs - win)
%for t = 21 : 120: (net.t_abs - win - 10000)
    %if ismember(t,net.mem_P_tref)
    %if ismember(t + 10000, net.mem_P_tref)
    if ismember(t, net.mem_P_tref)
        cpt_1 = cpt_1 + 1;
        mem_1 = mem_1 + mean(net.DYN_S{1}(:,t:t + win - 1)) * 2 * 1000;
        mem_Y_1 = mem_Y_1 + DYN_Y(:, t:t + win - 1);
        mem_I_1 = mem_I_1 + net.DYN_I{1}(:, t:t + win - 1);
    else
        cpt_0 = cpt_0 + 1;
        mem_0 = mem_0 + mean(net.DYN_S{1}(:,t:t + win - 1)) * 2 * 1000;
        mem_Y_0 = mem_Y_0 + DYN_Y(:, t:t + win - 1);
        mem_I_0 = mem_I_0 + net.DYN_I{1}(:, t:t + win - 1);
    end
end

figure(1)
clf
subplot(3,2,1);
imagesc((1:win)/2,1:100, -mem_I_1(1:100,:)/cpt_1,[-2.5 0])
colormap('gray')
axis('xy')
set(gca,'xtick',[])
ylabel('Neuron #')
title('REPEATED STIM.')
subplot(3,2,3);
imagesc((1:win)/2,1:100, -mem_Y_1(1:100,:)/cpt_1,[-2.5 0])
colormap('gray')
axis('xy')
set(gca,'xtick',[])
ylabel('Neuron #')
subplot(3,2,5)
plot((1:win)/2, mem_1 / cpt_1);
axis([1 win/2 0 60])
xlabel('Time (ms)')


subplot(3,2,2);
imagesc((1:win)/2,1:100, -mem_I_0(1:100,:)/cpt_0,[-2.5 0])
colormap('gray')
axis('xy')
set(gca,'xtick',[])
ylabel('Neuron #')
title('RANDOM STIM.')
subplot(3,2,4)
imagesc((1:win)/2,1:100, -mem_Y_0(1:100,:)/cpt_0,[-2.5 0])
colormap('gray')
axis('xy')
set(gca,'xtick',[])
ylabel('Neuron #')
subplot(3,2,6)
plot((1:win)/2, mem_0 / cpt_0);
axis([1 win/2 0 60])
xlabel('Time (ms)')







        
        
        