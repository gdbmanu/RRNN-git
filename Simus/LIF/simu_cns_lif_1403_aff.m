%load simu_cns_lif_1403_seed60_tau10_period60ms_RENEW066_noSTDP.mat

DYN_Y = lissage_lif(net,net.DYN_S{1}); 

cpt_0 = 0;
mem_0 = zeros(1,240);
mem_Y_0 = zeros(net.N(1), 240);
mem_I_0 = zeros(net.N(1), 240);

cpt_1 = 0;
mem_1 = zeros(1,240);
mem_Y_1 = zeros(net.N(1), 240);
mem_I_1 = zeros(net.N(1), 240);


%for t = 61 : 120: (net.t_abs - 240)
for t = 21 : 120: (net.t_abs - 240 - 10000)
    %if ismember(t,net.mem_P_tref)
    if ismember(t + 10000, net.mem_P_tref)
        cpt_1 = cpt_1 + 1;
        mem_1 = mem_1 + mean(net.DYN_S{1}(:,t:t + 239)) * 2 * 1000;
        mem_Y_1 = mem_Y_1 + DYN_Y(:, t:t + 239);
        mem_I_1 = mem_I_1 + net.DYN_I{1}(:, t:t + 239);
    else
        cpt_0 = cpt_0 + 1;
        mem_0 = mem_0 + mean(net.DYN_S{1}(:,t:t + 239)) * 2 * 1000;
        mem_Y_0 = mem_Y_0 + DYN_Y(:, t:t + 239);
        mem_I_0 = mem_I_0 + net.DYN_I{1}(:, t:t + 239);
    end
end

figure(1)
clf
subplot(3,2,[1 3]);
imagesc((1:240)/2,1:100, -mem_Y_1(1:100,:)/cpt_1,[-2.5 0])
colormap('gray')
axis('xy')
set(gca,'xtick',[])
ylabel('Neuron #')
title('REPEATED STIM.')
subplot(3,2,5)
plot((1:240)/2, mem_1 / cpt_1);
axis([1 120 0 60])
xlabel('Time (ms)')


subplot(3,2,[2 4]);
imagesc((1:240)/2,1:100, -mem_Y_0(1:100,:)/cpt_0,[-2.5 0])
colormap('gray')
axis('xy')
set(gca,'xtick',[])
ylabel('Neuron #')
title('RANDOM STIM.')
subplot(3,2,6)
plot((1:240)/2, mem_0 / cpt_0);
axis([1 120 0 60])
xlabel('Time (ms)')

figure(3)
clf








        
        
        