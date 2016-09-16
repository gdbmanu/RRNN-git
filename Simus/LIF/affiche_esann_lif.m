
clf

nbp_f = 400 + 10000 + 1000;
t_f = nbp_f * 0.5;

subplot(4,1,1)
DYN_Y = lissage_lif(net,net.DYN_S{1}); 
imagesc(0.5*(1:nbp_f),1:1000,-DYN_Y(:,1:nbp_f));
%hold on
for nbp = net.mem_P_tref
    text(nbp*0.5,1010,'*')
end
colormap(gray)
xlabel('Time(ms)')
axis('xy')
ylabel('Neuron #')
title('A C T I V I T Y')

subplot(4,1,2)
plot(0.5 * ([1,401:10400]),[net.dash.sum_J{1}{1}(1:200,1), net.dash.sum_J{1}{1}(1:200,:)]')
axis([1 t_f 0 0.4])
xlabel('Time(ms)')
title('W E I G H T S')

subplot(8,2,9)
%imagesc(0.5*(1:600),1:100,-net.mem_I1(:,1:600))
imagesc(0.5*(1:600),1:100,-net.DYN_I{1}(:,1:600))
hold on
for nbp = net.mem_P_tref
    plot([nbp*0.5:nbp*0.5+29],ones(30)*100,'r.')
end
axis('xy')
axis([0 300 .5 100.5])
ylabel('INPUT')
title('I N I T I A L')

subplot(8,2,10)
%imagesc(0.5*(9801:10400),1:100,-net.mem_I1(:,9801:10400))
imagesc(0.5*(9801:10400),1:100,-net.DYN_I{1}(:,9801:10400))
hold on
for nbp = net.mem_P_tref
    plot([nbp*0.5:nbp*0.5+29],ones(30)*100,'r.')
end
axis('xy')
axis([4900 5200 .5 100.5])
ylabel('Neuron #')
title('F I N A L')

subplot(8,2,[11 13])
plot(floor(find(net.DYN_S{1}(1:100,1:600)==1)/100)*0.5,mod(find(net.DYN_S{1}(1:100,1:600)==1),100),'.')
ylabel('ACTIVITY')

subplot(8,2,[12 14])
plot(floor(find(net.DYN_S{1}(1:100,9801:10400)==1)/100)*0.5 + 4900,mod(find(net.DYN_S{1}(1:100,9801:10400)==1),100),'.')
ylabel('Neuron #')

subplot(8,2,15)
plot(0.5*(1:600),mean(net.DYN_S{1}(:,1:600)) * 1000 / 0.5)
ylabel('Pop. freq. (Hz)')
xlabel('Time(ms)')
axis([0 300 0 100])

subplot(8,2,16)
plot(0.5*(9801:10400),mean(net.DYN_S{1}(:,9801:10400)) * 1000 / 0.5)
for nbp = net.mem_P_tref
    text(nbp*0.5,-50,'*')
end
xlabel('Time(ms)')
axis([4900 5200 0 100])


