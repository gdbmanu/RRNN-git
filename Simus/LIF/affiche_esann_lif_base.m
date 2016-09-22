
clf

subplot(5,1,1)
imagesc(0.5*(1:600),1:100,-net.mem_I1(:,1:600))
%imagesc(0.5*(1:600),1:100,-net.DYN_I{1}(:,1:600))
hold on
for nbp = net.mem_P_tref
    plot([nbp*0.5:nbp*0.5+59],ones(60)*100,'r.')
end
axis('xy')
axis([0 300 .5 100.5])
ylabel('INPUT')
colormap(gray)

subplot(5,1,[2 3 4])
plot(floor(find(net.DYN_S{1}(1:100,1:600)==1)/100)*0.5,mod(find(net.DYN_S{1}(1:100,1:600)==1),100),'.')
axis([0 300 0 100])
ylabel('ACTIVITY')

subplot(5,1,5)
plot(0.5*(1:600),mean(net.DYN_S{1}(:,1:600)) * 1000 / 0.5)
ylabel('Pop. freq. (Hz)')
xlabel('Time(ms)')
axis([0 300 0 100])



