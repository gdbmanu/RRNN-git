
map.M = {zeros(31,31), zeros(31,31)};
map.Q = {zeros(31,31), zeros(31,31)};
map.C = {zeros(31,31), zeros(31,31)};


load simu_ECM_1403_J10_tau10_unif_H1;
map.M{1}(1:2,:) = mem.M{1}(1:2,:);
map.M{2}(1:2,:) = mem.M{2}(1:2,:);
map.Q{1}(1:2,:) = mem.Q{1}(1:2,:);
map.Q{2}(1:2,:) = mem.Q{2}(1:2,:);
map.C{1}(1:2,:) = mem.C{1}(1:2,:);
map.C{2}(1:2,:) = mem.C{2}(1:2,:);

for i = 2:30
    
    if i < 10
        s = ['simu_ECM_1403_J10_tau10_unif_H1_I0',num2str(i)];
    else
        s = ['simu_ECM_1403_J10_tau10_unif_H1_I',num2str(i)];
    end
    
    load(s)
    map.M{1}(i+1,:) = mem.M{1}(i+1,:);
    map.M{2}(i+1,:) = mem.M{2}(i+1,:);
    map.Q{1}(i+1,:) = mem.Q{1}(i+1,:);
    map.Q{2}(i+1,:) = mem.Q{2}(i+1,:);
    map.C{1}(i+1,:) = mem.C{1}(i+1,:);
    map.C{2}(i+1,:) = mem.C{2}(i+1,:);
    
end

figure(1)
imagesc(0:0.1:2.3,0:0.1:2.3, map.M{1}(1:24,1:24) * 1000); axis('square'); axis('xy')
title('Excitatory population average activity')
colorbar('Location', 'EastOutside')
xlabel('\sigma_W')
ylabel('I')

figure(2)
imagesc(0:0.1:2.3,0:0.1:2.3, map.M{2}(1:24,1:24) * 1000); axis('square'); axis('xy')
title('Inhibitory population average activity')
colorbar('Location', 'EastOutside')
xlabel('\sigma_W')
ylabel('I')

figure(3)
imagesc(0:0.1:2.3,0:0.1:2.3, map.Q{1}(1:24,1:24)); axis('square'); axis('xy')
title('Excitatory population (Q)')
colorbar('Location', 'EastOutside')
xlabel('\sigma_W')
ylabel('I')

figure(4)
imagesc(0:0.1:2.3,0:0.1:2.3, map.Q{2}(1:24,1:24)); axis('square'); axis('xy')
title('Inhibitory population (Q)')
colorbar('Location', 'EastOutside')
xlabel('\sigma_W')
ylabel('I')

figure(5)
imagesc(0:0.1:2.3,0:0.1:2.3, map.C{1}(1:24,1:24)); axis('square'); axis('xy')
title('Excitatory population (C)')
colorbar('Location', 'EastOutside')
xlabel('\sigma_W')
ylabel('I')

figure(6)
imagesc(0:0.1:2.3,0:0.1:2.3, map.C{2}(1:24,1:24)); axis('square'); axis('xy')
title('Inhibitory population (C)')
colorbar('Location', 'EastOutside')
xlabel('\sigma_W')
ylabel('I')
