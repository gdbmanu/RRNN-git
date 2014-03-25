
%net.P11 = randn(net.N(1),1);%random('chi2',0.5,net.N(1),1); %randn(net.N(1),1); 
%net.P12 = randn(net.N(1),1);%random('chi2',0.5,net.N(1),1); %randn(net.N(1),1);
%net.P21 = randn(net.N(2),1); %random('chi2',0.5,net.N(2),1);%randn(net.N(2),1); 
%net.P22 = randn(net.N(2),1);%random('chi2',0.5,net.N(2),1);%randn(net.N(2), 1);

net.A1 =    randn(net.N(1),1); % 
net.phi1 =  rand(net.N(1),1) * 2 * pi; % 
net.A2 =    randn(net.N(2),1); % 
net.phi2 =  rand(net.N(2),1) * 2 * pi; % 

%net.mem_P = {net.P11,net.P12,net.P21,net.P22};
net.mem_P = {net.A1,net.phi1,net.A2,net.phi2};
% net.net.A1 = net.mem_P{1}; net.phi1 = net.mem_P{2}; net.net.A2 = net.mem_P{3}; net.phi2 = net.mem_P{4};

net.t_out = 0;
net.mem_t_out = [];
net.mem_I1 = [];
net.mem_I2 = [];
net.mem_P_tref = [];

net.FLAG_P_CHANGE = 1;
net.FLAG_P_NEW = 0;