
period_1 = net.ENV_PERIOD / net.delta_t; %50; %60; %6 ; %300; %* net.N(1) / sum(net.DYN_S{1}(:,net.t_abs)) / net.delta_t;
%period_2 = 20; %

net.t_out = net.t_abs; %

%net.t_out = net.t_out + sum(net.S{1}(:,20)) / net.N(1);%  / net.delta_t;

if ((2 * pi * mod(net.t_out/period_1, 1)) > pi) & ((2 * pi * mod(net.t_out/period_1, 1)) < 3 * pi/2) & (net.FLAG_P_CHANGE == 1) & (net.FLAG_P_PAIR == 1)
    %disp(net.t_out)
    if rand() < net.ENV_RENEWAL_RATE %1/3
        %net.FLAG_P_NEW = 1;
        %net.P12 = randn(net.N(1),1); net.P22 = randn(net.N(2), 1);
        net.A1 =  randn(net.N(1),1); % 0.35;%
        net.phi1 = rand(net.N(1),1) * 2 * pi; % 0;%
        
        %for i =51:100                    % correlations
        %    net.A1(i) = net.A1(51);
        %    net.phi1(i) = net.phi1(51);
        %end;
        

        %net.A1(1:300) = net.mem_P{1}(1:300);
        %net.phi1(1:300) = net.mem_P{2}(1:300);
        
        net.A2 =  randn(net.N(2),1); % 0.35;%
        net.phi2 = rand(net.N(2),1) * 2 * pi; % 0;%
    else
        %net.REWARD = 1;
        disp(['Pat : nbp = ' , num2str(net.t_abs)]);
        net.mem_P_tref = [net.mem_P_tref, net.t_abs];
        %net.P12 = net.mem_P{2};net.P22 = net.mem_P{4};
        net.A1 =  net.mem_P{1};
        net.phi1 = net.mem_P{2};
        net.A2 =  net.mem_P{3};
        net.phi2 = net.mem_P{4};
        net.FLAG_P_PAIR = 0;
    end
    net.FLAG_P_CHANGE = 0;
elseif ((2 * pi * mod(net.t_out/period_1, 1)) > pi) & ((2 * pi * mod(net.t_out/period_1, 1)) < 3 * pi/2) & (net.FLAG_P_CHANGE == 1)& (net.FLAG_P_PAIR == 0)
    %net.A1 =  randn(net.N(1),1); % 0.35;%
    %net.phi1 = rand(net.N(1),1) * 2 * pi; % 0;%
    net.A1 =  net.mem_P{5};
    net.phi1 = net.mem_P{6};
    net.A2 =  net.mem_P{7};
    net.phi2 = net.mem_P{8};
    net.FLAG_P_CHANGE = 0;
    net.FLAG_P_PAIR = 1;
end;

if ((2 * pi * mod(net.t_out/period_1, 1)) > 3 * pi/2) & (net.FLAG_P_CHANGE == 0)
    %if net.FLAG_P_NEW
    %    net.FLAG_P_NEW = 0;
    %    net.P11 = randn(net.N(1),1); net.P21 = randn(net.N(2),1);
    %else
    %    net.P11 = net.mem_P{1}; net.P21 = net.mem_P{3};
    %end
    net.FLAG_P_CHANGE = 1;
    if net.FLAG_P_PAIR == 1
        %net.REWARD = 0;
    end
end;


net.I{1} =   1 + net.ENV_FLAG_NON_FLAT * net.A1 .* (cos(2 * pi * net.t_out / period_1 + net.phi1));
%net.I{1} = net.I{1}.* (net.I{1} > 0);
%net.I{1} = net.P11 * cos(2 * pi * net.t_out / period_1) + net.P12 * sin(2 * pi * net.t_out / period_1);
%net.I{1} = (1 + sin(2 * pi * net.t_out / period_2)) * net.I{1};


net.I{2} =   1 + ENV_FLAG_NON_FLAT * net.A2 .* (cos(2 * pi * net.t_out / period_1 + net.phi2));
%net.I{2} = net.I{2}.* (net.I{2} > 0);
%net.I{2} = net.P21 * cos(2 * pi * net.t_out / period_1) + net.P22 * sin(2 * pi * net.t_out / period_1); 
%net.I{2} = (1 + sin(2 * pi * net.t_out / period_2)) * net.I{2};

net.mem_t_out = [net.mem_t_out, net.t_out];

%net.REWARD = 10 * (0.1 - sum(net.epsilon_i{1}(:,1)) / net.N(1) / 5);

net.mem_I1 = [net.mem_I1, net.I{1}];
net.mem_I2 = [net.mem_I2, net.I{2}];

%disp('hello')


