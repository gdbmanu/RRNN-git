net=init_param_rrnn_lif(0.5);

net.FLAG_REV=0;%1;
%net.V_moins = -1;

net=init_systeme_lif(net,1);
net=init_dyn_lif(net);
net=iter_dyn_lif(net,500,0);

%net.J{1}{1}(:,1:200)=abs(net.J{1}{1}(:,1:200));
%net.J{1}{1}(:,201:400)=-abs(net.J{1}{1}(:,1:200));
%net=iter_dyn_lif(net,200,0);
