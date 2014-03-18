function y =f_stat(net,u,sigma) %(u,f_data)

%f_data=[f_data;1];
%f_data=[0;net.data.mu;1];
indice_u=round(u/net.data.prec_u+net.data.shift_u);
indice_max=round(net.data.shift_u+net.data.u_max/net.data.prec_u);

axe_u_inf=indice_u<=0;
axe_u_mid=(1+indice_u).*(indice_u>0).*(indice_u<=indice_max);
axe_u_plus=(2+indice_max).*(indice_u>indice_max);

%axe_u_inf=u<-15;
%axe_u_mid=round(u*2+31).*(u>=-15).*(u<=25);
%axe_u_plus=82.*(u>25);

axe_u=round(axe_u_inf+axe_u_mid+axe_u_plus);

j=round(sigma/net.data.prec_sig+net.data.shift_sig);

y=net.data.tab_stat(axe_u,j)';
%y=net.data.mu(axe_u);
%y=f_data(axe_u);
%y=y';

