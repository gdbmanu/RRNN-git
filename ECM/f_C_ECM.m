
function y=f_C_ECM(h1,h2,net,p,sigma,H)

H = H +1;

membre_1=sqrt(abs(net.ECM.nu_h(p,1).*net.ECM.nu_h(p,H)-net.ECM.Delta_h_K(p,1).^2)./net.ECM.nu_h(p,H));
membre_2=net.ECM.Delta_h_K(p,1)./sqrt(net.ECM.nu_h(p,H));

% F_STAT
f1=f_stat(net,membre_1*h1+membre_2*h2+net.ECM.mu_h(p,1),sigma);%*net.delta_t/net.tau_r;
f2=f_stat(net,sqrt(net.ECM.nu_h(p,H))*h2+net.ECM.mu_h(p,H),sigma);%*net.delta_t/net.tau_r;   

% APPROX TANH
%f1=f_tanh_01(membre_1*h1+membre_2*h2+net.ECM.mu_h(p,1)-6/net.tau_r,0.25*net.tau_r);%*net.delta_t/net.tau_r;
%f2=f_tanh_01(sqrt(net.ECM.nu_h(p,net.nbp_m))*h2+net.ECM.mu_h(p,net.nbp_m)-6/net.tau_r,0.25*net.tau_r);%*net.delta_t/net.tau_r;

% F_FLORIAN
%f1=f_florian(membre_1*h1+membre_2*h2+net.ECM.mu_h(p,1),1,net.tau_m,net.tau_r)*net.delta_t/net.tau_r;
%f2=f_florian(sqrt(net.ECM.nu_h(p,net.nbp_m))*h2+net.ECM.mu_h(p,net.nbp_m),1,net.tau_m,net.tau_r)*net.delta_t/net.tau_r;

y=1/(2*pi)*exp(-h1.^2/2).*exp(-h2.^2/2).*f1.*f2;

