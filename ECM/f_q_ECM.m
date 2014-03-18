
function y=f_q_ECM(h,net,p,sigma)



% F_STAT
y=1/sqrt(2*pi)*exp(-h.^2/2).*(f_stat(net,(sqrt(net.ECM.nu_h(p,1)).*h)+net.ECM.mu_h(p,1),sigma(p))).^2;%*net.delta_t/net.tau_r).^2;

% APPROX TANH
%y=1/sqrt(2*pi)*exp(-h.^2/2).*(f_tanh_01(sqrt(net.ECM.nu_h(p,1)).*h+net.ECM.mu_h(p,1)-6/net.tau_r,0.25*net.tau_r)).^2;%*net.delta_t/net.tau_r).^2;


% F_FLORIAN
%y=1/sqrt(2*pi)*exp(-h.^2/2).*(f_florian((sqrt(net.ECM.nu_h(p,1)).*h)+net.ECM.mu_h(p,1),1,net.tau_m,net.tau_r)*net.delta_t/net.tau_r).^2;
