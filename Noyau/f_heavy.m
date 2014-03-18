% Fonction de transfert de type "Heaviside"
% Attention : la valeur en 0 est : 0 !!

function s=f_heavy(net,u,theta)

if net.delta_t<net.tau_r
    s=(u>=theta);
else
    %m=net.delta_t/net.tau_r;
    %s=((u>theta)&(u<=m*theta)).*(u./theta)+m*(u>m*theta);
    s=net.tau_m/net.tau_r*f_florian(u,theta,net.tau_m,net.tau_r);
    %s=((u>0)&(u<=m*theta)).*(u./theta)+m*(u>m*theta);
end;
        