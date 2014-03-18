
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%								                                 %
%	Iteration dans le temps     : Equations de champ moyen       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function net=iter_dyn_ECM_s(net,nbp,FLAG_CORR,prop_H);

nb_pop=length(net.J_barre);
h_max=3;

% pas de bruit

%net.sigma_W=0;

% Parametres neurone

delta_t=net.delta_t;%0.5;
tau_m=net.tau_m;%10;
tau_r=net.tau_r;%1;
gamma=1-delta_t/tau_m;

if isfield(net,'sigma_J_eff')
    mem_sigma=net.sigma_J;
    net.sigma_J=net.sigma_J_eff;
end;


BUF_MU=zeros(nb_pop,nbp);
BUF_NU=zeros(nb_pop,nbp);
BUF_DELTA=zeros(nb_pop,nbp);
BUF_MU_H=zeros(nb_pop,nbp);
BUF_NU_H=zeros(nb_pop,nbp);
BUF_DELTA_H=zeros(nb_pop,nbp);
BUF_M=zeros(nb_pop,nbp);
BUF_Q=zeros(nb_pop,nbp);

if FLAG_CORR
  BUF_C=zeros(nb_pop,nbp);
end;

sigma=1;
K = net.nbp_m;
H = prop_H * net.nbp_m;
%mem = max(H, net.nbp_d);
tau_ref = net.tau_min(1,1) + net.tau_moy(1,1);
mem = max(H+1, 2 * tau_ref - net.tau_min(1,1));
axe_tau = net.tau_min(1,1) : 2 * tau_ref - net.tau_min(1,1);

%if H == net.nbp_m
%    eps_kernel=ones(nb_pop,H);
%else
    rho = K / sum(gamma.^(0:(H-1)));
    eps_kernel = rho * ones(nb_pop,1) * (gamma.^(0:(H-1)));
%end


for l=1:K
    %eps_C(:,l)=1-sum(eps_kernel(:,1:l-1),2)./sum(eps_kernel,2); % donne la correlation selon l = k-k'
    eps_C(:,l)=sum(eps_kernel(:,l:H).^2,2)./sum(eps_kernel.^2,2);
end;


eps_I = zeros(net.nb_pop,1);
for p = 1:net.nb_pop
    eps_I(p) = net.I{p}/net.tau_m;
end

for t=1:nbp

    % CONSTANT DELAYS
	%mu = delta_t * (net.J_barre*net.ECM.m(:,net.nbp_d)+eps_I);   
	%nu = delta_t^2 * net.sigma_J.^2 * net.ECM.q(:,net.nbp_d);
    % UNIFORM DELAYS
	mu = delta_t * (net.J_barre*mean(net.ECM.m(:,axe_tau),2)+eps_I);   
	nu = delta_t^2 * net.sigma_J.^2 * mean(net.ECM.q(:,axe_tau), 2);
    
    if FLAG_CORR
        % CONSTANT DELAYS
        %Delta = min(delta_t^2 * net.sigma_J.^2 * net.ECM.C(:,net.nbp_d), nu);
        % UNIFORM DELAYS
        Delta = min(delta_t^2 * net.sigma_J.^2 * mean(net.ECM.C(:,axe_tau),2), nu);
    else
        Delta = nu;
    end;      
    
    sigma= min(sqrt(1.0*((nu-Delta).*((nu-Delta)>0))),0.99);

    net.ECM.mu=[mu,net.ECM.mu(:,1:mem)];
    net.ECM.nu=[nu,net.ECM.nu(:,1:mem)];
    net.ECM.Delta=[Delta,net.ECM.Delta(:,1:mem)];
       
    mu_h = sum(eps_kernel.*net.ECM.mu(:,1:H),2);
    
    %nu_h = sum(nbp_m*eps_kernel.*net.ECM.Delta(:,1:nbp_m) + 1/eps_kernel(0)*eps_kernel.*(net.ECM.nu(:,1:nbp_m)- net.ECM.Delta(:,1:nbp_m)),2);
    %nu_h = sum(eps_kernel.*sqrt(net.ECM.Delta(:,1:nbp_m)),2).^2+sum((eps_kernel.^2).*(net.ECM.nu(:,1:nbp_m)- net.ECM.Delta(:,1:nbp_m)),2);
    %nu_h = sum(eps_kernel.*sqrt(net.ECM.Delta(:,1:nbp_m)+(net.ECM.nu(:,1:nbp_m)- net.ECM.Delta(:,1:nbp_m))/2),2).^2+sum(eps_kernel.*(net.ECM.nu(:,1:nbp_m)- net.ECM.Delta(:,1:nbp_m))/2,2);
    %terme_1=net.ECM.Delta(:,1:nbp_m)+(1+gamma)*(net.ECM.nu(:,1:nbp_m)- net.ECM.Delta(:,1:nbp_m)).*eps_kernel/eps_kernel(1);
    %terme_2=(1+gamma)*(net.ECM.nu(:,1:nbp_m)- net.ECM.Delta(:,1:nbp_m)).*(1-(eps_kernel/eps_kernel(1)));
    %nu_h = sum(eps_kernel.*sqrt(terme_1),2).^2+sum(eps_kernel.*terme_2,2);

    %terme_1=(net.ECM.nu(:,1:nbp_m)-net.ECM.Delta(:,1:nbp_m));%.*eps_C;%+(1-eps_C)/nbp_m);
    %terme_2=(net.ECM.nu(:,1:nbp_m)-net.ECM.Delta(:,1:nbp_m)).*(1-eps_C);
    %var_h=sum(eps_C,2).*sum(eps_kernel.^2.*terme_1,2);%+sum((eps_kernel.^2).*terme_2,2);
    
    %nu_h = sum(eps_kernel.*(sqrt(net.ECM.Delta(:,1:nbp_m)+terme)),2).^2;
    
    
    %nu_h = sum(nbp_m/2*eps_kernel.^2.*net.ECM.nu(:,1:nbp_m),2);
    
    %nu_h = sum(eps_kernel.^2,2).*sum(eps_kernel.^2.*net.ECM.nu(:,1:nbp_m),2);
    
    nu_h = K * sum(eps_kernel.*net.ECM.nu(:,1:H),2);
    
    if FLAG_CORR

        Delta_h=K *  sum(eps_kernel.*net.ECM.Delta(:,1:H),2);
        
        %Delta_h=sum(eps_kernel.^2,2).*sum(eps_kernel.^2.*net.ECM.Delta(:,1:nbp_m),2);
        
        %Delta_h=sum(eps_kernel.*sqrt(net.ECM.Delta(:,1:nbp_m)),2).^2;
    else
        Delta_h=nu_h;
    end;
    

    %nu_h=var_h+Delta_h;
 
    
    net.ECM.mu_h=[mu_h,net.ECM.mu_h(:,1:mem)];
    net.ECM.nu_h=[nu_h+1e-20,net.ECM.nu_h(:,1:mem)];
    net.ECM.Delta_h=[Delta_h,net.ECM.Delta_h(:,1:mem)]; 
    net.ECM.Delta_h_K=[Delta_h,net.ECM.Delta_h_K(:,1:mem)]; 
    
    net.ECM.sigma=[sigma,net.ECM.sigma(:,1:mem)];

	BUF_MU(:,t)=mu;
	BUF_NU(:,t)=nu;
	BUF_DELTA(:,t)=nu;
    
	BUF_MU_H(:,t)=mu_h;
	BUF_NU_H(:,t)=nu_h;
	BUF_DELTA_H(:,t)=Delta_h;   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 		Integrales simples                            %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h_max = max(3, max(mu_h) + 2*sqrt(max(nu_h)))
    for p=1:nb_pop
        m(p)=quad(@f_m_ECM,-h_max,h_max,1e-6,[],net,p,sigma);
        q(p)=quad(@f_q_ECM,-h_max,h_max,1e-6,[],net,p,sigma);
    end;
    
    net.ECM.m=[m',net.ECM.m(:,1:mem)];
    net.ECM.q=[q',net.ECM.q(:,1:mem)];
    
	BUF_M(:,t)=m';
	BUF_Q(:,t)=q';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 		Integrale double                              %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if FLAG_CORR
        for p=1:nb_pop
            C(p)=dblquad(@f_C_ECM,-h_max,h_max,-h_max,h_max,1e-6,[],net,p,sigma(p),H);
        end;
    else
        C = q;
    end;
    net.ECM.C=[C',net.ECM.C(:,1:mem)];
    BUF_C(:,t)=C';
    
    disp(['****** t = ',num2str(t),'******']);
    disp(['mu = ',mat2str(mu)]);
    disp(['nu = ',mat2str(nu)]);
    disp(['Delta = ',mat2str(Delta)]);
    disp(['sigma = ',mat2str(sigma)]);
    disp(['mu_h = ',mat2str(mu_h)]);
    disp(['nu_h = ',mat2str(nu_h)]);
    disp(['Delta_h = ',mat2str(Delta_h)]);
    disp(['m = ',mat2str(m)]);
    disp(['q = ',mat2str(q)]);
    disp(['C = ',mat2str(C)]);
    
  	if rem(t,20)==0
	  disp(['t=',num2str(t)]);
    end;

end

net.ECM.DYN_MU=[net.ECM.DYN_MU BUF_MU];
net.ECM.DYN_NU=[net.ECM.DYN_NU BUF_NU];
net.ECM.DYN_DELTA=[net.ECM.DYN_DELTA BUF_DELTA];
net.ECM.DYN_MU_H=[net.ECM.DYN_MU_H BUF_MU_H];
net.ECM.DYN_NU_H=[net.ECM.DYN_NU_H BUF_NU_H];
net.ECM.DYN_DELTA_H=[net.ECM.DYN_DELTA_H BUF_DELTA_H];
net.ECM.DYN_M=[net.ECM.DYN_M BUF_M];
net.ECM.DYN_Q=[net.ECM.DYN_Q BUF_Q];

if FLAG_CORR
    net.ECM.DYN_C=[net.ECM.DYN_C BUF_C];
end;

if isfield(net,'sigma_J_eff')
    net.sigma_J=mem_sigma;
end;
