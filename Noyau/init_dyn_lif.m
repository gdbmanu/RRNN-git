% Initialisation de la dynamique
% Mod�le I&F (avril 2009)

function net=init_dyn_lif(net);

  if net.FLAG_STD == 1  
      net.J_eff = net.J;
  end;
  
  nb_pop=length(net.N);
    
  for p=1:nb_pop

      for q=1:nb_pop
          
          if net.connect(p,q)
            % Variable globale : net.delta_J (cellule sur p, q, tau)
            % memorise les modifications sur les poids
            % au cours de l'apprentissage
              
            net.trace{p}{q}=sparse(net.N(p),net.N(q));   

            if net.FLAG_REV   % pour reversal potentials
                net.rev{p}{q}=sparse(net.N(p),net.N(q));
                net.J_masque_plus{p}{q}=sparse(net.J{p}{q}>0);                   
                net.J_masque_moins{p}{q}=sparse(net.J{p}{q}<0);                   
            end;
            net.epsilon{p}{q}=zeros(net.N(p),net.N(q));
            net.epsilon_barre{p}{q}=0;
            
            %net.h{p}{q}=zeros(net.N(p),1); 
            net.h{p}{q}=zeros(net.N(p),net.tau_max(p)+1); 
            net.mem_h{p}{q}=zeros(net.N(p),1); 
            net.var_h{p}{q}=zeros(net.N(p),1); 
            net.m_h{p}{q}=zeros(net.N(p),1); 
            net.diff_h{p}{q}=zeros(net.N(p),1); 
                      
            net.F{p}{q}=zeros(net.N(p),net.tau_max(p)+1);
            
            % Pour f_bayes
            net.m_h_tire{p}{q}=5*ones(net.N(p),1);
            net.sigma_h_tire{p}{q}=ones(net.N(p),1);
            net.m_h_silent{p}{q}=zeros(net.N(p),1);
            net.sigma_h_silent{p}{q}=ones(net.N(p),1);
            
            net.bruit_h{p}{q}=zeros(net.N(p),1);
          
            % Variable globale : net.h (cellule sur p,q)    
            % Chaque vecteur, de taille net.N(p) memorise le champ local 
            % en provenance de la couche q
            % initialement nul
          
            net.DYN_H{p}{q}=[]; 
          
            % Variable globale : net.DYN_H (cellule sur p,q)
            % Chaque matrice, de taille net.N(p) x t, contient la dynamique
            % 
            % du champ local de la couche q vers la couche p
          end;
      end;
 
      % Pour f_bayes
      if net.FLAG_BAYES
        net.bayes{p}.p_tire=0.5*ones(net.N(p),1);
      
        net.bayes{p}.m_u_tire=5*ones(net.N(p),1);
        net.bayes{p}.var_u_tire=ones(net.N(p),1);
        net.bayes{p}.m_du_tire=ones(net.N(p),1);
        net.bayes{p}.var_du_tire=ones(net.N(p),1);
        net.bayes{p}.cov_tire=zeros(net.N(p),1);
      
        net.bayes{p}.m_u_silent=zeros(net.N(p),1);
        net.bayes{p}.var_u_silent=ones(net.N(p),1);
        net.bayes{p}.m_du_silent=-ones(net.N(p),1);
        net.bayes{p}.var_du_silent=ones(net.N(p),1);
        net.bayes{p}.cov_silent=zeros(net.N(p),1);
      
        net.mem_sum_h{p}=zeros(net.N(p),1);
        net.var_sum_h{p}=zeros(net.N(p),1);      
      end;
      
      % Pour reversal potentials
      net.norm_g{p}=ones(net.N(p),1);
     
      f0 = 20; %Hz
      net.S{p}=rand(net.N(p),net.tau_max(p)+1)< f0 / 1000 * net.delta_t;%1/mean(net.nbp_mb{p});% 1./mean(net.nbp_mb{1});%*(net.tau_max(p)+1);
      net.S_flat{p}=sparse(reshape(net.S{p}',1,net.N(p)*(net.tau_max(p)+1)));
      net.bool_refr{p}=zeros(net.N(p),net.tau_max(p)+1);
      net.U{p}=zeros(net.N(p),net.tau_max(p)+1);
      net.U_tire{p}=zeros(net.N(p),1);
      net.m_U{p}=zeros(net.N(p),1);   
      net.diff_U{p}=zeros(net.N(p),1);   
      net.V{p}=zeros(net.N(p),net.tau_max(p)+1);
      net.m_V{p}=zeros(net.N(p),1);   
      
      net.theta_ref{p}=net.theta{p};
      net.var_theta{p}=zeros(net.N(p),1);
      
      net.rnd_I{p}=0;
      net.buf_I{p}=0;
      
      % Variable globale : net.S (cellule sur p)
      % Chaque matrice, de taille net.N(p) x (net.tau_max(p)+1)
      % contient les tau_max valeurs passees du vecteur de l'activation de la couche
      % initialement : valeurs aleatoires uniformes entre 0 et 1 
      
      net.DYN_S{p}=[]; 
      
      % Variable globale : net.DYN_X (cellule sur p)
      % Chaque matrice, de taille net.N(p) x t, contient la dynamique 
      % du vecteur d'activation
      
      net.DYN_V{p}=[]; 
      
      net.DYN_U{p}=[]; 
      
      net.DYN_I{p}=[]; 
      
      % Variable globale : net.DYN_I (cellule sur p)
      % Chaque matrice, de taille net.N(p) x t, contient la dynamique 
      % du signal d'entree
      
      net.eta{p}=0*net.S{p};
      net.eta{p}(:,net.tau_max(p)+1)=net.S{p}(:,net.tau_max(p)+1).*net.theta{p};
      net.epsilon_i{p}=0*net.S{p};
      net.epsilon_i_flat{p}=reshape(net.epsilon_i{p}',1,net.N(p)*(net.tau_max(p)+1));
      net.mem_S{p}=0*net.S{p};%net.f_ref+0*net.x{p};
      
      % Variable globale : net.mem_x (cellule sur p)
      % Chaque vecteur, de taille net.N(p), contient la la valeur moyenne  
      % de l'activation
      if isfield(net,'FLAG_DASHBOARD')
          net.dash.dim{p}=[];
          net.dash.periode{p}=[];  
          for q=1:nb_pop
              %for i=1:1
                  net.dash.sum_J{p}{q}=[];
              %end;
            net.dash.mean_Delta_J{p}{q}=[];
            net.dash.std_Delta_J{p}{q}=[];
          end;
      end;

  end;
  
% Variable globale : net.t_abs
% temps absolu a partir de l'initialisation
  
net.t_abs=0;     
  
% Permet d'executer un script d'initialisation

if strcmp(net.script_init,'')==0  
    eval(net.script_init);
end;
 
disp('Initialisation de la dynamique');
