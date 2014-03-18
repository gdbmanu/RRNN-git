% Initialisation de la matrice des poids 
% Mod�le I&F (avril 2009)

function net=init_J_lif(net)
  
  rand('seed',net.num_reseau);
  
  nb_pop=length(net.N);
  
  % Structure de la matrice des poids : 
  % net.J est une CELLULE qui contient (nb_pop X nb_pop) cellules
  % pour p et q fixes, net.J{p}{q} est une cellule qui contient tau
  % matrices
  % ainsi net.J{p}{q}{tau} donne les valeurs des liens qui relient la population q 
  % � la population p avec un d�lai de transmission tau
  % Pour gagner de l'espace, ces matrices sont cod�es sous forme de matrices creuses
  % Cette fonction vise donc � initialiser net.J en fonction des parametres
  % fournis
    
  for p=1:nb_pop
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % On definit les constantes de membrane %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    % 
    
    %net.tau_mb{p}=10;
    %k=4;
    %net.tau_mb{p} = net.mb_inf+exp(k*(1:net.N(p))'/net.N(p))/exp(k)*(net.mb_sup-net.mb_inf);
    %net.tau_mb{p} = random('unif',net.mb_inf,net.mb_sup,net.N(p),1);
    %net.tau_mb{p} = random('exp',10/net.delta_t,net.N(p),1);
    
    %net.gamma{p}=(1-net.delta_t./net.tau_mb{p});
    if strcmp(net.FLAG_GRAD,'unif')
        net.gamma{p} = random('unif',0.5,0.99,net.N(p),1);    
    elseif strcmp(net.FLAG_GRAD,'linear')
        net.gamma{p} = 0.9 + (1:net.N(p))'/net.N(p) * 0.09;
    else
        net.gamma{p} = 1-net.delta_t/net.tau_m_ref(p);
    end;    
    net.tau_mb{p} = net.delta_t./(1-net.gamma{p});
    
    net.nbp_mb{p} = floor(net.tau_mb{p}/net.delta_t);    
    
    %tau_mb_ref=mean(net.tau_mb{p});
    %mu=spdiags(net.tau_mb{p}/tau_mb_ref,0,net.N(p),net.N(p));
    %for q=1:net.nb_pop
    %    if net.densite(p,q)>0
    %        net.J{p}{q}=mu*net.J{p}{q};
    %    end;
    %end;
      
  end;
  
  for p=1:nb_pop
      
    net.tau_max(p)=0;     
      
    for q=1:nb_pop

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TIRAGE DES POIDS D'UNE CLASSE DE LIENS (q vers p) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Parametre du reseau : net.distr (cellule de nb_pop chaines de caracteres)
        % Donne la nature de la distribution des poids arrivant sur la population p
        % Deux valeurs possibles 'unif' ou 'gauss'
        % (par defaut toute autre chaine de caractere donne une distribution gaussienne!)
        
        % Parametre du reseau : net.densite (matrice nb_pop x nb_pop)
        % donne la densite de la matrice des poids, comprise entre � (pas de liens) et 1 (connectivit� compl�te)
        
        % Variable temporaire : sp_J
        % Pour p et q fix�s, la matrice sp_J (variable temporaire) donne les valeurs des bases des poids
        % issus du tirage. 
        
        % Variable temporaire : J_masque
        % J_masque est le "masque" binaire
        % de la matrice des poids. D�signe les liens modifiables.
        
        
            
        if net.densite(p,q)>0
            
            % Si net.densite==0, la matrice n'est pas d�finie
            % Remarque : une matrice peut avoir une densit� non nulle et tous ses liens � z�ro

            

            masque_field=1;
            masque_field_plus=1;
            masque_grad=1;                                  
                                    
            if net.grad(p,q)==1
                %masque_grad=diag(1:net.N(p))*ones(net.N(p),net.N(q))/net.N(p)*2; 
                %masque_grad=ones(net.N(p),net.N(q))*diag(1:net.N(q))/net.N(q)*2;
                %masque_grad=ones(net.N(p),net.N(q))*log(diag(1:net.N(q))/net.N(q)+1)*2;
                k=4;masque_grad=ones(net.N(p),net.N(q))*(exp(k*diag(1:net.N(q))/net.N(q))-exp(0))/(exp(k)-2);%/(exp(1)-exp(0))*2;
            end;
            if net.grad(p,q)==-1
                %masque_grad=diag(1:net.N(p))*ones(net.N(p),net.N(q))/net.N(p)*2; 
                %masque_grad=1+ones(net.N(p),net.N(q))*diag(net.N(q):-1:1)/net.N(q)*2;
                %masque_grad=1+ones(net.N(p),net.N(q))*log(diag(net.N(q):-1:1)/net.N(q)+1)*2;
                %k=4;masque_grad=1+ones(net.N(p),net.N(q))*(exp(k*diag(net.N(q):-1:1)/net.N(q))-exp(0))/(exp(k)-2);%/(exp(1)-exp(0))*2;
                masque_grad=ones(net.N(p),net.N(q))*diag(1./net.gamma{q});
            end;
                                                                    %%%%%%%%%%%%%%%%%%%%%%%%
            if net.dim(p,q)~=0                                      %%%  Avec  topologie %%%
                                                                    %%%%%%%%%%%%%%%%%%%%%%%%
                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Ajout d'une structure topologique de type CHAMP NEURONAL %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                % Parametre du reseau : net.rayon (matrice nb_pop x nb_pop)
                % Le biais depend du "rayon" du champ : net.rayon(p,q)
                % Le rayon doit etre strictement superieur � 0
                % S'il est compris entre 0 et 1, il donne la proportion de neurones
                % qui appartiennent au voisinage d'un neurone, sur la population totale
                % Un rayon infini correspond � une absence de topologie

                % Parametre du reseau : net.dim (matrice nb_pop x nb_pop)
                % La dimension peut etre 0 (pas de topologie)
                % 1 (champ de dimension 1) ou 2 (champ de dimension 2)
                % et est donn� par la variable net.dim(p,q)
           
                % Parametre du reseau : net.rayon_plus (matrice nb_pop x nb_pop)
                % defini uniquement en cas de distribution gaussienne
                % donne le rayon du voisinage d'excitation
                
                % Variable temporaire : masque_field 
                % La topologie est m�moris�e dans la matrice masque_field 
                % qui d�finit le biais � introduire sur les poids en fonction de la distance entre 
                % les neurones.
            
                % Variable temporaire : masque_field_plus
                % La variable masque_field_plus sert �galement � d�finir le champ,
                % uniquement dans le cas de la distribution gaussienne
                
                % Variable temporaire : rayon (pour alleger l'ecriture)               
                
                rayon=net.rayon(p,q);
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                %%% D�finition de la distance entre les neurones de la population %%%                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                
                % Variable temporaire : dist (matrice de taille net.N(p) x net.N(q))
                % La distance est comprise entre 0 et pi
                % (Topologie circulaire ou torique)
                                                
                                                        %%%%%%%%%%%%%%%%%
                if net.dim(p,q)==1                      %% DIMENSION 1 %%
                                                        %%%%%%%%%%%%%%%%%
                                                        
                    disp(['Voisinage de dimension 1 : p=',num2str(p),' q=',num2str(q)]);
                    
                    X_i=diag(1:net.N(p))*ones(net.N(p),net.N(q))/net.N(p); 
                                                        % "position" sur l'axe des X selon les lignes de la matrice des poids
                                                        % normalise entre 0 et 1
                    
                    X_j=ones(net.N(p),net.N(q))*diag(1:net.N(q))/net.N(q); 
                                                        % "position" sur l'axe des X selon les colonnes de la matrice des poids
                                                        % normalise entre 0 et 1
                   
                    dist=abs(X_i-X_j)*2*pi;             % LINEAIRE
                    %dist=min(abs(X_i-X_j),1-abs(X_i-X_j))*2*pi;  % CIRCULAIRE
                                                        % distance comme difference des positions
                                                        % renormalise sur [0,pi]
                                            
                                                        %%%%%%%%%%%%%%%
                else                                    % DIMENSION 2 % 
                                                        %%%%%%%%%%%%%%%
                                                        % (on ne va pas au dela!!)
                                                        
                    disp(['Voisinage de dimension 2 : p=',num2str(p),' q=',num2str(q)]);
                    
                    dg=mod(0:net.N(p)-1,round(sqrt(net.N(p))));
                    X_i=diag(dg)*ones(net.N(p),net.N(q))/sqrt(net.N(p));
                                                        % position sur l'axe des X selon les lignes de la matrice des poids
                                                        % periodique selon sqrt(net.N(p))
                                                        % normalise entre 0 et 1
                    
                    dg=mod(0:net.N(q)-1,round(sqrt(net.N(q))));
                    X_j=ones(net.N(p),net.N(q))*diag(dg)/sqrt(net.N(q));
                                                        % position sur l'axe des Y selon les colonnes de la matrice des poids
                                                        % periodique selon sqrt(net.N(q))
                                                        % normalise entre 0 et 1
                    
                    dist_X=min(abs(X_i-X_j),1-abs(X_i-X_j))*2*pi;
                                                        % distance sur X comme difference des positions
                                                        % renormalise sur [0,pi]

                    dg=floor((0:net.N(p)-1)/round(sqrt(net.N(p))));
                    Y_i=diag(dg)*ones(net.N(p),net.N(q))/sqrt(net.N(p));
                                                        % position sur l'axe des Y selon les lignes de la matrice des poids
                                                        % periodique selon sqrt(net.N(p))
                                                        % normalise entre 0 et 1
                                                        
                    dg=floor((0:net.N(q)-1)/round(sqrt(net.N(q))));
                    Y_j=ones(net.N(p),net.N(q))*diag(dg)/sqrt(net.N(q));
                                                        % position sur l'axe des Y selon les colonnes de la matrice des poids
                                                        % periodique selon sqrt(net.N(q))
                                                        % normalise entre 0 et 1                    
                    
                    dist_Y=min(abs(Y_i-Y_j),1-abs(Y_i-Y_j))*2*pi;
                                                        % distance sur Y comme difference des positions
                                                        % renormalise sur [0,pi]
                    
                    dist=(dist_X.^2+dist_Y.^2).^(1/2);  % distance euclidienne
                    
                    rayon=2*sqrt(rayon/pi);             % le rayon est ajust� pour que la proportion de neurones
                                                        % voisins reste inchangee avec la dimension 2
                                                        % si on pose S1=pi*(pi r')^2
                                                        % et S2=4*pi^2
                                                        % on doit avoir S1/S2=r
                end;
                
                % gaussien
                %if q~=1
                    masque_field=(sqrt(2*pi)/rayon)^net.dim(p,q).*exp(-(dist/rayon).^2/2);
                % exponentiel
                    %masque_field=(pi./rayon).*exp(-dist./rayon);
                % power law
                %else 
                %    masque_field=(dist+100*(dist==0)).^-(1-rayon);
                %end;
                
                % masque_field est une matrice calcul�e selon une gaussienne (1D ou 2D)
                               
                % masque_field=masque_field/erf(pi/(sqrt(2)*rayon))^net.dim(p,q); 
                % Ajustement pour rayon grand (voir article Natural Computing p.15)
               
                % masque_field=sparse(masque_field.*(dist<(rayon*pi)));
                % On tronque les valeurs faibles pour gagner de l'espace m�moire
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Pour distribution gaussienne uniquement %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(net.distr{p},'gauss')
                    rayon=net.rayon_plus(p,q);
                   
                    % idem plus haut
                    masque_field_plus=(sqrt(2*pi)/rayon)^net.dim(p,q).*exp(-(dist/rayon).^2/2);
                    %masque_field_plus=masque_field_plus/erf(pi/(sqrt(2)*rayon))^net.dim(p,q); 
                    %masque_field_plus=sparse(masque_field_plus.*(dist<(rayon*pi)));
                   
                end;
            end; % MASQUE_FIELD
            
            FLAG_FIELD=1;
             
                                                %%%%%%%%%%%%%%%%%%%%%%
            if net.densite(p,q)>=1              %%% Matrice pleine %%%
                                                %%%%%%%%%%%%%%%%%%%%%%
                                                
                                                % (On considere uniquement densite=1)
                                                
                disp(['Matrice pleine p=',num2str(p),' q=',num2str(q)]);               
                
                J_masque=1;                     % Ici le masque est scalaire (tous les liens de la population sont modifiables)
                                                
                                                %%%%%%%%%%%%%%%%%%%%%%
            else                                %%% Matrice creuse %%%
                                                %%%%%%%%%%%%%%%%%%%%%%
                                                
                disp(['Matrice creuse p=',num2str(p),' q=',num2str(q)]);
                
                if net.dim(p,q)==0
                    masque=rand(net.N(p),net.N(q))<=net.densite(p,q);
                else
                    if (max(max(net.densite(p,q)*masque_field))<=1)
                        disp ('OK la densité des poids dépend du rayon');
                        masque=rand(net.N(p),net.N(q))<=net.densite(p,q)*masque_field;
                        FLAG_FIELD=0;
                    else
                        disp ('Erreur : la densité est trop forte!!');
                        masque=rand(net.N(p),net.N(q))<=net.densite(p,q);
                    end;
                end;    
                               
                %sp_J=sparse(rand(net.N(p),net.N(q))*masque);     
                
                J_masque=spones(masque);          % Ici le masque correspond � la structure binaire de la matrice creuse
                                                % (existence ou absence de lien)
            end;
            
            
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% INSTANCIATION DE LA MATRICE DES POIDS %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Cette instanciation est r�alis�e en deux etapes
            
            %%%%%%%%%%%
            % Etape 1 %
            %%%%%%%%%%%
            
            % on definit la matrice des poids moyens J
            % a l'aide de net.J_barre(p,q), net.sigma_J(p,q), J_masque et masque_field (et masque_field_plus)
           
            % Variable temporaire : K 
            % Nombre moyen de liens afferents � la population p (receptrice)
            
            K=net.N(q)*net.densite(p,q);
            
            % Variable globale : net.J_norm (matrice nb_pop x nb_pop)
            % normalisation de la valeur d'alpha (parametre d'apprentissage)
            % en fonction du nombre de liens afferents
            
	        net.J_norm(p,q)=net.norm_alpha(p,q)/K;   

            % Variable temporaire : mean_J
            % matrice contenant la partie constante de la distribution des poids
            
            mu_s_param=net.J_barre(p,q)/K; 
            
            % La multiplication par J_masque permet d'assurer que mean_J est une matrice creuse, si besoin est.
            % mean_J contient donc des liens constants ou nuls
            
            % Variable temporaire : sigma_J
            % matrice contenant la partie variable de la distribution des poids
            
            sigma_s_param=net.sigma_J(p,q)/sqrt(K);
            
            
            % Variable temporaire : J
            % Matrice des valeurs des liens (sans prendre en compte les delais) 
            % definition de la matrice a partir de mean_s et sigma_s 
            % ecart-type normalis� par le nombre de liens aff�rents K
           
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(net.distr{p},'unif')      %%% Distribution uniforme %%%
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                b_inf=mu_s_param-sigma_s_param*sqrt(3);
                b_sup=mu_s_param+sigma_s_param*sqrt(3);
                J = random('unif',b_inf,b_sup,net.N(p),net.N(q));
                J = J.* J_masque;               
                if FLAG_FIELD 
                    J = J.* masque_field; 
                end;
                % Ajout de la deformation liee � la topologie
                
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif strcmp(net.distr{p},'chi2')  %%%% Distribution chi 2 %%%%
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                
                k = 2*(mu_s_param/sigma_s_param)^2;
                J = mu_s_param/k*random('chi2',k,net.N(p),net.N(q));
                J = J.* J_masque;               
                if FLAG_FIELD 
                    J = J.* masque_field; 
                end;
                
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif strcmp(net.distr{p},'exp')   % Distribution exponentielle 
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                
                k = abs(mu_s_param);
                J = sign(mu_s_param)*random('exp',k,net.N(p),net.N(q));
                J = J.* J_masque;               
                if FLAG_FIELD 
                    J = J.* masque_field; 
                end;

                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            else                                %% Distribution gaussienne %%
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                
                J = random('norm',mu_s_param,sigma_s_param,net.N(p),net.N(q));
                J = J.* J_masque;
                
                if FLAG_FIELD & net.dim(p,q)~=0  % Ajout de la deformation liee � la topologie
                    
                    %!!mean_J=mean_J.*(1+masque_field_plus-masque_field);
                    %!!std_J= abs(std_J).*(masque_field_plus-masque_field*net.k); 
                    J=mu_s_param.*(1+masque_field_plus-masque_field).*J_masque+abs(J).*(masque_field_plus-masque_field*2);
                    
                end;
            end;
            
            J = J .* masque_grad;
            
            net.J{p}{q}=J;            
            
            net.J_masque{p}{q}=J_masque;
            
            % Pour accelerer, on repasse au format matrice pleine si le remplissage est superieur � 30 %
            if (nnz(net.J{p}{q})/prod(size(net.J{p}{q}))>0.3) % Taux de remplissage �lev� (plus de 30%)
                net.J{p}{q}=full(net.J{p}{q});
                net.J_masque{p}{q}=full(net.J_masque{p}{q});
                disp(['Matrice pleine couche ',num2str(q),' vers couche ',num2str(p)]);
            end;
            
          
          
            %%%%%%%%%%%
            % Etape 2 %
            %%%%%%%%%%%
            
            % On definit les delais
            
            % Variable temporaire : J_tau
            % Matrice des delais
            
            % Parametres : net.tau_min, net.tau_moy  (matrices nb_pop x nb_pop)
            % definissent les delais
            
            if net.tau_moy(p,q)==0
                J_tau=net.tau_min(p,q)*ones(net.N(p),net.N(q));
            else
                if net.dim(p,q)==0                   
                    if net.FLAG_UNIF
                    % Uniforme
                        J_tau=net.tau_min(p,q)+round(net.tau_moy(p,q)*rand(net.N(p),net.N(q))/net.delta_t);
                    else
                    % Poisson
                        J_tau=net.tau_min(p,q)+poissrnd(net.tau_moy(p,q),net.N(p),net.N(q));
                    end
                else
                    J_tau=net.tau_min(p,q)+round(net.tau_moy(p,q)*dist./rayon);
                end;
            end;
                        
            net.tau{p}{q}=J_tau;%.*J_masque;                        
            net.tau_max(p)=max(net.tau_max(p),max(max(J_tau)));
            
        end;  % if net.densite...
     end; % for q ...    
  end; % for p ...
  
  for p=1:net.nb_pop
      for q=1:net.nb_pop
          if net.densite(p,q)>0
             net.tau_eff{p}{q}=net.tau{p}{q}+(net.tau_max(q)+1)*ones(net.N(p),1)*(0:net.N(q)-1);
          end;  % if net.densite...
          if net.densite(q,p)>0
             net.tau_max(p)=max(net.tau_max(p),max(max(net.tau{q}{p})));
          end;
      end; % for q ...    
  end; % for p ...

  % net.connect : donne la connectivit� effective � l'initialisation
  
  for p=1:nb_pop
      for q=1:nb_pop
          if net.densite(p,q)>0
             net.J{p}{q}=net.J{p}{q}*net.connect(p,q);
          end;
      end;
  end;    
  
  
  % Variable globale : net.J_ref (cellule sur p, q, tau)
  % garde en memoire les valeurs initiales de la matrice des poids  
  
  net.J_ref=net.J;

  for p=1:nb_pop
      for q=1:nb_pop
          if net.densite(p,q)>0
            net.J_ref{p}{q}=net.J{p}{q};
            net.Delta_J_plus{p}{q} = sparse(net.N(p),net.N(q));              
            net.Delta_J_moins{p}{q} = sparse(net.N(p),net.N(q));
          end;
      end;
  end;

