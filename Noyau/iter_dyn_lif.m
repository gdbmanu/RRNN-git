
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                    %%%
%%%%                                                                    %
%%%%    Iteration de la dynamique d'un RRNN sur nbp pas de temps
%%%%    %%%%
%%%%                                                                    %%%%
%%%%    Mod�le I&F (avril 2009)
%%%%    %%%%
%%%%                                                                    %%%
%%%%                                                                    %
%%%%    E. Dauc�ne
%%%%    %%%%
%%%%                                                                    %%%
%%%%                                                                    %
%%%%    juin 2004
%%%%    %%%%
%%%%                                                                    %%%
%%%%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametre d'entree : nbp (nombre de pas de temps) %
% Ce parametre determine la duree de simulation     %
%                                                   %
% Parametre : FLAG_APP (booleen)                    %
% indique si l'apprentissage est activ�             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function net=iter_dyn_lif(net,nbp,FLAG_APP)


TEST=0;
for p=1:net.nb_pop
    net.m{p}=0;
end;

%N=sum(net.N);

% Pour f_bayes
tau_bayes=100;%net.tau_z;%100;
delta=1-net.delta_t/tau_bayes;

%m=net.theta_barre(1)/(1-exp(-net.tau_r/net.tau_m));
%s=1.6;%epsilon_max=m/2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Activation des entrees sur chaque neurone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parametre du reseau : net.flag_I
% booleen qui active ou desactive le signal sur chaque couche
%
% Variable temporaire : flag_in (cellule)
% flag_in{p} est un vecteur de booleens de taille N(p) qui active
% (ou desactive) l'entree sur chaque neurone

for p=1:net.nb_pop
    flag_in{p}=net.flag_I(p).*ones(net.N(p),1);
    if net.t_abs > 0 & isfield(net,'FLAG_CLEAR')==0
        I_in{p}=net.DYN_I{p}(:,net.t_abs);
    else
        I_in{p}=sparse(net.N(p),1);
    end;
end;

[buf,tau_I]=size(net.I{1});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation du parametre d'apprentissage (normalisation comprise) %
% sur chaque classe de liens                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parametre : net.alpha
% parametre d'apprentissage
%
% Variable globale : net.J_norm (matrice nb_pop x nb_pop)
% Normalisation de la valeur d'alpha (parametre d'apprentissage)
% en fonction du nombre de liens afferents
%
% Variable temporaire : alpha_eff (matrice nb_pop x nb_pop)
% d�finit la valeur effective du parametre d'apprentissage surchaque classe de liens
% en fonction du parametre global net.alpha

if FLAG_APP==1
    for p=1:net.nb_pop
        for q=1:net.nb_pop
            alpha_eff(p,q)=0;
            K=net.N(q)*net.densite(p,q);
            if K>0
                net.J_norm(p,q)=net.norm_alpha(p,q)/K;
                alpha_eff(p,q)=net.alpha*net.J_norm(p,q);
            end;
        end;
    end;
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On instancie ci-dessous differentes commandes qui dependent des parametres initiaux du reseau    %
% concernant le choix de la fonction de transfert ainsi que le choix de la methode d'apprentissage %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if FLAG_APP==1
    
    if net.FLAG_DELTA==0
        disp('Changement signe OK');
    elseif net.FLAG_DELTA==1
        disp('Pas de changement de signe');
    else
        disp('J_ref borne inf');
    end;
    
else
    disp('Pas d''apprentissage');
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocation de BUF_S (dynamique d'activation), BUF_I (dynamique de l'input) et BUF_H (dynamique des champs locaux) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables temporaires : BUF_S, BUF_I, BUF_H
% A chaque population p est associ� BUF_S{p} qui contient la dynamique pour cette population, de taille N(p) x nbp
% Cette matrice est completee � chaque pas de temps pendant l'iteration de la dynamique
% de meme pour BUF_I et BUF_H (dynamique associee � chaque classe de liens)


for p=1:net.nb_pop
    BUF_S{p}=zeros(net.N(p),nbp); % Allocation de la matrice memorisant les etats d'activation
    BUF_V{p}=zeros(net.N(p),nbp); % Allocation de la matrice memorisant les etats d'activation
    BUF_U{p}=zeros(net.N(p),nbp);
    BUF_I{p}=zeros(net.N(p),nbp); % Allocation de la matrice memorisant les entrees
    for q=1:net.nb_pop
        BUF_H{p}{q}=zeros(net.N(p),nbp); % Allocation de la matrice memorisant les champs locaux
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                    %%
%%     BOUCLE PRINCIPALE (1 boucle=1 pas de temps)    %%
%%                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variable temporaire : t
% Compteur de temps (entre 1 et nbp)
% l'activit� du r�seau est mise � jour � chaque pas de temps
%
% Variable temporaire : ARRET
% Permet d'arreter la boucle avant terme (break)

ARRET=0;
for t=1:nbp
    
    
    %t_buf=(t-300);
    %t_buf=t_buf.*(t_buf>0)+(t_buf<=0);
    
    %t_app=(t-round(net.nbp_z/3));
    %%t_app=t_app.*(t_app>0)+(t_app<=0);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                - 1 -                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mise � jour de la matrice d'activation source (in) pour l'apprentissage %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Variable temporaire : S_in (cellule)
    % Pour p,q,tau, S_in{p}{q} (matrice taille N(p)xN(q)) est mise � jour
    % cette mise � jour est r�alis�e AVANT la mise � jour de la dynamique
    % A partir des valeurs d'activation calcul�es au tour pr�c�dent
    
    
    if FLAG_APP==1
 
        for p=1:net.nb_pop
            net.epsilon_i_flat{p}=reshape(net.epsilon_i{p}',1,net.N(p)*(net.tau_max(p)+1));
        end;
       
        for p=1:net.nb_pop
            for q=1:net.nb_pop
                if alpha_eff(p,q)~=0
                    
                    % Calcul de m_in{p}{q}
                    m_in{p}{q} = sparse(net.masque_i{p}{q}, net.masque_j{p}{q},net.epsilon_i_flat{q}(net.tau_eff{p}{q}),net.N(p),net.N(q));%net.epsilon{p}{q};%
                    %m_in{p}{q} = sparse(net.J_masque{p}{q}.*net.epsilon_i_flat{q}(net.tau_eff{p}{q}));%net.epsilon{p}{q};%
                    %*spdiags(net.epsilon_i{q}(:,net.tau{p}{q}),0,net.N(q),net.N(q));
                    
                    % Calcul de S_in{p}{q}
                    if strcmp(net.app,'stdp')
                        S_in{p}{q} = sparse(net.masque_i{p}{q}, net.masque_j{p}{q},net.S_flat{q}(net.tau_eff{p}{q}),net.N(p),net.N(q));
                        %S_in{p}{q} = sparse(net.J_masque{p}{q}.*net.S_flat{q}(net.tau_eff{p}{q}));
                        %*spdiags(net.S{q}(:,net.tau{p}{q}),0,net.N(q),net.N(q));
                    end;
                end;
            end;
        end;
    end;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       - 2 -                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mise � jour de la dynamique pour chaque population %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for p=1:net.nb_pop
        
        %%%%%     Mise � jour de l'input      %%%%%
        %
        % Parametre du reseau : net.bruit_I
        % Taux de bruit de l'entree sur chaque couche
        %
        % Variable temporaire : I_in (cellule)
        % I_in{p} contient le motif present� � l'instant t
        % (possiblement issu de la dynamique de l'environnement simul� � t-1)
        % Ce vecteur peut etre corrompu par un terme de bruit
        
        if net.FLAG_REV
            r_plus_ref=(net.V_plus-net.m_V{p})/net.V_plus; %ones(net.N(p),1);% (1-net.bool_refr{p}(:,1)).*
            r_plus=spdiags(r_plus_ref.*(r_plus_ref>0),0,net.N(p),net.N(p));
            r_moins_ref=(net.m_V{p}-net.V_moins)/abs(net.V_moins); %ones(net.N(p),1);% (1-net.bool_refr{p}(:,1)).* %ones(net.N(p),1);%
            r_moins=spdiags(r_moins_ref.*(r_moins_ref>0),0,net.N(p),net.N(p));
        end;
       
        %%%%%%%%%
        %   I   %
        %%%%%%%%%
        
        I_courant=net.I{p}(:,rem(net.t_abs,tau_I)+1);
        if net.bruit_I(p)>0
            net.rnd_I{p}=net.delta*net.rnd_I{p}+randn(net.N(p),1)*net.bruit_I(p)*sqrt(2*net.delta_t/net.tau_b-(net.delta_t/net.tau_b)^2);
            buf_I_in=I_courant;
        else
            net.rnd_I{p}=0;
            buf_I_in=I_courant;
        end;
        
        if net.FLAG_REV
            net.buf_I{p}=net.gamma{p}.*net.norm_g{p}.*net.buf_I{p}+r_plus*((1-net.gamma{p}.*ones(net.N(p),1)).*buf_I_in);%buf_I_in*net.delta_t;%%-net.buf_I{p});
        else        
            net.buf_I{p}=net.gamma{p}.*net.buf_I{p}+(1-net.gamma{p}).*buf_I_in;%buf_I_in*net.delta_t;%%-net.buf_I{p});
            %net.buf_I{p}=net.gamma{p}.*net.norm_g{p}.*net.buf_I{p}+0.05.*buf_I_in;
        end
        
        I_in{p}=net.buf_I{p}+net.rnd_I{p};
        BUF_I{p}(:,t)=I_in{p};
        
        
        %%%%%      Mise � jour de l'activit�      %%%%%
        %
        % Variable globale : net.h (cellule)
        % net.h{p}{q} contient la valeur du champ local au temps t
        %
        % Variable globale : net.V (cellule)
        % net.V{p} contient la valeur du potentiel de membrane au temps t
        %
        % Variable globale : net.S (cellule)
        % net.S{p} contient la valeur de l'activation aux temps t, t-1,..., t-tau_max.
        % rappel : mis � jour avec com_trans, qui depend de la definition de la fonction de transfert
        % (voir plus haut)
        
        
        if net.FLAG_BAYES
            % Pour f_bayes
            tire=sum(net.S{p}(:,1:net.nbp_r),2);
            
            indices_tire=find(tire==1);
            net.bayes{p}.p_tire=delta*net.bayes{p}.p_tire+(1-delta)*tire;
            
            %u_tire=mean(net.U{p}(indices_tire,1:net.nbp_r),2);
            u_tire=net.U{p}(indices_tire,1);
            net.bayes{p}.m_u_tire(indices_tire)=delta*net.bayes{p}.m_u_tire(indices_tire)+(1-delta)*u_tire;
            net.bayes{p}.var_u_tire(indices_tire)=delta*net.bayes{p}.var_u_tire(indices_tire)+(1-delta)*((u_tire-net.bayes{p}.m_u_tire(indices_tire)).^2);
            
            du_tire=net.diff_U{p}(indices_tire);
            net.bayes{p}.m_du_tire(indices_tire)=delta*net.bayes{p}.m_du_tire(indices_tire)+(1-delta)*du_tire;
            net.bayes{p}.var_du_tire(indices_tire)=delta*net.bayes{p}.var_du_tire(indices_tire)+(1-delta)*((du_tire-net.bayes{p}.m_du_tire(indices_tire)).^2);
            
            u_cov=((u_tire-net.bayes{p}.m_u_tire(indices_tire)).*(du_tire-net.bayes{p}.m_du_tire(indices_tire)));
            net.bayes{p}.cov_tire(indices_tire)=delta*net.bayes{p}.cov_tire(indices_tire)+(1-delta)*u_cov;
            
            indices_silent=find(tire==0);
            
            %u_silent=mean(net.U{p}(indices_silent,1:net.nbp_r),2);
            u_silent=net.U{p}(indices_silent,1);
            net.bayes{p}.m_u_silent(indices_silent)=delta*net.bayes{p}.m_u_silent(indices_silent)+(1-delta)*u_silent;
            net.bayes{p}.var_u_silent(indices_silent)=delta*net.bayes{p}.var_u_silent(indices_silent)+(1-delta)*((u_silent-net.bayes{p}.m_u_silent(indices_silent)).^2);
            
            du_silent=net.diff_U{p}(indices_silent);
            net.bayes{p}.m_du_silent(indices_silent)=delta*net.bayes{p}.m_du_silent(indices_silent)+(1-delta)*du_silent;
            net.bayes{p}.var_du_silent(indices_silent)=delta*net.bayes{p}.var_du_silent(indices_silent)+(1-delta)*((du_silent-net.bayes{p}.m_du_silent(indices_silent)).^2);
            
            u_cov=((u_silent-net.bayes{p}.m_u_silent(indices_silent)).*(du_silent-net.bayes{p}.m_du_silent(indices_silent)));
            net.bayes{p}.cov_silent(indices_silent)=delta*net.bayes{p}.cov_silent(indices_silent)+(1-delta)*u_cov;
        end;
        
        
        %%%%%%%%%
        %   H   %
        %%%%%%%%%
        sum_h=0;
        for q=1:net.nb_pop
            if net.connect(p,q)
                %disp ([num2str(p),' ',num2str(q)])
                %if p==6
                %    1;
                %end;
                %tire = net.J_masque{p}{q}.*net.S_flat{q}(net.tau_eff{p}{q});
                tire = sparse(net.masque_i{p}{q}, net.masque_j{p}{q}, net.S_flat{q}(net.tau_eff{p}{q}),net.N(p),net.N(q));
                if net.FLAG_REV
                    epsilon_plus=r_plus*(net.J_masque_plus{p}{q}.*tire);%
                    epsilon_moins=r_moins*(net.J_masque_moins{p}{q}.*tire);
                    epsilon_new=epsilon_plus+epsilon_moins;
                    net.epsilon{p}{q}=spdiags(net.gamma{p}.*net.norm_g{p},0,net.N(p),net.N(p))*net.epsilon{p}{q}+epsilon_new;
                else
                    net.epsilon{p}{q}=spdiags(net.gamma{p}.*ones(net.N(p),1),0,net.N(p),net.N(p))*net.epsilon{p}{q}+tire;
                end;
                if net.FLAG_STD == 1
                    h=sum(net.J_eff{p}{q}.*net.epsilon{p}{q},2);
                else
                    h=sum(net.J{p}{q}.*net.epsilon{p}{q},2);
                end;
                sum_h=sum_h+h;
                BUF_H{p}{q}(:,t)=h;%sum(net.J{p}{q}.*tire,2); %
                net.h{p}{q}=[h,net.h{p}{q}(:,1:(net.tau_max(p)))];
                %
                if strcmp(net.app,'concurrent1')
                    net.m_h{p}{q}=delta*net.m_h{p}{q}+(1-delta)*net.h{p}{q}(:,1);
                else
                    net.m_h{p}{q}=net.gamma{p}.*net.m_h{p}{q}+(1-net.gamma{p}).*net.h{p}{q}(:,1);
                end;
                diff=(net.h{p}{q}(:,1)-net.h{p}{q}(:,2)).*net.nbp_mb{p};
                net.diff_h{p}{q}=net.gamma{p}.*net.diff_h{p}{q}+(1-net.gamma{p}).*diff;
            end;
        end;
        
        %%%%%%%%%
        %   U   %
        %%%%%%%%%
        
        net.U{p}=[sum_h+flag_in{p}.*I_in{p},net.U{p}(:,1:net.tau_max(p))];
        if strcmp(net.app,'concurrent1')
            net.m_U{p}=delta*net.m_U{p}+(1-delta)*net.U{p}(:,1);
        else
            net.m_U{p}=net.gamma{p}.*net.m_U{p}+(1-net.gamma{p}).*net.U{p}(:,1);
        end;
        %diff=(net.U{p}(:,1)-net.U{p}(:,2)).*net.nbp_mb{p};
        %net.diff_U{p}=net.gamma{p}.*net.diff_U{p}+(1-net.gamma{p}).*diff;
        BUF_U{p}(:,t)=net.U{p}(:,1);
        
        %%%%%%%%%
        %   V   %
        %%%%%%%%%
        
        V = net.U{p}(:,1)-net.eta{p}(:,1);
        if net.FLAG_REV
            net.norm_g{p}=ones(net.N(p),1);
            indices_plus=find(V>=1);%net.V_plus);%
            net.norm_g{p}(indices_plus)=1./V(indices_plus);%net.V_plus./V(indices_plus);%
            indices_moins=find(V<=net.V_moins);%-1);%
            net.norm_g{p}(indices_moins)=net.V_moins./V(indices_moins);%-1./V(indices_moins);%
            V = V.* net.norm_g{p};
        end;
        net.V{p}=[V,net.V{p}(:,1:net.tau_max(p))];
        net.m_V{p}=(1-net.delta_t/net.tau_r)*net.m_V{p}+net.delta_t/net.tau_r*V;
        BUF_V{p}(:,t)=net.V{p}(:,1);
    
    end;
    
    %%%%%%%%%
    %   S   %
    %%%%%%%%%
    for p=1:net.nb_pop
        if 1 %TEST==0
            net.S{p}=[f_heavy(net,net.V{p}(:,1),net.theta{p}).*(1-net.bool_refr{p}(:,1)),net.S{p}(:,1:(net.tau_max(p)))];
            net.S_flat{p}=sparse(reshape(net.S{p}',1,net.N(p)*(net.tau_max(p)+1)));
            net.bool_refr{p}=[sum(net.S{p}(:,1:net.nbp_r-1),2)>0,net.bool_refr{p}(:,1:(net.tau_max(p)))];
        else
            net.m{p}=net.m{p}+f_tanh_01(net.U{p}(:,1)-3,0.5)*net.delta_t/net.tau_r;
            net.S{p}=[net.m{p}>net.theta{p},net.S{p}(:,1:(net.tau_max(p)))];
            net.m{p}=net.m{p}-net.S{p}(:,1);
        end;
        if isfield(net,'FLAG_SFA')
            if net.FLAG_SFA   % spike frequency adaptation / time constant = 300 ms
                net.theta{p}=(1-net.delta_t/300)*net.theta{p}+100*(net.delta_t/300)*net.S{p}(:,1);
            end;
        end;
        BUF_S{p}(:,t)=net.S{p}(:,1);
    end;
    
    %%%%%%%%%
    %  ETA  %
    %%%%%%%%%
    for p=1:net.nb_pop
        net.epsilon_i{p}=[net.gamma{p}.*net.epsilon_i{p}(:,1)+net.S{p}(:,1),net.epsilon_i{p}(:,1:(net.tau_max(p)))];
        %if strcmp(net.app,'stdp')
        %    net.epsilon_i_flat{p}=reshape(net.epsilon_i{p}',1,net.N(p)*(net.tau_max(p)+1));
        %end;
        V=net.U{p}(:,1)-net.eta{p}(:,1);
        if net.FLAG_REV
            net.eta{p}=[net.gamma{p}.*(net.eta{p}(:,1)+net.S{p}(:,1).*V).*net.norm_g{p},net.eta{p}(:,1:(net.tau_max(p)))];
        else
            net.eta{p}=[net.gamma{p}.*(net.eta{p}(:,1)+net.S{p}(:,1).*V),net.eta{p}(:,1:(net.tau_max(p)))];
        end
    end;
    
    if (strcmp(net.app,'concurrent1')||strcmp(net.app,'concurrent2'))&&FLAG_APP==1  % ||strcmp(net.app,'stdp')
        for p=1:net.nb_pop
            for q=1:net.nb_pop
                if net.connect(p,q)
                    U=net.U{p}(:,1)-flag_in{p}.*net.rnd_I{p};
                    if 0 %TEST
                        F=f_tanh_01(U-3,0.5);
                    else
                        F=f_stat(net,U);
                    end;
                    net.F{p}{q}=[F,net.F{p}{q}(:,1:(net.tau_max(p)))];
                end;
            end;
        end;
    end;
    
    if strcmp(net.app,'test_fl2')
        for p=1:net.nb_pop
            for q=1:net.nb_pop
            end;
        end;
    end;
    
    
    %%% Short Term Depression (en test) %%%
    
    if net.FLAG_STD == 1
        for p=1:net.nb_pop
            for q=1:net.nb_pop
                %for tau=net.indices_tau{p}{q}
                % Calcul de m_in pour short term depression
                %m_in = net.J_masque{p}{q}{tau}*spdiags(net.epsilon{q}(:,tau),0,net.N(q),net.N(q));
                %net.J_eff{p}{q}{tau}=net.J{p}{q}{tau}.*(1-m_in/5);
                net.J_eff{p}{q}=net.J{p}{q}.*(1-net.epsilon{p}{q}*0.25).*((1-net.epsilon{p}{q}*0.25) > 0);
                %end;
            end;
        end;
    end;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         - FIN 2 -           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                - 3 -                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mise � jour de la matrice d'activation cible (out) pour l'apprentissage %
    % Et mise a jour de la matrice des poids.                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % (boucle sur p,q,tau...)
    %
    % Variable temporaire : a_out (matrice taille N(p)xN(q))
    % (voir plus haut : com_a_out lignes 120-130)
    % cette mise � jour est r�alis�e APRES la mise � jour de la dynamique
    % A partir des valeurs d'activation calcul�es plus haut
    %
    % Variable temporaire : delta_J
    % contient le produit S_out x S_in pond�r� par alpha_eff
    % i.e. la variation instantanee des poids
    %
    % Parametre du r�seau : net.renf (booleen)
    % active l'apprentissage par renforcement
    %
    % Parametre du reseau : net.lambda
    % en cas d'apprentissage par renforcement, donne le taux d'oubli de la trace (mem_delta_J)
    %
    % Variable globale : net.mem_delta_J (cellule)
    % net.mem_delta_J{p}{q}{tau} contient la trace des correlations les plus recentes dans le cadre
    % de l'apprentissage par renforcement.
    % remarque : net.mem_delta_J est initialisee dans init_dyn
    %
    % variable globale : net.REWARD
    % contient la valeur courante du signal de renforcement
    % (a priori limit� � -1, 0 ou 1)
    %
    % Variable temporaire : delta_J_rew
    % en cas d'apprentissage par renforcement,
    % contient la valeur de mem_delta_J pond�r�e par le renforcement net.REWARD
    
    
    
    
    if FLAG_APP==1
        if 0 %net.REWARD~=0
            disp(['REWARD=',num2str(net.REWARD)]);
        end;
        for p=1:net.nb_pop
            for q=1:net.nb_pop
                if alpha_eff(p,q)~=0
                    
                    %%%%%%%%%
                    % S_out %
                    %%%%%%%%%
                    
                    if strcmp(net.app,'hebb')||strcmp(net.app,'stdp')
                        %S_out = spdiags((1-net.F{p}{q}(:,1)).*net.S{p}(:,1),0,net.N(p),net.N(p));%*net.J_masque{p}{q};
                        S_out = spdiags(net.S{p}(:,1),0,net.N(p),net.N(p));%*net.J_masque{p}{q};
                    elseif strcmp(net.app,'hebb2')
                        S_out = spdiags(net.S{p}(:,1)-0.5*net.delta_t/net.tau_r,0,net.N(p),net.N(p));%*net.J_masque{p}{q};
                    elseif strcmp(net.app,'concurrent1')||strcmp(net.app,'concurrent2')||strcmp(net.app,'florian1')
                        S_out = spdiags((1-net.F{p}{q}(:,1))./(1-net.F{p}{q}(:,1)*net.delta_t/net.tau_r).*(net.S{p}(:,1)-net.F{p}{q}(:,1)*net.delta_t/net.tau_r),0,net.N(p),net.N(p));%*net.J_masque{p}{q};
                    end;
                    
                    %%%%%%%%%
                    % m_out %
                    %%%%%%%%%
                    
                    if strcmp(net.app,'stdp')||strcmp(net.app,'concurrent2')||strcmp(net.app,'test_fl2')
                        m_out = spdiags(net.epsilon_i{p}(:,1),0,net.N(p),net.N(p));%*net.J_masque{p}{q};
                    end;
                    
                    %%%%%%%%%%%
                    % Delta_J %
                    %%%%%%%%%%%
                    
                    if strcmp(net.app,'hebb')||strcmp(net.app,'hebb2')||strcmp(net.app,'concurrent1')||strcmp(net.app,'florian1')||strcmp(net.app,'test_fl')
                        %norm_alpha_out = sparse(spdiags((epsilon_max-net.epsilon{p}(:,1)).*((epsilon_max-net.epsilon{p}(:,1))>0),0,net.N(p),net.N(p)))*net.J_masque{p}{q};
                        delta_J=alpha_eff(p,q).*(S_out*m_in{p}{q});
                        %delta_J=alpha_eff(p,q).*norm_alpha_out.*(S_out.*m_in{p}{q});
                    elseif strcmp(net.app,'stdp')||strcmp(net.app,'concurrent2')||strcmp(net.app,'test_fl2')
                        %tt = spdiags(,0,net.N(p),net.N(p));
                        %delta_J=alpha_eff(p,q).*(tt*(S_out*m_in{p}{q} - m_out*S_in{p}{q}));
                        delta_J=alpha_eff(p,q).*(0.9 * S_out*m_in{p}{q} - m_out*S_in{p}{q});
                        %delta_J = delta_J .* (net.tau{p}{q}/20); % !!!??? rééquilibrage pour indépendance en tau !!!???
                    else
                        disp('R�gle non r�f�renc�e!!');
                        delta_J=sparse(net.N(p),net.N(q));
                    end;
                    
                    %%%%%%%%%
                    % Trace %
                    %%%%%%%%%
                    
                    if net.FLAG_RAZ
                        net.trace{p}{q}=net.trace{p}{q}+(1-net.lambda)*delta_J;
                    else
                        net.trace{p}{q}=net.lambda*net.trace{p}{q}+(1-net.lambda)*delta_J;
                    end;
                    
                    if net.REWARD~=0
                        J_scale=sparse(net.N(p),net.N(q));
                        norm_REWARD = net.REWARD/(1-net.lambda);
                        delta_J_rew=norm_REWARD*net.trace{p}{q};
                        
                        if net.FLAG_DELTA==0
                            J_test=1;
                        else
                            J_new=net.J{p}{q}+delta_J_rew;
                            if net.FLAG_DELTA==1
                                J_test=(J_new.*net.J_ref{p}{q})>0;
                            else
                                J_test=((J_new.*net.J_ref{p}{q})>0)&(net.J_ref{p}{q}>0);
                                %J_test=(abs(J_new)>abs(net.J_ref{p}{q}))&((J_new.*net.J_ref{p}{q})>0);
                            end;
                        end;
                        
                        delta_J_eff=delta_J_rew.*J_test;%.*net.masque_app;
                        
                        if net.FLAG_SCALING || net.FLAG_THETA
                            J_scale=J_scale+delta_J_eff;
                        end;
                        
                        if net.FLAG_SCALING
                            %delta_J_scale=spdiags(-sum(J_scale')',0,net.N(p),net.N(p))*net.J_masque{p}{q};
                            %net.J{p}{q}=net.J{p}{q}+delta_J_scale/(net.N(q)*net.densite(p,q));
                            %ampli=1;
                            
                            ampli=sum(abs(net.J{p}{q}'))'./sum(abs((net.J{p}{q}+delta_J_eff))')';
                            ampli=spdiags(ampli,0,net.N(p),net.N(p));
                        else
                            ampli=1;
                        end;
                        net.J{p}{q}=ampli*(net.J{p}{q}+delta_J_eff);
                        
                        
                        if net.FLAG_RAZ
                            net.trace{p}{q}=0;
                        end;
                        
                        
                        if net.FLAG_THETA
                            
                            Delta=sparse(net.N(p),net.N(q));
                            for tau=net.indices_tau{p}{q}
                                Delta=Delta+net.J{p}{q}-net.J_ref{p}{q};
                            end;
                            
                            %net.var_theta{p}=net.var_theta{p}+var(J_scale,0,2)./var(net.J_flat_ref{p}{q}+J_scale,0,2);
                            net.var_theta{p}=var(Delta,0,2)./var(net.J_flat_ref{p}{q}+Delta,0,2);
                            net.theta{p}=sqrt(net.theta_ref{p}+net.var_theta{p});
                            %net.theta{p}=net.theta_ref{p}+sqrt(net.var_theta{p});
                        end;
                        
                        
                    end; % if  net.REWARD~=0
                    
                end; % if alpha...
            end; % for q...
        end; % for p...
        
    end; % if FLAG_APP...
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         - FIN 3 -           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    net.t_abs=net.t_abs+1;     % nb de pdt depuis l'initialisation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Affichage de l'avancement (tous les 200 ms) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if net.FLAG_DASHBOARD
        for p=1:net.nb_pop
            for q=1:net.nb_pop
                if FLAG_APP==1 & alpha_eff(p,q)~=0
                    %for i=3:3
                    net.dash.sum_J{p}{q}= [net.dash.sum_J{p}{q},net.J{p}{q}(3,:)'];
                    %end;
                end;
            end;
        end;
    end;
    if rem(net.t_abs,(20/net.delta_t))==0
        disp(['t=',num2str(net.t_abs*net.delta_t)]);
        if net.FLAG_DASHBOARD
            for p=1:net.nb_pop
                % calcul de la dimension effective sur le signal d'input
                %t_min=max(1,t+1-200/net.delta_t);
                %C=cov(BUF_U{p}(:,t_min:t)');
                %E=eig(C);
                %pp=E/sum(E);
                %dim=exp(-sum(pp.*log(pp)));
                %net.dash.dim{p}=[net.dash.dim{p},[net.t_abs*net.delta_t;dim]];
                % Calcul de la periode
                %periode = [];
                %for i=1:net.N(p)
                %   periode(i)=calc_periode(BUF_U{p}(i,t_min:t),0)*net.delta_t;
                %end;
                %net.dash.periode{p}=[net.dash.periode{p},periode'];
                % Moyenne et ecart-type de delta_J
                if FLAG_APP==1
                    if alpha_eff(p,q)~=0
                        for q=1:net.nb_pop
                            net.dash.mean_Delta_J{p}{q}=[net.dash.mean_Delta_J{p}{q},[net.t_abs*net.delta_t;mean(nonzeros(net.J{p}{q}-net.J_ref{p}{q}))]];
                            net.dash.std_Delta_J{p}{q}=[net.dash.std_Delta_J{p}{q},[net.t_abs*net.delta_t;std(nonzeros(net.J{p}{q}-net.J_ref{p}{q}))]];
                        end;
                    end;
                end;
            end;
        end;
    end;
    
    if (ARRET==1)
        break;
    end;
    
    %%%%% Execution du script %%%%%
    %
    % Parametre du reseau : net.script_out
    % contient un script qui decrit l'evolution de l'environnement pendant un pas de temps
    % ce script permet de mettre � jour net.REWARD en cas d'apprentissage par renforcement
    %
    % Attention!!ce script etant realis� en parallele des etapes 1-2-3,
    % l'entree du processus du script doit reposer sur x(t-1)
    % c'est a dire concretement sur les net.x{p}(:,2)
    
    if strcmp(net.script_out,'')==0
        eval(net.script_out);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                    %%
%%            FIN DE LA BOUCLE PRINCIPALE             %%
%%                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Memorisation de la dynamique en sortie de boucle %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables globales : net.DYN_S, net.DYN_I, net.DYN_H
% recoivent le contenu de BUF_S, BUF_I, BUF_H

for p=1:net.nb_pop
    net.DYN_S{p}=[net.DYN_S{p},BUF_S{p}];
    net.DYN_V{p}=[net.DYN_V{p},BUF_V{p}];
    net.DYN_U{p}=[net.DYN_U{p},BUF_U{p}];
    net.DYN_I{p}=[net.DYN_I{p},BUF_I{p}];
    for q=1:net.nb_pop
        if net.connect(p,q)
            net.DYN_H{p}{q}=[net.DYN_H{p}{q},BUF_H{p}{q}];
        end;
    end;
end;
