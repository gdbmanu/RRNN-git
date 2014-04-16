% Initialisation des parametres modele P populations

function net=init_param_nips_s(delta_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 PARAMETRES GENERAUX                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


net.fichier_I={''};       % net.fichier_I d�crit le nom d'un fichier contenant 
                                % la description du signal � appliquer au r�seau
                        

net.script_out='iter_env_cns_lif';              % Pour ex�cuter un script au cours de la dynamique   

net.script_init='init_env_cns_lif';
                        
net.nb_pop = 2;                 % Nombre de populations 

net.N=[1000;                     % net.N : taille du r�seau
       800];
   
%net.N = [500;
%         100];
     
net.delta_t = delta_t;                      % R�solution  (en ms)

%net.tau_m = 10;                                 % Constante de temps du neurone (ms) 
%net.nbp_m = floor(net.tau_m/net.delta_t);  % Temps caract�ristique du neurone (en pas de temps)
%net.gamma=(1-net.delta_t/net.tau_m);            % Fuite du neurone
net.tau_m_ref = [10;
                 10];
%net.param_tau = [1   1;
%                 1   1];
             
net.tau_d = 10;%3; % (param/ fred) 6;%                                % D�lai moyen (ms)
net.nbp_d = floor(net.tau_d/net.delta_t);  % D�lai effectif (en pas de temps)

net.tau_r = 1 ;                         % Periode refractaire (ms)
net.nbp_r=floor(net.tau_r/net.delta_t);  % Periode refractaire (en pas de temps)
net.gamma_r=(1-net.delta_t/net.tau_r);

net.tau_z =  delta_t;% 10 ; %                       % Constante de temps de la trace (ms)
net.nbp_z=floor(net.tau_z/net.delta_t);
net.lambda=(1-net.delta_t/net.tau_z);   % Fuite de la trace

net.dim=zeros(net.nb_pop);              % net.dim : topologie des liens :
                                        % 0 : pas de topologie
                                        % 1 : topologie 1D
                                        % 2 : topologie 2D
                           
net.rayon_plus=zeros(net.nb_pop);
                
net.rayon=zeros(net.nb_pop);

net.grad = zeros(net.nb_pop);

net.densite=ones(net.nb_pop);          % densit� de connexion (entre 0 et 1)
                            
net.distr={'unif','unif'}; % net.distr : distribution de poids :
                                             % 'gauss' : distribution gaussienne
                                             % 'unif' : distribution uniforme
                                             % 'exp' : distribution exponentielle
                                             % 
                            
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LE SIGNAL                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.periode_I=zeros(net.nb_pop,1);   %periode de l'input
%net.periode_I = 1000 * ones(net.nb_pop,1);   %periode de l'input
              
net.rand_I=zeros(net.nb_pop,1);            

%net.sequence_I{1}=zeros(net.nb_pop,1);  %la sequence d'entree
%P11 = randn(net.N(1),1); P12 = randn(net.N(1),1);
%P21 = randn(net.N(2),1); P22 = randn(net.N(2), 1);
%net.sequence_I{1} = zeros(net.N(1), 1000);
%net.sequence_I{2} = zeros(net.N(2), 1000);
%for t = 1:1000 
%    net.sequence_I{1}(:,t) = P11 * cos(2 * pi * t / 1000) + P12 * sin(2 * pi * t / 1000); 
%    net.sequence_I{2}(:,t) = P21 * cos(2 * pi * t / 1000) + P22 * sin(2 * pi * t / 1000); 
%end;

net.flag_I = [1;                 % net.flag_I : bool�en qui active l'entree
              1];

net.bruit_I=[0;                  % net.bruit_I : taux de bruit sur l'entree
             0];


%net.num_I=0;     % seed alea

net.fichier_I = {'',''};
%net = init_sequence_s(net);  %genere la sequence du vecteur d'entree

I_ref = 1;%0;%1;%;
net.I{1} = I_ref;%0;%
net.I{2} = I_ref;%0;%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LES POIDS                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 0; % I/E ( Brunel / 4)

%d = sqrt(1 + g^2) / (g - I_ref) / 2; %4 %2.5 %   3 %    
% d = sqrt(1 + g^2) / alpha : std globale = alpha
% 1 - (1 - g) - I_ref = g - I_ref = distance au seuil
% std esperee : alpha * dist au seuil 

%d = sqrt(1 + g^2); % std globale = 1

J_ref= 10; %10; %3; % 20; %15; %10; %

sigma_J_ref = 1.5;  %2; % sigma_J_total 

d = J_ref * sqrt(1 + (1 + g/J_ref)^2) / sigma_J_ref;

%net.connex = J_ref* [1               -g;
%                     1               -g] ;   

tt = net.tau_m_ref(1) / net.tau_m_ref(2);
net.connex =[J_ref        -(J_ref + g) ;
             J_ref * tt   -(J_ref + g) * tt];     % diff excitation / inhibition
         
net.connect = net.connex ~= 0 ;          

net.J_barre = net.connex;

net.J_barre_eff = net.J_barre; % .* net.param_tau;

net.sigma_J_eff = abs(net.connex/d); % .* net.param_.tau;

net.N_aff=(net.N*ones(1,net.nb_pop))';                                              
                       
%%%% DENSITE CONTRAINTE : AJUSTEMENT DE LA DENSITE SELON J_BARRE ET SIGMA_J_EFF %%%%                                
     
d_eff = net.J_barre_eff./sqrt(net.sigma_J_eff.^2+(net.J_barre_eff.^2)./net.N_aff);

net.K = d_eff.^2;          % poids constants (loi uniforme)

%net.K = 4/3 * d_eff.^2;   % poids [0,m] (loi uniforme)

%net.K = 2*d_eff.^2;         % loi exponentielle

%%%% DENSITE NON CONTRAINTE %%%%

%net.K = [100 100;
%         100 100]; % 

%net.K = net.N_aff; %K; % loi chi2 ou gauss (pas d'ajustement)

%%%% PARAMETRAGE DE LA DENSITE ET AJUSTEMENT DE SIGMA_J %%%%

net.densite=net.K./net.N_aff;          % densit� de connexion (entre 0 et 1)

net.sigma_J = sqrt(abs(net.sigma_J_eff.^2-(1-net.densite).*(net.J_barre_eff.^2)./net.K));


%%%% AJUSTEMENTS LIES A LA DIMENSION %%%%
        
net.dim=[0 1;
         1 1];
     
net.rayon=[0 0.3;
           0.3 0.3];
        
net.grad=zeros(net.nb_pop);
net.FLAG_GRAD='';


        
%%% DELAIS

sigma_tau =  19; %8;%  %8; %8;% net.nbp_d - 1;%
net.tau_min= (20 - sigma_tau) * ones(net.nb_pop); %(net.nbp_d - sigma_tau) * ones(net.nb_pop);%net.nbp_d * ones(net.nb_pop);%net.nbp_d-4;%(param. fred)              % net.tau_min : delai de transmission minimal                

net.tau_moy= sigma_tau * ones(2);% [net.nbp_d-1    net.nbp_d-1;       %[0 0; % net.tau_moy : parametre pour la distribution de poisson
             %net.nbp_d-1    net.nbp_d-1];      % 0 0]; %
net.FLAG_UNIF = 1;                  % delais uniformes / Poisson
                          
                                           % (la valeur 0 indique des delais constants)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LES SEUILS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.theta_barre=ones(net.nb_pop,1);        % net.theta_barre : valeur moyenne du seuil

net.sigma_theta=zeros(net.nb_pop,1);          % net.sigma_theta : ecart-type du seuil

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       APPRENTISSAGE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.alpha=1;              % net.alpha : parametre d'apprentissage
                            % (regle l'intensit� de la modification des poids)


net.norm_alpha=[1 0;  %[1   -1 ;     % net.norm : normalisation du terme d'apprentissage           
                0 0]; % 0.1 0.1];     % (par population)
net.app='stdp';%'hebb';%
net.FLAG_DELTA=1;
net.FLAG_RAZ=0;
net.FLAG_DASHBOARD=1;
net.FLAG_STD = 0;
net.FLAG_REV = 1; 
    net.V_plus = 5;
    net.V_moins = -1;
net.FLAG_SCALING = 1;
net.FLAG_STDP_EPS = 0;
net.FLAG_BAYES = 0;
net.FLAG_THETA = 0;
net.FLAG_SFA = 1;


net.REWARD=1;                       

                   

