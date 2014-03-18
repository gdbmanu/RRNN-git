% Initialisation des parametres modele P populations

function net=init_param_rrnn_lif(delta_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 PARAMETRES GENERAUX                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


net.fichier_I={''};       % net.fichier_I d�crit le nom d'un fichier contenant 
                                % la description du signal � appliquer au r�seau
                        

net.script_out='iter_env_s';              % Pour ex�cuter un script au cours de la dynamique   

net.script_init='';
                        
net.nb_pop = 1;                 % Nombre de populations 

net.N=400;                     % net.N : taille du r�seau
        
net.delta_t = delta_t;                      % R�solution  (en ms)

net.mb_inf = 1;                                 % Constante de temps du neurone (ms) 
net.mb_sup = 40;                                 % Constante de temps du neurone (ms) 


net.tau_d = 6;% (param/ fred) 6;%                                % D�lai moyen (ms)
net.nbp_d = floor(net.tau_d/net.delta_t);  % D�lai effectif (en pas de temps)

net.tau_r = 1 ;                         % Periode refractaire (ms)
net.nbp_r=floor(net.tau_r/net.delta_t);  % Periode refractaire (en pas de temps)
net.gamma_r=(1-net.delta_t/net.tau_r);

net.tau_z =  delta_t;% 100; %                       % Constante de temps de la trace (ms)
net.nbp_z=floor(net.tau_z/net.delta_t);
net.lambda=(1-net.delta_t/net.tau_z);   % Fuite de la trace

                            
net.distr={'unif'}; % net.distr : distribution de poids :
                                             % 'gauss' : distribution gaussienne
                                             % 'unif' : distribution uniforme
                            
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LE SIGNAL                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.periode_I=1;   %periode de l'input

net.rand_I=0;            

net.sequence_I{1}=0;  %la sequence d'entree

net.flag_I = 1;                 % net.flag_I : bool�en qui active l'entree

net.bruit_I=0;                  % net.bruit_I : taux de bruit sur l'entree

%net.num_I=0;     % seed alea

%net = init_sequence_s(net);  %genere la sequence du vecteur d'entree
net.I{1}=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LES POIDS                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.connect = 1;

J_ref=1.3;% param g : 0.9-->0.99 %1.5;

net.J_barre = J_ref;

d=5;

net.sigma_J_eff = J_ref/d;  %2 * net.tau_m/10; % (param. fred) 1;

net.N_aff=(net.N*ones(1,net.nb_pop))';                                              
                          
%%%% DENSITE CONTRAINTE : AJUSTEMENT DE LA DENSITE SELON J_BARRE ET SIGMA_J_EFF %%%%                                
     
d_eff = net.J_barre./sqrt(net.sigma_J_eff.^2+(net.J_barre.^2)./net.N_aff);

net.K = d_eff.^2;          % poids constants (loi uniforme)

%net.K = 4/3 * d_eff.^2;   % poids [0,m] (loi uniforme)

%net.K = 2*d_eff.^2;         % loi exponentielle

%%%% DENSITE NON CONTRAINTE %%%%

%net.K = 50;%net.N_aff; %K; % loi chi2 ou gauss (pas d'ajustement)

%%%% PARAMETRAGE DE LA DENSITE ET AJUSTEMENT DE SIGMA_J %%%%

net.densite=net.K./net.N_aff;          % densit� de connexion (entre 0 et 1)

net.sigma_J = sqrt(abs(net.sigma_J_eff.^2-(1-net.densite).*(net.J_barre.^2)./net.K));

%%%% AJUSTEMENTS LIES A LA DIMENSION %%%%
        
net.dim=1; %1;
     
net.rayon=0.2; %0.2;
        
net.grad=-1;

net.FLAG_GRAD='unif';

        
%%% DELAIS


net.tau_min=1;%net.nbp_d;%net.nbp_d-4;%(param. fred)              % net.tau_min : delai de transmission minimal
                 

net.tau_moy=net.nbp_d-1;%;%4;%0;% %(param. fred)              % net.tau_moy : parametre pour la distribution de poisson

                                           % (la valeur 0 indique des delais constants)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LES SEUILS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.theta_barre=1;        % net.theta_barre : valeur moyenne du seuil

net.sigma_theta=0;          % net.sigma_theta : ecart-type du seuil

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       APPRENTISSAGE                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

net.alpha=0.1;              % net.alpha : parametre d'apprentissage
                            % (regle l'intensit� de la modification des poids)


net.norm_alpha=1;     % net.norm : normalisation du terme d'apprentissage           
                       % (par population)
net.app='stdp';
net.FLAG_DELTA=0;
net.FLAG_RAZ=0;
net.FLAG_DASHBOARD=1;
net.FLAG_STD=0;
net.FLAG_REV=0;
net.FLAG_SCALING=0;
net.FLAG_STDP_EPS=0;
net.FLAG_BAYES=0;
net.FLAG_THETA=0;
%net.FLAG_GRAD='linear';




net.REWARD=1;                       

                   

