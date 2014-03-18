% Initialisation des seuils

function net=init_theta(net)
   
  % net.theta_barre : VECTEUR de taille nb_pop donnant le seuil moyen sur chaque couche 
  % net.sigma_theta : VECTEUR de taille nb_pop donnant l'ecart-type sur chaque couche
  % net.N : VECTEUR de taille nb_pop donnant les tailles des populations 
  % net.theta : CELLULE contenant p vecteurs de taille net.N(p)
 
  randn('seed',net.num_theta);
  nb_pop=length(net.N);

  for p=1:nb_pop
    net.theta{p}=net.theta_barre(p)-net.sigma_theta(p)*randn(net.N(p),1);
  end;