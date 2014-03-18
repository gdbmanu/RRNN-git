
% Initialisation des poids et des seuils 
% Mod�le I&F (octobre 2006) 

function net=init_systeme_lif(net,num_reseau)

  % Remarque : num_reseau sert sert d'amorce pour le generateur pseudo-al�atoire de Matlab 
  % Il fixe donc la distribution des poids et des seuils
  % et d�termine ainsi l'"identit�" du r�seau
  
  net.num_reseau=num_reseau;
  net.num_theta=num_reseau;
  
  % Appel � la fonction init_theta (initialisation des seuils)
  net=init_theta(net);
  % Appel � la fonction init_J (initialisation des poids)  
  net=init_J_lif(net);
  
  % Execution d'un script d'initialisation si celui-ci est d�fini
  if strcmp(net.script_init,'')==0  
      eval(net.script_init);
  end;
