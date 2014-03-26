% Initialisation de la dynamique
% Mod�le I&F (octobre 2006)

function net=clear_dyn_s(net);
      
  for p=1:net.nb_pop
      for q=1:net.nb_pop
          net.DYN_H{p}{q}=[]; 
      end;
      net.DYN_S{p}=[]; 
      net.DYN_U{p}=[]; 
      net.DYN_V{p}=[]; 
      net.DYN_I{p}=[]; 
  end;
  
net.FLAG_CLEAR=1;
