
% Lissage du diagramme d'activité DYN_X selon la résolution 

function DYN_Y = lissage_lif(net,DYN_X)

[N,nbp]=size(DYN_X);
l=floor(0.5*mean(net.nbp_mb{1})/net.delta_t);
DYN_Y=zeros(N,nbp);

for i=1:N
    trace=conv(DYN_X(i,:),ones(1,l));
    DYN_Y(i,:)=trace(l/2+1:nbp+l/2);
end;