function [ edge_medio, modgrad ] = edge_grad_implicito( u )
%UNTITLED Summary of this function goes here
%   Data la matrice dell'immagine u restituisce:
%   edge_medio= valore medio del gradiente di X,
%   grad= matrice del modulo del gradiente ottenuto per differenze centrali sia nella
%           direzione x che nella direzione y.

[m,n]=size(u);
xi=[2:m-1];
yi=[2:n-1];

uxp=u(xi+1,yi);
uxm=u(xi-1,yi);
uyp=u(xi,yi+1);
uym=u(xi,yi-1);
uc=u(xi,yi);

duxp=uxp-uc;
duxm=uc-uxm;
duyp=uyp-uc;
duym=uc-uym;

duxc=(duxp+duxm)/2;
duyc=(duyp+duym)/2;

modsqrgrad=duxc.^2+duyc.^2;
modgrad=sqrt(modsqrgrad);

edge_medio=mean(mean(modgrad));

end

