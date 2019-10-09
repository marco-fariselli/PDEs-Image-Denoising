function [gn,ge,gs,gw]=my_g(u,b)

% g(u)= 1/ (1+b |grad(u)|^2)

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

gn=1./(1+b*(duxc.^2+duyp.^2));
ge=1./(1+b*(duxp.^2+duyc.^2));
gs=1./(1+b*(duxc.^2+duym.^2));
gw=1./(1+b*(duxm.^2+duyc.^2));

% gn=exp(-sqrt(b*(duxc.^2+duyp.^2)));
% ge=exp(-sqrt(b*(duxp.^2+duyc.^2)));
% gs=exp(-sqrt(b*(duxc.^2+duym.^2)));
% gw=exp(-sqrt(b*(duxm.^2+duyc.^2)));

return