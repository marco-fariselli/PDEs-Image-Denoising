function [ A ] = semi_implicit_matrix( u, dt, gn, ge, gs, gw )
% input: 
%   -u : immagine sotto forma di matrice
%   -b : coefficiente interno a g (b=1/K^2)
%   -g = [gn, ge, gs, gw] : matrice della diffusione
%   -dt: passo temporale

[rows,cols]=size(u);
[m,n]=size(gn);
%metto le matrici g in dei vettori ordinati per righe
gn=reshape(gn', 1, m*n);
ge=reshape(ge', 1, m*n);
gs=reshape(gs', 1, m*n);
gw=reshape(gw', 1, m*n);

diag0=1+dt.*(gn+ge+gw+gs);

diage=-dt.*ge;
diagw=-dt.*gw;
diags=-dt.*gs;
diagn=-dt.*gn;

diage=[0,diage];
diagw=[diagw,0];
diags=[zeros(1,rows),diags];
diagn=[diagn,zeros(1,rows)];

for i=1:cols-1
    diage(i*rows+1)=0;
    diagw(i*rows+1)=0;
end

B=zeros(length(diag0),5);
B(:,1)=diag0;
B(:,2)=diage(1:length(diage)-1);
B(:,3)=diagw(2:length(diagw));
B(:,4)=diags(1:length(diags)-rows);
B(:,5)=diagn(rows+1:length(diagn));

d=[0, 1, -1, rows, -rows];

A=spdiags(B, d, rows*cols, rows*cols);
