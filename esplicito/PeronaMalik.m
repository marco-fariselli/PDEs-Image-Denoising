%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image toolbox
% MATLAB file
% 
% PeronaMalikExample.m
% Perona Malik anisotropic filter
%
% function B = PeronaMalik( A, b, tfinal, outfilename )
%
%
%   u' = div[ g(u) grad(u) ]
%   du/dn |_{d Omega} = 0
%   u(t=0) = u0
%
%   g(u)= 1/ (1+b |grad(u)|^2)
%
% input:  A: image file name
%         b: parameter in the edge detector (b=1/k^2)
%         tfinal: artificial time (tempo a cui si vuole fermare il
%         processo)
%         outfilename: output file name
% output: B: denoised image
% example: B = gaussian( A, 0.1, 10 );
%
% created:       07/28/2009
% author:        syleung@gmail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B, peacksnr] = PeronaMalik( A, b, tfinal, outfilename )

fileout = 5;    % numero di passi che si intendono fare

dtout = tfinal/fileout;
time=0;
[m,n,p]=size(A);
xi=2:m+1;
yi=2:n+1;
B=zeros(m+2,n+2);
B(xi,yi) = A(:,:,1);
j=1; %serve per creare il vettore dei peacksnr

while time < tfinal
tic    
    % du/dn=0 on the boundary: vengono imposte le condizioni di Neumann al
    % contorno tramite l'allargamento della matrice dell'immagine di una
    % riga sopra, una sotto, una colonna a destra e una a sinistra, e
    % copiando la cornice più esterna dell'immagine originale anche nella
    % nuova cornice creata. In questo modo i pixel nelle prime due e ultime
    % due righe/colonne avranno lo stesso valore --> derivata numerica = 0.
    
    B(1,yi)= B(2,yi);
    B(end,yi)= B(end-1,yi);
    B(xi,1)= B(xi,2);
    B(xi,end)= B(xi,end-1);
        
    % g(|grad u|^2)
    
    C = gaussian(B,2,5);
    % viene prima filtrata in modo più superficiale e uniforme l'immagine
    % originale tramite un filtro gaussiano con matrice di convoluzione
    % settabile (in questo caso 5x5). (si veda function gaussian)
    
    [gn,ge,gs,gw]=g(C,b);
    % vengono create le varie g riferite ai gradienti della matrice C
    % (destro-w, sinistro-e, superiore-n, inferiore-s). (function g in
    % fondo a questo file)
    
    Bn = B(xi,yi+1);
    Be = B(xi+1,yi);
    Bs = B(xi,yi-1);
    Bw = B(xi-1,yi);
    Bc = B(xi,yi);
    
    % u^k+1 = u^k + dt * div[ g(u) grad(u) ]
    
    temp = gn.*Bn + ge.*Be + gs.*Bs + gw.*Bw ...
        - (gn+ge+gs+gw).*Bc;
    
    dt = 0.25*max(max( max(max(max(gn,ge),gs),gw) ));   
    %essendo g una funzione sicuramente minore uguale di 1, il dt massimo
    %non supererà 0.25 coerentemente a quanto detto per il metodo esplicito 
    dt = min(dt,tfinal-time);    
    %viene ulteriormente ridotto nel caso in cui, dopo tot passi, ci si trovi a
    %un passo inferiore di 0.25 dal tfinal.
    time=time+dt;
  %  [time dt tfinal]
    
    B(xi,yi)=B(xi,yi)+dt*temp;

    % check for intermediate output
    
    if (floor(time/dtout) ~= floor((time-dt)/dtout) )
        figure
        filename=[outfilename num2str(floor(time/dtout)) '.jpg'];
        imwrite(uint8(B(xi,yi)),filename)
        peacksnr(j)=psnr(B(2:m+1,2:n+1)/255, A/255);
        imshow(uint8(B(xi,yi))), title(filename), xlabel(strcat('psnr = ',string(peacksnr(j))));
        j=j+1;
        toc
    end  
end

B=B(2:end-1,2:end-1);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gn,ge,gs,gw]=g(u,b)

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

duxc=(duxp+duxm)/2; %matrici delle derivate centrali rispetto a x e y
duyc=(duyp+duym)/2;

gn=1./(1+b*(duxc.^2+duyp.^2));
ge=1./(1+b*(duxp.^2+duyc.^2));
gs=1./(1+b*(duxc.^2+duym.^2));
gw=1./(1+b*(duxm.^2+duyc.^2));

return