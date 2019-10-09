close all
clear all

M=menu('Seleziona il tipo di immagine', 'fotografica (cameraman)', 'medica (radiografia mano)');

switch M
    case 1
        I_original=imread('Cameraman256.png'); % load image
        imgstrg='cameraman semi implicito ';
    case 2
        I_original=imread('mano.jpg'); % load image
        imgstrg='mano semi implicito ';
end

I_original=(I_original(:,:,1));
I_noise=imnoise(uint8(I_original),'gaussian',0,0.01);
I_noise=double(I_noise); % cut a piece, convert to double
I_original=double(I_original);

[rows,cols]=size(I_noise);

I_iter=I_noise;
numpassi=10;
peack_snr=zeros(1,numpassi);
peack_snr_iniziale=psnr(I_noise,I_original,255);
edge_mean=zeros(1,numpassi);

for h=1:numpassi
  tic  
%creo la matrice dell'immagine con le condizioni di Neumann
    xi=2:rows+1;
    yi=2:cols+1;
    I_neumann_pre_filt=zeros(rows+2,cols+2);
    I_neumann_pre_filt(xi,yi) = I_iter(:,:,1);

    I_neumann_pre_filt(1,yi)= I_neumann_pre_filt(3,yi);
    I_neumann_pre_filt(end,yi)= I_neumann_pre_filt(end-2,yi);
    I_neumann_pre_filt(xi,1)= I_neumann_pre_filt(xi,3);
    I_neumann_pre_filt(xi,end)= I_neumann_pre_filt(xi,end-2);
    
   % faccio un filtraggio blando con un filtro gaussiano
   I_neumann_post_filt=imgaussfilt(I_neumann_pre_filt,0.001);
    
    %calcolo il valore medio del gradiente e lo uso per determinare k
    [edge_mean(h), modgrad] = edge_grad_implicito(I_neumann_post_filt);
    
    switch M
        case 1        
            %per immagine normale
            K=0.75*edge_mean(h);
            dt=0.75;
        case 2 
            %per immagine medica
            K=3*edge_mean(h);
            dt=0.75;
    end
    
    b=1/K^2;  %b=1/k^2
   
    [gn,ge,gs,gw]=my_g(I_neumann_post_filt,b);
    mat_term_noti=I_iter;

    %% Impongo le Neumann nel vettore dei termini noti
    
    %Aggiungo i 4 vettori che andranno poi nei termini noti per imporre le Neumann

    for i=1:rows
        mat_term_noti(i,1)=I_iter(i,1)+dt*gw(i,1)*I_iter(i,1);
        mat_term_noti(i,cols)=I_iter(i,cols)+dt*ge(i,cols)*I_iter(i,cols);    
    end
    for i=1:cols
        mat_term_noti(1,i)=I_iter(1,i)+dt*gn(1,i)*I_iter(1,i);
        mat_term_noti(rows,i)=I_iter(rows,i)+dt*gs(rows,i)*I_iter(rows,i);    
    end

%% Sistemo la matrice in un vettore (row-wise)

    vett_term_noti=reshape(mat_term_noti',cols*rows,1);

%% Costruzione della matrice e risoluzione del sistema lineare

%costruisco la matrice A e risolvo il sistema A* vett_new_image=vett_term_noti
%utilizzando il metodo del gradiente coniugato precondizionato già implementato da Matlab

    A  = semi_implicit_matrix( I_noise, dt, gn, ge, gs, gw );

    tol=1e-7;
    maxit=50;
    [vett_new_image,flags,relres,iter]=pcg(A,vett_term_noti,tol,maxit); % senza precondizionatore

%% Riordino il vettore creando una matrice

%costruisco l'immagine nuova risistemando il vettore new_image in una
%matrice della stessa dimensione di quella di partenza

    new_image=reshape(vett_new_image,rows,cols);
    new_image=new_image';
    
    I_iter=new_image;
    peack_snr(h)=psnr(I_iter,I_original,255);

    figure
    imshow(uint8(new_image)),title(strcat(imgstrg,strcat(' ',string(h),' iterazione'))), xlabel(strcat('psnr=',string(peack_snr(h))));
toc
end

figure
imshow(uint8(I_original)),title(strcat(imgstrg,' originale'));
figure 
imshow(uint8(I_noise)),title(strcat(imgstrg,' con rumore')),xlabel(strcat('psnr=',string(peack_snr_iniziale)));
figure
plot(peack_snr), xlabel('passo'), ylabel('peack SNR');
