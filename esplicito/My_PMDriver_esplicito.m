clear all
close all
format compact
format long

M=menu('Seleziona il tipo di immagine', 'fotografica (cameraman)', 'medica (radiografia mano)');

switch M
    case 1
        I_original=imread('Cameraman256.png'); % load image
        imgstrg='cameraman esplicito ';
    case 2
        I_original=imread('mano.jpg'); % load image
        imgstrg='mano esplicito ';
end


I_original=I_original(:,:,1);
I_noise=imnoise(uint8(I_original),'gaussian',0,0.01);  
%aggiungo del rumore all'immagine originale facendo finta
%che sia questa l'immagine a mia disposizione:
%tramite delle elaborazioni dovrò cercare di tornare il più possibile a
%quella di partenza.

I_noise=double(I_noise); % cut a piece, convert to double
I_original=double(I_original);

[ edge_mean2, modgrad2]=edge_grad_2(I_noise);

 switch M
%K consigliato in base a prove sperimentali:
        case 1        
            %per immagine normale
            K=edge_mean2/20;
        case 2 
            %per immagine medica
            K=edge_mean2/15;
 end
    
    b=1/K^2;  %b=1/k^2
%b consigliato in base a prove sperimentali:

[B2, peacksnr2] = PeronaMalik( I_noise, b, 10, strcat(imgstrg,'2') );
%guardando il peacksnr delle immagini filtrate sembrerebbe migliore quello
%di PM21, poichè, anche se di poco, più grande rispetto alle altre immagini.
%All'occhio umano risulta però migliore l'immagine PM25.
