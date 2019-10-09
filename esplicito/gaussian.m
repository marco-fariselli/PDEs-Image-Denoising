%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image toolbox
% MATLAB file
% 
% gaussian.m
% Gaussian image denoising filer
%
% function B = gaussian( A, sd, R)
%
% input:  A: image file name; 
%         sd: standard deviation of the Gaussian kernal;
%         R: size of the local mask
% output: B: denoised image
% example: B = gaussian( A, 2, 5 );
%
% created:       07/23/2009
% author:        syleung@gmail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = gaussian(A,sd,R)

B = zeros(size(A));

for i = 1:size(A,1)
for j = 1:size(A,2)
    zsum = 0;
    for m = -R:R
        if i+m >= 1 && i+m <= size(A,1)
            for n = -R:R
                if j+n >= 1 && j+n <= size(A,2)
                    z = exp(-(m^2 + n^2)/(2*sd^2));
                    zsum = zsum + z;
                    B(i,j) = B(i,j) + z*A(i+m,j+n);
                end
            end
        end
    end
    B(i,j) = B(i,j)/zsum;
end
end
