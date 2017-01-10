function image = image2DIFFT(coefficient)
% do 2D inverse discreat Fourier transform on a image
% result is normalized so that energy is preserved 
[R, C]=size(coefficient);
image=ifft(coefficient,[],1)*sqrt(R);
image=ifft(image,[],2)*sqrt(C);