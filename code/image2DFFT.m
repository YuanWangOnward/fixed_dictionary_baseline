function coefficient = image2DFFT(image)
% do 2D discreat Fourier transform on a image
% result is normalized so that energy is preserved 
[R, C]=size(image);
coefficient=fft(image,[],1)/sqrt(R);
coefficient=fft(coefficient,[],2)/sqrt(C);