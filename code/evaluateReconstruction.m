function [rMSE, PSNR]=evaluateReconstruction(DRecon,DRefer)
% evaluate reconstruction by relative mean square error and peak signal
% noise ratio

% input
% DRecon: the reconstructed data
% DRefer: the reference data

% output
% rMSE: relative mean square error
% peak signal noise ratio

n=length(DRecon(:));
rMSE=norm(DRecon(:)-DRefer(:))/norm(DRefer(:));
MSE=sum((DRecon(:)-DRefer(:)).*conj(DRecon(:)-DRefer(:)))/n;
% PSNR=10*log10(max(abs(DRefer(:)))^2/MSE);   
PSNR=10*log10(max(abs(DRefer(:)))^2/MSE);   
