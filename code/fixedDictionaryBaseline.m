function [imageImageRepresentation]=fixedDictionaryBaseline(measure,dictionary,parameters)
% Fixed-Dictionay Baseline MRI compressive sensing reconstruction algorithm
%
% Input
% measure:    k-space measurement, in 2D array, of the same size as the  
%             target image, un-sampled location is zero
% dictionary: sparsify dictionary
% parameters: a structure containing the follow attributes 
%   (mandatory)
%   -lambda:             the weighting factor of sparse penalty
%   -mu:                 the weighting factor for measurement
%   -projectMatrix:      reshaping and partition operator matrices, 
%                        please see createProjectionMatrix.m
%   -backProjectMatrix:  reshaping and partition operator matrices, 
%                        please see createProjectionMatrix.m
%   -stopTolerance:      if relative change of sparse coefficients in 
%                        adjacent iterations is smaller than this value, 
%                        reconstruction stops
%   -maxIterationNumber: maximal iterations allowed
%   (optional)
%   -L:                  constant used by FISTA, control convergence&speed
%   -ifPlot:             whether show intermediate image in each iteration
%
% Ouput
% imageImageRepresentation: reconstructed image
%
% -------------------------------------------------------------------------
% Written by Yuan Wang, New York University, 2016, yw1225@nyu.edu
% -------------------------------------------------------------------------


%% check parameters
if(~isfield(parameters,'L')) parameters.L=1; end
if(~isfield(parameters,'ifPlot')) parameters.ifPlot=0; end

%% initialization
measure=double(measure);
xkm1=zeros(size(parameters.projectMatrix));   % series x, x_k-1
yk=zeros(size(parameters.projectMatrix));     % series y, y_k
tk=1;   % series t, t_k
nIte=0;

%% reconstruction
% First step, find the optimal sparse coefficient
while(nIte<parameters.maxIterationNumber)
    
    nIte=nIte+1;
    if(mod(nIte,10)==0)
        display(['nIte=',num2str(nIte)]);
    end
    
    %%%% updata sparse coefficient sequence one xk %%%%
    xk = updateSparseCoefficient(yk,measure,dictionary,parameters);
    
    %%%% update tkp1 %%%%
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    %%%% update sparse coefficient sequence two ykp1 %%%%
    ykp1=xk+((tk-1)/tkp1)*(xk-xkm1);
    
    if(parameters.ifPlot)
        imageImageRepresentation = reconstructImage(ykp1,dictionary,measure,parameters.mu,parameters.backProjectMatrix);
        imagesc(abs(imageImageRepresentation(:,:,1)))
        colormap(gray);
        title(['reconstruction after ', num2str(nIte),' iterations'])
        pause(0.01); %Without it, matlab will not plot in the recon process
    end
    
    % check stop condition
    if (norm(ykp1(:)-yk(:))/norm(yk(:))<parameters.stopTolerance)
        ST = dbstack;
        currentFunctionName = ST.name;
        display([currentFunctionName,'() ends because it has reached stop tolerance'])
        break;
    end
    
    %%%% prepare for next loop %%%%
    tk=tkp1;
    xkm1=xk;
    yk=ykp1;  
end

% Second step, reconstruct image
imageImageRepresentation = reconstructImage(ykp1,dictionary,measure,parameters.mu,parameters.backProjectMatrix);

end
