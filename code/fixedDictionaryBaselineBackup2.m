function [imageImageRepresentation,timeArray]=fixedDictionaryBaseline(measure,dictionary,parameters,iOriginal)
% Fixed-Dictionay Baseline MRI compressive sensing reconstruction algorithm
%
% input
% measure: k-space measurement, in multiway array structure of the same
% size as target image, un-sampled location is zero
% B: block  bases
% lambda: the weighting factor of sparse penalty
% blockSize: size of each block
% stepSize: step size used when partition the image into blocks
% tolerance: if image changing becomes smaller than it, stop
% PM: projection matrix, segment image into block matrix
% BPM: back projection matrix, rearrange block matrix into image
%
% ouput
% image: reconstructed image
% coreCell: for calculate cost function
% timeArray: record time consumption
%
% parameters is a structure including atributes:
% lambda
% mu
% projectMatrix
% backProjectMatrix
% tolerance
% maxIterationNumber
% L
% ifPlot

% parameters.mu = parameters.mu;
% parameters.projectMatrix = parameters.projectMatrix;
% parameters.backProjectMatrix = parameters.backProjectMatrix;
% parameters.stopTolerance = parameters.stopTolerance;
% parameters.maxIterationNumber = parameters.maxIterationNumber;
% L = parameters.L;    % used for FISTA


%% configure
%ifSpeedTest=true;
%timeArray=zeros(parameters.maxIterationNumber,1);
%startTime=cputime;

%% initialization
% cost=zeros(parameters.maxIterationNumber,1);
% sparseCost=zeros(parameters.maxIterationNumber,1);
% coreCell=zeros(parameters.maxIterationNumber,1);

measure=double(measure);
xkm1=zeros(size(parameters.projectMatrix));   % series x, x_k-1
yk=zeros(size(parameters.projectMatrix));     % series y, y_k
tk=1;   % series t, t_k
nIte=0;

%% reconstruction
% First step, FISTA based, find the optimal sparse coefficient
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
    
    if(parameter.ifPlot)
        imageImageRepresentation = reconstructImage(ykp1,dictionary,measure,parameters.mu,parameters.backProjectMatrix);
        imagesc(abs(imageImageRepresentation(:,:,1)))
        colormap(gray);
        title(['reconstruction after ', num2str(nIte),' iterations'])
    end
    
    
    % calculate statistics
    %     if(ifSpeedTest==false)
    %         dk=(measure-kSpaceCoefficientNext).*downSamplingPattern;
    %         kRMSE(nIte)=(norm(dk(downSamplingPattern(:)>0))/norm(measure(downSamplingPattern(:)>0)));
    %         cost(nIte)=norm(dk(downSamplingPattern(:)>0)).^2+lambda*(sum(abs(coreNext(:))>0));
    %         sparseCost(nIte)=sum(abs(coreNext(:))>0);
    %     end
    
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
    
    
    %     timeArray(nIte)=cputime-startTime;
    %     kSparse=double(ttm(tensor(iTemp),F));
    %     kSparse(abs(measure)>0)=(lambda*kSparse(abs(measure)>0)+measure(abs(measure)>0))/(lambda+1);
    %     iTemp=double(ttm(tensor(kSparse),F,'t'));
    %     coreCell(nIte)=norm(abs(iTemp(:))-abs(iOriginal(:)))/norm(iOriginal(:));
    
end

% Second step, reconstruct target image
% reconstruct image by balance k-space measurement and sparse representation
imageImageRepresentation = reconstructImage(ykp1,dictionary,measure,parameters.mu,parameters.backProjectMatrix);

end
