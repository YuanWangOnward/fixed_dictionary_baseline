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
% projectionMatrix
% backProjectionMatrix
% tolerance
% maxIterationNumber
% L

lambda = parameters.lambda;
mu = parameters.mu;
projectionMatrix = parameters.projectionMatrix;
backProjectionMatrix = parameters.backProjectionMatrix;
tolerance = parameters.tolerance;
maxIterationNumber = parameters.maxIterationNumber;
L = parameters.L;    % used for FISTA 


%% configure
ifPlot=false; % whether plot images during the reconstruction
ifDisplay=false; % whether display debug information during the reconstrution
ifSpeedTest=true;
timeArray=zeros(maxIterationNumber,1);
startTime=cputime;

%% set parameters
imageSize = size(measure);
%% initialization
measure=double(measure);
downSamplingPattern=abs(measure)>0;

cost=zeros(maxIterationNumber,1);
sparseCost=zeros(maxIterationNumber,1);
coreCell=zeros(maxIterationNumber,1);


xkm1=zeros(size(projectionMatrix));   % series x, x_k-1
yk=zeros(size(projectionMatrix));     % series y, y_k
tk=1;   % series t, t_k
nIte=0;

%% reconstruction
% First step, FISTA based, find the optimal sparse coefficient
while(nIte<maxIterationNumber)
    
    nIte=nIte+1;
    if(mod(nIte,10)==0)
        display(['nIte=',num2str(nIte)]);
    end
    
    %%%% updata sparse coefficient sequence one xk %%%%
    xk = updateSparseCoefficient(yk,measure,dictionary,parameters);
%     imageBlockRepresentation=B*yk;
%     imageImageRepresentation = blocks2image(imageBlockRepresentation, BPM,imageSize);
%     kSpaceCoefficientCurrent = image2DFFT(imageImageRepresentation);
%     % data consistancy
%     dk=(measure-kSpaceCoefficientCurrent).*downSamplingPattern;
%     kSpaceCoefficientNext=kSpaceCoefficientCurrent+1/L*dk;
%     % find sparse coefficients
%     imageImageRepresentation = image2DIFFT(kSpaceCoefficientNext);
%     imageBlockRepresentation = image2blocks(imageImageRepresentation,PM);
%     xk=B'*imageBlockRepresentation;
%     % soft thresholding
%     xk=SoftThresh(xk,lambda/2/L);
    
    %%%% update tkp1 %%%%
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    %%%% update sparse coefficient sequence two ykp1 %%%%
    ykp1=xk+((tk-1)/tkp1)*(xk-xkm1);
    
    
    if(ifPlot)
        % plot k-space before and after data consistancy
        % plot only the first slice of
        figure(1)
        subplot(1,2,1)
        imagesc(log(abs(kSpaceCoefficientCurrent(:,:,1))))
        title('before data consistancy')
        subplot(1,2,2)
        imagesc(log(abs(kSpaceCoefficientNext(:,:,1))))
        title('after data consistancy')
        shg
        %pause(0.1)
    end
    
    if(ifDisplay)
        display(['norm measure=    ',num2str(norm(measure(downSamplingPattern(:)>0)))]);
        display(['norm kCurrent=  ', num2str(norm(kSpaceCoefficientCurrent(downSamplingPattern(:)>0)))]);
        display(['norm kNext=  ', num2str(norm(kSpaceCoefficientNext(downSamplingPattern(:)>0)))]);
        display(['norm dk=  ', num2str(norm(dk(downSamplingPattern(:)>0)))]);
        display(['rMSE (after data consistancy)=  ', num2str(norm(dk(downSamplingPattern(:)>0))/norm(measure(downSamplingPattern(:)>0)))]);
    end
    
    if(ifPlot)
        % show coefficinet core before and after soft thresholding
        % only show the first slice
        figure(2)
        subplot(1,2,1)
        imagesc(log(abs(yk(:,:,1))))
        title('beginning core')
        shg
        subplot(1,2,2)
        imagesc(log(abs(ykp1(:,:,1))))
        title('modified core')
        shg
    end
    
    if(ifPlot)
        % reconstruct image and k-space data
        imageBlockRepresentation=dictionary*ykp1;
        %         image=blockMatrix2subjectArray(temp,size(measure),blockSize);
        imageImageRepresentation = blocks2image(imageBlockRepresentation, backProjectionMatrix,imageSize);
        kSpaceCoefficientNext=double(ttm(tensor(imageImageRepresentation),F));
        
        % plot reconstructed image
        % plot only the first frame
        figure(3)
        %         subplot(1,2,1)
        imagesc(abs(imageImageRepresentation(:,:,1)))
        colormap(gray);
        title(['reconstruction after ', num2str(nIte),' iterations'])
        %         subplot(1,2,2)
        %         imagesc(abs(image(:,:,1))-abs(iOriginal(:,:,1)))
        %         colormap(gray);
        %         title(['reconstruction error after ', num2str(nIte),' iterations'])
    end
    
    
    % calculate statistics
    if(ifSpeedTest==false)
        dk=(measure-kSpaceCoefficientNext).*downSamplingPattern;
        kRMSE(nIte)=(norm(dk(downSamplingPattern(:)>0))/norm(measure(downSamplingPattern(:)>0)));
        cost(nIte)=norm(dk(downSamplingPattern(:)>0)).^2+lambda*(sum(abs(coreNext(:))>0));
        sparseCost(nIte)=sum(abs(coreNext(:))>0);
    end
    
    % check stop condition
    if (norm(ykp1(:)-yk(:))/norm(yk(:))<tolerance)
        display('BP ends because it has reached error tolerance')
        break;
    end
    
    
    %%%% prepare for next loop %%%%
    tk=tkp1;
    xkm1=xk;
    yk=ykp1;
    
    
    timeArray(nIte)=cputime-startTime;
    %     kSparse=double(ttm(tensor(iTemp),F));
    %     kSparse(abs(measure)>0)=(lambda*kSparse(abs(measure)>0)+measure(abs(measure)>0))/(lambda+1);
    %     iTemp=double(ttm(tensor(kSparse),F,'t'));
    %     coreCell(nIte)=norm(abs(iTemp(:))-abs(iOriginal(:)))/norm(iOriginal(:));
    
end

% Second step, reconstruct target image
% reconstruct image by balance k-space measurement and sparse representation
imageImageRepresentation = reconstructImage(ykp1,dictionary,measure,mu,backProjectionMatrix);
% imageBlockRepresentation = B*ykp1;
% imageImageRepresentation = blocks2image(imageBlockRepresentation, BPM,imageSize);
% 
% kSpaceCoefficientCurrent = double(ttm(tensor(imageImageRepresentation),F));
% kSpaceCoefficientCurrent(abs(measure)>0)=(kSpaceCoefficientCurrent(abs(measure)>0)+mu*measure(abs(measure)>0))/(1+mu);
% imageImageRepresentation=double(ttm(tensor(kSpaceCoefficientCurrent),F,'t'));

end
