function [image,coreCell,timeArray]=fixedDictionaryBaseline(measure,B,lambda,mu,PM,BPM, tolerance,NIte,iOriginal)
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


%% configure
ifPlot=true; % whether plot images during the reconstruction
ifDisplay=false; % whether display debug information during the reconstrution
ifSpeedTest=true;
timeArray=zeros(NIte,1);
startTime=cputime;
%% set parameters
L=1;    % parameter used for ISTA to contral step size, the larger the smaller the step

%% initialization
measure=double(measure);
% kMask=SP;
kMask=abs(measure)>0;

cost=zeros(NIte,1);
sparseCost=zeros(NIte,1);
coreCell=zeros(NIte,1);


% find bases to transfor data between k-space and image domain
F={};
for ii=1:ndims(measure)
    F{ii}=fft(eye(size(measure,ii)),[],1)/sqrt(size(measure,ii));
end

xkm1=zeros(size(PM));   % series x, x_k-1
yk=zeros(size(PM));     % series y, y_k

tk=1;   % series t, t_k
nIte=0;

%% reconstruction
% nIte=0;
while(nIte<NIte)
    
    nIte=nIte+1;
    if(mod(nIte,1)==0)
    display(['nIte=',num2str(nIte)]);
    end
    
     %%%% updata xk %%%%
    %----------------- pL ----------------------
    % for pL, see FISTA by BECK2009
%     iTemp=blockMatrix2subjectArray(B*coreNext,size(measure),blockSize);

    iBlockMatrix=B*yk;
    iTemp=iBlockMatrix(BPM);
    iTemp=mean(iTemp,1);
    iTemp=reshape(iTemp,size(measure));
    kCurrent=double(ttm(tensor(iTemp),F));
    % data consistancy
    dk=(measure-kCurrent).*kMask;
    kNext=kCurrent+1/L*dk;
    % find coefficients
    iTemp=double(ttm(tensor(kNext),F,'t'));
%     iTemp=subjectArray2BlockMatrix(iTemp,blockSize);
    iBlockMatrix=iTemp(PM);
    xk=B'*iBlockMatrix;
    %---------------------------------------------
    % soft thresholding
    xk=SoftThresh(xk,lambda/2/L);
    
    %%%% update tkp1 %%%%
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    %%%% update ykp1 %%%%
    ykp1=xk+((tk-1)/tkp1)*(xk-xkm1);
    
    
    if(ifPlot)
        % plot k-space before and after data consistancy
        % plot only the first slice of
        figure(1)
        subplot(1,2,1)
        imagesc(log(abs(kCurrent(:,:,1))))
        title('before data consistancy')
        subplot(1,2,2)
        imagesc(log(abs(kNext(:,:,1))))
        title('after data consistancy')
        shg
        %pause(0.1)
    end
    
    if(ifDisplay)
        display(['norm measure=    ',num2str(norm(measure(kMask(:)>0)))]);
        display(['norm kCurrent=  ', num2str(norm(kCurrent(kMask(:)>0)))]);
        display(['norm kNext=  ', num2str(norm(kNext(kMask(:)>0)))]);
        display(['norm dk=  ', num2str(norm(dk(kMask(:)>0)))]);
        display(['rMSE (after data consistancy)=  ', num2str(norm(dk(kMask(:)>0))/norm(measure(kMask(:)>0)))]);
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
        iBlockMatrix=B*ykp1;
        %         image=blockMatrix2subjectArray(temp,size(measure),blockSize);
        image=iBlockMatrix(BPM);
        image=mean(image,1);
        image=reshape(image,size(measure));
        kNext=double(ttm(tensor(image),F));
        
        % plot reconstructed image
        % plot only the first frame
        figure(3)
        %         subplot(1,2,1)
        imagesc(abs(image(:,:,1)))
        colormap(gray);
        title(['reconstruction after ', num2str(nIte),' iterations'])
        %         subplot(1,2,2)
        %         imagesc(abs(image(:,:,1))-abs(iOriginal(:,:,1)))
        %         colormap(gray);
        %         title(['reconstruction error after ', num2str(nIte),' iterations'])
    end
    
    
    % calculate statistics
    if(ifSpeedTest==false)
        dk=(measure-kNext).*kMask;
        kRMSE(nIte)=(norm(dk(kMask(:)>0))/norm(measure(kMask(:)>0)));
        cost(nIte)=norm(dk(kMask(:)>0)).^2+lambda*(sum(abs(coreNext(:))>0));
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
    kSparse=double(ttm(tensor(iTemp),F));
    kSparse(abs(measure)>0)=(lambda*kSparse(abs(measure)>0)+measure(abs(measure)>0))/(lambda+1);
    iTemp=double(ttm(tensor(kSparse),F,'t'));
    coreCell(nIte)=norm(abs(iTemp(:))-abs(iOriginal(:)))/norm(iOriginal(:));
    
end

% find output image by balance k-space measurement and sparse
% representation Ravishankar2011
iBlockMatrix=B*ykp1;
image=iBlockMatrix(BPM);
image=mean(image,1);
image=reshape(image,size(measure));

kSparse=double(ttm(tensor(image),F));
kSparse(abs(measure)>0)=(kSparse(abs(measure)>0)+mu*measure(abs(measure)>0))/(1+mu);
image=double(ttm(tensor(kSparse),F,'t'));








end
