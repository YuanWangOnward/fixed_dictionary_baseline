function [PM, BPM]=createProjectionMatrix(frameSize,blockSize,stepSize)
% Create operator matrices used to transfer image from/to natural image  
% representation to/from block representation.
%
% Natural image representation treats an image as a 2D matrix which can be
% shown and watch by people directly.
% Block representation is also a 2D matrix, but each column is a 
% vectorized block cut from the image. So the row number of the matrix is
% determined by the block size; the column number of the matrix is
% determined by image size and block stride. Notably, blocks can have
% overlapping in the image.
%
% Input
% frameSize: image size, example: [512, 512]
% blockSize: the size of block, example: [8, 8]
% stepSize: the stride used to move the block in the image, example: [1, 1]
%
% Output
% PM: project matrix
% BPM: back project matrix
%
% The underlying conceptual idea is easy,
% but the coding feels a little hacky.
%
% Usage:
% From image representation to block representation:
% blockRepresentation=imageRepresentation(PM);
% From block representation to image representation:
% temp=blockRepresentation(BPM);
% temp=mean(temp,1);
% imageRepresentation=reshape(temp,size(frameSize));
%
% -------------------------------------------------------------------------
% Written by Yuan Wang, New York University, 2016, yw1225@nyu.edu
% -------------------------------------------------------------------------


PM=zeros(prod(blockSize),prod(frameSize)/prod(stepSize));
BPM=zeros(prod(blockSize)/prod(stepSize),prod(frameSize));
indexArray=zeros(1,prod(frameSize));

count=0;
for r=1:stepSize(1):frameSize(1)
    for c=1:stepSize(2):frameSize(2)
        rRange=r:(r+blockSize(1)-1);
        cRange=c:(c+blockSize(2)-1);
        
        rRange(rRange<1)=rRange(rRange<1)+frameSize(1);
        cRange(cRange<1)=cRange(cRange<1)+frameSize(1);
        
        rRange(rRange>frameSize(1))=rRange(rRange>frameSize(1))-frameSize(1);
        cRange(cRange>frameSize(1))=cRange(cRange>frameSize(1))-frameSize(1);
        
        iBPM=repmat((cRange-1)*frameSize(1),blockSize(1),1)+repmat(rRange',1,blockSize(2));
        iPM=count*prod(blockSize)+[1:prod(blockSize)]';
        
        indexArray(iBPM(:))=indexArray(iBPM(:))+1;
        indexBPM=(iBPM(:)-1)*size(BPM,1)+indexArray(iBPM(:))';
        BPM(indexBPM)=iPM;
        
        count=count+1;
        PM(:,count)=iBPM(:);
    end
end

end
