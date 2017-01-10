function image = blocks2image(blockMatrix, backProjectionMatrix,imageSize)
% Transfer image from block representation to natual image representation.
%
% Natual image representation represents an image as a 2D matrix
%
% Block representation of an image is also a 2D matrix, each column is a 
% vectorized block cut from the image. So the row number of the matrix is
% determined by the block size; the column number of the matrix is
% determined by image size and block stride. Notably, blocks can have
% overlapping in the image

image=blockMatrix(backProjectionMatrix);
image=mean(image,1);
image=reshape(image,imageSize);
end