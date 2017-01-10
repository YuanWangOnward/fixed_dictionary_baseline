function blockMatrix = image2blocks(image,projectionMatrix)
% Transfer image from natual image representation to block  representation.
%
% Natual image representation represents an image as a 2D matrix
%
% Block representation of an image is also a 2D matrix, each column is a 
% vectorized block cut from the image. So the row number of the matrix is
% determined by the block size; the column number of the matrix is
% determined by image size and block stride. Notably, blocks can have
% overlapping in the image

blockMatrix=image(projectionMatrix);

end