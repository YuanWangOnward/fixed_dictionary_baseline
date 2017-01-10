function imageReconstructed = reconstructImage(sparseCoefficient,dictionary,measurement,mu,backProjectionMatrix)
% With given sparse coefficient, reconstruct image by balance k-space 
% measurement and sparse representation
% It's a subproblem of MRI CS reconstruction

imageSize = size(measurement);

imageBlockRepresentation = dictionary*sparseCoefficient;
imageImageRepresentation = blocks2image(imageBlockRepresentation, backProjectionMatrix,imageSize);
kSpaceCoefficient = image2DFFT(imageImageRepresentation);

% modify k-space coefficient with measurement
samplePattern = abs(measurement)>0;
kSpaceCoefficient(samplePattern)=(kSpaceCoefficient(samplePattern)+mu*measurement(samplePattern))/(1+mu);
imageReconstructed = image2DIFFT(kSpaceCoefficient);

