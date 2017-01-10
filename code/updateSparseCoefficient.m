function sparseCoefficientOut = updateSparseCoefficient(sparseCoefficientIn,measurement,dictionary,parameters)
% Update sparse coefficient by doing measurement consistancy and soft thresholding 
% It's a subproblem of MRI CS reconstruction
% 
% parameters should have attributes:
% projectionMatrix
% backProjectionMatrix
% L
% lambda


imageSize = size(measurement);
downSamplingPattern = abs(measurement)>0;

% from sparse coefficient to k-space coefficient
imageBlockRepresentation = dictionary*sparseCoefficientIn;
imageImageRepresentation = blocks2image(imageBlockRepresentation,parameters.backProjectMatrix,imageSize);
kSpaceCoefficientCurrent = image2DFFT(imageImageRepresentation);
% data consistancy
dk=(measurement-kSpaceCoefficientCurrent).*downSamplingPattern;
kSpaceCoefficientNext=kSpaceCoefficientCurrent+1/parameters.L*dk;
% from k-space coefficient to sparse coefficient
imageImageRepresentation = image2DIFFT(kSpaceCoefficientNext);
imageBlockRepresentation = image2blocks(imageImageRepresentation,parameters.projectMatrix);
sparseCoefficientTemp=dictionary'*imageBlockRepresentation;
% soft thresholding
sparseCoefficientOut=SoftThresh(sparseCoefficientTemp,parameters.lambda/2/parameters.L);

end