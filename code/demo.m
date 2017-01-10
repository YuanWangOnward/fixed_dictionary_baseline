% This demo shows how to use the Fixed-dictionary Baseline(FDB) MRI 
% compressive sensing reconstruction algorithm
% -------------------------------------------------------------------------
% Written by Yuan Wang, New York University, 2016, yw1225@nyu.edu
% -------------------------------------------------------------------------

%% IMPORTANT !!! 
% modify as needed
DIRECTORY_CONTAINING_THIS_FILE = 'set proper path here';
try
    cd(DIRECTORY_CONTAINING_THIS_FILE)
catch expection
    fprintf(2,[expection.message,'\n\n']);
    error('Have you change matlab working directory to the one containing this demo code, please?')
end

%% read in and normalize image
originalImagePath = fullfile('..','data','originalImage.mat');
try
    originalImage = load(originalImagePath);
catch expection
    fprintf(2,[expection.message,'\n\n']);
    error('Have you change matlab working directory to the one containing this demo code, please?')
end
originalImage = originalImage.originalImage;
maxValue=max(originalImage(:));
minValue=min(originalImage(:));
originalImage=(originalImage-minValue)/(maxValue-minValue)*255;
imagesc(originalImage)
title('original image')
colormap(gray)
shg

%% simulate down sampling
samplePatternPath = fullfile('..','data','radialSampling5.mat');
try
    samplePattern=load(samplePatternPath);
catch expection
    fprintf(2,[expection.message,'\n\n']);
    error('Have you change matlab working directory to the one containing this demo code, please?')
end
samplePattern = fftshift(samplePattern.samplePattern);
measure = makeAScan(originalImage,samplePattern);
imagesc(fftshift(log(abs(measure))))
title('k-space measurement (log)')
colormap(gray)
shg

%% reconstruction
% create dictionary
% the tight frame dictionary is too large to be represented as a matrix explicitly
dictionary = kron(haarmtx(parameters.blockSize(1))',haarmtx(parameters.blockSize(2))');    
% other settings
parameters = struct;
parameters.blockSize = [8, 8];
parameters.stepSize = [1, 1];
[parameters.projectMatrix, parameters.backProjectMatrix]=...
    createProjectionMatrix(size(measure),parameters.blockSize,parameters.stepSize);
parameters.stopTolerance = 0.000125;
parameters.maxIterationNumber = 100;
parameters.lambda = 0.2;
parameters.mu = 1;
parameters.L = 1;
parameters.ifPlot = 0;
% reconstruction
reconstructedImage = fixedDictionaryBaseline(measure,dictionary,parameters); 
reconstructedImage = abs(reconstructedImage);

%% evaluate and show results                     
[rMSE, PSNR]=evaluateReconstruction(abs(reconstructedImage),abs(originalImage));

figure(1)
subplot(1,2,1)
imagesc(originalImage)
title('original image')
axis square 
subplot(1,2,2)
imagesc(reconstructedImage)
title('reconstructed image')
colormap(gray)
axis square 

figure(2)
imagesc(abs(originalImage-reconstructedImage))
title(['error map, PSNR=',num2str(PSNR)])
colormap default
axis square 
shg









