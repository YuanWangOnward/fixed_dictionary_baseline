function measurement = makeAScan(originalImage,samplePattern)
% simulate measurement in the k-space with given smaple pattern

coefficient = image2DFFT(originalImage);
measurement=samplePattern.*coefficient;


