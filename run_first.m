function [imageHistograms,imageBitDepth] = run_first(inputDirectory,outputDirectory,imageExtension)
% RUN_FIRST Run this function first to perform some necessary basic image processing
% 
% [imageHistograms] = run_first(inputDirectory,outputDirectory,imageExtension)
% takes three strings as input and returns status of image histograms where
% number of bins is equal to 2^imageBitDepth and bin locations are 0 to 
% (2^imageBitDepth)-1.
%
% The purpose of this function is to perform flip over horizontal axis of 
% the image to correct for issues in the TIFF library used to capture images
% (see SoftTriggerOSC.cpp). The function will also copy the input files to 
% a working directory to prevent overwriting of original data and minimizing
% use of the original data disk. The return values are the histogram counts 
% and the image bit depth.
%
% Examples:
% [imageHistograms,imageBitDepth] = run_first('inputDir','workingDir','tif')
% figure; plot(0:2^imageBitDepth-1,imageHistograms(:,1))
% 
% See also:
% OTHER_USEFUL_FUNCTIONS

assert(isdir(inputDirectory),'Not a valid input directory.');
assert(isdir(outputDirectory),'Not a valid output directory.');
assert(numel(imformats(imageExtension))==1,'Not a valid image file format.');

imageOutputList = dir([outputDirectory '*.' imageExtension]);
continueOrNot = 'y';
if ~isempty(imageOutputList)
    continueOrNot = input(['Output directory contains ' upper(imageExtension)...
        ' files. Continue with processing [y/n]?\n'],'s');
end

if ~strcmp(continueOrNot,'y')
    disp('Cancelling the process.');
    imageHistograms = nan;
    imageBitDepth = nan;
else
    disp('Processing files...');
    imageInputList = dir([inputDirectory '*.' imageExtension]);
    imageInfo = imfinfo([inputDirectory imageInputList(1).name],imageExtension);
    imageBitDepth = imageInfo.BitDepth;
    
    nImages = numel(imageInputList);
    imageHistograms = nan(2^imageBitDepth,nImages);
    for iImage = 1:nImages
        imageName = imageInputList(iImage).name;
        imageTemp = imread([inputDirectory imageName]);
        
        imageTemp = flipud(imageTemp);
        
        imageHistogram = imhist(imageTemp,2^imageBitDepth);
        imageHistograms(:,iImage) = imageHistogram;
        
        imwrite(imageTemp,[outputDirectory imageName],imageExtension);
    end
    disp('Completed the process.');
end

end