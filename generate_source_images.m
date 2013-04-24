function [] = generate_source_images(nRowsPixelRegion,nRegionsRowIncidentImage)
% GENERATE_SOURCE_IMAGES Generate binary pixel region images for use with
% OSMAR instrument data acquisition process.
% 
% generate_source_images takes as input two integer parameters
% 
% 'nRowsPixelRegion' indicates the number of rows that a square pixel region
% should have, so that the number of pixels is such a region is then
% nRowsPixelRegion^2
% 
% 'nRegionsRowIncidentImage' indicates the number of rows (of pixel regions)
% that should be present in the incident image summation, so that there
% will be nRegionsRowIncidentImage^2 and the total number of generated
% image files for an acquisition set is 2*(nRegionsRowIncidentImage^2)
% including dark images
% 
% This function will generate binary images containing pixel regions of
% specified size. Every odd frame is a dark image, which is used to
% distinctly separate the acquired images.
% 
% In addition, the routine will write an image 'img_pixels_full.png' that
% contains all the pixel regions present in the acquisition set.
% 
% Images are currently generated as 480x480 binary PNG files, to be used
% with MicroVision PicoP Projector.
% 
% See also:
% OTHER_RELEVANT_FUNCTIONS

% image dimension size
imageDimension = 480;

% index corresponding to pixel region size
regionSizeIndex = floor(nRowsPixelRegion/2);

% index corresponding to number of pixel regions
numRegionsIndex = (imageDimension-24*nRegionsRowIncidentImage+24)/(2*24);
regionPositionsY = 24*numRegionsIndex:24:imageDimension-(24*numRegionsIndex);
regionPositionsX = regionPositionsY;

nRegionsY = numel(regionPositionsY); nRegionsX = numel(regionPositionsX);

countImageWrite = 0;
imagePixelsFull = zeros(imageDimension,imageDimension);
for iRegionY = 1:nRegionsY
    for iRegionX = 1:nRegionsX
        
        % write dark frame
        imagePixels = zeros(imageDimension,imageDimension);
        
        countImageWrite = countImageWrite + 1;
        
        if countImageWrite < 10
            imageName = ['img_pixels_000' num2str(countImageWrite) '.png'];
        elseif countImageWrite >= 10 && countImageWrite < 100
            imageName = ['img_pixels_00' num2str(countImageWrite) '.png'];
        elseif countImageWrite >= 100 && countImageWrite < 1000
            imageName = ['img_pixels_0' num2str(countImageWrite) '.png'];
        else imageName = ['img_pixels_' num2str(countImageWrite) '.png'];
        end
        imwrite(imagePixels,[pwd '/generated_source_images/' imageName],'png');
        
        % write pixel frame
        imagePixels = zeros(imageDimension,imageDimension);
        
        imagePixels(regionPositionsY(iRegionX)-regionSizeIndex:regionPositionsY(iRegionX)+regionSizeIndex,...
            regionPositionsX(iRegionY)-regionSizeIndex:regionPositionsX(iRegionY)+regionSizeIndex) = 1;
        imagePixels = imagePixels == 1;     % convert to logical (binary)

        countImageWrite = countImageWrite + 1;
        
        if countImageWrite < 10
            imageName = ['img_pixels_000' num2str(countImageWrite) '.png'];
        elseif countImageWrite >= 10 && countImageWrite < 100
            imageName = ['img_pixels_00' num2str(countImageWrite) '.png'];
        elseif countImageWrite >= 100 && countImageWrite < 1000
            imageName = ['img_pixels_0' num2str(countImageWrite) '.png'];
        else imageName = ['img_pixels_' num2str(countImageWrite) '.png'];
        end
        imwrite(imagePixels,[pwd '/generated_source_images/' imageName],'png');
        
        imagePixelsFull = imagePixelsFull + imagePixels;
        
    end
end

imwrite(imagePixelsFull,[pwd '/generated_source_images/img_pixels_full.png'],'png');

end