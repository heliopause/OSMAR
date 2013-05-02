function [newPixelValue] = multiple_exposures_calculate_new(pixelValues,weightingValues,medianPowerRatios)
% function [newPixelValue] = multiple_exposures_calculate_new(pixelValues,weightingValues,medianPowerRatios)
% Function to merge multiple exposure images
%
% Step 1 - Calculate mean power ratios
%          Load two power measurement sets, determine where they are
%          correlated, shift based on correlation, get ratio, get mean and
%          variance of ratio.
% Step 2 - Determine weighting values for given pixel position
%          Load all (usually three) images corresponding to single incident
%          angle, store pixel values for each position, calculate weighting
%          value based on metric.
% Step 3 - Calculate new pixel intensities
%          Using mean power ratios and weight values, calculate the new
%          pixel intensities for each image position based on the
%          algorithm.

% 		1. 	determine mean power ratio (mPR) for a given image set as described above
% 		2. 	for each image, iterate through every pixel position
% 		3.	determine the weighting value (wi) for a given pixel intensity based on its
% 			  value BEFORE adjusting by mPR
% 		4.	divide pixel intensity (zi) by mPR and multiply by wi
% 		5.	sum up (pixel intensity * weighting value / mean power ratio) and divide by
% 			  sum of weighting values (w1 + w2 + w3 + ... + wi)
% 		EX:
% 			for a single pixel position...
% 			sum (z1 * (z1-Zmin OR zmax-z1) + z2 * (z2-Zmin OR zmax-z2) / mPR2
% 				+ z3 * (z3-Zmin OR zmax-z3) / mPR3) / sum (w1 + w2 + w3)
% 			where wi = (zi-Zmin OR zmax-zi) depending on zi <= or > .5*(Zmin+Zmax)
% 				  mPR1 = 1

% pixelValues and weightingValues are now (h,w,N)
% medianPowerRatios is still (1,N)

nDims = size(pixelValues,3);
newPixelValueNumer = zeros([size(pixelValues,1) size(pixelValues,2)]);
newPixelValueDenom = newPixelValueNumer;
for iDim = 1:nDims
    newPixelValueNumer = newPixelValueNumer + pixelValues(:,:,iDim).*weightingValues(:,:,iDim)/medianPowerRatios(iDim);
    newPixelValueDenom = newPixelValueDenom + weightingValues(:,:,iDim);
end

newPixelValue = newPixelValueNumer./newPixelValueDenom;

% newPixelValueNum = sum(pixelValues(1)*weightingValues(1) + pixelValues(2)*weightingValues(2)/medianPowerRatios(2) + pixelValues(3)*weightingValues(3)/medianPowerRatios(3));
% newPixelValueDen = sum([weightingValues(1) weightingValues(2) weightingValues(3)]);
% newPixelValue = newPixelValueNum/newPixelValueDen;

end