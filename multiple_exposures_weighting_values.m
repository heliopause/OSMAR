function [weightingValue] = multiple_exposures_weighting_values(pixelValue,trueBitDepth)
% Function to get weighting values

% Step 2 - Determine weighting values for given pixel position
%          Load all (usually three) images corresponding to single incident
%          angle, store pixel values for each position, calculate weighting
%          value based on metric.

%   	  	the weighting should ignore saturated (or nearly saturated) values!
%   	  		so maybe the weighting function should emphasize the values nearest the
%   	  		  middle of the intensity range (e.g. around 2000 rather than higher or
%   	  		  lower
%   	  		this is what the Debevec and Malik weighting function does anyway
%   	  			w(z) = z-Zmin 	for z <= 1/2 * (Zmin+Zmax)
%   	  				   Zmax-z	for	z >  1/2 * (Zmin+Zmax)
%   	  		note w(z) is increasingly small for values near saturation

% 		2. 	for each image, iterate through every pixel position
% 		3.	determine the weighting value (wi) for a given pixel intensity based on its
% 			  value BEFORE adjusting by mPR

minPixelValue = 0;
maxPixelValue = 2^trueBitDepth - 1;

middleValue = (1/2) * (minPixelValue + maxPixelValue);

lowVec = pixelValue <= middleValue;
highVec = pixelValue > middleValue;
weightingValueLow = pixelValue - minPixelValue;
weightingValueHigh = maxPixelValue - pixelValue;

weightingValue = weightingValueLow.*lowVec + weightingValueHigh.*highVec;

end