% as part of undistortion, zero pad image to square (680x680 or 
% largest dim X largest dim should be enough)
%
% possibly make square based on calculated actual center point
% this would be best, so that images for different wavelengths would be 
% centered on the same point, which would make it easier to combine them,
% since would only need to stretch/compress
function [imageUndistort] = calibration_geometric_undistort_image_single(imageTemp,setColor)

mainProgramDirectory = pwd;
if strncmp(setColor,'r',1)
    addpath([mainProgramDirectory '/calibration_data/geometric/red/']);
    Calib_Results;
    rmpath([mainProgramDirectory '/calibration_data/geometric/red/']);
elseif strncmp(setColor,'g',1)
    addpath([mainProgramDirectory '/calibration_data/geometric/grn/']);
    Calib_Results;
    rmpath([mainProgramDirectory '/calibration_data/geometric/grn/']);
elseif strncmp(setColor,'b',1)
    addpath([mainProgramDirectory '/calibration_data/geometric/blu/']);
    Calib_Results;
    rmpath([mainProgramDirectory '/calibration_data/geometric/blu/']);
end

[imageY,imageX] = size(imageTemp);
nZeroRowsHalf = (imageX - imageY)/2;

% make image square by zero padding
% also shift the center point
imageTemp = padarray(imageTemp,nZeroRowsHalf,0,'both');
if iImage == 1
    cc
    cc(2) = cc(2) + nZeroRowsHalf;
    cc
end

% moved into loop since cc was not redefined otherwise
KK = [fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2) ; 0 0 1];

% undistort the image
% rect() is part of camera calib toolbox
imageUndistort = rect(imageTemp,eye(3),fc,cc,kc,alpha_c,KK);

end
