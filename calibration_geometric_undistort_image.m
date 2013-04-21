% as part of undistortion, zero pad image to square (680x680 or 
% largest dim X largest dim should be enough)
%
% possibly make square based on calculated actual center point
% this would be best, so that images for different wavelengths would be 
% centered on the same point, which would make it easier to combine them,
% since would only need to stretch/compress
function [] = calibration_geometric_undistort_image(imageInputDirectory,setColor)

mainProgramDirectory = '/Users/justin/Documents/School/Scripps/Jaffe Lab/MURI project/BRDF project/programs/instrument_revision/OSMAR/';
if strncmp(setColor,'r',1)
    addpath([mainProgramDirectory 'calibration_data/red/']);
    Calib_Results;
    rmpath([mainProgramDirectory 'calibration_data/red/']);
elseif strncmp(setColor,'g',1)
    addpath([mainProgramDirectory 'calibration_data/grn/']);
    Calib_Results;
    rmpath([mainProgramDirectory 'calibration_data/grn/']);
elseif strncmp(setColor,'b',1)
    addpath([mainProgramDirectory 'calibration_data/blu/']);
    Calib_Results;
    rmpath([mainProgramDirectory 'calibration_data/blu/']);
end

imageList = dir([imageInputDirectory '*.tiff']);  	% Get list of images
nImages = numel(imageList);

imageOutputDirectory = [imageInputDirectory 'undistort/'];
if isdir(imageOutputDirectory)
    yesOrNo = input('Output directory exists. Overwrite (y/n)? ','s');
    if strcmp(yesOrNo,'y')
        delete([imageOutputDirectory '*.*']);
        rmdir(imageOutputDirectory);
    else
        disp('Aborting image undistortion process.');
        return;
    end
    mkdir(imageOutputDirectory);
else
    mkdir(imageOutputDirectory);
end

for iImage = 1:nImages
    imageName = imageList(iImage).name;
    imageTemp = double(imread([imageInputDirectory imageName]));
    
    [imageY,imageX] = size(imageTemp);
    nZeroRows = (imageX - imageY)/2;

    % make image square by zero padding
    imageTemp = padarray(imageTemp,nZeroRows,0,'both');
    if iImage == 1
        cc
        cc(2) = cc(2) + nZeroRows;          % shift center point after zero padding
        cc
    end
    
    % moved into loop since cc was not redefined otherwise
    KK = [fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2) ; 0 0 1];
    
    % undistort the image
    % rect() is part of camera calib toolbox
    imageUndistort = rect(imageTemp,eye(3),fc,cc,kc,alpha_c,KK);

    if iImage < 10
        imageNameNew = ['img_undistort_00' num2str(iImage) '.tiff'];
    elseif iImage >= 10 && iImage < 100
        imageNameNew = ['img_undistort_0' num2str(iImage) '.tiff'];
    else imageNameNew = ['img_undistort_' num2str(iImage) '.tiff'];
    end

    % save rectified image
    imwrite(uint16(round(imageUndistort)),[imageOutputDirectory imageNameNew],'tiff');
    disp(['Writing ' imageNameNew ' ...']);
end
