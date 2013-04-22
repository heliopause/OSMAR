function [] = calibration_geometric_perform(calibrationImageDirectory,setColor)
% CALIBRATION_GEOMETRIC_PERFORM Perform camera calibration
% 
% This function uses the camera calibration toolbox found on
% http://www.vision.caltech.edu/bouguetj/calib_doc/ to get parameters
% needed to perform image rectification (geometric undistortion).
% 
% Procedure...
% 
% 1. Change to desired directory
% 2. Load calibration images
% 3. Load correspondence points
% 4. Perform calibration operation
%       NOTE: some iterations may be necessary here if initial calibration
%       is unsuccessful
% 5. Display extrinsic parameters
% 6. Save calibration parameters
% 7. Undistort calibration images
% 
% See also:
% CALIBRATION_GEOMETRIC_GET_POINTS, CALIBRATION_GEOMETRIC_SET_POINTS,
% CALIBRATION_GEOMETRIC_UNDISTORT_IMAGE

% calibrationImageDirectory = '/Users/justin/Documents/School/Scripps/Jaffe Lab/MURI project/BRDF project/programs/instrument_revision/OSMAR/calibration_data/geometric/';
% setColor = 'grn';
originalDirectory = pwd;
cd([calibrationImageDirectory setColor]);

% load calibration images
kc = zeros(5,1);
data_calib;

% load previously determined correspondence points
calibration_geometric_set_points;

% perform camera calibration and save results
go_calib_optim;
go_calib_optim;

% iterate again if initially unsuccessful
if ~exist('ex','var')
    center_optim = 0;
    go_calib_optim;
    add_suppress;
    center_optim = 1;
    go_calib_optim;
    go_calib_optim;
end

% display extrinsic parameters
ext_calib;

% save calibration results
saving_calib;

% undistort calibration images
undistort_image;

cd(originalDirectory);

end