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
%       NOTE: sometimes need to set center_optim=0, go_calib_optim, and
%       then center_optim=1, go_calib_optim again
% 5. Display extrinsic parameters
% 6. Save calibration parameters
% 7. Undistort calibration images
% 
% See also:
% CALIBRATION_GEOMETRIC_GET_POINTS, CALIBRATION_GEOMETRIC_SET_POINTS,
% CALIBRATION_GEOMETRIC_UNDISTORT_IMAGE

% calibrationImageDirectory = '/Users/justin/Documents/School/Scripps/Jaffe Lab/MURI project/BRDF project/programs/instrument_revision/OSMAR/calibration_data/geometric/';
% setColor = 'grn';
cd([calibrationImageDirectory setColor]);

% load calibration images
kc = zeros(5,1);
data_calib;

% load previously determined correspondence points
calibration_geometric_set_points;

% perform camera calibration and save results
go_calib_optim;

% should write something in here to take care of re-iterations until there
% is no error (basically, can check for existence of variable 'ex') and if
% it is not there, center_optim=0, reiterate
ext_calib;
saving_calib;

% undistort calibration images
undistort_image;

end