% function [] = calibration_geometric_set_points
% CALIBRATION_GEOMETRIC_SET_POINTS Set camera correspondence points for
% pinhole occluder with 18 rotations. These will be used by camera
% calibration toolbox routines.
%
% This function uses points determined by CALIBRATION_GEOMETRIC_GET_POINTS
% 
% See also:
% CALIBRATION_GEOMETRIC_GET_POINTS, CALIBRATION_GEOMETRIC_PERFORM
% CALIBRATION_GEOMETRIC_UNDISTORT_IMAGE

% load pinhole positions (actual and image)
load('pinhole_positions_image.mat');
load('pinhole_positions_actual.mat');

% from 'go_calib_optim.m'
%INPUT: x_1,x_2,x_3,...: Feature locations on the images     - image
%       X_1,X_2,X_3,...: Corresponding grid coordinates      - actual

% define actual pinhole position coordinates (replicate)
X_1 = pinholePositionsActual;
X_1 = [X_1'; ones(1,size(pinholePositionsActual,1))];

X_2 = X_1;
X_3 = X_1;
X_4 = X_1;
X_5 = X_1;
X_6 = X_1;
X_7 = X_1;
X_8 = X_1;
X_9 = X_1;
X_10 = X_1;
X_11 = X_1;
X_12 = X_1;
X_13 = X_1;
X_14 = X_1;
X_15 = X_1;
X_16 = X_1;
X_17 = X_1;
X_18 = X_1;

% define image pinhole position coordinates

x_1 = pinholePositionsImage(:,:,1)';
x_2 = pinholePositionsImage(:,:,2)';
x_3 = pinholePositionsImage(:,:,3)';
x_4 = pinholePositionsImage(:,:,4)';
x_5 = pinholePositionsImage(:,:,5)';
x_6 = pinholePositionsImage(:,:,6)';
x_7 = pinholePositionsImage(:,:,7)';
x_8 = pinholePositionsImage(:,:,8)';
x_9 = pinholePositionsImage(:,:,9)';
x_10 = pinholePositionsImage(:,:,10)';
x_11 = pinholePositionsImage(:,:,11)';
x_12 = pinholePositionsImage(:,:,12)';
x_13 = pinholePositionsImage(:,:,13)';
x_14 = pinholePositionsImage(:,:,14)';
x_15 = pinholePositionsImage(:,:,15)';
x_16 = pinholePositionsImage(:,:,16)';
x_17 = pinholePositionsImage(:,:,17)';
x_18 = pinholePositionsImage(:,:,18)';

% get center point (toolbox should be able to find this on its own)
% 
% center_points = pinholePositionsImage(1,:,:);
% center_points = reshape(center_points,[2 18]);
% figure; plot(center_points(1,:),center_points(2,:),'b*');
% [x_center,y_center,R_circle,a_circle] = circfit(center_points(1,:),center_points(2,:));
% 
% center_optim = 0;
% cc = [x_center ; y_center]';
