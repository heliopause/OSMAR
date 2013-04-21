% script for testing OSMAR functions

close all; clear all; clc;
addpath(pwd);

% % load config/parameter file OR specify input arguments
% 
% inputDirectory = '/Users/justin/Desktop/Hatchetfish test data/Spc16a-5 (red-air, mid-behind-head, 2x2, .4, 17ms, 2-10, g0)/';
% outputDirectory = '/Users/justin/Desktop/Hatchetfish test data/test/';
% % inputDirectory = 'TEST DATA/inputDirectory/';
% % outputDirectory = 'TEST DATA/outputDirectory/';
% imageExtension = 'tiff';
% 
% % test run_first with directory input
% [imageHistograms,imageBitDepth] = run_first(inputDirectory,outputDirectory,imageExtension);
% 
% % test other scripts
% dataDirectory = outputDirectory;

% % test calibration_geometric scripts
% imageInputDirectory = '/Volumes/Calibration Data (v2)/07-13-12/pinhole_occluder/';
% imageOutputDirectory = '/Users/justin/Documents/School/Scripps/Jaffe Lab/MURI project/BRDF project/programs/instrument_revision/OSMAR/calibration_data/';
% setColor = 'grn';
% calibration_geometric(imageInputDirectory,imageOutputDirectory,setColor);

% % perform geometric calibration procedure
% calibrationImageDirectory = [pwd '/calibration_data/'];
% setColor = 'grn';
% calibration_geometric_perform(calibrationImageDirectory,setColor);

% test image undistort routine
imageInputDirectory = '/Users/justin/Documents/School/Scripps/Jaffe Lab/MURI project/BRDF project/programs/instrument_revision/OSMAR/TEST DATA/inputDirectory/';
setColor = 'grn';
calibration_geometric_undistort_image(imageInputDirectory,setColor);
