% script for testing OSMAR functions

close all; clear all; clc;
addpath(pwd);

% % load config/parameter file OR specify input arguments
inputDirectory = [pwd '/TEST DATA/inputDirectory/S16a-8 (air, blu, mid-behind-head, 2x2, .4, 17ms, 2-10, g0)/'];
% % outputDirectory = '/Users/justin/Desktop/Hatchetfish test data/test/';
% % inputDirectory = 'TEST DATA/inputDirectory/';
outputDirectory = [pwd '/TEST DATA/outputDirectory/'];
imageExtension = 'tiff';

% test run_first with directory input
[imageHistograms,imageBitDepth] = run_first(inputDirectory,outputDirectory,imageExtension);

% % test other scripts
% dataDirectory = outputDirectory;

% % test calibration_geometric scripts
% imageInputDirectory = [pwd '/calibration_data/geometric/pinhole_occluder/'];
% imageOutputDirectory = [pwd '/calibration_data/geometric/'];
% setColor = 'blu';
% calibration_geometric_get_points(imageInputDirectory,imageOutputDirectory,setColor);

% % perform geometric calibration procedure
% calibrationImageDirectory = [pwd '/calibration_data/geometric/'];
% setColor = 'blu';
% calibration_geometric_perform(calibrationImageDirectory,setColor);

% test image undistort routine
imageInputDirectory = outputDirectory;
setColor = 'blu';
calibration_geometric_undistort_image(imageInputDirectory,setColor);

% ~~~~~~~~~ also need to get correspondence points for blue ~~~~~~~~~~

