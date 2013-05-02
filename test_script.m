% script for testing OSMAR functions

close all; clear all; clc;
addpath(pwd);

% % load config/parameter file OR specify input arguments
% inputDirectory = [pwd '/TEST DATA/inputDirectory/S16a-8 (air, blu, mid-behind-head, 2x2, .4, 17ms, 2-10, g0)/'];
% % outputDirectory = '/Users/justin/Desktop/Hatchetfish test data/test/';
% % inputDirectory = 'TEST DATA/inputDirectory/';
% outputDirectory = [pwd '/TEST DATA/outputDirectory/'];

% inputDirectory = [pwd '/calibration_data/geometric/original_images/mirror_images/E02 pixels_15x15_blu_2/'];
% outputDirectory = [pwd '/calibration_data/geometric/blu/E02 pixels_15x15_blu_2/'];
% imageExtension = 'tiff';

% % test run_first with directory input
% [imageHistograms,imageBitDepth] = run_first(inputDirectory,outputDirectory,imageExtension);

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

% % test image undistort routine
% imageInputDirectory = outputDirectory;
% setColor = 'blu';
% calibration_geometric_undistort_image(imageInputDirectory,setColor);

% % test source image generation routine
% nRowsPixelRegion = 15;
% nRegionsRowIncidentImage = 13;
% generate_source_images(nRowsPixelRegion,nRegionsRowIncidentImage);

% test multiple exposure algorithm
powerDir = [pwd '/TEST DATA/inputDirectory/power/'];
powerFiles2 = {'S16a-1.txt','S16a-2.txt'};
powerFiles3 = {'S16a-1.txt','S16a-3.txt'};

[medianPowerRatio2] = multiple_exposures_mean_power_ratio(powerDir,powerFiles2);
[medianPowerRatio3] = multiple_exposures_mean_power_ratio(powerDir,powerFiles3);
medianPowerRatios = [1 medianPowerRatio2 medianPowerRatio3];

bitRate = 12;
pixelValues = 2000*medianPowerRatios;

[weightingValue1] = multiple_exposures_weighting_values(pixelValues(1),bitRate);
[weightingValue2] = multiple_exposures_weighting_values(pixelValues(2),bitRate);
[weightingValue3] = multiple_exposures_weighting_values(pixelValues(3),bitRate);
weightingValues = [weightingValue1 weightingValue2 weightingValue3];

[newPixelValue] = multiple_exposures_calculate_new(pixelValues,weightingValues,medianPowerRatios)

