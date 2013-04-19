% script for testing OSMAR functions

close all; clear all; clc;

% load config/parameter file OR specify input arguments
% test run_first with directory input

inputDirectory = '/Users/justin/Desktop/Hatchetfish test data/Spc16a-5 (red-air, mid-behind-head, 2x2, .4, 17ms, 2-10, g0)/';
outputDirectory = '/Users/justin/Desktop/test/';
% inputDirectory = 'TEST DATA/inputDirectory/';
% outputDirectory = 'TEST DATA/outputDirectory/';
imageExtension = 'tiff';

[imageHistograms,imageBitDepth] = run_first(inputDirectory,outputDirectory,imageExtension);
