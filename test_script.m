% script for testing OSMAR functions

close all; clear all; clc;

% load config/parameter file OR specify input arguments
% test run_first with directory input

inputDirectory = 'TEST DATA/inputDirectory/';
outputDirectory = 'TEST DATA/outputDirectory/';
imageExtension = 'tif';

[imageHistograms,imageBitDepth] = run_first(inputDirectory,outputDirectory,imageExtension);
