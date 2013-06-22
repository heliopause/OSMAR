% Script to calculate BRDF values from processed files

clear all; close all; clc;
addpath(pwd);
addpath('/Users/justin/Documents/MATLAB/circle_fit/');

% select set to process
whichSetToProcess = 'D10';
expectedBRDF = str2double(whichSetToProcess(2:3))/(100*pi);
% select color to process
whichColor = 'blu';
% where to save BRDF mappings
mainSaveDirBRDF = [pwd '/TEST DATA/outputDirectory/'];
saveDirBRDF = [mainSaveDirBRDF whichSetToProcess '/' whichColor '/'];
if ~exist(saveDirBRDF,'dir')
    mkdir(saveDirBRDF)
end
% specify reference set
whichSetReference = 'D40';

% load files to process and reference files
imageDirToProcess = [pwd '/TEST DATA/inputDirectory/' whichSetToProcess '/final_' whichColor '/'];
imageDirReference = [pwd '/TEST DATA/inputDirectory/' whichSetReference '/final_' whichColor '/'];
fileList = dir([imageDirReference '*.mat']);
nFiles = numel(fileList);

% load angle mappings
angleMappingDir = [pwd '/angle_mappings/'];
incidentMapping = load([angleMappingDir 'incident_points_actual_' whichColor '.mat']);
incidentMappingAngles = incidentMapping.imageIncidentAnglesActual;
incidentMappingPoints = incidentMapping.imageIncidentPointsActual;
viewingMappingTheta = load([angleMappingDir 'viewing_angle_mapping_' whichColor '_theta.mat']);
viewingMappingTheta = viewingMappingTheta.viewingAngleMappingTheta;
viewingMappingPhi = load([angleMappingDir 'viewing_angle_mapping_' whichColor '_phi.mat']);
viewingMappingPhi = viewingMappingPhi.viewingAngleMappingPhi;

% calculate circles for theta = 15, 30, 45
thetaVector = linspace(0,2*pi,360)';
[ct15_i,ct15_j] = find(round(viewingMappingTheta) == 15);
[xc15,yc15,Re15,a15] = circfit(ct15_i,ct15_j);
[ct30_i,ct30_j] = find(round(viewingMappingTheta) == 30);
[xc30,yc30,Re30,a30] = circfit(ct30_i,ct30_j);
[ct45_i,ct45_j] = find(round(viewingMappingTheta) == 45);
[xc45,yc45,Re45,a45] = circfit(ct45_i,ct45_j);
% get indices for values outside of 45 degrees
outerPointCircle = find(viewingMappingTheta > 45);

medianBRDFValue = nan(nFiles,1);
meanBRDFValue = nan(nFiles,1);
for iFile = 1:nFiles
    imageToProcess = load([imageDirToProcess fileList(iFile).name]);
    imageToProcess = imageToProcess.imageFinal;
    imageReference = load([imageDirReference fileList(iFile).name]);
    imageReference = imageReference.imageFinal;
    
    % calculate BRDF mapping using reference value
    referenceBRDFValue = str2double(whichSetReference(2:3))/(100*pi);
    mappingBRDF = (imageToProcess./imageReference)*referenceBRDFValue;
    mappingBRDF(outerPointCircle) = nan;
    
    % calculate basic statistics
    medianBRDFValue(iFile) = nanmedian(nanmedian(mappingBRDF));
    meanBRDFValue(iFile) = nanmean(nanmean(mappingBRDF));
    
    figure(1000); imshow(mappingBRDF,[0 expectedBRDF*2]); impixelinfo; hold on;
    
    % plot theta circles
    circleColor = [0 .3 .1];
    plot(yc15,xc15,'.','Color',circleColor);
    plot(yc15+Re15*sin(thetaVector),xc15+Re15*cos(thetaVector),'-','Color',circleColor,'LineWidth',2);
    plot(yc30+Re30*sin(thetaVector),xc30+Re30*cos(thetaVector),'-','Color',circleColor,'LineWidth',2);
    plot(yc15+Re45*sin(thetaVector),xc45+Re45*cos(thetaVector),'-','Color',circleColor,'LineWidth',5);

    % plot incident points
    plot(incidentMappingPoints(iFile,2),incidentMappingPoints(iFile,1),'b*','MarkerSize',20);

%     % plot max point
%     [max_i,max_j] = find(BRDF_data_temp == max_val(ii));
%     plot(max_j,max_i,'m*','MarkerSize',20);

%     % create scanline (have to manually select!!!!!!)
%     [scan_x,scan_y,scan_val,scan_xi,scan_yi] = improfile('bicubic');
%     % plot scanline
%     plot(scan_x,scan_y,'--b','LineWidth',2); % b for (a) and r for (b)
%     % create theta vector
%     theta_vec = improfile(viewing_angle_mapping_theta,scan_xi,scan_yi,length(scan_val));
%     % make first half of values negative ('incident side')
%     theta_vec_min = find(theta_vec == nanmin(theta_vec));
%     for kk = 1:theta_vec_min
%         theta_vec(kk) = -theta_vec(kk);
%     end
    
%     % re-plot incident and max points to be on top of line
%     plot(incident_pts_corrected(ii,2),incident_pts_corrected(ii,1),'*','MarkerSize',20,'Color',[.9 .9 .9]);
%     plot(max_j,max_i,'m*','MarkerSize',20);
    
    currentFileName = fileList(iFile).name;
    newFileName = ['BRDFmap' currentFileName(6:end-4)];
    save([saveDirBRDF newFileName '.mat'],'mappingBRDF');

    disp(['Processed BRDF mapping ' num2str(iFile) ' of ' num2str(nFiles) '.']);
    pause(0.5);
%     waitforbuttonpress;
end

expectedBRDF
meanMedianBRDF = nanmean(medianBRDFValue)
meanMeanBRDF = nanmean(meanBRDFValue)
varianceMedianBRDF = nanvar(medianBRDFValue)
varianceMeanBRDF = nanvar(meanBRDFValue)

figure(500); plot(1:nFiles,medianBRDFValue,'b-'); hold on;
plot(1:nFiles,meanBRDFValue,'r-'); hold on;
title('mean and median plots for BRDF values');
legend('median','mean','Location','NorthWest');

