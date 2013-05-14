% script to perform normal estimation
%
% first try will be to simply find the point of maximum intensity and draw
% a line from it to the initial incident angle (from the mappings
% calculated from mirror samples). This may work better for some samples,
% e.g. Pterygio rather than hatchetfish.
%
% load .mat files generated from 'process_samples' and iterate through
% while determining max intensity point
% also load incident angle mappings and plot corresponding initial incident
% angle

close all; clear all; clc;

whichColor = 'blu';
inputDirectory = [pwd '/TEST DATA/inputDirectory/S16/final_' whichColor '/'];
angleMappingDirectory = [pwd '/angle_mappings/'];
geoCalibDirectory = [pwd '/calibration_data/geometric/'];

imageInputList = dir([inputDirectory '*.mat']);
nImages = length(imageInputList);

% load center points
centerPoints = load([geoCalibDirectory 'center_points_shifted.mat']);
centerPoints = centerPoints.cc_red;

% load incident angle mappings
incidentAnglesActual = load([angleMappingDirectory 'incident_points_actual_' whichColor '.mat']);
incidentAnglesActualPoints = incidentAnglesActual.imageIncidentPointsActual;
incidentAnglesActualAngles = incidentAnglesActual.imageIncidentAnglesActual;
incidentAnglesReflected = load([angleMappingDirectory 'incident_points_reflected_' whichColor '.mat']);
incidentAnglesReflectedPoints = incidentAnglesReflected.imageIncidentPointsReflected;

% load viewing angle mappings
viewingAngleMappingTheta = load([angleMappingDirectory 'viewing_angle_mapping_' whichColor '_theta.mat']);
viewingAngleMappingTheta = viewingAngleMappingTheta.viewingAngleMappingTheta;
viewingAngleMappingPhi = load([angleMappingDirectory 'viewing_angle_mapping_' whichColor '_phi.mat']);
viewingAngleMappingPhi = viewingAngleMappingPhi.viewingAngleMappingPhi;

% get image parameters
imageTemp = load([inputDirectory imageInputList(1).name]);
imageTemp = imageTemp.imageFinal;
imageHeight = size(imageTemp,1);
imageWidth = size(imageTemp,2);
nPixels = imageHeight*imageWidth;

diffBetweenActualAndEstimated = nan(nImages,1);
angleDiffBetweenActualAndEstimated = nan(nImages,1);
for iImage = 1:nImages
% for iImage = 1
    imageName = imageInputList(iImage).name;
    imageTemp = load([inputDirectory imageName]);
    imageTemp = imageTemp.imageFinal;
    
    maxIntensityVal = max(max(imageTemp));
%     maxIntensityVec = (imageTemp == maxIntensityVal);
    [maxValX,maxValY] = find(imageTemp == maxIntensityVal);
    
%     figure(1000); imshow(imageTemp,[0 maxIntensityVal]); impixelinfo;
%     hold on; plot(maxValY,maxValX,'m*','MarkerSize',20);
%     plot(centerPoints(1),centerPoints(2),'r*','MarkerSize',20);
%     plot(incidentAnglesActualPoints(iImage,2),incidentAnglesActualPoints(iImage,1),'g*','MarkerSize',20);
%     plot(incidentAnglesReflectedPoints(iImage,1),incidentAnglesReflectedPoints(iImage,2),'b*','MarkerSize',20);
    
    % convert from spherical to cartesian coordinates
    calibReflectedTheta = incidentAnglesActualAngles(iImage,1);
    calibReflectedPhi = incidentAnglesActualAngles(iImage,2)+180;
    cartCalibReflectedPoint = [sind(calibReflectedTheta)*cosd(calibReflectedPhi) sind(calibReflectedTheta)*sind(calibReflectedPhi) cosd(calibReflectedTheta)];
    measuredReflectedTheta = viewingAngleMappingTheta(maxValX,maxValY);
    measuredReflectedPhi = viewingAngleMappingPhi(maxValX,maxValY);
    cartMeasuredReflectedPoint = [sind(measuredReflectedTheta)*cosd(measuredReflectedPhi) sind(measuredReflectedTheta)*sind(measuredReflectedPhi) cosd(measuredReflectedTheta)];
    
    diffBetweenActualAndEstimated(iImage) = pdist([cartCalibReflectedPoint;cartMeasuredReflectedPoint]);
    disp(diffBetweenActualAndEstimated(iImage));
    
    angleDiffBetweenActualAndEstimated(iImage) = acosd(dot(cartCalibReflectedPoint,cartMeasuredReflectedPoint)/(norm(cartCalibReflectedPoint)*norm(cartMeasuredReflectedPoint)));
    disp(angleDiffBetweenActualAndEstimated(iImage));
    
%     pause(0.5);
end

figure; plot(1:nImages,diffBetweenActualAndEstimated,'.');
mode(diffBetweenActualAndEstimated)

figure; plot(1:nImages,angleDiffBetweenActualAndEstimated,'.');
mode(angleDiffBetweenActualAndEstimated)


