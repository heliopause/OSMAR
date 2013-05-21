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

% load center points
centerPoints = load([geoCalibDirectory 'center_points_shifted.mat']);
centerPoints = centerPoints.cc_blu;                                                 % !!!!!!!!!!!!!!!!!!!!

imageInputList = dir([inputDirectory '*.mat']);
nImages = length(imageInputList);

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

diffBetweenActualAndEstimated = nan(nImages,1);
angleDiffBetweenActualAndEstimated = nan(nImages,1);
normalVectorSpherical = nan(nImages,2);
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
    
    % convert from spherical to cartesian coordinates (radius = 1)
    calibReflectedTheta = incidentAnglesActualAngles(iImage,1);
    calibReflectedPhi = incidentAnglesActualAngles(iImage,2)+180;
    cartCalibReflectedPoint = [sind(calibReflectedTheta)*cosd(calibReflectedPhi) sind(calibReflectedTheta)*sind(calibReflectedPhi) cosd(calibReflectedTheta)];
    measuredReflectedTheta = viewingAngleMappingTheta(maxValX,maxValY);
    measuredReflectedPhi = viewingAngleMappingPhi(maxValX,maxValY);
    cartMeasuredReflectedPoint = [sind(measuredReflectedTheta)*cosd(measuredReflectedPhi) sind(measuredReflectedTheta)*sind(measuredReflectedPhi) cosd(measuredReflectedTheta)];
    
    % get bisecting vector between incident and measured reflected rays (surface normal)
    calibIncidentTheta = incidentAnglesActualAngles(iImage,1);
    calibIncidentPhi = incidentAnglesActualAngles(iImage,2);
    cartCalibIncidentPoint = [sind(calibIncidentTheta)*cosd(calibIncidentPhi) sind(calibIncidentTheta)*sind(calibIncidentPhi) cosd(calibIncidentTheta)];
    normalVectorCart = 0.5*(cartCalibIncidentPoint + cartMeasuredReflectedPoint);
    % --begin test vectors
%     testIncidentTheta = incidentAnglesActualAngles(iImage,1);
%     testIncidentPhi = incidentAnglesActualAngles(iImage,2);
%     testReflectedTheta = testIncidentTheta;
%     testReflectedPhi = testIncidentPhi + 180;
%     cartTestIncidentPoint = [sind(testIncidentTheta)*cosd(testIncidentPhi) sind(testIncidentTheta)*sind(testIncidentPhi) cosd(testIncidentTheta)];
%     cartTestReflectedPoint = [sind(testReflectedTheta)*cosd(testReflectedPhi) sind(testReflectedTheta)*sind(testReflectedPhi) cosd(testReflectedTheta)];
%     normalVector = 0.5*(cartTestIncidentPoint + cartTestReflectedPoint);
    % --end test vectors
    normalVectorX = normalVectorCart(1);
    normalVectorY = normalVectorCart(2);
    normalVectorZ = normalVectorCart(3);
    normalVectorRho = sqrt(normalVectorX^2 + normalVectorY^2 + normalVectorZ^2);
    normalVectorS = sqrt(normalVectorX^2 + normalVectorY^2);
    normalVectorTheta = acosd(normalVectorZ/normalVectorRho);
    normalVectorPhi = asind(normalVectorY/normalVectorS);
    if normalVectorX < 0
        normalVectorPhi = 180 - normalVectorPhi;
    end
    if normalVectorPhi < 0
        normalVectorPhi = normalVectorPhi + 360;
    end
    normalVectorSpherical(iImage,:) = [normalVectorTheta normalVectorPhi];
%     disp(normalVectorSpherical(iImage,:));

%     diffBetweenActualAndEstimated(iImage) = pdist([cartCalibReflectedPoint;cartMeasuredReflectedPoint]);
%     disp(diffBetweenActualAndEstimated(iImage));
    
    angleDiffBetweenActualAndEstimated(iImage) = acosd(dot(cartCalibReflectedPoint,cartMeasuredReflectedPoint)/(norm(cartCalibReflectedPoint)*norm(cartMeasuredReflectedPoint)));
%     disp(angleDiffBetweenActualAndEstimated(iImage));
    
%     pause(0.5);
end

% figure; plot(1:nImages,diffBetweenActualAndEstimated,'.');
% mode(diffBetweenActualAndEstimated)

figure; plot(1:nImages,angleDiffBetweenActualAndEstimated,'.');
modeAngleDiffBetweenActualAndEstimated = mode(angleDiffBetweenActualAndEstimated)

% mode should be most important, since I am assuming that the estimated
% reflected point is not necessarily accurate much of the time
% form histogram of values from which mode will be drawn
figure; plot(1:nImages,normalVectorSpherical(:,1),'.');
binWidthTheta = 1;
edgesTheta = 0:binWidthTheta:45;
[thetaN,thetaBin] = histc(normalVectorSpherical(:,1),edgesTheta);
modeTheta = mode(thetaBin);
modeThetaValue = edgesTheta([modeTheta,modeTheta+1])
figure; hist(normalVectorSpherical(:,1),edgesTheta+binWidthTheta/2);

modeNormalVectorTheta = mode(normalVectorSpherical(:,1))

figure; plot(1:nImages,normalVectorSpherical(:,2),'.');
binWidthPhi = 1;
edgesPhi = 0:binWidthPhi:360;
[phiN,phiBin] = histc(normalVectorSpherical(:,2),edgesPhi);
modePhi = mode(phiBin);
modePhiValue = edgesPhi([modePhi,modePhi+1])
figure; hist(normalVectorSpherical(:,2),edgesPhi+binWidthPhi/2);

modeNormalVectorPhi = mode(normalVectorSpherical(:,2))

