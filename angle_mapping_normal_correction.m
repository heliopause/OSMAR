% script to use results of normal estimation to modify angle mappings
% 
% this will be done once for each data set
% 
% Procedure:
% == viewing angles ==
% 1. run 'normal_estimation.m' to get (theta,phi)_tilt
% 2. divide theta mapping into two parts by drawing a line perpendicular to
%    the phi_tilt value
% 3. make theta values on tilt side (based on theta_tilt) positive
% 4. make theta values on other side negative
% 5. subtract theta_tilt from all values
% 6. convert all theta values to positive
% 7. symmetry holds in azimuthal (phi) direction so could adjust phi
%    mapping so that phi = 0 -> phi = phi_tilt by adding phi_tilt to all
%    phi values
% == incident angles ==
% 1. choose whether to add or subtract theta_tilt based on which side
%    theta_incident falls
% 2. search new viewing angle mapping for new theta_incident value to get
%    updated index value

close all; clear all; clc;

% specify image set location and wavelength in 'normal_estimation'
% store theta and tilt values
normal_estimation;
thetaTilt = mean(modeThetaValueFromHistogram);
phiTilt = mean(modePhiValueFromHistogram);
% test values
% thetaTilt = 30;

% get tilt axis vector
maxThetaValue = max(max(viewingAngleMappingTheta));
maxPhiValue = max(max(viewingAngleMappingPhi));
% tiltAxisStart = phiTilt + 90;
% if tiltAxisStart > maxPhiValue;
%     tiltAxisStart = tiltAxisStart - maxPhiValue;
% end
tiltAxisEnd = phiTilt - 90;
if tiltAxisEnd < 0
    tiltAxisEnd = tiltAxisEnd + 360;
end

% rename for clarity
thetaOriginal = viewingAngleMappingTheta;
phiOriginal = viewingAngleMappingPhi;

[X,Y,Z] = sph2cart(phiOriginal*pi/180,(90-thetaOriginal)*pi/180,1);

[xRotAxis,yRotAxis,zRotAxis] = sph2cart(tiltAxisEnd*pi/180,0,1);
rotateArbAxis = makehgtform('axisrotate',[xRotAxis yRotAxis zRotAxis],thetaTilt*pi/180);
rotateArbAxis = rotateArbAxis(1:3,1:3);

% need to get rows of 3d cartesian points, then can multiply rotateArbAxis*XYZcoords
newPoints = (rotateArbAxis*[X(:) Y(:) Z(:)]');
% newPoints = (rotateZ*[X(:) Y(:) Z(:)]')*180/pi;
Xp = newPoints(1,:); Yp = newPoints(2,:); Zp = newPoints(3,:);

% convert back to spherical (should have all R = 1)
[phiCorrectedVec,thetaCorrectedVec,RsbOne] = cart2sph(Xp,Yp,Zp);

% reshape angle mappings
viewingAngleMappingThetaCorrected = 90-(reshape(thetaCorrectedVec,680,680)*180/pi);
viewingAngleMappingPhiCorrected = (reshape(phiCorrectedVec,680,680)*180/pi);
% define new phi = 0 as phiTilt
negativeIdx = viewingAngleMappingPhiCorrected < 0;
viewingAngleMappingPhiCorrected(negativeIdx) = viewingAngleMappingPhiCorrected(negativeIdx)+360;

% plot old and new angle mappings for comparison
figure; subplot(1,2,1); imagesc(viewingAngleMappingTheta); colorbar; axis square; impixelinfo;
title('original theta mapping');
subplot(1,2,2); imagesc(viewingAngleMappingThetaCorrected); colorbar; axis square; impixelinfo;
title('corrected theta mapping');
figure; subplot(1,2,1); imagesc(viewingAngleMappingPhi); colorbar; axis square; impixelinfo;
title('original phi mapping');
subplot(1,2,2); imagesc(viewingAngleMappingPhiCorrected); colorbar; axis square; impixelinfo;
title('corrected phi mapping');

%% incident angles

% get incident angle coordinates and multiply with rotation matrix
incidentAnglesActualTheta = incidentAnglesActualAngles(:,1);
incidentAnglesActualPhi = incidentAnglesActualAngles(:,2);
[incActX,incActY,incActZ] = sph2cart(incidentAnglesActualPhi*pi/180,(90-incidentAnglesActualTheta)*pi/180,1);
incidentAnglesActualAnglesNew = (rotateArbAxis*[incActX(:) incActY(:) incActZ(:)]');
incActXp = incidentAnglesActualAnglesNew(1,:);
incActYp = incidentAnglesActualAnglesNew(2,:);
incActZp = incidentAnglesActualAnglesNew(3,:);
[incPhiCorrectedVec,incThetaCorrectedVec,incRsbOne] = cart2sph(incActXp,incActYp,incActZp);

incidentAnglesActualThetaCorrected = 90-(incThetaCorrectedVec*180/pi);
incidentAnglesActualPhiCorrected = incPhiCorrectedVec*180/pi;
incNegativeIdx = incidentAnglesActualPhiCorrected < 0;
incidentAnglesActualPhiCorrected(incNegativeIdx) = incidentAnglesActualPhiCorrected(incNegativeIdx)+360;

% store corrected incident angles
incidentAnglesActualAnglesCorrected = [incidentAnglesActualThetaCorrected' incidentAnglesActualPhiCorrected'];

% % index into viewing angle mapping to get new indices
% % find (x,y) points (indices) for incident point
% nPoints = size(incidentAnglesActualAnglesCorrected,1);
% correctedIndexI = nan(nPoints,1);
% correctedIndexJ = nan(nPoints,1);
% for iPoint = 1:nPoints
%     [newIndexI,newIndexJ] = find(floor(viewingAngleMappingThetaCorrected) == floor(incidentAnglesActualThetaCorrected(iPoint))...
%         & floor(viewingAngleMappingPhi) == floor(incidentAnglesActualPhiCorrected(iPoint)));
%     % relax condition for match if floor() does not yield result
%     if isempty(newIndexI) || isempty(newIndexJ)
%         [newIndexI,newIndexJ] = find(round(viewingAngleMappingThetaCorrected) == round(incidentAnglesActualThetaCorrected(iPoint))...
%             & round(viewingAngleMappingPhi) == round(incidentAnglesActualPhiCorrected(iPoint)));
%     end
% %     if isempty(newIndexI) || isempty(newIndexJ)
% %         userInputFlag = 1;
% %         [newIndexJ,newIndexI] = myginput(1,'crosshair');
% %         newIndexI = round(newIndexI); newIndexJ = round(newIndexJ);
% %     end
%     % take the centermost point as the corrected incident point indices
%     correctedIndexI(iPoint) = median(newIndexI);
%     correctedIndexJ(iPoint) = median(newIndexJ);
% end


