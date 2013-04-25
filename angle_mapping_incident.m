% Script to generate incident angle mapping from flipped and undistorted
% E02 mirror images for a given acquisition set.

% Procedure
% 1. Load flipped and undistorted E02 mirror images.
% 2. Find centers of pixel regions using edge detection/region props.
% 3. Load shifted (to account for image size increase) center points.
% 4. Reflect the points about the center (theta angle about normal).
% 5. Use viewing angle mapping to search for corresponding incident angle.
% 6. Save incident angle positions corresponding to given acquisition set.

close all; clear all; clc;

fileSubDirectory = 'blu/E02 pixels_15x15_blu_2/';
% set noise threshold to decrease number of edges detected
% for 06/03/2012 E02 mirror images:
%   red -> 1000, grn -> 1400, blu -> 
noiseThreshold = 800;

geoCalibDirectory = [pwd '/calibration_data/geometric/'];
angleMappingDirectory = [pwd '/angle_mappings/'];

inputDirectory = [geoCalibDirectory fileSubDirectory];
wavelengthString = fileSubDirectory(1:3);

imageList = dir([inputDirectory 'img_undistort_*.tiff']);
nImages = size(imageList,1);

if ~exist([angleMappingDirectory 'incident_points_reflected_' wavelengthString '.mat'],'file')
    
    imageSum = 0;
    imageIncidentPointsReflected = nan(nImages/2,2);
    % iterate through but skip dark images
    for iImage = 2:2:nImages
        imageTemp = im2double(imread([inputDirectory imageList(iImage).name]),'indexed');
        imageTempOriginal = imageTemp;
        imageTemp(imageTemp < noiseThreshold) = 0;    	% drop some noise
        
        % get region edges
        imageEdges = edge(imageTemp,'canny');
        % find pinhole centroids
        imageRegionProps = regionprops(imageEdges);
        nPoints = size(imageRegionProps,1);             % ideally this is 1 value
        noiseThresholdTemp = noiseThreshold;
        while nPoints == 0
            noiseThreshold = noiseThreshold/2;          % reduce noise threshold
            imageTemp = im2double(imread([inputDirectory imageList(iImage).name]),'indexed');
            imageTemp(imageTemp < noiseThreshold) = 0;
            
            % get region edges
            imageEdges = edge(imageTemp,'canny');
            % find pinhole centroids
            imageRegionProps = regionprops(imageEdges);
            nPoints = size(imageRegionProps,1);
        end
        noiseThreshold = noiseThresholdTemp;            % reset value

        % sum images for display later
        imageSum = imageSum + imageTemp;
        
        % extract pinhole centroid points from structure
        ctrGP = 0;
        imats = nan(nPoints,2);
        imageAreaPoints = nan(nPoints,2);
        for iPoint = 1:nPoints
            ctrGP = ctrGP + 1;
            imageCentroidPoints(ctrGP,:) = imageRegionProps(ctrGP).Centroid;
            imageAreaPoints(ctrGP,:) = imageRegionProps(ctrGP).Area;
        end
        
        % display all found centroid points
        imageWithEdgeOverlay = imoverlay(mat2gray(imageTempOriginal), imageEdges, [0 1 0]);
        figure(100); imshow(imageWithEdgeOverlay); title(['reflectance image ' num2str(iImage) ' for ' wavelengthString]);
        xlabel(['there are ' num2str(nPoints) ' points found']); impixelinfo;
        hold on; plot(imageCentroidPoints(:,1),imageCentroidPoints(:,2),'r*','MarkerSize',12); hold off;
        pause(.05);
        % check found centroid points to find edge points
        if nPoints == 1
            imageIncidentPointsReflected(iImage/2,:) = imageCentroidPoints(1,:);
        elseif nPoints > 1
            for iPoint = 1:nPoints
                figure(100); imshow(imageWithEdgeOverlay); title(['reflectance image ' num2str(iImage)]); impixelinfo;
                hold on; plot(imageCentroidPoints(iPoint,1),imageCentroidPoints(iPoint,2),'r*','MarkerSize',12); hold off;
                yesOrNo = mvdlg('Is the point correct? [y/n]','Check point',[.8 .8 .2 .15]);
                if strncmp(yesOrNo,'y',1) == 1
                    % point is good so store it
                    imageIncidentPointsReflected(iImage/2,:) = imageCentroidPoints(iPoint,:);
                    continue;
                else
                    correctOrSkip = mvdlg('Correct OR skip? [c/s]','Incorrect point',[.8 .8 .2 .15]);
                    if strncmp(correctOrSkip,'c',1) == 1
                        figure(100); imshow(imageWithEdgeOverlay,[0 2^12-1]); title(['reflectance image ' num2str(iImage)]); impixelinfo;
                        hold on; plot(imageIncidentPointsReflected(iImage/2,1),imageIncidentPointsReflected(iImage/2,2),'g*','MarkerSize',12);
                        plot(imageCentroidPoints(iPoint,1),imageCentroidPoints(iPoint,2),'r*','MarkerSize',12); hold off; impixelinfo;
                        % get correct point (using myginput to show crosshair)
                        [newX,newY] = myginput(1,'crosshair');
                        imageIncidentPointsReflected(iImage/2,:) = [newX newY];
                        figure(100); imshow(imageWithEdgeOverlay); title(['reflectance image ' num2str(iImage)]); impixelinfo;
                        hold on; plot(imageIncidentPointsReflected(iImage/2,1),imageIncidentPointsReflected(iImage/2,2),'c*','MarkerSize',12); hold off;
                        waitforbuttonpress;
                    end
                end
            end
        end
        
    end
    % save incident point values
    save([angleMappingDirectory 'incident_points_reflected_' wavelengthString '.mat'],'imageIncidentPointsReflected','imageSum');

else load([angleMappingDirectory 'incident_points_reflected_' wavelengthString '.mat']);
end

% display all points
figure(201); imshow(imageSum,[1 2^12]); title(['total reflectance image for ' wavelengthString]); impixelinfo;
hold on; plot(imageIncidentPointsReflected(:,1),imageIncidentPointsReflected(:,2),'r*','MarkerSize',12); hold off;
waitforbuttonpress;

% load viewing angle mapping
load([angleMappingDirectory 'viewing_angle_mapping_' wavelengthString '_theta.mat']);
load([angleMappingDirectory 'viewing_angle_mapping_' wavelengthString '_phi.mat']);

% set center point
load([geoCalibDirectory 'center_points_shifted.mat']);
if strncmp(wavelengthString,'r',1)
    centerPoint = cc_red;
elseif strncmp(wavelengthString,'g',1)
    centerPoint = cc_grn;
elseif strncmp(wavelengthString,'b',1)
    centerPoint = cc_blu;
end

% reflect the angles about the normal
ctrMP = 0;
userInputFlag = 0;
for iImage = 2:2:nImages
    iImageTemp = iImage;
    
    imageTemp = im2double(imread([inputDirectory imageList(iImage).name]),'indexed');
    imageTemp = medfilt2(imageTemp);        % a bit of smoothing

    % skip points where no centroid was found
    if isnan(imageIncidentPointsReflected(iImage/2,:))
        continue
    end
    % iterate through found points and use with viewing angle mapping
    incidentPointReflectedX = imageIncidentPointsReflected(iImage/2,1);
    incidentPointReflectedY = imageIncidentPointsReflected(iImage/2,2);
    incidentPointReflectedXRounded = round(incidentPointReflectedX);
    incidentPointReflectedYRounded = round(incidentPointReflectedY);
    incidentPointReflectedTheta = viewingAngleMappingTheta(incidentPointReflectedYRounded,incidentPointReflectedXRounded);
    incidentPointReflectedPhi = viewingAngleMappingPhi(incidentPointReflectedYRounded,incidentPointReflectedXRounded);
    
    % show image and plot center point and mirror point
    figure(300); imshow(imageTemp,[1 2^12]); impixelinfo;
    hold on; plot(centerPoint(1),centerPoint(2),'g+','MarkerSize',12);
    plot(incidentPointReflectedX,incidentPointReflectedY,'r*','MarkerSize',12);
    title(['reflectance image ' num2str(iImage)]);
%     waitforbuttonpress;
    
    % reflect point about the normal
    incidentPointActualTheta = incidentPointReflectedTheta;
    incidentPointActualPhi = incidentPointReflectedPhi + 180;
    if incidentPointActualPhi >= 360
        incidentPointActualPhi = incidentPointActualPhi - 360;
    end
    % find (x,y) points (indices) for incident point
    [newIndexI,newIndexJ] = find(floor(viewingAngleMappingTheta) == floor(incidentPointActualTheta) & floor(viewingAngleMappingPhi) == floor(incidentPointActualPhi));
    % relax condition for match if floor() does not yield result
    if isempty(newIndexI) || isempty(newIndexJ)
        [newIndexI,newIndexJ] = find(round(viewingAngleMappingTheta) == round(incidentPointActualTheta) & round(viewingAngleMappingPhi) == round(incidentPointActualPhi));
    end
    if isempty(newIndexI) || isempty(newIndexJ)
        userInputFlag = 1;
        [newIndexJ,newIndexI] = myginput(1,'crosshair');
        newIndexI = round(newIndexI); newIndexJ = round(newIndexJ);
    end
    
    % use indices to get angle values
    nPotentialVals = length(newIndexI);
    thetaVals = nan(1,nPotentialVals);
    phiVals = nan(1,nPotentialVals);
    for iPotentialVal = 1:nPotentialVals;
        thetaVals(iPotentialVal) = viewingAngleMappingTheta(newIndexI(iPotentialVal),newIndexJ(iPotentialVal));
        phiVals(iPotentialVal) = viewingAngleMappingPhi(newIndexI(iPotentialVal),newIndexJ(iPotentialVal));
    end
    
    % divide by actual angle value for comparison
    thetaRatio = thetaVals/incidentPointActualTheta;
    phiRatio = phiVals/incidentPointActualPhi;
    
    % value of ratio vector closest to 1 is the best match
    [incidentPointActualThetaMatch,incidentPointActualThetaMatchPosition] = min(abs(thetaRatio(:)-1));
    [incidentPointActualPhiMatch,incidentPointActualPhiMatchPosition] = min(abs(phiRatio(:)-1));
   
    % assign new (x,y) point as reflection point
    if incidentPointActualThetaMatchPosition < incidentPointActualPhiMatchPosition
        imageIncidentPointActualX = newIndexI(incidentPointActualThetaMatchPosition);
        imageIncidentPointActualY = newIndexJ(incidentPointActualThetaMatchPosition);
    else
        imageIncidentPointActualX = newIndexI(incidentPointActualPhiMatchPosition);
        imageIncidentPointActualY = newIndexJ(incidentPointActualPhiMatchPosition);
    end

    % save actual incident angles and points
    ctrMP = ctrMP + 1;
    imageIncidentAnglesActual(ctrMP,1) = incidentPointActualTheta;
    imageIncidentAnglesActual(ctrMP,2) = incidentPointActualPhi;
    imageIncidentPointsActual(ctrMP,1) = imageIncidentPointActualX;
    imageIncidentPointsActual(ctrMP,2) = imageIncidentPointActualY;
    
    % plot incident point
    plot(imageIncidentPointActualY,imageIncidentPointActualX,'cx','MarkerSize',12); hold off;
    if userInputFlag == 1
        waitforbuttonpress;
        userInputFlag = 0;
    end
    
    pause(.01); 
end

save([angleMappingDirectory 'incident_points_actual_' wavelengthString '.mat'],'imageIncidentPointsActual','imageIncidentAnglesActual');

 