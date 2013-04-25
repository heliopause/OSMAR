% Script to generate viewing angle mapping from flipped and undistorted
% images.

% Procedure
% 1. Get edge points from rectified calibration image
% 2. Calculate distance between center point and edge points
% 3. Calculate mean distance of above values
% 4. Use this distance from the center as point for maximum viewing angle
% 5. Generate viewing angle mapping for each wavelength

close all; clear all; clc;

geoCalibDirectory = [pwd '/calibration_data/geometric/'];
load([geoCalibDirectory 'center_points.mat']); 	% load center points from calibration
                                                % created using 'save_center_points.m'
proper_angle_calc;                              % previously written function
maxTrueDist = 12;                 % [mm]

mappingSaveDirectory = [pwd '/angle_mappings/'];

imageList = dir([geoCalibDirectory 'pinhole_occluder__rect01_*.tif']);
nImages = size(imageList,1);
imageLimit = 150;

nEdgePinholes = 16;
imageEdgePoints = nan(nEdgePinholes,2,nImages);
imageEdgePointDistances = nan(nEdgePinholes,nImages);
for iImage = 1:nImages
    imageName = imageList(iImage).name;
    wavelengthString = imageName(end-6:end-4);
    imageTemp = im2double(imread([geoCalibDirectory imageName]),'indexed');
    figure(1); imshow(imageTemp,[1 imageLimit]); impixelinfo;
    
    imageTemp = medfilt2(imageTemp);        % a bit of smoothing

    % wavelength specific changes
    if strncmp(wavelengthString,'r',1)
        imageTemp(imageTemp < 40) = 0;      % drop some noise
        imageTemp(173:245,351:435) = 0;   % get rid of reflection junk
        centerPoint = cc_red;          % set center point
        thetaMax = theta_red(3);       % set theta max
        axialDist = z(1);              % set axial (z) distance
    elseif strncmp(wavelengthString,'g',1)
        imageTemp(imageTemp < 40) = 0;
        imageTemp(178:250,335:425) = 0;
        centerPoint = cc_grn;
        thetaMax = theta_grn(3);
        axialDist = z(2);
    elseif strncmp(wavelengthString,'b',1)
        imageTemp(imageTemp < 10) = 0;
        imageTemp(180:250,339:419) = 0;
        centerPoint = cc_blu;
        thetaMax = theta_blu(3);
        axialDist = z(3);
    end
    
    % skip iteration if already done
    if ~exist([mappingSaveDirectory 'im_edge_points_' wavelengthString '.mat'],'file')
        % get region edges
        imageEdges = edge(imageTemp,'canny');
        % find pinhole centroids
        imageRegionProps = regionprops(imageEdges);
        nPointsFound = size(imageRegionProps,1);         % ideally this 1 + 8 + 12-1 + 16 = 36
        % extract pinhole centroid points from structure
        ctrGP = 0;
        imageCentroidPoints = nan(nPointsFound,2);
        for iPointFound = 1:nPointsFound
            ctrGP = ctrGP + 1;
            imageCentroidPoints(ctrGP,:) = imageRegionProps(ctrGP).Centroid;
        end
        % display all found centroid points
        figure(100); imshow(imageEdges); title(['rectified image for ' wavelengthString]);
        xlabel(['there are ' num2str(nPointsFound) ' points found']); impixelinfo;
        hold on; plot(imageCentroidPoints(:,1),imageCentroidPoints(:,2),'r*','MarkerSize',12); hold off;
        waitforbuttonpress;
        % check found centroid points to find edge points
        ctrDPC = 0;                         % number of correct pinholes
        for iPointFound = 1:nPointsFound
            figure(100); imshow(imageEdges); title(['rectified image for ' wavelengthString]); impixelinfo;
            hold on; plot(imageEdgePoints(:,1,iImage),imageEdgePoints(:,2,iImage),'g+','MarkerSize',20);
            plot(imageCentroidPoints(iPointFound,1),imageCentroidPoints(iPointFound,2),'r*','MarkerSize',12); hold off;
            yesOrNo = mvdlg('Is the point an edge point? [y/n]','Check point',[.8 .8 .2 .15]);
            if strncmp(yesOrNo,'y',1) == 1
                ctrDPC = ctrDPC + 1;
                % point is good so store it
                imageEdgePoints(ctrDPC,:,iImage) = imageCentroidPoints(iPointFound,:);
                % calculate distance between center point and each edge point
                positionPoint = [centerPoint(1),centerPoint(2);imageCentroidPoints(iPointFound,1),imageCentroidPoints(iPointFound,2)];
                calcDist = pdist(positionPoint,'euclidean');
                imageEdgePointDistances(ctrDPC,iImage) = calcDist;
                continue;
            else
                correctOrSkip = mvdlg('Skipping point','Incorrect point',[.8 .8 .2 .15]);
            end
        end
        % save verified edge points
        imageEdgePointsSingle = imageEdgePoints(:,:,iImage);
        imageEdgePointsSingleDistances = imageEdgePointDistances(:,iImage);
        save([mappingSaveDirectory 'im_edge_points_' wavelengthString '.mat'],'imageEdgePointsSingle','imageEdgePointsSingleDistances');
    else
        load([mappingSaveDirectory 'im_edge_points_' wavelengthString '.mat']);
    end
    
    % calculate mean distance value for current wavelength
    meanDist = mean(imageEdgePointsSingleDistances);
    
    % need to convert from [mm] to [pxl] using 'max_true_dist'
    mmToPixels = meanDist/maxTrueDist;
    
    % create matrix of pixel distances from 'center_point'
    [nImageRows,nImageCols] = size(imageTemp);
    
    % changes to make square image (680x680)
    centerPoint(2) = centerPoint(2) + (nImageCols - nImageRows)/2;
    nImageRows = nImageCols;

    distFromCenter = nan([nImageRows nImageCols]);
    for iImageRow = 1:nImageRows
        for iImageCol = 1:nImageCols
            distFromCenter(iImageRow,iImageCol) = hypot(iImageRow-centerPoint(2),iImageCol-centerPoint(1));
        end
    end
    % convert from pixel to metric distances
    distFromCenter = distFromCenter/mmToPixels;
    
    % calculate (theta) viewing angle mapping using 'axial_dist' [degrees]
    viewingAngleMappingTheta = atand(distFromCenter/axialDist);
    figure(101); imagesc(viewingAngleMappingTheta);
    impixelinfo; colorbar; axis equal;
    title(['viewing angle mapping (theta) for ' wavelengthString ' wavelength']);
    
    % save mapping
    save([mappingSaveDirectory 'viewing_angle_mapping_' wavelengthString '_theta.mat'],'viewingAngleMappingTheta');
    waitforbuttonpress;
    
    % calculate (phi) viewing angle mapping
    viewingAngleMappingPhi = nan([nImageRows nImageCols]);
    for iImageRow = 1:nImageRows
        for iImageCol = 1:nImageCols
            % calculate angle of zero vector and point vector and compare
            angleOne = atan2(centerPoint(1) - centerPoint(1),0 - centerPoint(2));
            angleTwo = atan2(centerPoint(2) - iImageRow,centerPoint(1) - iImageCol);
            viewingAngleMappingPhi(iImageRow,iImageCol) = (angleOne - angleTwo)*180/pi;
        end
    end
    figure(102); imagesc(viewingAngleMappingPhi);
    impixelinfo; colorbar; axis equal;
    title(['viewing angle mapping (phi) for ' wavelengthString ' wavelength']);

    % save mapping
    save([mappingSaveDirectory 'viewing_angle_mapping_' wavelengthString '_phi.mat'],'viewingAngleMappingPhi');
    waitforbuttonpress;
end
