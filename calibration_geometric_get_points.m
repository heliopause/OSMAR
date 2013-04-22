function [] = calibration_geometric_get_points(imageInputDirectory,imageOutputDirectory,setColor)
% CALIBRATION_GEOMETRIC_GET_POINTS Get correspondence points for geometric calibration
% (undistortion) based on input images with mirror as sample
% 
% [outputOne,outputTwo] = calibration_geometric(imageInputDirectory,imageOutputDirectory,setColor)
% takes input and output directories, and specified color as input strings
%
% The purpose of this function is to:
% 
%   Step 1 - For a given wavelength, read in each folder of position images,
%            take the mean image from each folder, flip it UD, and write to
%            specified directory.
%   Step 2 - Get pixel positions for center of each pinhole in each image.
%            Some manual input is necessary here due to reflections.
%   Step 3 - Input center positions for actual pinholes and store in a
%            matrix. This is based on physical pinhole dimensions.
%   Step 4 - Sort estimated pinhole positions so that the points are
%            ordered the same as the physical pinhole centers. The data 
%            from Steps 3 and 4 will be used to perform geometric calibration
%            using the camera calibration toolbox.
%
% Examples:
% calibration_geometric_set_points('inputDir','workingDir','clr')
% 
% See also:
% CALIBRATION_GEOMETRIC_SET_POINTS, CALIBRATION_GEOMETRIC_PERFORM
% CALIBRATION_GEOMETRIC_UNDISTORT_IMAGE

%% Step 0 - Manipulate inputs
colorSubDirectory = [imageInputDirectory 'occluder_' setColor '/'];
positionSubDirectoryBegin = 'position_';
nRotations = numel(dir([colorSubDirectory positionSubDirectoryBegin '*']));

%% Step 1 - Obtain rotation images
imageTemp = imread([colorSubDirectory positionSubDirectoryBegin '01/snap0001.tiff']);
imagesPinholeOccluder = zeros([size(imageTemp) nRotations]);
for iRotation = 1:nRotations
    if iRotation < 10
        positionSubDirectory = [positionSubDirectoryBegin '0' num2str(iRotation) '/'];
    else positionSubDirectory = [positionSubDirectoryBegin num2str(iRotation) '/'];
    end
    positionDirectory = [colorSubDirectory positionSubDirectory];
    
    % process each position folder (20 images in each - 10 dark)
    imageInputList = dir([positionDirectory 'snap*.tiff']);
    nImages = numel(imageInputList);
    imageTemp = zeros(size(imageTemp));
    for iImage = 1:nImages
        % skip dark images
        if mod(iImage,2) == 0
            imageTemp = imageTemp + im2double(imread([positionDirectory imageInputList(iImage).name]),'indexed');
        else continue
        end
    end
    imageTemp = imageTemp/(nImages/2);                  % mean image
    imageTemp(imageTemp < 180) = 0;                     % drop some noise
    imageTemp = flipud(imageTemp);                      % correct for TIFF library issue
    imagesPinholeOccluder(:,:,iRotation) = imageTemp;
    disp(['Processed ' num2str(iRotation) ' of ' num2str(nRotations) ' position images.']);
end

% Write averaged and flipped pinhole occluder images to disk
imageWriteDirectory = [imageOutputDirectory setColor '/'];
if ~isdir(imageWriteDirectory)
    mkdir(imageWriteDirectory);
end
for iRotation = 1:nRotations
    if iRotation < 10
        imageName = ['pinhole_occluder_0' num2str(iRotation) '.tif'];
    else imageName = ['pinhole_occluder_' num2str(iRotation) '.tif'];
    end
    imageTemp = imagesPinholeOccluder(:,:,iRotation);
    imageTemp = imageTemp/(2^12-1);
    imwrite(im2uint8(imageTemp),[imageWriteDirectory imageName],'TIFF');
end

%% Step 2 - Obtain pinhole centroid positions
% This would work almost perfectly if not for the reflections. Maybe redo
% with external light source (other projector)?

rotationCentroidFiles = dir([imageWriteDirectory 'im_centroid_pts_position_*.mat']);
rotationRestart = size(rotationCentroidFiles,1) + 1;

nPinholes = 36;                  % actual number of pinholes (points)
imageCentroidPointsCorrected = nan(nPinholes,2,nRotations);
figure(100);
for iRotation = rotationRestart:nRotations
    % detect edges of pinhole images (converts to BW image)
    imageTemp = medfilt2(imagesPinholeOccluder(:,:,iRotation));
    imageTempOriginal = imageTemp;
    % get rid of reflection junk that confuses edge detection
    if strncmp(setColor,'r',1)          % (vert, horiz)
        imageTemp(173:245,351:435) = 0;
    elseif strncmp(setColor,'g',1)
        imageTemp(178:250,335:425) = 0;
    elseif strncmp(setColor,'b',1)
        imageTemp(180:250,339:419) = 0;
    end
    % get region edges
    imageEdges = edge(imageTemp,'canny');
    
    % find pinhole centroids
    imageRegionProps = regionprops(imageEdges);
    nPointsFound = size(imageRegionProps,1);         % ideally this 1 + 8 + 12-1 + 16 = 36
    ctrCPF = 0;                                     % num pinhole centroids found
    imageCentroidPoints = nan(nPointsFound,2);
    for iPointFound = 1:nPointsFound
        ctrCPF = ctrCPF + 1;
        imageCentroidPoints(ctrCPF,:) = imageRegionProps(ctrCPF).Centroid;
    end
    
    % display all found centroid points for single rotation
    imageWithEdgeOverlay = imoverlay(mat2gray(imageTempOriginal), imageEdges, [1 1 1]);
    figure(100); imshow(imageWithEdgeOverlay); title(['rotation image ' num2str(iRotation)]);
    xlabel(['there are ' num2str(nPointsFound) ' points found']); impixelinfo;
    hold on; plot(imageCentroidPoints(:,1),imageCentroidPoints(:,2),'r*','MarkerSize',12); hold off;
    waitforbuttonpress;
    % check found centroid points to see if correct -- if not, fix them
    ctrPC = 0;                         % number of correct pinholes
    for iPointFound = 1:nPointsFound
        figure(100); imshow(imageWithEdgeOverlay); title(['rotation image ' num2str(iRotation)]); impixelinfo;
        hold on; plot(imageCentroidPointsCorrected(:,1,iRotation),imageCentroidPointsCorrected(:,2,iRotation),'g+','MarkerSize',20);
        plot(imageCentroidPoints(iPointFound,1),imageCentroidPoints(iPointFound,2),'r*','MarkerSize',12); hold off;
        yesOrNo = mvdlg('Is the point correct? [y/n]','Check point',[.8 .8 .2 .15]);
        if strncmp(yesOrNo,'y',1) == 1
            ctrPC = ctrPC + 1;
            % point is good so store it
            imageCentroidPointsCorrected(ctrPC,:,iRotation) = imageCentroidPoints(iPointFound,:);
            continue;
        else
            correctOrSkip = mvdlg('Correct OR skip? [c/s]','Incorrect point',[.8 .8 .2 .15]);
            if strncmp(correctOrSkip,'c',1) == 1
                ctrPC = ctrPC + 1;
                figure(100); imshow(imagesPinholeOccluder(:,:,iRotation),[0 2^12-1]); title(['rotation image ' num2str(iRotation)]); impixelinfo;
                hold on; plot(imageCentroidPointsCorrected(:,1,iRotation),imageCentroidPointsCorrected(:,2,iRotation),'g*','MarkerSize',12);
                plot(imageCentroidPoints(iPointFound,1),imageCentroidPoints(iPointFound,2),'r*','MarkerSize',12); hold off;
                % get correct point (using myginput to show crosshair)
                [newX,newY] = myginput(1,'crosshair');
                imageCentroidPointsCorrected(ctrPC,:,iRotation) = [newX newY];
                figure(100); imshow(imageEdges); title(['rotation image ' num2str(iRotation)]); impixelinfo;
                hold on; plot(imageCentroidPointsCorrected(ctrPC,1,iRotation),imageCentroidPointsCorrected(ctrPC,2,iRotation),'c*','MarkerSize',12); hold off;
                waitforbuttonpress;
            end
        end
    end
    numCentroidPointsCorrected = sum(~isnan(imageCentroidPointsCorrected(:,1,iRotation)));
    % display all checked and corrected points
    figure(100); imshow(imageEdges); impixelinfo; hold on;
    plot(imageCentroidPointsCorrected(:,1,iRotation),imageCentroidPointsCorrected(:,2,iRotation),'r*','MarkerSize',12); hold off;
    title(['found ' num2str(numCentroidPointsCorrected) ' corrected points']);
    waitforbuttonpress;

    % get pinhole ID numbers
    pinholeNumberOrdered = nan(nPinholes,nRotations);
    figure(100); imshow(imageEdges); title(['rotation image ' num2str(iRotation)]); impixelinfo;
    for iPinhole = 1:nPinholes
        hold on; plot(imageCentroidPointsCorrected(iPinhole,1,iRotation),imageCentroidPointsCorrected(iPinhole,2,iRotation),'r*','MarkerSize',12);
        getPinholeNumber = mvdlg('Which point is this?','Point numbering',[.8 .8 .2 .15]);
        pinholeNumberOrdered(iPinhole,iRotation) = str2double(getPinholeNumber{1});
        plot(imageCentroidPointsCorrected(iPinhole,1,iRotation),imageCentroidPointsCorrected(iPinhole,2,iRotation),'g*','MarkerSize',12);
        text(imageCentroidPointsCorrected(iPinhole,1,iRotation)-10,imageCentroidPointsCorrected(iPinhole,2,iRotation)-10,getPinholeNumber,...
        'VerticalAlignment','bottom','HorizontalAlignment','right','Color','g');
    end
    
    figure(100); imshow(imageEdges); title(['rotation image ' num2str(iRotation)]); impixelinfo;
    hold on; plot(imageCentroidPointsCorrected(:,1,iRotation),imageCentroidPointsCorrected(:,2,iRotation),'r*','MarkerSize',12); hold off;
    pinholeImageLabels = cellstr(num2str(pinholeNumberOrdered(:,iRotation)));
    text(imageCentroidPointsCorrected(:,1,iRotation)-10,imageCentroidPointsCorrected(:,2,iRotation)-10,pinholeImageLabels,...
        'VerticalAlignment','bottom','HorizontalAlignment','right','Color','r');
    
    saveOrNot = mvdlg('Save the corrected points? [y/n]','Save points',[.8 .8 .2 .15]);
    if iRotation < 10
        centroidPointsFileName = ['im_centroid_pts_position_0' num2str(iRotation) '.mat'];
    else
        centroidPointsFileName = ['im_centroid_pts_position_' num2str(iRotation) '.mat'];
    end
    if strncmp(saveOrNot,'y',1) == 1
        if exist(centroidPointsFileName,'file')
            str = input([centroidPointsFileName ' already exists, overwrite? [y/n] '],'s');
            if strncmp(str,'y',1)
                imageCentroidPointsCorrectedSinglePosition = imageCentroidPointsCorrected(:,:,iRotation);
                imageCentroidPointsCorrectedSinglePositionLabels = pinholeNumberOrdered(:,iRotation);
                save([imageWriteDirectory centroidPointsFileName],'imageCentroidPointsCorrectedSinglePosition','imageCentroidPointsCorrectedSinglePositionLabels');
                disp([centroidPointsFileName ' saved.']);
                clear imageCentroidPointsCorrectedSinglePosition;
                clear imageCentroidPointsCorrectedSinglePositionLabels;
            else disp('Not writing file.');
            end
        else
            imageCentroidPointsCorrectedSinglePosition = imageCentroidPointsCorrected(:,:,iRotation);
            imageCentroidPointsCorrectedSinglePositionLabels = pinholeNumberOrdered(:,iRotation);
            save([imageWriteDirectory centroidPointsFileName],'imageCentroidPointsCorrectedSinglePosition','imageCentroidPointsCorrectedSinglePositionLabels');
            disp([centroidPointsFileName ' saved.']);
            clear imageCentroidPointsCorrectedSinglePosition;
            clear imageCentroidPointsCorrectedSinglePositionLabels;
        end
    end

end
% % save pinhole centroid points to file
% allCentroidPointsFileName = 'im_centroid_pts_corrected.mat';
% if exist(allCentroidPointsFileName,'file')
%     str = input([allCentroidPointsFileName ' already exists, overwrite? [y/n] '],'s');
%     if strncmp(str,'y',1)
%         save([imageWriteDirectory allCentroidPointsFileName],'imageCentroidPointsCorrected','pinholeNumberOrdered');
%         disp([allCentroidPointsFileName ' saved.']);
%     else disp('Not writing file.');
%     end
% else 
%     save([imageWriteDirectory allCentroidPointsFileName],'imageCentroidPointsCorrected','pinholeNumberOrdered');
%     disp([allCentroidPointsFileName ' saved.']);
% end

%% Step 3 - Create pinhole ground truth mapping
% This is the same for any wavelength since it is based on actual pinhole
% occluder.

ringRadii = [4 8 12];                                % ring radii [mm]
pinholePositionsActual = nan(nPinholes,2);

% Note there are three rings and each has points evenly angularly spaced.
pinholePositionsActual(1,:) = [0 0];
pinholePositionsActual(2:9,:) = ringRadii(1)*[cos((2:9)*2*pi/8) ; sin((2:9)*2*pi/8)]';
pinholePositionsActual(10:15,:) = ringRadii(2)*[cos((3:8)*2*pi/12) ; sin((3:8)*2*pi/12)]';
pinholePositionsActual(16:20,:) = ringRadii(2)*[cos((10:14)*2*pi/12) ; sin((10:14)*2*pi/12)]';
pinholePositionsActual(21:36,:) = ringRadii(3)*[cos((4:19)*2*pi/16) ; sin((4:19)*2*pi/16)]';

figure; plot(pinholePositionsActual(:,1),pinholePositionsActual(:,2),'.');
axis square; xlabel('[mm]'); ylabel('[mm]');

pinholeLabels = cellstr(num2str((1:36)'));
text(pinholePositionsActual(:,1),pinholePositionsActual(:,2),pinholeLabels,...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
title('Pinhole occluder position and numbering');

save([imageWriteDirectory 'pinhole_positions_actual.mat'],'pinholePositionsActual');

%% Step 4 - Preparation for radial distortion correction (geometric calibration)
% Need to have correspondence between pinhole centroids on the rotation
% images and the actual pinhole locations on the occluder. Make sure things
% are in the right order based on the position of the blocked pinhole and
% the mapping from Step 3.

centroidPointsList = dir([imageWriteDirectory 'im_centroid_pts_position_*.mat']);
nCentroidPoints = size(centroidPointsList,1);
% load(centroidPointsFileName);

% order the found centroid points same as actual pinhole positions
imageCentroidPointsCalibrated = nan(nPinholes,2,nCentroidPoints);
for iCentroidPoint = 1:nCentroidPoints
    load(centroidPointsList(iCentroidPoint).name);
    centroidPointsTemp = [imageCentroidPointsCorrectedSinglePosition imageCentroidPointsCorrectedSinglePositionLabels];
    centroidPointsTemp = sortrows(centroidPointsTemp,3);
    imageCentroidPointsCalibrated(:,:,iCentroidPoint) = centroidPointsTemp(:,1:2);
end

% check the new ordered points
for iCentroidPoint = 1:nCentroidPoints
    figure(101); imshow(imagesPinholeOccluder(:,:,iCentroidPoint),[0 2^12-1]); title(['rotation image ' num2str(iCentroidPoint)]); impixelinfo;
    hold on; plot(imageCentroidPointsCalibrated(:,1,iCentroidPoint),imageCentroidPointsCalibrated(:,2,iCentroidPoint),'b*','MarkerSize',12);
    pinholeLabels = cellstr(num2str((1:36)'));
    text(imageCentroidPointsCalibrated(:,1,iCentroidPoint)-10,imageCentroidPointsCalibrated(:,2,iCentroidPoint)-10,pinholeLabels,...
        'VerticalAlignment','bottom','HorizontalAlignment','right','Color','r');
    pause(.5);
%     waitforbuttonpress;
end

pinholePositionsImage = imageCentroidPointsCalibrated;
save([imageWriteDirectory 'pinhole_positions_image.mat'],'pinholePositionsImage');

end
