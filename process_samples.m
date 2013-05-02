% This is a full processing script that will:
% 
% 	1)   load data set parameters
%   2)   perform 'run_first' on image set to flip and copy images
%   3)   subtract dark image
% 	4)   correct for gain
%   5)   perform HDR adjustment
%   6)   undistort normalized images
%   7)   normalize based on measurement parameters ???????
%           divide out integration time ratio (normalize to 17020us)
%           divide out power reading (normalize to 1uW=1e-6W)
%
% Note: Image should be loaded once and saved once (no iterim saves). The
%       final results should be saved as .mat files.

% Want to be able to specify the main image directory, then have the
% program process each color individually. So, for each main image
% directory, get list of image subdirs (using loaded NFO file), figure out
% which ones correspond to which color, and process completely.
% So there would be three main loops.

close all; clear all; clc;
addpath(pwd);

%% Step 1: Load data set parameters
% -------------------------------------------------------------------

% Acquisition parameters (same for a given instrument, AVT GC1380 here)
% match = 'init*.tiff';                   % acquired reflectance images

% % Reference values (prevent small or large calculated values)
% integration_time_reference = 17020;     % [us]
% power_value_reference = 1e-6;         % [W]

mainFileDir = [pwd '/TEST DATA/inputDirectory/S16/'];
trueBitDepth = 12;                   	% number of recorded bits per pixel
integrationTimeMinUS = 17020;        	% microseconds

% parameters (medium, wavelength, sample location, binning, QC intensity,
%             integration time, QC timing, gain value)
nfoFileList = dir([mainFileDir '*.nfo']);
nfoFileID = fopen([mainFileDir nfoFileList(1).name]);
subImageDirs = textscan(nfoFileID,'%s','Delimiter','\n');
nImageDirs = size(subImageDirs{1},1);

for iImageDir = 1:nImageDirs

    subImageDir = subImageDirs{1}{iImageDir};
    
    imageSetName{1}{iImageDir} = subImageDir(1:6);
    
    paramBeginPos = find(subImageDir == '(',1,'last');
    paramEndPos = find(subImageDir == ')',1,'last');
    paramStr = subImageDir(paramBeginPos+1:paramEndPos-1);
    
    % Separate individual parameters from parameter string
    paramInfo = textscan(paramStr,'%s %s %s %s %s %s %s %s','Delimiter',',');
     
    % Set wavelength color string
    wavelengthStr{1}{iImageDir} = char(paramInfo{2});

    % Set integration time
    integrationTimeStrMS = char(paramInfo{6});          % milliseconds [ms]
    integrationTimeMS = str2double(integrationTimeStrMS(1:end-2));
    integrationTime{1}{iImageDir} = integrationTimeMinUS*(integrationTimeMS/round(integrationTimeMinUS/1000));
    
    % Set gain value
    gainValueTemp = char(paramInfo{8});
    if numel(gainValueTemp) == 2
        gainValue{1}{iImageDir} = str2double(gainValueTemp(2));
    elseif numel(gainValueTemp) == 3
        gainValue{1}{iImageDir} = str2double(gainValueTemp(2:3));
    end

end

%% Step 2: Perform 'run_first' on image sets
% This will flip and copy all images to specified working directory.

imageExtension = 'tiff';
for iImageDir = 1:nImageDirs
    inputDirectory = [mainFileDir subImageDirs{1}{iImageDir} '/'];
    outputDirectory = [inputDirectory '/init/'];
    if ~exist(outputDirectory,'dir')
        mkdir(outputDirectory);
    end
    run_first(inputDirectory,outputDirectory,imageExtension);
    disp(['Processed directory ' '~/' subImageDirs{1}{iImageDir} '/']);
end

%% Steps 3 through ?
% average all dark images for a given set and subtract that from each image
% make sure to change negative values to 0. For image '2' use dark image
% '3' since sometimes '1' is weird.

ctrRedWD = 0; ctrGrnWD = 0; ctrBluWD = 0;
for iImageDir = 1:nImageDirs
    initDirectory = [mainFileDir subImageDirs{1}{iImageDir} '/init/'];
    
    imageInputList = dir([initDirectory '*.' imageExtension]);
    nImages = length(imageInputList);
    
% -------------------------------------------------------------------
%   Step 3a: Calculate dark image and store MAT file.
% -------------------------------------------------------------------
    % get image parameters
    imageTemp = imread([initDirectory imageInputList(1).name]);
    imageHeight = size(imageTemp,1);
    imageWidth = size(imageTemp,2);
    nPixels = imageHeight*imageWidth;

    darkImageDir = [mainFileDir subImageDirs{1}{iImageDir} '/dark/'];
    if ~exist(darkImageDir,'dir')
        mkdir(darkImageDir);
        imageDark = zeros([imageHeight imageWidth nImages/2]);
        % loop over dark (odd) images
        for iImage = 1:2:nImages-1
            imageName = imageInputList(iImage).name;
            imageTemp = imread([initDirectory imageName]);
            imageDark(:,:,iImage) = imageTemp;
        end
        % calculate average dark frame
        imageDark = median(imageDark,3);
        disp('Saved dark image file.');
        save([darkImageDir 'imageDark.mat'],'imageDark');
    else
        disp('Already saved dark image file.');
    end

    currentWavelength = wavelengthStr{1}{iImageDir};
    
    if strncmp('r',currentWavelength,1)
        ctrRedWD = ctrRedWD + 1;
        redImageDirs(ctrRedWD) = iImageDir;
    elseif strncmp('g',currentWavelength,1)
        ctrGrnWD = ctrGrnWD + 1;
        grnImageDirs(ctrGrnWD) = iImageDir;
    elseif strncmp('b',currentWavelength,1)
        ctrBluWD = ctrBluWD + 1;
        bluImageDirs(ctrBluWD) = iImageDir;
    end
end
clear imageTemp imageDark;

% number of directories for each of R G B wavelengths
nWavelengthDirs = [numel(redImageDirs) numel(grnImageDirs) numel(bluImageDirs)];
% indices for wavelength directories (order R G B)
allWavelengthDirs = [redImageDirs grnImageDirs bluImageDirs];

% iterate over number of wavelengths
nWavelengths = 3;
for iWavelength = 1:nWavelengths
    
    % iterate over number of directories for each wavelength (order R G B)
    % may be simpler to just use a counter
    if iWavelength == 1
        whichWavelengthDirs = allWavelengthDirs(1:nWavelengthDirs(iWavelength));
    elseif iWavelength == 2
        whichWavelengthDirs = allWavelengthDirs(nWavelengthDirs(iWavelength-1)+(1:nWavelengthDirs(iWavelength)));
    elseif iWavelength == 3
        whichWavelengthDirs = allWavelengthDirs(nWavelengthDirs(iWavelength-2)+nWavelengthDirs(iWavelength-1)+(1:nWavelengthDirs(iWavelength)));
    end
    
    nWavelengthSubDirs = nWavelengthDirs(iWavelength);
    % create and populate nWavelengthSubDirs-dimensional arrays
    clear wavelengthSubDir wavelengthDir wavelengthImageInputList;
    clear wavelengthDarkDir imageDark powerFileNames setColor;
    for iWavelengthSubDir = 1:nWavelengthSubDirs
        wavelengthSubDir(:,:,iWavelengthSubDir) = subImageDirs{1}{allWavelengthDirs(whichWavelengthDirs(iWavelengthSubDir))};
        wavelengthDir(:,:,iWavelengthSubDir) = [mainFileDir wavelengthSubDir(:,:,iWavelengthSubDir) '/init/'];
        wavelengthImageInputList(:,:,iWavelengthSubDir) = dir([wavelengthDir(:,:,iWavelengthSubDir) '*.' imageExtension]);
        % load dark image
        wavelengthDarkDir(:,:,iWavelengthSubDir) = [mainFileDir wavelengthSubDir(:,:,iWavelengthSubDir) '/dark/'];
        imageDarkStruct = load([wavelengthDarkDir(:,:,iWavelengthSubDir) 'imageDark.mat'],'imageDark');
        imageDark(:,:,iWavelengthSubDir) = imageDarkStruct.imageDark;
        wavelengthSetNameList = imageSetName{1}{allWavelengthDirs(whichWavelengthDirs(iWavelengthSubDir))};
        powerFileNames(:,:,iWavelengthSubDir) = [wavelengthSetNameList '.txt'];
        setColor(:,:,iWavelengthSubDir) = wavelengthStr{1}{allWavelengthDirs(whichWavelengthDirs(iWavelengthSubDir))};
    end
 
    % -------------------------------------------------------------------
    %   Step 5a: Get mean power ratio values.
    % -------------------------------------------------------------------
    powerDir = [mainFileDir 'power/'];
    powerFileBase = powerFileNames(:,:,1);
    medianPowerRatios(1) = 1;
    for iWavelengthSubDir = 2:nWavelengthSubDirs
        powerFileComp = powerFileNames(:,:,iWavelengthSubDir);
        medianPowerRatios(iWavelengthSubDir) = multiple_exposures_mean_power_ratio(powerDir,powerFileBase,powerFileComp);
    end
    
    % loop over data (even) images
    nImages = size(wavelengthImageInputList,1);
    for iImage = 2:2:nImages
        for iWavelengthSubDir = 1:nWavelengthSubDirs
            imageName(:,:,iWavelengthSubDir) = wavelengthImageInputList(iImage,iWavelengthSubDir).name;
            imageTemp(:,:,iWavelengthSubDir) = double(imread([wavelengthDir(:,:,iWavelengthSubDir) imageName(:,:,iWavelengthSubDir)]));
            
            % -------------------------------------------------------------------
            %   Step 3b: Subtract dark image.
            % -------------------------------------------------------------------
            % subtract dark image
            imageTemp(:,:,iWavelengthSubDir) = imageTemp(:,:,iWavelengthSubDir) - imageDark(:,:,iWavelengthSubDir);
            % correct for negative values
            imageTempTemp = imageTemp(:,:,iWavelengthSubDir);
            imageTempTemp(imageTempTemp < 0) = 0;
            imageTemp(:,:,iWavelengthSubDir) = imageTempTemp;
            
            % -------------------------------------------------------------------
            %   Step 4: Correct for gain response.
            % -------------------------------------------------------------------
            currentGainValue = gainValue{1}{allWavelengthDirs(whichWavelengthDirs(iWavelengthSubDir))};
            if currentGainValue ~= 0
                imageTemp(:,:,iWavelengthSubDir) = process_samples_gain_correction(imageTemp(:,:,iWavelengthSubDir),currentGainValue);
                disp('Gain response corrected.');
            end
        end
        
        imageTempHDR = nan([imageHeight imageWidth]);
        % iterate over all pixels (slow!)
        for iPixel = 1:nPixels
            % -------------------------------------------------------------------
            %   Step 5b: Get weighting values.
            % -------------------------------------------------------------------
            pixelValues = nan(1,nWavelengthSubDirs);
            weightingValues = nan(1,nWavelengthSubDirs);
            for iWavelengthSubDir = 1:nWavelengthSubDirs
                imageTempTemp = imageTemp(:,:,iWavelengthSubDir);
                pixelValues(iWavelengthSubDir) = imageTempTemp(iPixel);
                weightingValues(iWavelengthSubDir) = multiple_exposures_weighting_values(pixelValues(iWavelengthSubDir),trueBitDepth);
            end
            % -------------------------------------------------------------------
            %   Step 5c: Get new pixel values.
            % -------------------------------------------------------------------
            newPixelValue = multiple_exposures_calculate_new(pixelValues,weightingValues,medianPowerRatios);
            imageTempHDR(iPixel) = newPixelValue;
        end
        
        % -------------------------------------------------------------------
        % Step 6: Undistort normalized image.
        % -------------------------------------------------------------------
        % Note that image will be square (imageHeight X imageHeight) after this.
        imageFinal = calibration_geometric_undistort_image_single(imageTempHDR,setColor(:,:,iWavelengthSubDir));

        figure; imshow(imageFinal,[0 2^trueBitDepth-1]); impixelinfo;
        waitforbuttonpress;
%         pause(0.5);
        disp(['Processed image ' num2str(iImage) ' of ' num2str(nImages)]);
    end
    
    disp('Processed one wavelength...');
end

%%

% % % -------------------------------------------------------------------
% % % Step 4: Normalize by integration time and power value.
% % % -------------------------------------------------------------------
% % %     incident_angle_pts = round(im_incident_points(im_ctr,:));
% % %     incident_angle_theta = viewing_angle_mapping_theta(incident_angle_pts(2),incident_angle_pts(1));
% % %     cosine_adjustment(im_ctr) = cosd(incident_angle_theta);
% % %     power_value_adjusted(im_ctr) = power_values(ii/2)*cosine_adjustment(im_ctr);
% % 
% %     im_temp = im_temp/(integration_time_us/integration_time_reference);
% % %     im_temp = im_temp/((power_values(ii/2)/cosine_adjustment(im_ctr))/power_value_reference);
% %     im_temp = im_temp/(power_values(ii/2)/power_value_reference);
% % %     im_temp = im_temp/(power_value_adjusted(im_ctr)/power_value_reference);
% % % -------------------------------------------------------------------
% % % Save and display corrected image.
% % % -------------------------------------------------------------------
% %     im_processed(:,:,im_ctr) = im_temp;
% % %     figure(200); imagesc(im_processed(:,:,im_ctr)); impixelinfo; colorbar; axis square;
% % %     title(['image number ' num2str(ii)]);
% % % %     title(['image number ' num2str(ii) ' shown and has median value ' num2str(median_val(im_ctr))]);
% % %     waitforbuttonpress;
% %     pause(.05);
% % end
% % 
% % save_fname = ['im_processed' save_power_fname(13:end)];
% % if ~exist(['im_processed_' powr_file(1:7) '.mat'],'file')
% %     disp('Saving file...');
% %     save(save_fname,'im_processed');
% %     disp('File saved.');
% % end
% % 
% % % num_ims_actual = num_ims/2;
% % % figure; plot(1:num_ims_actual,median_val,'b-'); title('averaged values');
% % % hold on; plot(1:num_ims_actual,mean_val,'r-');
% % % plot(1:num_ims_actual,median_val,'b.');
% % % hold off; legend('median','mean');
% % 
% % end
