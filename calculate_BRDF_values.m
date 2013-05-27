% Script to calculate BRDF value from processed files

clear all; close all; clc;
addpath(pwd);

% files to process
imageDirToProcess = [pwd '/TEST DATA/inputDirectory/D10/final_grn/'];

% reference files
imageDirReference = [pwd '/TEST DATA/inputDirectory/D40/final_grn/'];

% load reference file
ref_file_dir = '/Volumes/Calibration Data (v2)/processed files/im_processed_v3/HiC/DRS/';
im_ref_file = 'im_processed_06-08-12_DRS40-2_grn.mat';
% ref_file_dir = '/Volumes/Calibration Data (v2)/processed files/im_processed_v2/grn/';
% im_ref_file = 'im_processed_08-08-12_DRS40-2_grn.mat';
load([ref_file_dir im_ref_file]);
im_processed_ref = im_processed;
wlen_str = im_ref_file(31:33);

file_list = dir([file_dir 'im_processed_*']);
num_files = size(file_list,1);
for mm = 1:num_files
% for mm = 2
% load data file                                        % CHANGE THIS ONLY!
im_file = file_list(mm).name;
% im_file = 'im_processed_08-09-12_DRS20-1_red.mat';
load([file_dir im_file]);

% load angle mappings
incident_angles_corrected_file = ['mappings_with_flip/incident_pts_corrected_with_flip_' wlen_str '.mat'];
load(incident_angles_corrected_file);
% load('viewing angle mapping/center_points');
% incident_angles_file = ['incident_points_' wlen_str '.mat'];
viewing_angles_files = ['mappings_with_flip/viewing_angle_mapping_with_flip_' wlen_str];
% load(incident_angles_file);
load([viewing_angles_files '_theta.mat']);
load([viewing_angles_files '_phi.mat']);

% % FLIP VIEWING ANGLE MAPPING AS WELL... DAMMIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% viewing_angle_mapping_theta = flipud(viewing_angle_mapping_theta);
% viewing_angle_mapping_phi = flipud(viewing_angle_mapping_phi);

% calculate circles for theta = 15, 30, 45
th = linspace(0,2*pi,360)';
[ct15_i,ct15_j] = find(round(viewing_angle_mapping_theta) == 15);
[xc15,yc15,Re15,a15] = circfit(ct15_i,ct15_j);
[ct30_i,ct30_j] = find(round(viewing_angle_mapping_theta) == 30);
[xc30,yc30,Re30,a30] = circfit(ct30_i,ct30_j);
[ct45_i,ct45_j] = find(round(viewing_angle_mapping_theta) == 45);
[xc45,yc45,Re45,a45] = circfit(ct45_i,ct45_j);

% get indices for values outside of 45 degrees
outer_pts = find(viewing_angle_mapping_theta > 45);

% calculate BRDF
BRDF_ref = (str2double(im_ref_file(26:27))/100)/pi;

num_ims = size(im_processed_ref,3);
median_val = nan(1,num_ims);
mean_val = nan(1,num_ims);
min_val = nan(1,num_ims);
max_val = nan(1,num_ims);

% if mm == 1 
%     which_ii = 162;
% %     which_ii = 30;
% elseif mm == 2
%     which_ii = 132;
% %     which_ii = 49;
% end

for ii = 1:num_ims
% for ii = 30
% for ii = 20:160
% for ii = which_ii
%     incident_angle_pts = round(im_incident_points(ii,:));
%     incident_angle_theta = viewing_angle_mapping_theta(incident_angle_pts(2),incident_angle_pts(1));
%     incident_cosine_adjustment = cosd(incident_angle_theta);

    % skip for now
    if isnan(incident_pts_corrected(ii,1))
        continue;
    end
    
    % which incident angle is this?
    incident_angle_theta = viewing_angle_mapping_theta(incident_pts_corrected(ii,1),incident_pts_corrected(ii,2))
    incident_angle_phi = viewing_angle_mapping_phi(incident_pts_corrected(ii,1),incident_pts_corrected(ii,2))
    
    im_data = im_processed(:,:,ii);
    im_ref = im_processed_ref(:,:,ii);

%     % filter/smooth images before dividing
%     % but how to choose filter intelligently?           !!!
% %     H_gauss = fspecial('gaussian',[11 11],5);
%     H_gauss = [1 1 1; 1 1 1; 1 1 1];
%     for jj = 1:10
%         im_data = imfilter(im_data,H_gauss);
%         im_ref = imfilter(im_ref,H_gauss);
%     end
% %     im_ref = 1;

    % divide by actual reference image.. will cause artifacts
    BRDF_data_temp = (im_data./im_ref)*BRDF_ref;
    
%     area_selection = 200:500;                % AVERAGE OVER A SELECTED AREA
%     % divide by median of reference 'image'
%     ref_med = median(median(im_processed_ref(area_selection,area_selection,ii)));
%     BRDF_data_temp = (im_data./ref_med)*BRDF_ref;
%     % angle corrections if dividing by single value
%     BRDF_data_temp = BRDF_data_temp./cosd(viewing_angle_mapping_theta)/incident_cosine_adjustment;

%     % calculate basic statistics
%     median_val(ii) = median(median(BRDF_data_temp(area_selection,area_selection)));
%     mean_val(ii) = mean(mean(BRDF_data_temp(area_selection,area_selection)));
%     min_val(ii) = min(min(BRDF_data_temp(area_selection,area_selection)));
%     max_val(ii) = max(max(BRDF_data_temp(area_selection,area_selection)));

%     % histogram equalization (contrast enhancement for display purposes)
%     BRDF_data_temp = histeq(BRDF_data_temp,256);
    
    % zero out values outside of 45 degrees
    BRDF_data_temp(outer_pts) = nan;

    % calculate basic statistics
    median_val(ii) = nanmedian(nanmedian(BRDF_data_temp));
    disp(median_val(ii));
    mean_val(ii) = nanmean(nanmean(BRDF_data_temp));
    min_val(ii) = nanmin(nanmin(BRDF_data_temp));
    max_val(ii) = nanmax(nanmax(BRDF_data_temp));

    % plot calculated BRDF
%     hh = figure(100); imshow(BRDF_data_temp(:,:),[0 .8]); hold on;
    hh = figure(100); imshow(BRDF_data_temp(:,:),[0 2]); hold on;
%     impixelinfo; %colorbar;
%     title(['image number ' num2str(ii)]); %axis square;
    
    % plot theta circles
    which_color = [0 .3 .1];
    plot(yc15,xc15,'.','Color',which_color);
    plot(yc15+Re15*sin(th),xc15+Re15*cos(th),'-','Color',which_color,'LineWidth',2);
    plot(yc30+Re30*sin(th),xc30+Re30*cos(th),'-','Color',which_color,'LineWidth',2);
    plot(yc15+Re45*sin(th),xc45+Re45*cos(th),'-','Color',which_color,'LineWidth',5);

    % plot incident points
    plot(incident_pts_corrected(ii,2),incident_pts_corrected(ii,1),'w*','MarkerSize',20);

    % plot max point
    [max_i,max_j] = find(BRDF_data_temp == max_val(ii));
    plot(max_j,max_i,'m*','MarkerSize',20);

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
%     % save scanline
%     % this doesn't work.. so immediately copy/paste the last 6 lines
%     theta_vec_all(mm,:) = theta_vec;
%     scan_val_all(mm,:) = scan_val;
% %     
%     % re-plot incident and max points to be on top of line
%     plot(incident_pts_corrected(ii,2),incident_pts_corrected(ii,1),'*','MarkerSize',20,'Color',[.9 .9 .9]);
%     plot(max_j,max_i,'m*','MarkerSize',20);
    
%     % flipud for full image !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     % need to do this until I correct the flip in the raw processing algorithm
%     set(gca,'ydir','normal');
    
    % save figure
    if ii < 10
        num_im = ['00' num2str(ii)];
    elseif ii < 100
        num_im = ['0' num2str(ii)];
    else
        num_im = num2str(ii);
    end
    export_fig([im_file(14:21) '_img_' num_im],'-png','-native',hh);
    
%     waitforbuttonpress; pause(.001);
%     pause(.5);
    disp(ii);
end

% figure(200); plot(1:num_ims,median_val,'b-'); hold on;
% plot(1:num_ims,mean_val,'r-');
% plot(1:num_ims,min_val,'g-');
% plot(1:num_ims,max_val,'m-');
% plot(1:num_ims,median_val,'b*');
% plot(1:num_ims,mean_val,'r*');
% plot(1:num_ims,min_val,'g*');
% plot(1:num_ims,max_val,'m*'); hold off;
% legend('median','mean','min','max');
% ylim([0 .35]);

BRDF_expected = (str2double(im_file(26:27))/100)/pi
median_of_median_val = nanmedian(median_val)
mean_of_mean_val = nanmean(mean_val)
median_of_max_val = nanmedian(max_val)
median_of_min_val = nanmedian(min_val)

% save_fname = ['BRDF_' im_file(23:29) '_' wlen_str];
% % save_fname = ['BRDF_' im_file(14:21) '_' wlen_str];
% save(save_fname,'median_val','mean_val','min_val','max_val');
% disp('File saved.');

end

