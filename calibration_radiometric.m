% Routine to obtain radiometric response from images of same scene but
% different (known) exposures.
% 
% Uses images created from 'avg_images.m' for radiometric calibrations.
% 
% Algorithm is mostly following Grossberg & Nayar (2003). Intensity
% mapping functions are found using image pairs (ordered dark to light)
% and comparison intensities of corresponding pixels.
% 
% This routine works by:
% 
%   Determination of comparagram (count intensities of corresponding
%   pixels, assuming perfectly registered). Fit values to a polynomial
%   to obtain intensity mapping function.
%
%   Linear search performed to obtain inverse intensity mapping functions.
% 
%   Exposure ratios are known since exposures for each image are known. 
%   Inverse response function can then be calculated from the intensity
%   mapping functions (fitting B1 and B2 values to a polynomial).
% 
% NOTE TO ANY FUTURE USERS: Sorry for the poor coding style. I didn't feel
% the need to rewrite this since it only really needs to be run once for a
% given camera. The AVT ProSilica cameras generally have a highly linear
% response anyway (since they are designed for machine vision).

close all; clear all; clc;

%% Load averaged images and generate comparagrams
loc_ims = [pwd '/calibration_data/radiometric/camera_response_6exp/'];
im_list = dir([loc_ims '*.tiff']);               % get list of images
im_list = im_list(1:end-1);                                                             % IGNORE LAST ONE (SATURATION)
n = numel(im_list);                                     % number of images
im_temp = imread([loc_ims im_list(1).name]);
[ht wd] = size(im_temp);                                % image dimensions

bitspp = 12;                                            % bits per pixel
L = 2^bitspp;                                           % number of intensity levels
Lm1 = L - 1;
intens = 0:Lm1;                                         % actual intensity values

cgram = zeros(L,L,n-1);                                 % comparagram
im_avg_exp = nan(1,n);
for ii = 1:n-1
    im_name = im_list(ii).name;
    im_temp = imread([loc_ims im_name]);                % load average image
    im_temp = im_temp(:);
    
    im_name_next = im_list(ii+1).name;
    im_temp_next = imread([loc_ims im_name_next]);      % load next average image
    im_temp_next = im_temp_next(:);

    for jj = 1:(ht*wd)                                  % create comparagram
        B1_val = im_temp(jj);
        B2_val = im_temp_next(jj);
        cgram(B2_val+1,B1_val+1,ii) = cgram(B2_val+1,B1_val+1,ii) + 1;
    end
    
    im_avg_exp(ii) = str2double(im_name(end-12:end-7)); % calc exposure ratios
    im_avg_exp(ii+1) = str2double(im_name_next(end-12:end-7));

end

%% Obtain intensity mapping functions from comparagrams
extract_Bval = nan(L,1,n-1);
tau_imf = extract_Bval;
tau_imf_inv = tau_imf;
for ii = 1:n-1                                          % extract values from comparagram
    [B2_idx, B1_idx] = find(cgram(:,:,ii) > 0);         % ignore zero counts
    for jj = 1:length(B1_idx)
        extract_Bval(B1_idx(jj),1,ii) = B2_idx(jj);     % assign values
    end
end

figure;
for ii = 1:n-1                                          % fit polynomial to comparagram
    B_temp = extract_Bval(:,1,ii);
    plot(intens,B_temp,'.'); hold all;
    axis([0 Lm1 0 Lm1]);
    
    B_fit = polyfit(intens(~isnan(B_temp))',B_temp(~isnan(B_temp)),1);      % maybe try higher order here..
    tau_imf(:,1,ii) = polyval(B_fit,intens);

   tau_imf_inv_temp = nan(1,2^bitspp);
   for jj = 1:2^bitspp
        bindx = find(tau_imf(:,1,ii) == intens(jj));
        if ~isempty(bindx)
            tau_imf_inv_temp(jj) = bindx(end);
        else
            continue;
        end
    end
    inv_nan = ~isnan(tau_imf_inv_temp);
    if sum(inv_nan) ~= 0
        tau_imf_inv(:,1,ii) = interp1(intens(inv_nan),tau_imf_inv_temp(inv_nan),intens);
    else
        tau_imf_inv(:,1,ii) = tau_imf_inv_temp;
    end
end

title('Intensity mapping functions and fits'); axis square;
xlabel('Intensity in darker image');
ylabel('Intensity in lighter image');
% legend('17ms to 34ms','34ms to 51ms','51ms to 68ms','68ms to 85ms',...
%     '85ms to 102ms','Location','SouthEast');

for ii = 1:n-1
    plot(intens,tau_imf(:,1,ii),'k');
end

%% Determine inverse response function of camera
% For each pair of images, generate a set of pairs of intensity values B1
% and B2 = tau(B1). One pair is (B1,tau(B1)) and other pair is
% (invtau(B2)),B2).
%
% Need to solve the equation g(tau(B))=kg(B) for g()
% Get constraints for B and tau(B) (B1, B2) pairs and solve for
% coefficients of g, assuming polynomial. k is exposure ratio

% exp_ratio = [im_avg_exp(2)/im_avg_exp(1) im_avg_exp(3)/im_avg_exp(2)...
%     im_avg_exp(4)/im_avg_exp(3) im_avg_exp(5)/im_avg_exp(4) im_avg_exp(6)/im_avg_exp(5)];

exp_ratio = nan(1,n-1);
for ii = 1:n-1
    exp_ratio(ii) = im_avg_exp(ii+1)/im_avg_exp(ii);
end

B1tB1 = nan(L,2,n-1);                                   % point pairs
iB2B2 = B1tB1;
for ii = 1:n-1
    for jj = 1:L
        B1 = jj;
        tB1 = tau_imf(B1,1,ii);
        B2 = B1;
        iB2 = tau_imf_inv(B2,1,ii);                     % NaN if no matching value
        B1tB1(jj,:,ii) = [B1 tB1] / Lm1;
        iB2B2(jj,:,ii) = [iB2 B2] / Lm1;
    end
    
    B1tB1(:,1,ii) = B1tB1(:,1,ii)*exp_ratio(ii);        % multiply by exposure ratio
    iB2B2(:,1,ii) = iB2B2(:,1,ii)*exp_ratio(ii);
end

B1_constr = [B1tB1(:,1,1) ; iB2B2(:,1,1)];              % separate dark and light points
B2_constr = [B1tB1(:,2,1) ; iB2B2(:,2,1)];
for aa = 2:n-1
    B1_constr = [B1_constr ; B1tB1(:,1,aa) ; iB2B2(:,1,aa)];
    B2_constr = [B2_constr ; B1tB1(:,2,aa) ; iB2B2(:,2,aa)];
end

B1nan = isnan(B1_constr);                                % ignore NaNs
B1_constr = B1_constr(~B1nan);
B2_constr = B2_constr(~B1nan);
B2nan = isnan(B2_constr);                                % ignore NaNs
B1_constr = B1_constr(~B2nan);
B2_constr = B2_constr(~B2nan);

B1_constr = B1_constr/max(B1_constr);
B2_constr = B2_constr/max(B2_constr);

% polynomial order DOESN'T HAVE TO BE ORDER 6! (make sure it's not ill-conditioned though!)
poly_ord = 2;
irf_poly_coeffs = polyfit(B2_constr,B1_constr,poly_ord);
irf_poly_values = polyval(irf_poly_coeffs,B2_constr);

figure; plot(B2_constr,irf_poly_values,'ko','LineWidth',2); hold on;
plot(B2_constr,B1_constr,'b.');
axis([0 1 0 1]); axis square;
title('Inverse response curve');
xlabel('Image intensity'); ylabel('Image irradiance');

