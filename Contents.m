% MATLAB Support Package for OSMAR instrument
% Version 2.0 (R2012a), J. M. Haag,  07-Mar-2013
%
% This is a set of MATLAB code that allows for calibration and processing
% of data gathered from the OSMAR (Optical Scatterometer for Measurement of
% Angular Reflectance) instrument.
%
% FILES:
%
% Contents                      - this file (needed for matlab help)
% readme.txt                    - file containing additional information
% 
% run_first                     - run this file first on any image set
%
% calibration_data              - folder containing calibration data
% calibration_geometric_*       - perform geometric camera calibration
% calibration_radiometric       - perform radiometric camera calibration
% 
% angle_mapping_viewing         - obtain viewing angle mapping
% angle_mapping_incident        - obtain incident angle mapping
% angle_mappings                - folder containing angle mappings
% 
% generate_source_images        - create incident angle projector images
% generated_source_images       - folder containing source images
% 
% process_samples               - process data from sample measurements
% 
% calculate_BRDF_values         - calculate BRDF values for processed samples
% calculate_BRDF_statistics     - calculate BRDF statistics and error
% calculated_values             - folder containing calculated values
% 
% display_BRDF_values           - display BRDF values
% display_BRDF_statistics       - display BRDF statistics and error
% display_plots                 - folder containing plots of results
%
%
% HISTORY:
%
% Ver  1.0 - --- 2012  - Initial version
% Ver  2.0 - Mar 2013  - Beginning of complete rewrite and refactoring
