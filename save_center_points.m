% Script to generate MAT file of center points from geometric calibration

geoCalibDirectory = [pwd '/calibration_data/geometric/'];

cc_red = load([geoCalibDirectory 'red/Calib_Results.mat'],'cc');
cc_red = cc_red.cc;

cc_grn = load([geoCalibDirectory 'grn/Calib_Results.mat'],'cc');
cc_grn = cc_grn.cc;

cc_blu = load([geoCalibDirectory 'blu/Calib_Results.mat'],'cc');
cc_blu = cc_blu.cc;

save([geoCalibDirectory '/center_points.mat'],'cc_red','cc_grn','cc_blu');

