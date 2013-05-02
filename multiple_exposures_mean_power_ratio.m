function [medianPowerRatio] = multiple_exposures_mean_power_ratio(powerDir,powerFileBase,powerFileComp)
% Function to get mean power ratios

% Step 1 - Calculate average power ratios
%          Load two power measurement sets, determine where they are
%          correlated, shift based on correlation, get ratio, get mean and
%          variance of ratio.

% 		1. 	determine mean power ratio (mPR) for a given image set as described above

% close all; clear all; clc;

% powerDir = '/Users/justin/Documents/School/Scripps/Jaffe Lab/MURI project/BRDF project/programs/instrument_revision/OSMAR/TEST DATA/inputDirectory/power/';
% powerFiles = {'S16a-4.txt','S16a-5.txt'};
% powerFileBase = powerFiles{1};
% powerFileComp = powerFiles{2};

baseFile = fopen([powerDir powerFileBase]);
compFile = fopen([powerDir powerFileComp]);

baseFileData = textscan(baseFile,'%s%s%s');
fclose(baseFile);
compFileData = textscan(compFile,'%s%s%s');
fclose(compFile);

ratioStr = ['ratio of ' powerFileComp ' to ' powerFileBase];
nRatioStr = numel(ratioStr);
hLineSep = repmat('=',1,nRatioStr);
disp(hLineSep);
disp(ratioStr);
disp(hLineSep);

% extract the power values
basePowerDataRaw = baseFileData{3};
basePowerDataRaw = basePowerDataRaw(3:end);
basePowerDataRaw = basePowerDataRaw(1:2:end);
compPowerDataRaw = compFileData{3};
compPowerDataRaw = compPowerDataRaw(3:end);
compPowerDataRaw = compPowerDataRaw(1:2:end);

% convert strings to doubles (probably a way to do this directly from structure)
nMeasurements = min([numel(basePowerDataRaw) numel(compPowerDataRaw)]);
measurementVector = 1:nMeasurements;
basePowerData = nan(1,nMeasurements);
compPowerData = nan(1,nMeasurements);
for iMeasurement = 1:nMeasurements
    basePowerData(iMeasurement) = str2double(basePowerDataRaw{iMeasurement});
    compPowerData(iMeasurement) = str2double(compPowerDataRaw{iMeasurement});
end
basePowerData = basePowerData(end:-1:1);                    % reverse order (see time stamps)
compPowerData = compPowerData(end:-1:1);                    % reverse order (see time stamps)

% figure; plot(measurementVector,basePowerData,'-b'); hold on;
% plot(measurementVector,compPowerData,'-r');
% legend('basePowerData','compPowerData');

[powerCorrelation,powerLags] = xcorr(basePowerData,compPowerData,'coeff');
% correlationValue = powerCorrelation(find(powerLags == 0))
[correlationValue,delayValue] = max(powerCorrelation)
% figure; subplot(3,1,1); plot(basePowerData);
% subplot(3,1,2); plot(compPowerData);
% subplot(3,1,3); plot(powerCorrelation(measurementVector));

% find peaks in power data, ignore extremely low values
baseThreshold = .01*max(basePowerData);
[basePowerPeaks,basePowerPeakLocations] = findpeaks(basePowerData,'MINPEAKHEIGHT',baseThreshold);
nBasePeaks = numel(basePowerPeakLocations);
compThreshold = .01*max(compPowerData);
[compPowerPeaks,compPowerPeakLocations] = findpeaks(compPowerData,'MINPEAKHEIGHT',compThreshold);
nCompPeaks = numel(compPowerPeakLocations);

defPeakVec = 180;
diffBaseDef = defPeakVec-nBasePeaks;
basePowerPeaks = [basePowerPeaks zeros(1,diffBaseDef)];
diffCompDef = defPeakVec-nCompPeaks;
compPowerPeaks = [compPowerPeaks zeros(1,diffCompDef)];

% figure; plot(basePowerPeaks,'-b'); hold on; plot(compPowerPeaks,'-r');
[peaksCorrelation,peakLags] = xcorr(basePowerPeaks,compPowerPeaks,'coeff');
[peaksCorrelationValue,peaksDelayValue] = max(peaksCorrelation)
% figure; subplot(3,1,1); plot(basePowerPeaks);
% subplot(3,1,2); plot(compPowerPeaks);
% subplot(3,1,3); plot(peaksCorrelation(1:numel(basePowerPeaks)));

delayDiff = peaksDelayValue - defPeakVec;
if delayDiff > 0
    basePowerPeaksAdjust = [basePowerPeaks zeros(1,delayDiff)];
    compPowerPeaksAdjust = [zeros(1,delayDiff) compPowerPeaks];
elseif delayDiff < 0
    basePowerPeaksAdjust = [zeros(1,abs(delayDiff)) basePowerPeaks];
    compPowerPeaksAdjust = [compPowerPeaks zeros(1,abs(delayDiff))];    
elseif delayDiff == 0
    basePowerPeaksAdjust = basePowerPeaks;
    compPowerPeaksAdjust = compPowerPeaks;
end

% peaksDelayValue is 180 if the two measurement sets are perfectly
% correlated, so need to shift based on the difference in delay value
% before taking the ratio
% figure; plot(basePowerPeaksAdjust,'-b'); hold on; plot(compPowerPeaksAdjust,'-r');
[peaksCorrelationAdjust,peakLagsAdjust] = xcorr(basePowerPeaksAdjust,compPowerPeaksAdjust,'coeff');
[peaksCorrelationValueAdjust,peaksDelayValueAdjust] = max(peaksCorrelationAdjust)
% figure; subplot(3,1,1); plot(basePowerPeaksAdjust);
% subplot(3,1,2); plot(compPowerPeaksAdjust);
% subplot(3,1,3); plot(peaksCorrelationAdjust(1:numel(basePowerPeaksAdjust)));

powerRatio = basePowerPeaksAdjust ./ compPowerPeaksAdjust;
powerRatio(isinf(powerRatio)) = nan;
figure; plot(powerRatio,'-b'); hold on;
powerRatioFilt = medfilt1(powerRatio,10);
plot(powerRatioFilt,'-r');
title(ratioStr);

meanPowerRatio = nanmean(powerRatio)
medianPowerRatio = nanmedian(powerRatio)
variancePowerRatio = nanvar(powerRatio);
stddevPowerRatio = nanstd(powerRatio);

meanPowerRatioFilt = nanmean(powerRatioFilt)
medianPowerRatioFilt = nanmedian(powerRatioFilt)

hLineEnd = repmat('~',1,nRatioStr);
disp(hLineEnd);

end
