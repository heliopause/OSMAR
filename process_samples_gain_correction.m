function [imageTemp] = process_samples_gain_correction(imageTemp,currentGainValue)
% Function to perform gain correction

%           gain (in dB) = 20log(So/Si)             So -> recorded intensity
%           Si = So/10^(gain/20)                    Si -> original intensity

    imageTemp = imageTemp./10^(currentGainValue/20);

end