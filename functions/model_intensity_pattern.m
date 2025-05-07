
function [intensity, mask, image] = model_intensity_pattern(squaresize, z, lambda, SNR)
addpath('./functions')

imageorig = abs(cell2mat(struct2cell(load("demo_logo.mat"))));
imageorig = imbinarize(imageorig(:,:,1));

% Apply padding for oversampling
padding = 801;
image = padarray(imageorig , [padding, padding], 'both');
mask = zeros(size(image));
mask(padding+200:end-padding-200,padding+80:end-padding-80) = 1; %more tailored for eli logo

% Simulate intensity pattern
[~,~,farfield] = propFR(image,lambda,z,squaresize);
intensity = abs(farfield).^2;
    
% Apply for reduction of dynamic range: 16 bit (65535 counts)
%scaled_vector = rescale(intensity, 0, 65535);
%intensity  = double(uint16(scaled_vector)); 

% Introduce noise
signal_power = mean(intensity(:).^2);
noise_power = signal_power / SNR; % asuuming Gaussian noise model
noise = sqrt(noise_power) * randn(size(intensity)); 
intensity = intensity + noise;
   
end

