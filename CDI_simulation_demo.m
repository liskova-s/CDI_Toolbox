% Demo of shrinkwrap reconstruction of object with size ~ 780 um from its
% diffraction pattern intensity. Shrinkwrap with HIO update core performed with Fresnel
% propagation propFR and inverse ipropFR.

clear;
clc;
addpath('./functions/')

lambda = (792/17)*1e-9;     % source wavelength [m]
squaresize = 1e-6;          % size of one object domain pixel [m]
z = 10;                     % object- detector distance [m]
SNR = 1e5;                  % signal-to-noise ratio noise level to add to simulated diffraction pattern

% load prepared/measured intensity pattern and object mask
% intensity_pattern = abs(cell2mat(struct2cell(load("intensity_demo.mat")))); 
% reference = abs(cell2mat(struct2cell(load("reference_demo.mat"))));
% object_mask = abs(cell2mat(struct2cell(load("mask_demo.mat"))));


[intensity_pattern, object_mask, reference] = model_intensity_pattern(squaresize, z, lambda, SNR);

% Parameter specification: 
% HIO + shrinkwrap 
beta = 0.8;                         %HIO update coefficient
N_iter = 150;                       % number of iterations
shrinkwrap_iteration_start = 20;    % starting shrinkwrap support updates after 
shrinkwrap_iteration_end = 150;     % last shrinkwrap support update
shrinkwrap_sigma = 50;              % sigma for initial blur while forming support update
shrinkwrap_freq = 10;               % updating support every {shrinkwrap_freq}-th iteration
shrinkwrap_threshold = 0.18;        % threshold value for shrinkwrap support forming
end_sigma = 7;                     % sigma for last blur while forming support update
delta_sigma = (shrinkwrap_sigma - end_sigma)/((shrinkwrap_iteration_end-shrinkwrap_iteration_start)/shrinkwrap_freq - 1);

rng(4, "twister");                  % fixing random seed for result replicability

% Pre-initialization for speed
error_object_dom = zeros([N_iter,1]);
error_fourier_MSE = zeros([N_iter,1]);
count = 1;
shrinkwrap_count = 1;

% Forming initial input
measured_amp = sqrt(double(intensity_pattern));     % get amplitude from intensity
random_phase = 2 * pi * rand(size(measured_amp));   % initialize random phase

G0 = measured_amp.*exp(1i * random_phase);
[~, ~, input_g] = ipropFR(G0, lambda, z, squaresize); 

for k = 1:N_iter
    fprintf("Iteration %d/%d\n", k, N_iter);
    
    % 1) FT of input
    [~,~,G_k]= propFR(input_g,lambda,z,squaresize); 

    % 2) Replacing amplitude with original (measured) one
    G_kk = measured_amp.*exp(1i.*angle(G_k));

    % 3) IFFT back to object domain
    [~,~,output_g] = ipropFR(G_kk,lambda,z,squaresize);
    output_g = real(output_g);  % Ensure the output is real
    input_g = real(input_g); 

    error_object_dom(k) = error_object_domain(output_g,object_mask);
        
    %SHRINKWRAP support update
    if mod(k,shrinkwrap_freq)==0 && k>=shrinkwrap_iteration_start && k<shrinkwrap_iteration_end
        N = size(object_mask,1);
        [X, Y] = meshgrid(-N/2+1:N/2, -N/2+1:N/2);
        sigma = shrinkwrap_sigma - shrinkwrap_count*delta_sigma;
        G = 1 * exp(-(X.^2 + Y.^2) / (2 * sigma^2)); % gaussian kernel 
        ob_mask = abs(ifftshift(ifft2(fft2(fftshift(G)).*fft2(fftshift(abs(input_g)))))); %convolution with G kernel -> blur
        ob_mask = rescale(abs(ob_mask), 0, 1);
        temp_mask = zeros(size(ob_mask));
        temp_mask(ob_mask>shrinkwrap_threshold) = 1; % threshold blurred estimate
        object_mask = object_mask.*temp_mask; % update object mask
        shrinkwrap_count = shrinkwrap_count + 1;
    end
  
     % error in fourier domain for current prediction
    [~,~,output_in_fd]= propFR(output_g,lambda,z,squaresize); 
    error_fourier_MSE(k)= error_MSE_intensity(abs(output_in_fd).^2, intensity_pattern);
               
    % constraints in the object plane
    % 4) positivity + object constraints: form new input
    output_g = real(output_g);  % Take the real part before applying the constraint
    
    % HIO update: form new input
    mask_violation = (object_mask == 0)|(output_g <= 0); % positivity + object mask constraints
    new_g = output_g;
    
    % Intensity constraint
    temp = abs(output_g);
    temp(temp>1)=1;
    multi_output_g(:,:)=temp;
    
    new_g(mask_violation) = input_g(mask_violation) - beta * output_g(mask_violation);
    input_g = new_g;  % Update the input for the next iteration
    
    
    % visualisation
    figure(1)
    visual = abs(output_g);
    visual(visual>1)=1;
    imagesc(rot90(rot90(visual)))
    axis equal; axis image;
    colorbar()
    title(sprintf("Iteration number %d",k))
    
    figure(2)
    subplot 121
    imagesc(rot90(rot90(object_mask))+rot90(rot90(reference)))
    axis equal; axis image;
    title("Object mask + reference")
    subplot 122
    imagesc(rot90(rot90(object_mask))+reference)
    axis equal; axis image;
    title("Object mask + reference")
    pause(0.01);
end
 
figure(4)
plot(error_fourier_MSE)
xlabel("Iterations")
ylabel("E_{FD} Fourier domain MSE")

figure(5)
plot(error_object_dom)
xlabel("Iterations")
ylabel("Fienup's object domain error")


