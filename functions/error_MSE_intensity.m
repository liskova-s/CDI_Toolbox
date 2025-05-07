function int_error = error_MSE_intensity(reconstructed,measured)
    % MSE of intensity in fourier domain
     %R = (abs(reconstructed)./max(max(abs(reconstructed)))).^2;
     %M = (abs(measured)./max(max(measured))).^2;
     int_error = sum(sum((reconstructed - measured).^2))/sum(sum(measured.^2));
end