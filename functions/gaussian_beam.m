function U = gaussian_beam(lambda,P,w0,z0,z,matrix_size,pixel_size)
    k = 2*pi/lambda;
    w = w0*(1+(z/z0)^2)^(1/2);
    R = z*(1+(z0/z)^2);
    I0 = 2*P/(pi*w^2);
    A0 = sqrt(2*I0/(pi*w0^2));
    x = ((1:matrix_size) - (matrix_size+1)/2) * pixel_size;  % centered x-coordinates
    y = ((1:matrix_size) - (matrix_size+1)/2) * pixel_size;  % centered y-coordinates
    [X, Y] = meshgrid(x, y);
    coord = cat(3, X, Y); 
    rho2 = coord(:,:,1).^2+coord(:,:,2).^2;
    U = A0.*w0./w.*exp(-rho2./w.^2).*exp(-1i.*k.*z-1i.*k.*rho2./(2.*R)+1i.*atan(z/z0));
end