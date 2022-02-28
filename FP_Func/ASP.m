function [ prop_image ] = ASP(image, d, wl, psize, NA)
%This function simulates the complex image propagation using angular spectrum propagation method.
% d:     propagation distance (Unit: meter)
% wl:    wavelength (Unit: meter)
% psize: pixel size of input image (Unit: meter)
% NA:    system NA
% pupil: existing pupil function (aberrations) if pupil = 1, it is assumed
%        no aberrations 
%*******************Example parameters***************
% radius = 4e-6;
% delta_n = 0.05;
% wl = 632e-9;
% psize = 0.275e-6;
% m = 400; 
% n = 400;
% dx = 0; dy = 0;
% d = -10e-6;
% image = SimulatedBead(radius, delta_n, wl, psize, m, n, dx, dy);

[n2, m2] = size(image);% # of pixels. n2 in y direction, m2 in x direction
m = m2*2;
n = n2*2;
image_new = padarray(image,[n2/2,m2/2],0,'both'); % padd by extending the OUTREACH with inner value
% figure, imshow(image_new,[])
k0 = 2*pi/wl;
kmax = pi/(psize);
dkx = kmax/m;
dky = kmax/n;
[kx, ky] = meshgrid(linspace(-kmax,kmax-dkx,m),linspace(-kmax,kmax-dky,n));
mask = ones(n,m); % (kx.^2+ky.^2<=(k0*NA)^2).*
Fimage = fftshift(fft2(ifftshift(image_new)));


Amask = mask.*(exp(1i*real(sqrt(k0^2-kx.^2-ky.^2))*d));

% size(kx)
Fprop_image = Fimage.*Amask;
prop_image_temp = fftshift(ifft2(ifftshift(Fprop_image)));
prop_image = prop_image_temp(n2/2:n2*3/2-1,m2/2:m2*3/2-1);
% figure, imshow(abs(prop_image).^2,[])
end