%% edited by helen, adding spherical aberration on April 27
function [fmask_ab] = mask_aberration_4(nx,ny,NA,spsize,wl,para, ShowPupil)
%function [fmask,F_pupil_ab] = mask_aberration(nx,ny,NA,spsize,wl,defo,ax,ay,comax,comay, sph,tiltx,tilty, ShowPupil)
%generating an aberrated mask F_pupil_ab, and a perfect low-pass filter
%F_mask
%inputs:
% nx,ny: mask size(imagesize), x=1:nx;y=1:ny
% NA pupil size
% spsize: the pixel size of image
% wl: wavelength
% including 
% edited by Helen
defo = para(1)*1e-6;
ax = para(2);
ay = para(3);
sph = para(4);
kmax=pi/spsize;
ky2=-kmax:kmax/((ny-1)/2):kmax;
kx2=-kmax:kmax/((nx-1)/2):kmax; %odd nx
[kxm, kym]=meshgrid(kx2,ky2);
k0 = 2*pi/wl;
kzm=sqrt(k0^2-kxm.^2-kym.^2);
H2=exp(1j.*defo.*real(kzm)).*exp(-abs(defo).*abs(imag(kzm)));
fmask=double(((kxm)).^2+(kym).^2<=(2*pi*NA/wl)^2);%old method, perfect reconstruction
NAfilx=round(NA/wl*spsize*ny);
NAfily=round(NA/wl*spsize*nx);
zn=0*gzn(max(ny,nx),2*max(round(NAfily),round(NAfilx)),0,2)+...
    ax*gzn(max(ny,nx),2*max(round(NAfily),round(NAfilx)),2,2)+...
    ay*gzn(max(ny,nx),2*max(round(NAfily),round(NAfilx)),-2,2)+...
    sph*gzn(max(ny,nx),2*max(round(NAfily),round(NAfilx)),0,4);
 
zn=imresize(zn,[ny,nx]);
fmask_ab=fmask.*exp(pi*1j.*zn).*H2;
if ShowPupil==1
figure, imshow(angle(fmask_ab),[]);
end

end

