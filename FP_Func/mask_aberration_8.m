%% edited by helen, adding spherical aberration on April 27
function [fmask,fmask_ab] = mask_aberration(n,m,NA,psize,wl,defo,ax,ay,comax,comay, sph,tiltx,tilty, ShowPupil)
%function [fmask,F_pupil_ab] = mask_aberration(n,m,NA,psize,wl,defo,ax,ay,comax,comay, sph,tiltx,tilty, ShowPupil)
%generating an aberrated mask F_pupil_ab, and a perfect low-pass filter
%F_mask
%inputs:
% n,m: mask size(imagesize), x=1:n;y=1:m
% NA pupil size
% psize: the pixel size of image
% wl: wavelength
% including 
% edited by Helen
kmax=pi/psize;
kx2=-kmax:kmax/((n-1)/2):kmax;
ky2=-kmax:kmax/((m-1)/2):kmax; %odd N
[kxm, kym]=meshgrid(kx2,ky2);
fmask=double(((kxm)).^2+(kym).^2<=(2*pi*NA/wl)^2);%old method, perfect reconstruction
NAfilx=round(NA/wl*psize*n);
NAfily=round(NA/wl*psize*m);
zn=defo*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),0,2)+...
    ax*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),2,2)+...
    ay*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),-2,2)+...
    comax*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),1,3)+...
    comay*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),-1,3)+...
    sph*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),0,4)+...
    tilty*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),-1,1)+...
    tiltx*gzn(max(m,n),2*max(round(NAfily),round(NAfilx)),1,1);

zn=imresize(zn,[m,n]);
fmask_ab=fmask.*exp(pi*1j.*zn);
if ShowPupil==1
figure, imshow(angle(fmask_ab),[]);colormap jet
end

end

