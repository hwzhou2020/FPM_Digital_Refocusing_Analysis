function imfil = filim(im, NA, wlength, spsize)
%   filim computes the low res seq,  im is squre image
%   x is along horizonal direction; y is along vertical direction
[m n]=size(im);[M N]=meshgrid(1:m,1:n);
NAfilx=NA*(1/wlength)*n*spsize;NAfily=NA*(1/wlength)*m*spsize;
fmask2=(((M-(n+1)/2)/(2*NAfily)).^2+((N-(m+1)/2)/(2*NAfilx)).^2<=1);
imFTtemp=fftshift(fft2(im));
imfil=real(ifft2(ifftshift(imFTtemp.*fmask2)));
end