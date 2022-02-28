function fmask = gOTFmask(im,NA,wl,psize)
% Low pass filter the image 'im'
[m1, n1]=size(im);factor=2;
mup=m1*factor;nup=n1*factor;fmask=ones(mup,nup);
NA=(32/32)*NA;
NAfilxup=NA*(1/wl)*n1*psize*factor; NAfilyup=NA*(1/wl)*m1*psize*factor; 
afil=fspecial('gaussian',round(2.5*NAfilxup/32),0.7*round(NAfilxup/32)); %fspecial('disk',round(NAfilxup/32));
[M1, N1]=meshgrid(1:mup,1:nup);fmask =fmask.*(((N1-(mup+1)/2)/NAfilyup).^2+((M1-(nup+1)/2)/NAfilxup).^2<=1);
fmask=imfilter(fmask,afil);
%x=-((nup-0.5)*wl)/(2*psize*nup):wl/(psize*nup):((nup-0.5)*wl)/(2*psize*nup);%plot(x,fmask(mup/2,:));
fmask=imresize(fmask,[m1 n1]);
end