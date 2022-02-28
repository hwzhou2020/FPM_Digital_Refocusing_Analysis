%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file to implement Fourier Ptychography reconstruction algorithm
% ref
% Lei Tian, et.al, Biomedical Optics Express 5, 2376-2389 (2014).
%
% last modified on 10/07/2015
% by Lei Tian, lei_tian@alum.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To do list for the user: (marked by 'TODO#')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) specify where the data located in 'filedir'
% 2) specify where you want to store the results in 'out_dir'
% 3) specify a threshold value, above which the estimated background value is signal rather than noise.
% 4) make sure the LED index (used for taking the images) are properly defined in 'lit'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Personal use notes Haowen Zhou
% 1) The sample ROI need to be a square region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
close all;
%% Z scan
zscan = 0:6;
Oz = zeros(1024,1024,length(zscan),4);

case_correct = 0;

for dfc = 1:2

for iz = 1:length(zscan)
%% Reconstruction library locates here
addpath(['./FP_Func']);
FilePath = 'E:\2021_FPM_Haowen\Data_08172021_238';
DataName = ['\08172021_Siemens_df',num2str(zscan(iz)+2),'_20x_g']; % +2 Siemens 
DataDir = [FilePath, DataName,'_preprocess.mat'];
filename = DataDir;
[filepath,name,ext] = fileparts(filename);
eval(['load ' filename]);

out_dir = ['Test'];
mkdir(out_dir);

%% Function operators
%Operator
F  = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) ifftshift(ifft2(fftshift(x)));
logamp = @(x) log10(abs(x)+1);

%% Set all necessary parameters (Unit: m)
if DataDir(end-15) == 'r'
    lambda = 632.3e-9;
elseif DataDir(end-15) == 'g'
    lambda = 522.6e-9;    
elseif DataDir(end-15) == 'b'
    lambda = 471.2e-9;
else
    disp('ERROR: Wavelength Undefined!');
    pause(100)
end

NA_ape = 0.40; % objective NA in FPMA
magnif = 20;
spsize = (4.6e-6)/magnif; % pixel size on sample plane
ds_led = 4e-3; % LED pitch
h = h*1e-6; % 5.4e-2% distance from LED to sample ## loaded from data
if h > 0.1 || h <0.01
    disp('error h !!!!!!!!!!!!!');
    pause(100);
end
USE_LED_IDX = 2;  % bright field only = 1; all led = 2
% physical offset of ROI
frame_c = [3126 4072];
% frame_c = [3175,3891]; % 6004*7920

LED_PC = 'none';%'cal_bf'; % 'none' no calibration; 'cal_bf' bright field correction;
                   % 'cal_nRO' no remove out. Method Laura group code.

EPRY = 1; % 0: not using; 1: use
if case_correct == 1
    defocus = -zscan(iz)*1e-6
    defocus_post_correct = 0;
else
    defocus = 0;
    defocus_post_correct = -zscan(iz)*1e-6
end


%% Prepare input for AlterMin
% prepare I
I = double(Isum);
[Nmy,Nmx,Nimg] = size(I);

Np = Nmy; % ASSUME Nmy = Nmx
% sampling size at Fourier plane set by the image size (FoV) = 1/FoV, FoV = Np*psize;
if mod(Np,2) == 1
    du = 1/(spsize*(Np-1));
else
    du = 1/(Np*spsize);
end

% prepare Ns (pixel index of subaperture center in f space)
% set up LED coordinates
Nled = Nimg;
ledc_x = x0; ledc_y = y0; % center LED position
led_idx = led_idx';
hhled = (led_idx(1,:)-ledc_x);
vvled = (led_idx(2,:)-ledc_y);
% physical offset of ROI
x = ROI_c(1)-frame_c(1);y = ROI_c(2)-frame_c(2);
img_center = [x,y]*spsize;
% corresponding angles for each LED
dd = sqrt((hhled*ds_led+img_center(1)).^2+(vvled*ds_led+img_center(2)).^2+h.^2);
sin_thetav = (hhled*ds_led+img_center(1))./dd;
sin_thetah = (vvled*ds_led+img_center(2))./dd;

% correct_k = 0; % only for NUScell
% if correct_k == 1
%     correct_k_NUScell();
%     sin_thetav = sin_thetav_corr;
%     sin_thetah = sin_thetah_corr;
% end


NAillu = sqrt(sin_thetav.^2+sin_thetah.^2);
vled = sin_thetav/lambda;
uled = sin_thetah/lambda;
freqUV = [vled;-uled];

%%
% spatial freq index for each plane wave relative to the center
switch LED_PC
    case 'none'
        idx_u = round(uled/du);
        idx_v = round(vled/du);
        Ns = [];
        Ns(:,:,1) = -idx_u;
        Ns(:,:,2) = idx_v; % [list of spatial freq values,m = 1:Nimg,v or h]
    case 'cal_bf'
        % LED Position Correction
        [uvled_cal,uvled_cal_nRO,NA_cal] = LEDPositionCorrection(I,freqUV',lambda,spsize*magnif,NA_ape,magnif);

        idx_u = round(-uvled_cal(:,2)'*1e6/du);
        idx_v = round(uvled_cal(:,1)'*1e6/du);
        Ns = [];
        Ns(:,:,1) = -idx_u;
        Ns(:,:,2) = idx_v; % [list of spatial freq values,m = 1:Nimg,v or h]
    case 'cal_nRO'
        % LED Position Correction
        [uvled_cal,uvled_cal_nRO,NA_cal] = LEDPositionCorrection(I,freqUV',lambda,spsize*magnif,NA_ape,magnif);

        idx_u = round(-uvled_cal_nRO(:,2)'*1e6/du);
        idx_v = round(uvled_cal_nRO(:,1)'*1e6/du);
        Ns = [];
        Ns(:,:,1) = idx_u;
        Ns(:,:,2) = -idx_v; % [list of spatial freq values,m = 1:Nimg,v or h]
end

% prepare No (size of reconstructed image)
NA_max = max(NAillu(:))+NA_ape; % maximum synthesized NA 
um_p = NA_max/lambda;
N_obj = round(2*um_p/du)*2;
N_obj = 2*Np;%ceil(N_obj/Np)*Np;
No = [N_obj,N_obj];

% generate w_NA -- CTF by subaperture NA & support of subaperture NA
% optional for this small section
m = 1:Np;
[mm,nn] = meshgrid(m-round((Np+1)/2));
ridx = sqrt(mm.^2+nn.^2);
um_m = NA_ape/lambda;
um_idx = um_m/du;
w_NA = double(ridx<um_idx);
rawimage_NA = double(ridx<(um_idx*2)); % lpf for low-res raw images

% generate H0
k0 = 2*pi/lambda;
kmax = pi/(spsize);
dkx = kmax/(Np/2);
dky = kmax/(Np/2);
[kx, ky] = meshgrid(linspace(-kmax,kmax-dkx,Np),linspace(-kmax,kmax-dky,Np));
mask = (kx.^2+ky.^2<=(k0*NA_ape)^2).*ones(Np,Np);

H0 = mask.*(exp(1i*real(sqrt(k0^2-kx.^2-ky.^2))*defocus));
% figure;
% subplot(121),imshow(mask,[]);
% subplot(122),imshow(angle(H0),[]);

s_NA = zeros(No); % updated during the iteration

%% initialize phase with the DPC result

%% reconstruction algorithm options: opts
opts.tol = 1e-15;
opts.maxIter = 50;
opts.minIter = 20;
opts.monotone = 1;

opts.display = '0';%'iter';
% 'full', display every subroutine
% 'iter', display only results from outer loop
% 0, no display
opts.mode = 'real'; % 'real'
opts.saveIterResult = 0;
opts.out_dir = out_dir;

upsamp = @(x) padarray(x,double([(N_obj-Np)/2,(N_obj-Np)/2]));
opts.O0 = F(sqrt(I(:,:,round(Nimg/2))));  % SOMETIMES NEED TO CHECK
opts.O0 = upsamp(opts.O0);
opts.P0 = w_NA;
opts.Os = s_NA; % support constraint for O0
opts.Ps = w_NA; % support constraint for P0
opts.AS = 1;    % Adaptive stepsize
opts.eta = 0.01; % threshold to stop adaptive stepsize
opts.OP_alpha = 0.05;% 0.1 Simens 0.01 BloodSmear
opts.OP_beta = 0.01; % 1 for Simens
opts.H0 = H0; % known portion of the aberration function

% LED intensity correction
opts.flag_inten_corr = 1;
opts.iter_inten_corr = 50; % when to start intensity correction

% LED position correction
opts.poscalibrate = '0'; % '0','sa','ga'
opts.calbratetol = 1e-1; % only for 'sa'
opts.iter_pos_corr = 25; % when to start intensity correction

% Use EPRY or not
opts.EPRY = EPRY;  % 0: not using; 1: use

opts.F = F;
opts.Ft = Ft;
opts.logamp = logamp;

%% algorithm starts
% lpf low-res raw images with 2NA_obj
lpf_rawimage = 1;
if lpf_rawimage == 1
    for i = 1:Nimg
        temp = I(:,:,i);
        I(:,:,i) = abs(Ft(F(temp).*rawimage_NA)).^2;
    end
end

%sort I and Ns according to the NAillu
[NAillu_reorder,order] = sort(NAillu);
I_reorder = I(:,:,order);
Ns_reorder = Ns(:,order,:);

%choose I and Ns to be used (for example, only Bright Field)
BF_only = 1;
if BF_only == USE_LED_IDX
    Nused = sum(squeeze(NAillu_reorder <= NA_ape));
else
    Nused = Nled;
end
idx_used = 1:Nused;
I_used = I_reorder(:,:,idx_used);
Ns_used = Ns_reorder(:,idx_used,:);

opts.scale = ones(length(idx_used),1); % LED brightness map, but never used

%show used low-res raw images before reconstruction
last_check = 0;
if last_check == 1
    for i = 1:Nused
        figure(13),imshow(I_used(:,:,i),[]);
        title([num2str(i) '-' num2str(NAillu_reorder(i))]);
        pause(0.1);
    end
end

%% real start
[O,P,err_pc,Ns_cal] = AlterMin(I_used,No,round(Ns_used),opts);
fprintf('processing completes\n');

if case_correct == 1
    Oz(:,:,iz) = O;
end
%% display results
I_kohler = sum(I_used,3);
% figure;imshow(I_kohler,[]);title('Köhler illumination');

% figure;
% subplot(221);
% imshow(abs(O), []);title('Sample amplitude');
% subplot(222);
% imshow(angle(O), []);title('Sample phase');
% subplot(223);
% imshow(abs(P), []);title('Probe amplitude');
% subplot(224);
% imshow(angle(P),[]);title('Probe phase');

% figure(1110+iz), imagesc(abs(O));
% axis off, axis square;
% colormap('gray'), colorbar;

% figure(111), imagesc(angle(P));
% axis off, axis square;
% colormap('jet'), colorbar;
% 
% figure;
% plot(1:length(err_pc),err_pc,'r-','LineWidth',1);
% set(gca,'YScale','log');
% axis square;xlabel('Iterations'),ylabel('Error');

%% show LED position before and after calibration
% figure;
% for k = 1:Nused
%     % [list of spatial freq values,m = 1:Nimg,v or h]
%     plot(Ns_used(:,k,1),Ns_used(:,k,2),'r*','MarkerSize',16);
%     hold on;
%     plot(Ns_cal(:,k,1),Ns_cal(:,k,2),'b*');
% end
% max(Ns_used(:)-Ns_cal(:))

%% digital refocusing
addpath(['./func_Autofocus']);

z = defocus_post_correct;
%% direct prop
prop_image = ASP_r(O, z, lambda, Np(1)*spsize/No(1), 1);
intensity  = abs(prop_image);
Odir_prop = prop_image;

% figure(116), imagesc(mat2gray(abs(FPM_prop)));
% axis off, axis normal;axis square;
% colormap('gray'), colorbar;

%% Syn prop
Spec_zdf = F(O);
H_syn = zeros(size(O));
H_ones = zeros(size(O));
kx0 = uled(:,order);
ky0 = vled(:,order);
idx_kx0 = round(kx0/du);
idx_ky0 = round(ky0/du);
for i = 1:Nused
    maski = ((kx).^2+(ky).^2<=(k0*NA_ape)^2).*ones(Np,Np);
    Hi = maski.*(exp(1i*real(sqrt(k0^2-kx.^2-ky.^2)) ...
        *z));
    Ho = maski;
    aperture_cx = round(size(O,1)/2) - idx_kx0(i);
    aperture_cy = round(size(O,2)/2) - idx_ky0(i);
    secsize = Nmx/2;
    secx = (aperture_cx - secsize) : (aperture_cx + secsize - 1);
    secy = (aperture_cy - secsize) : (aperture_cy + secsize - 1);
    H_syn(secx,secy) = H_syn(secx,secy) + Hi;
    H_ones(secx,secy) = H_ones(secx,secy) + Ho;
end
H_uni = H_ones./(H_ones + eps);
H_syn = H_syn./(H_ones + eps);

Spec_z = Spec_zdf.*H_syn;
Osyn_prop = Ft(Spec_z);

%%
if case_correct == 0
    Oz(:,:,iz,2) = Odir_prop;
    Oz(:,:,iz,3) = Osyn_prop;
    Oz(:,:,iz,4) = O;
end
end
case_correct = 1;
end
%% save results

if opts.EPRY == 0
    fn = [name,'_',num2str(Nused),'_zscan_',num2str(zscan(1)),'_',num2str(zscan(length(zscan))),'.mat'];
    save([out_dir,'\',fn],'Oz','lambda','zscan','h','NA_ape','magnif');
elseif opts.EPRY == 1
    fn = [name,'_',num2str(Nused),'_zscan_',num2str(zscan(1)),'_',num2str(zscan(length(zscan))),'_EPRY.mat'];
    save([out_dir,'\',fn],'Oz','lambda','zscan','h','NA_ape','magnif');
else
    disp('ERROR: opts.EPRY not correct. Values should be 0 or 1.');
end



% if opts.EPRY == 0
%     fn = [name,'_',num2str(Nused),'_drf',num2str(defocus*1e6),'_prop',num2str(defocus_post_correct*1e6),'.mat'];
%     save([out_dir,'\',fn],'Oz','Pz','err_pc','Ns_cal','FPM_prop','lambda','defocus','defocus_post_correct');
% elseif opts.EPRY == 1
%     fn = [name,'_',num2str(Nused),'_drf',num2str(defocus*1e6),'_prop',num2str(defocus_post_correct*1e6),'_EPRY','.mat'];
%     save([out_dir,'\',fn],'O','P','err_pc','Ns_cal','FPM_prop','lambda','defocus','defocus_post_correct');
% else
%     disp('ERROR: opts.EPRY not correct. Values should be 0 or 1.');
% end
