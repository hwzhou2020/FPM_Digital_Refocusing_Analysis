%%%%% Illumination Calibration Script %%%%%
function [freqUV_cal,freqUV_noRemoveOut,NA_cal] = LEDPositionCorrection(I,freqUV,lambda,pixel_pitch,NA,mag)
addpath('Utilities');

%%
% *** Input Values *** %
% Data     I
% lambda
% pixel_pitch 6.5
% NA = 0.4;
% mag = 20;

% *** Output Values *** %
% na_calib = freqUV_cal.*lambda;
% na_nRO = freqUV_noRemoveOut.*lambda;
% na_cal = NA_cal;

% *** Imaging system parameters *** %
% Define:
%   I: image stack, NxNxM (for M images)
%   NA: objective NA
%   mag: total magnification
%   dpix_c: pixel size (um)
%   lambda: imaging wavelength (um)
%   freqUV: expected spatial frequency of illumination (1/um), Mx2, ('x' in
%           column 1, 'y' in column 2)


pixel_pitch = pixel_pitch*1e6;
lambda = lambda*1e6;
freqUV = freqUV/1e6;

warning_count = 1;
warning_list = [];
for z_defocus = 0 % unit: um

    % *** Controls for removing outliers *** %
    removeOut = 1; %Whether you want to remove outliers and extend to darkfield (generally a good idea)
    method='rigidScale'; %Option for transformation from original guess to new
    %guess, used to remove outliers that do not fit that transformation.
    %Possibilities are: 'rigid' (rotation and shift),'rigidScale' (rotation, shift, & scale),
    %'affine'(rotation, shift, scale, & shear), 'projective'
    alpha=3; %Multiplier by stdev to get outliers (usually 2-3)
    scale=0.75; %Value between 0.01 and 0.99; what the weight of suspected outliers is multiplied by at each iteration
    tol=0.05; %Level where we consider something an outlier (should be low, ~0.05)

    % *** Controls for circle-finding *** %
    rScan=[5 0.5]; %[+/- to scan, delta(r)] %Radius scan boundaries and granularity
    %Row #1 = radius estimation, row #2 = center estimation
    thScan=[10 1; 5  0.5]; %Theta scan boundaries and granularity
    dScan=[20 1; 5 0.5]; %Distance scan boundaries and granularity
    calRad=1; %Whether or not to calibrate the radius (recommended to do so)

    %% Define system
    %Set up space
    imSz=size(I);
    N=imSz(1:2);
    numImg=imSz(3);
    imSz(3)=[];

    tic;

    %Define the coordinates in frequency space in terms of pixels and k-space
    %(u,v); convert pixel coord to polar (centD, theta); get the conversion
    %factor from k-space to pixels (con)
    [freqXY, con, radP, xI, yI, uI, vI, XYmid ] = calCoord(freqUV,imSz,pixel_pitch,mag, NA, lambda );

    %Convert to distance (pixels), theta (degrees) of circle center from center frequency
    freqDTh = cart2Pol(freqXY, XYmid);

    %Get Fourier space amplitude images that will be processed
    sigmaG=3; %Amount of blurring to use in a Gaussian filter (generally = 2 gives good results)
    [FIdiv, FIdivG, FI, w_2NA] = calFI(I, xI, yI, XYmid, radP, sigmaG);

    %% Identify darkfield
%     [ DFI ] = calDF( FI, XYmid ); %%%%%%%%%%%%

    NA_illum = sqrt((freqUV(:,1)*lambda).^2 + (freqUV(:,2)*lambda).^2);
    [IDX] = find(NA_illum > NA);
    DFI = zeros(length(freqUV(:,1)),1);
    DFI(IDX) = 1;
    
    %Print number of brightfield and darkfield images
    fprintf('%i BF, %i DF\n',sum(~DFI),sum(DFI))

    %% Find circles
    %Find circles for the brightfield images only
    [freqDTh3, rad_cal] = calCircEdge(FIdivG(:,:,~DFI), I(:,:,~DFI), radP, freqDTh(~DFI,:), XYmid, xI, yI, sigmaG, rScan, thScan, dScan, calRad, con, lambda);

    freqDTh2=freqDTh;
    freqDTh2(~DFI,:)=freqDTh3; %Replace the brightfield values only

    freqXY_noRemoveOut=pol2Cart(freqDTh2,XYmid); %Save the result before removing outliers to compare

    if removeOut
        [ freqDTh2 ] = removeOutliers(freqDTh, freqDTh2, XYmid, method,alpha, scale, tol, DFI);
        %When include DFI, it doesn't use the darkfield to fit but it does
        %replaces darkfield with fit values (if you don't include the DFI
        %variable, it does not extend to darkfield)
    end

    freqXY2=pol2Cart(freqDTh2,XYmid);

    %% Format output & save
    %Convert back to k-space values
    freqUV_cal=(freqXY2-repmat(XYmid,[numImg 1]))./con;
    freqUV_noRemoveOut=(freqXY_noRemoveOut-repmat(XYmid,[numImg 1]))./con;
    NA_cal=(rad_cal./con).*lambda;

    t_cal=toc
    fprintf('Time elapsed: %.3f seconds\n',t_cal)

    % Return values
    na_design = freqUV.*lambda;
    na_calib = freqUV_cal.*lambda;
    na_nRO = freqUV_noRemoveOut.*lambda;
    na_cal = NA_cal;

    %Check error
%     na_err = sqrt((freqUV(:,1)-freqUV_cal(:,1)).^2+(freqUV(:,2)-freqUV_cal(:,2)).^2);
%     if max(na_err)>(NA/2)
%         warning_list{warning_count} = outCaliFile;
%         warning_count = warning_count+1;
%     end

    %% Display results
%     figure();
%     plot(freqUV(:,1),freqUV(:,2),'*',freqUV_cal(:,1),freqUV_cal(:,2),'*');
%     title('Uncalibrated and calibrated angles, k-space');
%     legend('Original','SelfCal');
%     axis square;
% 
%     sliderDisplayImVC2(log10(abs(FI(:,:,~DFI))+1), cat(3,[freqXY(~DFI,:) radP.*ones(sum(~DFI),1)],[freqXY2(~DFI,:) rad_cal.*ones(sum(~DFI),1)]),{'title(''Brightfield Only: (red) Uncalibrated, (green) Calibrated'')'});
%     sliderDisplayImVC2(log10(abs(FI)+1), cat(3,[freqXY radP.*ones(numImg,1)],[freqXY2 rad_cal.*ones(numImg,1)]),{'title(''All Images: (red) Uncalibrated, (green) Calibrated'')'});

end