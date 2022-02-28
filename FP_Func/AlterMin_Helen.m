function [O, P, err, scale, cenLED] = AlterMin_Helen( imf, No, cenLED, opts)
%AlterMin Implements alternating minimization sequentially on a stack of
%measurement I (n1 x n2 x nz). It consists of 2 loop. The main loop update
%the reconstruction results O and P. the inner loop applies projectors/minimizers
%P1 and P2 on each image I and steps through the entire dataset
%   Outputs:
%   O: reconsturcted high-resolution complex object
%   P: reconstructed complex pupil function
%   err: errors at each iteration
%   scale: LED brightness map
%   Ns: estimated LED positions using local search algorithms
%
%   Inputs:
% Measurements data
%   I: intensity measurements by different LEDs
% Reconstruction parameters
%   No = [Ny_obj,Nx_obj]: size of the reconstructed image
% Illumination coding parameters
%   Ns = [Nsy,Nsx]: centers of corresponding lpf (low pass filter) regions for
%   the illumination pattern

% Iteration parameters: opts
%   tol: maximum change of error allowed in two consecutive iterations
    %   maxIter: maximum iterations 
    %   minIter: minimum iterations
    %   monotone (1, default): if monotone, error has to monotonically dropping
    %   when iters>minIter
%   display: display results (0: no (default) 1: yes)
    %   saveIterResult: save results at each step as images (0: no (default) 1: yes)
    %   out_dir: saving directory
%   O0, P0: initial guesses for O and P
    %   OP_alpha: regularization parameter for O
    %   OP_beta: regularization parameter for P
%   scale: LED brightness map
%   H0: known portion of the aberration function, 
        % e.g. sample with a known defocus induce a quardratic aberration
        % function can be defined here
%   poscalibrate: flag for LED position correction using
    % '0': no correction
    % 'sa': simulated annealing method
        % calbratetol: parameter in controlling error tolence in sa
    % 'ga': genetic algorithm
    % caution: takes consierably much longer time to compute a single iteration
%   F, Ft: operators of Fourier transform and inverse

% Last modified on 10/07/2017
% by Lei Tian, lei_tian@alum.mit.edu
% edited by Hangwen 08/01/2017
F = @(x) fftshift(fft2(ifftshift(x)));
IF = @(x) fftshift(ifft2(ifftshift(x)));

%% derived constants
% size of measurement
[ny,nx,Nimg] = size(imf);
Np = [ny,nx];
% M = No(1)/Np(1);
% r0 defines # of LEDs lit up in each pattern
[r0,~,~] = size(cenLED);
cen0 = round((No+1)/2); % center
row = @(x) x(:).';

%% options for the algorithms 
% if nargin<4
%     % default values (if opts are not input we need to define the opts)
%     opts.tol = 1;
%     opts.MaxIter = 50;
%     opts.MinIter = 3;
%     opts.disp = 'None';
%     opts.save = 'No';
%     opts.O0 = F((imf(:,:,1)))/r0;
%     opts.O0 = padarray(opts.O0,(No-Np)/2);
%     opts.P0 = ones(Np);
%     opts.alpha = 1;
%     opts.beta = 1;
%     opts.scale = ones(Nimg,1);
%     opts.aber = zeros(1,4);
%     opts.monotone = 1;
%     opts.F = @(x) fftshift(fft2(ifftshift(x)));
%     opts.IF = @(x) fftshift(ifft2(ifftshift(x)));
%     opts.epry = 0;
%     opts.H0 = ones(Np);
% % if opts exists, but some of its parameters are not defined, then define
% % them separately.
% else
    
    if ~isfield(opts,'tol')
        opts.tol = 1;
    end
    if ~isfield(opts,'MaxIter')
        opts.MaxIter = 50;
    end
    if ~isfield(opts,'MinIter')
        opts.MinIter = 3;
    end
    if ~isfield(opts,'disp')
        opts.disp = 'None';
    end
    if ~isfield(opts,'save')
        opts.save = 'No';
    end
    if ~isfield(opts,'O0')
        opts.O0 = F((imf(:,:,1)))/r0;
        opts.O0 = padarray(opts.O0,(No-Np)/2);
    end
    if ~isfield(opts,'P0')
        opts.P0 = ones(Np);
    end
    if ~isfield(opts,'alpha')
        opts.alpha = 1;
    end
    if ~isfield(opts,'beta')
        opts.beta = 1;
    end
    if ~isfield(opts,'Ps')
        opts.Ps = 1;
    end
    if ~isfield(opts,'iters')
        opts.iters = 10;
    end
    if ~isfield(opts,'scale')
        opts.scale = ones(led,1);
    end

    if ~isfield(opts,'aber')
        opts.aber = zeros(1,4);
    end
    if ~isfield(opts,'F')
        opts.F = @(x) fftshift(fft2(ifftshift(x)));
    end
    if ~isfield(opts,'IF')
        opts.IF = @(x) fftshift(ifft2(ifftshift(x)));
    end
    if ~isfield(opts,'epry')
        opts.epry = 0;
    end
    if ~isfield(opts,'monotone')
        opts.monotone = 1;
    end
    if ~isfield(opts,'NA')
        opts.NA = 1;
    end
    if ~isfield(opts,'wl')
        opts.wl = 1e-6;
    end
    if ~isfield(opts,'spsize')
        opts.spsize = 0.275e-6; % sensor pixel size
    end

% end
% % known aberration
aber = opts.aber;
if aber~= zeros(1,4)
    NA = opts.NA;
    spsize = opts.spsize;
    wl = opts.wl;
    opts.H0 = mask_aberration_4(ny,nx,NA,spsize,wl,aber,1);
else
    opts.H0 = opts.P0;
end

H0 = opts.H0; % known aberration
F = opts.F;
IF = opts.IF;

%% Operators
Ps = opts.P0; % mask in low resolution Fourier space

% operator to crop region of O from proper location at the O plane (for
% example, selecting a spacific tile in the original full frame image; or
% selecting a certian region in the Fourier domain during each updation)
downsamp = @(x,cen) x(cen(1)+round((No(1)+1)/2)-floor(ny/2):cen(1)+round((No(1)+1)/2)-floor(ny/2)+ny-1,...
    cen(2)+round((No(2)+1)/2)-floor(nx/2):cen(2)+round((No(2)+1)/2)-floor(nx/2)+nx-1);

T0 = clock;

fprintf('| iter |  rmse    |\n'); % rmse root mean square error
for j=1:20, fprintf('-'); end
fprintf('\n');



%% initialization in FT domain
% P = opts.P0;%opts.P0 = 0;
P = opts.P0.*H0;% Edited by Helen0808   
O = opts.O0; %opts.O0 = 0;
err1 = inf;
err2 = 100;
err = [];
iter = 0;
scale = opts.scale;
scale = reshape(scale,r0,Nimg);
ImshowOP(O,P,opts.disp)

if ~strcmp(opts.save,'No')
    export_fig(f1,[opts.save,'/R_',num2str(iter),'.png'],'-m4');
end


%% main algorithm starts here
% stopping criteria: when relative change in error falls below some value,
% can change this value to speed up the process by using a larger value but
% will trading off the reconstruction accuracy
% error is defined by the difference b/w the measurement and the estimated
% images in each iteration
fprintf('| %2d   | %.2e |\n',iter,err1);

sp0 = max(row(abs(cenLED(:,1,:)-cenLED(:,2,:)))); % used for position correction

while abs(err1-err2)>opts.tol && iter<opts.MaxIter 
    err1 = err2;
    err2 = 0;
    iter = iter+1 ;  
    for m = 1:Nimg % iteration over frame
        % initilize psi for correponing image, ROI determined by cen
        Psi0 = zeros(Np(1),Np(2),r0);
        Psi_scale = zeros(Np(1),Np(2),r0);
        cen = zeros(2,r0);
        scale0 = zeros(r0,1);
        for p = 1:r0 % iteration over led in a single frame
            cen(:,p) = row(cenLED(p,m,:));
            scale0(p) = scale(p,m);
            Psi0(:,:,p) = downsamp(O,cen(:,p)').*P; % H0 deleted by Helen 0808
%             Psi0(:,:,p) = downsamp(O,cen(:,p)).*P.*H0; 
            Psi_scale(:,:,p) = sqrt(scale0(p))*Psi0(:,:,p); % scaled (LED intensity)
        end
        % measured intensity
        I_mea = abs(imf(:,:,m)).^2;
        % compute field in real space
        psi0 = IF(Psi_scale);
        % estimated intensity w/o correction term
        I_est = sum(abs(psi0).^2,3);
        Psi = Proj_Fourier_v2(psi0, I_mea, I_est, scale0, F);
        % for non-multiplexing case, this is a simple update and then
        % Fourier transform
        % projection 2
        dPsi = Psi-Psi0; % this factor is added because the M-factor mismatch between enlarged image and captured image
        Omax = max(max(abs(downsamp(O,cen(:,1)'))));% maximum value in the subsection
        if r0 == 1
            P2 = @(O,P,dpsi,Omax,cen)...
                GDUpdate_Multiplication_rank1(O,P,dpsi,Omax,cen,Ps,...
                opts.alpha,opts.beta,opts.epry);
%         else
%             P2 = @(O,P,dpsi,Omax,cen)...
%                 GDUpdate_Multiplication_rank_r(O,P,dpsi,Omax,cen,Ps,...
%                 opts.OP_alpha,opts.OP_beta);
        end
%         [O,P] = P2(O,P,dPsi.*repmat(H0,[1,1,r0]),Omax,cen);
        [O,P] = P2(O,P,dPsi.*repmat(opts.P0,[1,1,r0]),Omax,cen); % edited by Helen 0808
        %%
%         ImshowOP(O,P,opts.disp)
    % compute the  error of current image frame
    err2 = err2+sqrt(sum(sum((I_mea-I_est).^2)));        
    end
    ImshowOP(O,P,opts.disp)
    % add to the error array with iteration
    err = [err,err2];
    
    fprintf('| %2d   | %.2e |\n',iter,err2);
    
    if opts.save ~= 'No'
        export_fig(f1,[opts.save,'\R_',num2str(iter),'.png'],'-m4');
        %     saveas(f2,[opts.out_dir,'\Ph_',num2str(iter),'.png']);
    end
    
    if opts.monotone && iter>opts.MinIter
        if err2>err1
            fprintf('not monotone decreasing')
            break;
        end
    end
%     err1
%     err2
% abs(err1-err2)> opts.tol   
%  iter<opts.MaxIter  
% judge = abs(err1-err2)>opts.tol && iter<opts.MaxIter     
end

%% show the image result
    o = IF(O);
    if opts.epry == 0
        figure(11), subplot(1,2,1); imshow(abs(o).^2,[])
        title(' the reconstructed sample intensity')
        subplot(1,2,2); imshow(angle(o),[])
        title(' the reconstructed sample phase')
    else
        figure(11), subplot(2,2,1); imshow(abs(o).^2,[])
        title(' the reconstructed sample intensity')
        subplot(2,2,2); imshow(angle(o),[])
        title(' the reconstructed sample phase')  
        subplot(2,2,3); imshow(abs(P),[])
        title(' recovered pupil function')
        subplot(2,2,4);imshow(angle(P),[])
        title(' recovered pupil function (phase)')
    end

fprintf('elapsed time: %.0f seconds\n',etime(clock,T0));

end

