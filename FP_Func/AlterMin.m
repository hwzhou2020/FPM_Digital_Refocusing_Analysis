function [O, P, err, Ns] = AlterMin(I, No, Ns, opts)
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
%   Ns = [Nsy,Nsx]: centers of corresponding lpf regions for
%   the illumination pattern

% Iteration parameters: opts
%   tol: maximum change of error allowed in two consecutive iterations
    %   maxIter: maximum iterations 
    %   minIter: minimum iterations
    %   monotone (1, default): if monotone, error has to monotonically dropping when iters>minIter
%   display: display results (0: no (default) 1: yes)
    %   saveIterResult: save results at each step as images (0: no (default) 1: yes)
    %   mode: display in 'real' space or 'fourier' space.
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
    % caution: takes considerably much longer time to compute a single iteration
%   Use EPRY
    % 1: use
    % 0: not use
%   F, Ft: operators of Fourier transform and inverse

% Last modified on 10/07/2017
% by Lei Tian, lei_tian@alum.mit.edu

%% derived constants
% size of measurement
[Nmy,Nmx,Nimg] = size(I);
Np = [Nmy,Nmx];
MAGimg = No(1)/Np(1);
[r0,~,~] = size(Ns); % r0 defines # of LEDs lit up in each pattern
cen0 = round((No+1)/2);
row = @(x) x(:).';

%% options for the algorithms
if nargin<4
    % default values
    opts.tol = 1;
    opts.maxIter = 50;
    opts.minIter = 3;
    opts.monotone = 1;
    opts.display = 0;
    opts.saveIterResult = 0;
    opts.out_dir = [];
    opts.O0 = Ft(sqrt(I(:,:,1)))/r0;
    opts.O0 = padarray(opts.O0,(No-Np)/2);
    opts.P0 = ones(Np);
    opts.OP_alpha = 1;
    opts.OP_beta = 1;
    opts.mode = 'real';
    opts.scale = ones(Nled,1);
    opts.H0 = ones(Np);
    opts.poscalibrate = 0;
    opts.calbratetol = 1e-1;
    opts.F = @(x) fftshift(fft2(x));
    opts.Ft = @(x) ifft2(ifftshift(x));
    opts.EPRY = 0;
else
    if ~isfield(opts,'tol')
        opts.tol = 1;
    end
    if ~isfield(opts,'maxIter')
        opts.maxIter = 50;
    end
    if ~isfield(opts,'minIter')
        opts.minIter = 3;
    end
    if ~isfield(opts,'monotone')
        opts.monotone = 1;
    end
    if ~isfield(opts,'display')
        opts.display = 0;
    end
    if ~isfield(opts,'saveIterResult')
        opts.saveIterResult = 0;
    end
    if ~isfield(opts,'out_dir')
        opts.out_dir = ['IterResults'];
        if opts.saveIterResult
            mkdir(opts.out_dir);
        end
    end
    if ~isfield(opts,'O0')
        opts.O0 = Ft(sqrt(I(:,:,1)))/r0;
        opts.O0 = padarray(opts.O0,(No-Np)/2);
    end
    if ~isfield(opts,'P0')
        opts.P0 = ones(Np);
    end
    if ~isfield(opts,'OP_alpha')
        opts.OP_alpha = 1;
    end
    if ~isfield(opts,'OP_beta')
        opts.OP_beta = 1;
    end
    if ~isfield(opts,'mode')
        opts.mode = 'real';
    end
    if ~isfield(opts,'Ps')
        opts.Ps = 1;
    end
    if ~isfield(opts,'scale')
        opts.scale = ones(Nled,1);
    end

    if ~isfield(opts,'H0')
        opts.H0 = ones(Np);
    end
    if ~isfield(opts,'poscalibrate')
        opts.poscalibrate = 0;
    end
    if ~isfield(opts,'F')
        opts.F = @(x) fftshift(fft2(ifftshift(x)));
    end
    if ~isfield(opts,'Ft')
        opts.Ft = @(x) fftshift(ifft2(ifftshift(x)));
    end
    if ~isfield(opts,'calbratetol')
        opts.calbratetol = 1e-1;
    end
    if ~isfield(opts,'EPRY')
        opts.EPRY = 0;
    end
end

H0 = opts.H0;
F = opts.F;
Ft = opts.Ft;
logamp = opts.logamp;
EPRY = opts.EPRY;

%% Operators
Ps = opts.Ps;
Os = opts.Os;
% operator to crop region of O from proper location at the O plane
downsamp = @(x,cen) x((cen(1)-floor(Np(1)/2)):(cen(1)-floor(Np(1)/2)+Np(1)-1),...
                      (cen(2)-floor(Np(2)/2)):(cen(2)-floor(Np(2)/2)+Np(2)-1));

T0 = clock;

fprintf('| iter |  rmse    |\n');
for j=1:20, fprintf('-'); end
fprintf('\n');

%% initialization in FT domain
P = opts.P0; opts.P0 = 0;
O = opts.O0; opts.O0 = 0;
err1 = inf;
err2 = 50;
err = [];
iter = 0;
flag_inten_corr = opts.flag_inten_corr;
iter_inten_corr = opts.iter_inten_corr;
iter_pos_corr = opts.iter_pos_corr;
AS = opts.AS;
eta = opts.eta;

if opts.display
    f1 = figure(88);
    if strcmp(opts.mode,'fourier')
        subplot(221); imagesc(logamp(O));axis image; colormap gray; colorbar;
        title('ampl(O)');
        subplot(222); imagesc(angle(O)); axis image; colormap gray; colorbar;
        title('phase(O)');
    elseif strcmp(opts.mode,'real')
        o = Ft(O);
        subplot(221); imagesc(abs(o));axis image; colormap gray; colorbar;
        title('ampl(o)');
        subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
        title('phase(o)');
    end
    subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
    title('ampl(P)');
    subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
    title('phase(P)');
    drawnow;
end

if opts.saveIterResult
    export_fig(f1,[opts.out_dir,'\R_',num2str(iter),'.png'],'-m4');
end

%% main algorithm starts here
% stopping criteria: when relative change in error falls below some value,
% can change this value to speed up the process by using a larger value but
% will trading off the reconstruction accuracy
% error is defined by the difference b/w the measurement and the estimated
% images in each iteration
fprintf('| %2d   | %.2e |\n',iter,err1);

sp0 = max(row(abs(Ns(:,1,:)-Ns(:,2,:))));
tt = zeros(1,Nimg);

while abs(err1-err2)>opts.tol&&iter<opts.maxIter
    err1 = err2;
    err2 = 0;
    iter = iter+1;
    Os1 = Os;
    for m = 1:Nimg
        % initilize psi for correponing image, ROI determined by cen
        Psi0 = zeros(Np(1),Np(2),r0);
        cen = zeros(2,r0);
        for p = 1:r0
            cen(:,p) = cen0-row(Ns(p,m,:));
            Psi0(:,:,p) = downsamp(O,cen(:,p)).*P.*H0;%./(MAGimg^2)
        end
        % compute field in real space
        psi0 = Ft(Psi0);
        % estimated intensity w/o correction term
        I_est = sum(abs(psi0).^2,3);
        % measured intensity
        I_mea = I(:,:,m);
        
        %% raw image intensity correction
        if (flag_inten_corr == 1) && (iter >= iter_inten_corr)
            tt(1,m)=(mean(mean(abs(I_est)))/mean(mean(abs(I_mea))));  
            I(:,:,m)=I(:,:,m)*tt(1,m);
        end
        
        % projection 1: amplitude constraint in real space
        Psi = Proj_Fourier_v2(psi0, I_mea, I_est, 1, F);
        % projection 2: support constraint in frequency space & GD update
        dPsi = (Psi-Psi0);%.*(MAGimg^2);
        Omax = abs(O(cen0(1),cen0(2)));
        if r0 == 1
            P2 = @(O,P,dpsi,Omax,cen)...
                GDUpdate_Multiplication_rank1(O,P,dpsi,Omax,cen,Ps,...
                opts.OP_alpha,opts.OP_beta,EPRY);
        else
            P2 = @(O,P,dpsi,Omax,cen)...
                GDUpdate_Multiplication_rank_r(O,P,dpsi,Omax,cen,Ps,...
                opts.OP_alpha,opts.OP_beta);
        end
        dPsi_temp = dPsi./repmat(H0,[1,1,r0]);
        dPsi_temp(isnan(dPsi_temp)|isinf(dPsi_temp)) = 0;
        [O,P] = P2(O,P,dPsi_temp,Omax,cen);

        %% position correction
        if iter >= iter_pos_corr
            poscost = @(ss) sum(sum((abs(Ft(downsamp(O,ss).*P.*H0)).^2-I_mea).^2));
            if strcmp(opts.poscalibrate,'sa')
                optsanneal = saoptimset('Display','off','TolFun',opts.calbratetol);
                cen_correct = round(simulannealbnd(poscost,...
                    cen(:,1),cen(:,1)-sp0/3,cen(:,1)+sp0/3,optsanneal));
                Ns(:,m,:) = cen0-cen_correct';
            elseif strcmp(opts.poscalibrate,'ga')
                optsanneal = saoptimset('Display','off');
                cen_correct = ga(poscost,2,[],[],[],[],...
                    cen(:,1)-sp0/3,cen(:,1)+sp0/3,[],[1,2],optsanneal);
                Ns(:,m,:) = cen0-cen_correct;
            end
        end
        
        % compute the total difference to determine stopping criterion
        err2 = err2+sqrt(sum(sum((I_mea-I_est).^2)));

        Np1 = Np(:);
        n1 = cen-floor(Np1/2);
        n2 = n1+Np1-1;
        Os1(n1(1):n2(1),n1(2):n2(2)) = Ps+Os1(n1(1):n2(1),n1(2):n2(2)).*(ones(size(Ps))-Ps);
 
        if strcmp(opts.display,'full')
            f1 = figure(88);
            if strcmp(opts.mode,'fourier')
                subplot(221); imagesc(logamp(O));axis image; colormap gray; colorbar;
                title('ampl(O)');
                subplot(222); imagesc(angle(O)); axis image; colormap gray; colorbar;
                title('phase(O)');
            elseif strcmp(opts.mode,'real')
                o = Ft(O);
                subplot(221); imagesc(abs(o));axis image; colormap gray; colorbar;
                title('ampl(o)');
                subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
                title('phase(o)');
            end
            subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
            title(['ampl(P) ite = ' num2str(iter)]);
            subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
            title(['phase(P) raw img = ' num2str(m)]);
            drawnow;
        end

    end
   
    O = O.*Os1;

    %% compute error
    % record the error and can check the convergence later.
    err = [err,err2];
    
    if strcmp(opts.display,'iter')
        f1 = figure(88);
        if strcmp(opts.mode,'fourier')
            subplot(221); imagesc(logamp(O));axis image; colormap gray; colorbar;
            title('ampl(O)');
            subplot(222); imagesc(angle(O)); axis image; colormap gray; colorbar;
            title('phase(O)');
        elseif strcmp(opts.mode,'real')
            o = Ft(O);
            subplot(221); imagesc(abs(o));axis image; colormap gray; colorbar;
            title('ampl(o)');
            subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
            title('phase(o)');
        end
        subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
        title(['ampl(P) ite = ' num2str(iter)]);
        subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
        title(['phase(P) raw img = ' num2str(m)]);
        drawnow;
    end
    
    fprintf('| %2d   | %.2e |\n',iter,err2);
    
    % adaptive stepsize
    if  AS==1 && iter>1 && (err1-err2)/err1<0.005 
        % Reduce the stepsize when no sufficient progress is made
        opts.OP_alpha = opts.OP_alpha/2;
        opts.OP_beta = opts.OP_beta/2;
        % Stop the iteration when alpha is less than eta(convergenced)
        if (opts.OP_alpha < eta)
            opts.OP_alpha = 0;
        end
    end
    fprintf('\n');
    fprintf('alpha = %4.2f',opts.OP_alpha);
    fprintf('\n');

    if opts.saveIterResult
        export_fig(f1,[opts.out_dir,'\R_',num2str(iter),'.png'],'-m4');
        %     saveas(f2,[opts.out_dir,'\Ph_',num2str(iter),'.png']);
    end
    
    if opts.monotone&&iter>opts.minIter
        if err2>err1
            O = O_bef; % use the result from last iteration
            break;
        end
    end
    O_bef = O;
    
    if opts.OP_alpha < 0.01
        break;
    end
end

if strcmp(opts.mode,'real')
    O = Ft(O); % always output the reconstruction in real space
end


fprintf('elapsed time: %.0f seconds\n',etime(clock,T0));

end