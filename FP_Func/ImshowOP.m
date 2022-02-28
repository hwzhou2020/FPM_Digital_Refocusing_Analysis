function ImshowOP( O,P,disp )
IF = @(x) fftshift(ifft2(ifftshift(x)));
F = @(x) fftshift(fft2(ifftshift(x)));

        if ~strcmp(disp,'None')

            if strcmp(disp,'FS')
                o = O;
                f1 = figure(88);
                subplot(221); imagesc(log(abs(o))); axis image; colormap gray; colorbar;
                title('Fourier amplitude Log(O)');
                subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
                title('Fourier Phase(O)');
                subplot(223); imagesc((abs(P))); axis image; colormap gray; colorbar;
                title('ampl(P)');
                subplot(224); imagesc(angle(P)); axis image; colorbar;
                title('phase(P)');
                drawnow;  
            elseif strcmp(disp,'RS')
                o = IF(O);            
                f1 = figure(88);
                subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
                title('ampl(o)');
                subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
                title('phase(o)');
                subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
                title('ampl(P)');
                subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
                title('phase(P)');
                drawnow;

            end
        end

end

