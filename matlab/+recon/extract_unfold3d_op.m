function T = extract_unfold3d_op(sz,nvol,vol_per_cyc,lam,fwhm)
% function to extract spectral peaks of UNFOLD aliasing artifacts
% by David Frey
%
% inputs:
% sz - input 3D matrix size
% nvol - total number of volumes in timeseries
% vol_per_cyc - number of volumes that get repeated in each cycle
% lam - regularization parameter
% fwhm - full width half max of gaussian filter (odd number of freq.
%   components)
%
% outputs:
% T - linear operator which extracts the spectral peak signals
%

    % set default notch filter width
    if nargin < 5 || isempty(fwhm)
        fwhm = 3;
    end

    % create notch filter
    if mod(nvol,2) % odd number of samples
        H = freq_filt(-(nvol-1)/2:(nvol-1)/2, fwhm, nvol/vol_per_cyc);
    else
        H = freq_filt(-nvol/2:nvol/2-1, fwhm, nvol/vol_per_cyc);
    end

    % define forward extraction operator
    function xFH = extract_spec_fwd(x)
    
        % reshape input and take temporal fft
        x = reshape(x,[],nvol);
        xF = fftshift(fft(x,[],2),2) / sqrt(nvol);
    
        % filter xF
        xFH = sqrt(lam) * xF .* H(:).';
    
    end
    
    % define adjoint extraction operator
    function yHtFt = extract_spec_adj(y)
    
        % reshape output
        y = reshape(y,[],nvol);
    
        % create spectral subsampling operator
        yHt = sqrt(lam) * y .* H(:)';
        
        % take inverse temporal fft
        yHtFt = ifft(fftshift(yHt,2),[],2) * sqrt(nvol);
    
    end

    % create the fatrix operator
    T = fatrix2('idim', [sz,nvol], ...
        'odim', [sz,nvol], ...
        'forw', @(~,x) extract_spec_fwd(x), ...
        'back', @(~,y) extract_spec_adj(y));

end

function Hf = freq_filt(f, fwhm, f_alias)

    % create gaussian attenuation filter at aliasing frequencies
    g = @(x) exp(-(pi*x).^2);
    Hf = g((mod(f - floor(f_alias/2),f_alias) - ceil(f_alias/2))/fwhm);
    
    % DC term - no attenuation
    Hf(abs(f) <= f_alias/2) = 0;

end