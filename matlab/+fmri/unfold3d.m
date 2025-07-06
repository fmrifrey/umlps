function [x_unfold,H] = unfold3d(x,vol_per_cyc,fwhm)
% function to perform UNFOLD filtering via notch filters at aliasing
% frequency
% by David Frey
%
% inputs:
% x - timeseries data (N x N x N x nvol)
% vol_per_cyc - number of volumes that get repeated in each cycle
% fwhm - full width half max of gaussian filter (odd number of freq.
%   components)
%
% outputs:
% x_unfold - clean data
% H - UNFOLD notch filter
%

    % set default notch filter width
    if nargin < 3 || isempty(fwhm)
        fwhm = 3;
    end

    % get length of timeseries
    nvol = size(x,4);

    % create notch filter
    if mod(nvol,2) % odd number of samples
        H = freq_filt(-(nvol-1)/2:(nvol-1)/2, fwhm, nvol/vol_per_cyc);
    else
        H = freq_filt(-nvol/2:nvol/2-1, fwhm, nvol/vol_per_cyc);
    end
    
    % apply filter
    X = fftshift(fft(x,[],4),4);
    x_unfold = ifft(fftshift(X.*permute(H(:),[2:4,1]),4),[],4);

end

function Hf = freq_filt(f, fwhm, f_alias)

    % create gaussian attenuation filter at aliasing frequencies
    g = @(x) exp(-(pi*x).^2);
    Hf = 1-g((mod(f - floor(f_alias/2),f_alias) - ceil(f_alias/2))/fwhm);
    
    % DC term - no attenuation
    Hf(abs(f) <= f_alias/2) = 1;

end