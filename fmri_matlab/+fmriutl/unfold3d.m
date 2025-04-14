function [x_unfold,H] = unfold3d(x,vol_per_cyc)
% function to perform UNFOLD filtering via notch filters at aliasing
% frequency
% by David Frey
%
% inputs:
% x - timeseries data (N x N x N x nvol)
% vol_per_cyc - number of volumes that get repeated in each cycle
%
% outputs:
% x_unfold - clean data
% H - UNFOLD notch filter
%

    % get length of timeseries
    nvol = size(x,4);

    % make notch filter
    H = ones(nvol,1);
    H(1:nvol/vol_per_cyc:end) = 0; % null the aliasing freqs
    H(nvol/2+1) = 1; % restore DC
    
    % apply filter
    X = fftshift(fft(x,[],4),4);
    x_unfold = ifft(fftshift(X.*permute(H(:),[2:4,1]),4),[],4);

end

