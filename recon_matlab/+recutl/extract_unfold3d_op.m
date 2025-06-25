function T = extract_unfold3d_op(sz,nvol,vol_per_cyc,lam)
% function to extract spectral peaks of UNFOLD aliasing artifacts
% by David Frey
%
% inputs:
% sz - input 3D matrix size
% nvol - total number of volumes in timeseries
% vol_per_cyc - number of volumes that get repeated in each cycle
% lam - regularization parameter
%
% outputs:
% T - linear operator which extracts the spectral peak signals
%

    % define forward extraction operator
    function xFPsi = extract_spec_fwd(x)
    
        % reshape input and take temporal fft
        x = reshape(x,[],nvol);
        xF = fftshift(fft(x,[],2),2) / sqrt(nvol);
    
        % create spectral subsampling operator
        k = 1:nvol/vol_per_cyc:nvol;
        xFPsi = zeros(size(xF));
        xFPsi(:,k) = sqrt(lam) * xF(:,k);
        xFPsi(:,nvol/2+1) = 0; % ignore DC component
    
    end
    
    % define adjoint extraction operator
    function yPsitFt = extract_spec_adj(y)
    
        % reshape output
        y = reshape(y,[],nvol);
    
        % create spectral subsampling operator
        k = 1:nvol/vol_per_cyc:nvol;
        yPsit = zeros(size(y));
        yPsit(:,k) = sqrt(lam) * y(:,k);
        yPsit(:,nvol/2+1) = 0;
        
        % take inverse temporal fft
        yPsitFt = ifft(fftshift(yPsit,2),[],2) * sqrt(nvol);
    
    end

    % create the fatrix operator
    T = fatrix2('idim', [sz,nvol], ...
        'odim', [sz,nvol], ...
        'forw', @(~,x) extract_spec_fwd(x), ...
        'back', @(~,y) extract_spec_adj(y));

end