function T = specdiff_unfold3d_op(sz,nvol,lam)
% function to return finite differencing operator for fourier spectrum
% by David Frey
%
% inputs:
% sz - input 3D matrix size
% nvol - total number of volumes in timeseries
% lam - regularization parameter
%
% outputs:
% T - finite differencing linear operator in fourier dimension
%

    % define forward extraction operator
    function xFD = diff_spec_fwd(x)
    
        % reshape input and take temporal fft
        x = reshape(x,[],nvol);
        xF = fftshift(fft(x,[],2),2) / sqrt(nvol);
    
        % create spectral finite differencing operator
        xFD = sqrt(lam) * xF(:,2:end) - xF(:,1:end-1);
    
    end
    
    % define adjoint extraction operator
    function yDtFt = diff_spec_adj(y)
    
        % reshape output
        y = reshape(y,[],nvol-1);
    
        % adjoint finite differencing
        yDt = sqrt(lam) * [-y(:,1), y(:,1:end-1) - y(:,2:end), y(:,end)];
        
        % take inverse temporal fft
        yDtFt = ifft(fftshift(yDt,2),[],2) * sqrt(nvol);
    
    end

    % create the fatrix operator
    T = fatrix2('idim', [sz,nvol], ...
        'odim', [sz,nvol-1], ...
        'forw', @(~,x) diff_spec_fwd(x), ...
        'back', @(~,y) diff_spec_adj(y));

end