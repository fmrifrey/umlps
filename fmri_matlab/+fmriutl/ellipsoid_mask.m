function msk = ellipsoid_mask(sz,off,sma)
% function to create an ellipsoid shaped mask
% by David Frey
%
% inputs:
% sz - size of 3d image [x,y,z]
% off - center offset (normalized units) [x,y,z]
% sma - semi-major axis (normalized units) [x,y,z]
%
% outputs:
% msk - boolean mask array
%
    [X,Y,Z] = ndgrid(linspace(-1,1,sz(1)), ...
        linspace(-1,1,sz(2)), ...
        linspace(-1,1,sz(3)));
    msk = ((X - off(1))/sma(1)).^2 + ...
        ((Y - off(2))/sma(2)).^2 + ...
        ((Z - off(3))/sma(3)).^2 <= 1;
end

