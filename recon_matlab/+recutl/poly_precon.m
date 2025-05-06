function P = poly_precon(A,d,AHA_maxeig,AHA_mineig)

    % get coefficients of the polynomial
    r_coeffs_numer = recutl.chebyshev_affine(d+1, -2/(AHA_maxeig - AHA_mineig), ...
        (AHA_maxeig + AHA_mineig)/(AHA_maxeig - AHA_mineig));
    r_coeffs_denom = polyval(recutl.chebyshev_affine(d+1, 0, ...
        (AHA_maxeig + AHA_mineig)/(AHA_maxeig - AHA_mineig)), 1);
    pz_coeffs = -r_coeffs_numer./r_coeffs_denom;
    p_coeffs = pz_coeffs(1:end-1); % 1 - last coeff = 0!

    % define forward operation of polynomial p(A^HA)
    function Px = poly_precon_fwd(c,A,x)
        if isscalar(c) % if one coefficient, just constant
            Px = c(end)*x;
        else % recurse with higher order coefficients
            Px = c(end)*x + reshape(A'*(A*poly_precon_fwd(c(1:end-1),A,x)), A.idim);
        end
    end
    
    % form preconditioner forward operator
    P = fatrix2('idim', A.idim, ...
        'odim', A.idim, ...
        'forw', @(~,x) poly_precon_fwd(p_coeffs,A,x) ...
        );

end
