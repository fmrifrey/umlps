function coeffs = chebyshev_affine(n, a, b)
    
    % base cases
    T0 = 1; % T_0(x) = 1
    T1 = [a b]; % T_1(x) = a*z + b
    
    if n == 0
        coeffs = T0;
        return
    elseif n == 1
        coeffs = T1;
        return
    end
    
    % recursive building
    for k = 2:n
        % multiply (a*z + b) * T_{k-1}
        T1_shift = conv([a b], T1);
        
        % 2*(a*z+b)*T_{k-1} - T_{k-2}
        Tn = 2*T1_shift - [zeros(1, length(T1_shift) - length(T0)) T0];
        
        % update T0, T1
        T0 = T1;
        T1 = Tn;
    end
    
    coeffs = Tn;
end