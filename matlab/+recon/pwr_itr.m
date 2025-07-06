function [v1,lam1] = pwr_itr(A,niter,tol)

    % initialize
    if isa(A,'fatrix2')
        v1 = rand(A.idim);
    else
        v1 = rand(size(A,2),1);
    end
    Av1 = A*v1;
    lam1 = v1(:)' * Av1(:);
    v1_kp1 = Av1/norm(Av1(:));
    err = norm(v1_kp1(:) - v1(:));

    % loop through iterations
    for k = 1:niter
        fprintf('power iter %d, lam_1 = %.3g, norm error = %.3g\n', k, lam1, err);

        % save last iteration
        v1 = v1_kp1;
        Av1 = A*v1;

        % compute corresponding eigenvalue from Raleigh quotient
        lam1 = v1(:)' * Av1(:);

        % compute the next power iteration step
        v1_kp1 = Av1/norm(Av1(:));

        % calculate error from last iter and determine if convergence criteria is met
        err = norm(v1_kp1(:) - v1(:));
        if err < tol
            fprintf('convergence critera met! terminating early...\n');
            break
        end

    end

end