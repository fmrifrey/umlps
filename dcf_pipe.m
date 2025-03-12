function W = dcf_pipe(F,niter)
% returns density compensation weighting Fatrix for non-uniform nufft F
% using Jim Pipe's method

    if nargin < 2 || isempty(niter)
        niter = 10;
    end

    wi = ones(size(F,1),1);
    for itr = 1:niter
        wd = real( F.arg.st.interp_table(F.arg.st, ...
            F.arg.st.interp_table_adj(F.arg.st, wi) ) );
        wi = wi ./ wd;
    end
    W = Gdiag(wi / sum(abs(wi)));

end