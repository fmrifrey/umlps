function W = dcf_pipe(A)
    wi = ones(size(A,1),1);
    for itr = 1:10 % 10 iterations
        wd = real( A.arg.st.interp_table(A.arg.st, ...
            A.arg.st.interp_table_adj(A.arg.st, wi) ) );
        wi = wi ./ wd;
    end
    W = Gdiag(wi / sum(abs(wi)));
end