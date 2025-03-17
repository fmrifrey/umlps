function [Ws_in,Ws_out] = setup_dcfs(Fs_in,Fs_out,varargin)

    arg.mode = 'pipe';
    arg.niter = 3;
    arg.usepar = false;

    arg = vararg_pair(arg,varargin);

    nvol = length(Fs_in);

    switch arg.mode
        case 'pipe'
            dcf = @(Fv) pipe_dcf(Fv,arg.niter);
        case 'cheap'
            dcf = @(Fv) Gdiag(vecnorm(Fv.arg.arg{1},2,2)/pi);
        otherwise
            error('invalid dcf mode: %s',arg.mode);
    end

    Ws_in = cell(nvol,1);
    Ws_out = cell(nvol,1);
    if arg.usepar
        parfor v = 1:nvol
            Ws_in{v} = dcf(Fs_in{v});
            Ws_out{v} = dcf(Fs_out{v});
        end
    else
        for v = 1:nvol
            Ws_in{v} = dcf(Fs_in{v});
            Ws_out{v} = dcf(Fs_out{v});
        end
    end

end

function W = pipe_dcf(F,niter)

    wi = ones(size(F,1),1);
    for itr = 1:niter
        wd = real( F.arg.st.interp_table(F.arg.st, ...
            F.arg.st.interp_table_adj(F.arg.st, wi) ) );
        wi = wi ./ wd;
    end
    W = Gdiag(wi / sum(abs(wi)));

end