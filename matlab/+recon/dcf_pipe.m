function W = dcf_pipe(F,niter)
% returns density compensation weighting Fatrix for non-uniform nufft F
% using Jim Pipe's method
% reference:
% Pipe JG, Menon P. Sampling density compensation in MRI: rationale and an
% iterative numerical solution. Magn Reson Med. 1999 Jan;41(1):179-86.
% doi: 10.1002/(sici)1522-2594(199901)41:1<179::aid-mrm25>3.0.co;2-v.
% by David Frey (djfrey@umich.edu)
%
% inputs:
% F - NUFFT operator(s)
% niter - number of iterations
% usepar - option to parallelize computation over volumes
%
% outputs:
% W - density compensation weighting operator (diagonal fatrix)
%

    if nargin < 2 || isempty(niter)
        niter = 5;
    end
    
    % convert to cell if single vol
    nvol = length(F);
    if ~iscell(F)
        F = {F};
    end
    W = cell(nvol,1);

    if nvol > 1 && usepar % recurse with parfor
        parfor i = 1:length(F)
            W{i} = recon.dcf_pipe(F(i), niter);
        end
    elseif nvol > 1 % recurse
        for i = 1:length(F)
            W{i} = recon.dcf_pipe(F(i), niter);
        end
    else % compute density compensation
        wi = ones(size(F,1),1);
        for itr = 1:niter
            wd = real( F.arg.st.interp_table(F.arg.st, ...
                F.arg.st.interp_table_adj(F.arg.st, wi) ) );
            wi = wi ./ wd;
        end
        W{1} = Gdiag(wi / sum(abs(wi)));
    end

    % convert back if single vol
    if nvol == 1
        W = W{1};
    end
    
end