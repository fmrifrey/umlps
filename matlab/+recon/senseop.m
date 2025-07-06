function AS = senseop(A,smaps)
% adds sensitivity encoding to forward operator A
% by David Frey
%
% inputs:
% A - forward fatrix operator
% smaps - sensitivity maps
%
% outputs:
% AS - forward fatrix operator including sensitivity encoding
%

    % get number of coils
    nc = size(smaps,ndims(smaps));

    % get image mask
    msk = A.imask;

    % mask out the smaps
    smaps = reshape(smaps,[],nc);
    smaps_msk = smaps(msk(:),:);

    % create cell array of sensitivity operators
    AS_cell = cell(nc,1);
    for i = 1:nc
        AS_cell{i} = A*Gdiag(smaps_msk(:,i), 'mask', msk);
    end

    % concatenate into block column matrix
    AS = cat(1,AS_cell{:});

end
