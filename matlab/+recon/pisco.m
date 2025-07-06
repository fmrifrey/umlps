function smaps = estimate_smaps_pisco(kdata,varargin)

    % set defaults
    arg.Nacs = 32; % ACS region size
    arg.kernel_shape = 'sphere'; % kernel shape
    arg.tau = 3; % neighborhood radius
    arg.rank = 500; % overestimated rank of C matrix
    arg.pad = true; % option to pad the convolutional matrix

    % get data size
    sz = size(kdata);
    if ndims(kdata) == 3 % 2D case - translate to 3D
        kdata = reshape(kdata,sz(1),sz(2),1,sz(3));
        sz = size(kdata);
    end

    % extract acs data
    even_RL = @(i) not(rem(i,2));
    acs_idx = -floor(arg.Nacs/2):floor(arg.Nacs/2) - even_RL(arg.Nacs/2);
    acs_idx_x = acs_idx + ceil(sz(1)/2) + even_RL(sz(1)/2);
    acs_idx_y = acs_idx + ceil(sz(2)/2) + even_RL(sz(2)/2);
    if sz(3) > 1
        acs_idx_z = acs_idx + ceil(sz(2)/2) + even_RL(sz(2)/2);
    else
        acs_idx_z = 1;
    end
    kdata_cal = kdata(acs_idx_x, acs_idx_y, acs_idx_z);

    % construct the kernel
    if sz(3) > 1
        [kern_idx_x, kern_idx_y, kern_idx_z] = ...
            meshgrid(-arg.tau:arg.tau, -arg.tau:arg.tau, -arg.tau:arg.tau);
    else
        [kern_idx_x, kern_idx_y] = ...
            meshgrid(-arg.tau:arg.tau, -arg.tau:arg.tau);
        kern_idx_z = zeros(size(kern_idx_x));
    end
    kern_idx_x = kern_idx_x(:);
    kern_idx_y = kern_idx_y(:);
    kern_idx_z = kern_idx_z(:);
    switch lower(arg.kernel_shape)
        case 'sphere'
            kern_msk = kern_idx_x.^2 + kern_idx_y.^2 + kern_idx_z.^2 <= arg.tau^2;
            kern_idx_x = kern_idx_x(kern_msk);
            kern_idx_y = kern_idx_y(kern_msk);
            kern_idx_z = kern_idx_z(kern_msk);
        case 'rect'
            % indicies are already good
        otherwise
            error('invalid kernel shape')
    end

    % construct the C matrix
    if arg.pad
        Nacsp = 2^ceil(log2(arg,Nacs+2*arg.tau));
    else
        Nacsp = arg.Nacs;
    end
    if sz(3) > 1
        inds = sub2ind(arg.Nacs*ones(1,3), ...
            floor(Nacsp/2) + 1 - kern_idx_x' + kern_idx_x, ...
            floor(Nacsp/2) + 1 - kern_idx_y' + kern_idx_y, ...
            floor(Nacsp/2) + 1 - kern_idx_z' + kern_idx_z);
        
        [kern_idx_x, kern_idx_y, kern_idx_z] = ...
            meshgrid(-arg.tau:arg.tau, -arg.tau:arg.tau, -arg.tau:arg.tau);
    else
        [kern_idx_x, kern_idx_y] = ...
            meshgrid(-arg.tau:arg.tau, -arg.tau:arg.tau);
        kern_idx_z = zeros(size(kern_idx_x));
    end
    kern_idx_x = kern_idx_x(:);
    kern_idx_y = kern_idx_y(:);
    kern_idx_z = kern_idx_z(:);

end