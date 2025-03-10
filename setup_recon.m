function [A,W,b] = setup_recon(varargin)
% Set up the reconstruction model for SPI3D data:
% b = WAx + noise
%
% by David Frey (djfrey@umich.edu)
%
% Inputs:
% safile - scanarchive data file name
% (also reads arguments from seq_args.mat file in current directory)
%
% Outputs:
% A - NUFFT operator (block diagonal if multiple scales are used)
% W - density compensation weighting matrix
% b - formatted measurements
%
% Note: each scale is treated as a seperate frame of the reconstruction
%
    
    % set default arguments
    arg.safile = 'data.h5'; % scanarchive data file name

    % parse arguments
    arg = vararg_pair(arg,varargin);

    % load in sequence arguments
    seq_args = load('seq_args.mat');
    
    % load first shot and get data size
    archive = GERecon('Archive.Load', arg.safile);
    shot = GERecon('Archive.Next', archive);
    [ndat,nc] = size(shot.Data);
    
    % load data
    d = zeros(ndat, nc, seq_args.nrep*seq_args.nint*seq_args.nprj);
    d(:, :, 1) = shot.Data;
    for l = 2:seq_args.nrep*seq_args.nint*seq_args.nprj
        shot = GERecon('Archive.Next', archive);
        d(:, :, l) = shot.Data;
    end
    
    % compress coils
    d = permute(d,[1,3,2]); % n x nrot x nc
    if nc > 1
        d = ir_mri_coil_compress(d,'ncoil',1);
    end

    [~,~,~,k_in,k_out] = gen_lps_waveforms( ...
        'fov', seq_args.fov, ... % fov (cm)
        'N', seq_args.N, ... % nominal matrix size
        'nspokes', seq_args.nspokes, ... % number of lps spokes
        'nseg', seq_args.nseg, ... % number of samples/segment
        'nrf', seq_args.nrf, ... % number of samples/rf pulse
        'fa', seq_args.fa, ... % rf flip angle (deg)
        'gmax', seq_args.gmax, ... % max gradient amplitude (G/cm)
        'smax', seq_args.smax, ... % max slew rate (G/cm/s)
        'dt', seq_args.dt ... % raster time (s))
        );
    
    % rotate the kspace trajectory
    k = zeros(3,ndat,seq_args.nint*seq_args.nprj);
    for iint = 1:seq_args.nint
        for iprj = 1:seq_args.nprj
            R = rot_3dtga(iprj,iint);
            k(:,:,(iint-1)*seq_args.nprj + iprj) = R * k0;
        end
    end
    kspace = permute(k,[2,3,1]);
    
    % make sure sizes match
    d = d(1:ndat,:);
    kspace = kspace(1:ndat,:,:);
    
    % create nufft object
    omega = 2*pi*seq_args.fov/seq_args.N*reshape(kspace,[],3);
    omega_msk = vecnorm(omega,2,2) < pi;
    omega = omega(omega_msk,:);
    nufft_args = {seq_args.N*ones(1,3), 6*ones(1,3), 2*seq_args.N*ones(1,3), seq_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};
    A = Gnufft(true(seq_args.N*ones(1,3)),[omega,nufft_args]);
    
    % calculate density compensation using the Pipe method
    wi = ones(size(A,1),1);
    for itr = 1:10 % 10 iterations
        wd = real( A.arg.st.interp_table(A.arg.st, ...
            A.arg.st.interp_table_adj(A.arg.st, wi) ) );
        wi = wi ./ wd;
    end
    W = Gdiag(wi / sum(abs(wi)));

    % repeat for each scale
    A = kronI(seq_args.nscl,A);
    W = kronI(seq_args.nscl,W);

    % mask out data
    d_r = reshape(d,[],seq_args.nscl);
    b = d_r(omega_msk,:);

end