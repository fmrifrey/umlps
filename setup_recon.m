function [A,W,F,b] = setup_recon(varargin)
% Set up the reconstruction model for lps data:
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
    arg.cutoff = 0.8;
    arg.rolloff = 0.1;
    arg.dcf = 'cheap';
    arg.dcf_niter = 3;
    arg.prjs2use = [];
    arg.ints2use = [];
    arg.reps2use = [];
    arg.rpv = [];

    % parse arguments
    arg = vararg_pair(arg,varargin);

    % load in sequence arguments
    seq_args = load('seq_args.mat');
    
    % load first shot and get data size
    archive = GERecon('Archive.Load', arg.safile);
    shot = GERecon('Archive.Next', archive);
    [ndat,nc] = size(shot.Data);
    
    % load data
    kdata = zeros(ndat, nc, seq_args.nrep*seq_args.nint*seq_args.nprj);
    kdata(:, :, 1) = shot.Data;
    for l = 2:seq_args.nrep*seq_args.nint*seq_args.nprj
        shot = GERecon('Archive.Next', archive);
        kdata(:, :, l) = shot.Data;
    end
    kdata = reshape(kdata,[ndat,nc,seq_args.nrep,seq_args.nint,seq_args.nprj]);
    kdata = permute(kdata,[1,3:5,2]); % n x nrep x nint x nprj x nc
    
    % compress to 1 coil
    if nc > 1
        nc = 1;
        kdata = ir_mri_coil_compress(kdata,'ncoil',nc);
    end

    % generate kspace trajectory
    [~,~,~,k_in0,k_out0] = gen_lps_waveforms( ...
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

    % convert to 3D
    k_in0 = padarray(k_in0,[0,1],0,'post');
    k_out0 = padarray(k_out0,[0,1],0,'post');
    
    % rotate the kspace trajectory
    k_in = zeros(ndat,3,1,seq_args.nint,seq_args.nprj);
    k_out = zeros(ndat,3,1,seq_args.nint,seq_args.nprj);
    for iint = 1:seq_args.nint
        for iprj = 1:seq_args.nprj
            R = rot_3dtga(iprj,iint);
            k_in(:,:,1,iint,iprj) = k_in0 * R';
            k_out(:,:,1,iint,iprj) = k_out0 * R';
        end
    end
    k_in = repmat(k_in,[1,1,seq_args.nrep,1,1]); % [n x 3 x nrep x nint x nprj]
    k_out = repmat(k_out,[1,1,seq_args.nrep,1,1]);
    k_in = permute(k_in,[1,3:5,2]); % [n x nrep x nshots x nrots x 3]
    k_out = permute(k_out,[1,3:5,2]);
    
    % determine loop indicies to use
    if isempty(arg.ints2use)
        arg.ints2use = 1:seq_args.nint; % use all unique in-plane rotations
        nint = seq_args.nint;
    else
        nint = length(arg.ints2use);
    end
    if isempty(arg.prjs2use)
        arg.prjs2use = 1:seq_args.nprj; % use all unique thru-plane rotations
        nprj = seq_args.nprj;
    else
        nprj = length(arg.prjs2use);
    end
    if isempty(arg.reps2use)
        arg.reps2use = 1:seq_args.nrep; % use all repetitions
        nrep = seq_args.nrep;
    else
        nrep = length(arg.reps2use);
    end

    % index desired loop indicies
    kdata = kdata(:,arg.reps2use,arg.ints2use,arg.prjs2use,:);
    k_in = k_in(:,arg.reps2use,arg.ints2use,arg.prjs2use,:);
    k_out = k_out(:,arg.reps2use,arg.ints2use,arg.prjs2use,:);

    % split data and trajectory into 3D volumes
    if isempty(arg.rpv)
        arg.rpv = nint*nprj; % each rep is a vol
    end
    if mod(nrep*nint*nprj, arg.rpv)
        error('total planes (%d) must be divisible by rpv (%d)', ...
            nrep*nint*nprj, arg.rpv)
    else
        nvol = nrep*nint*nprj / arg.rpv;
    end
    kdata = reshape(kdata,[],nvol,nc);
    k_in = reshape(k_in,[],nvol,3);
    k_out = reshape(k_out,[],nvol,3);

    % only recon the first vol for now...
    kdata = squeeze(kdata(:,1,:));
    k_in = squeeze(k_in(:,1,:));
    k_out = squeeze(k_out(:,1,:));

    % convert trajectory to spatial frequencies
    omega_in = 2*pi*seq_args.fov/seq_args.N * k_in;
    omega_out = 2*pi*seq_args.fov/seq_args.N * k_out;

    % create kspace filter
    w_cut = pi*arg.cutoff;
    kfilt = @(w) (w < w_cut) .* 1./(1 + exp(2*pi*(w - w_cut)/(pi*arg.rolloff)));
    filt_in = kfilt(vecnorm(omega_in,2,2));
    filt_out = kfilt(vecnorm(omega_out,2,2));
    F_in = Gdiag(filt_in);
    F_out = Gdiag(filt_out);
    
    % create NUFFT objects
    nufft_args = {seq_args.N*ones(1,3), 6*ones(1,3), 2*seq_args.N*ones(1,3), ...
        seq_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};
    A_in = Gnufft(true(seq_args.N*ones(1,3)),[omega_in,nufft_args]);
    A_out = Gnufft(true(seq_args.N*ones(1,3)),[omega_out,nufft_args]); 

    % calculate dcf
    switch lower(arg.dcf)
        case 'cheap'
            W_in = Gdiag(vecnorm(omega_in,2,2)/pi);
            W_out = Gdiag(vecnorm(omega_out,2,2)/pi);
        case 'pipe'
            W_in = dcf_pipe(A_in,arg.dcf_niter);
            W_out = dcf_pipe(A_out,arg.dcf_niter);
        otherwise
            error('invalid dcf mode: %s',arg.dcf);
    end

    % make y
    b = reshape(kdata,[],nvol);

    % repeat for each volume
    if nvol > 1
        A_in = kronI(nvol,A_in);
        W_in = kronI(nvol,W_in);
        F_in = kronI(nvol,F_in);
        A_out = kronI(nvol,A_out);
        W_out = kronI(nvol,W_out);
        F_out = kronI(nvol,F_out);
    end

    % return
    A = {A_in,A_out};
    W = {W_in,W_out};
    F = {F_in,F_out};

end