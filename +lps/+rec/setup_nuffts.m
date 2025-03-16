function [Fv_in,Fv_out,b] = setup_nuffts(varargin)
% by David Frey (djfrey@umich.edu)
%
% Inputs:
% safile - scanarchive data file name
% prjs2use - projection indicies to use (array of indicies)
% ints2use - interleaf indicies to use (array of indicies)
% reps2use - repetition indicies to use (array of indicies)
% rpv - rotations per volume (number of prjs/ints/reps to use per vol)
% stride - number of rotations overlapping in each volume
% (also reads arguments from seq_args.mat file in current directory)
%
% Outputs:
% Fv_in - cell array of nufft operators per volume (echo in)
% Fv_out - cell array of nufft operators per volume (echo out)
% b - formatted kspace measurements
%
% Note: each scale is treated as a seperate frame of the reconstruction
%
    
    % set default arguments
    arg.safile = 'data.h5'; % scanarchive data file name
    arg.prjs2use = [];
    arg.ints2use = [];
    arg.reps2use = [];
    arg.rpv = [];
    arg.stride = 0;

    gam = 4258; % GMR of H+ (Hz/G)

    % parse arguments
    arg = vararg_pair(arg,varargin);

    % load in sequence arguments
    seq_args = load('seq_args.mat');
    
    % load first shot and get data size
    archive = GERecon('Archive.Load', arg.safile);
    shot = GERecon('Archive.Next', archive);
    [ndat,nc] = size(shot.Data);
    
    % load data
    fprintf('loading data...\n');
    kdata = zeros(ndat, nc, seq_args.nrep*seq_args.nint*seq_args.nprj);
    kdata(:, :, 1) = shot.Data;
    for l = 2:seq_args.nrep*seq_args.nint*seq_args.nprj
        shot = GERecon('Archive.Next', archive);
        kdata(:, :, l) = shot.Data;
    end
    kdata = reshape(kdata,[ndat,nc,seq_args.nrep,seq_args.nint,seq_args.nprj]);
    kdata = permute(kdata,[1,5:-1:3,2]); % n x nprj x nint x nrep x nc

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
    
    % rotate the kspace trajectory & excitation grads
    k_in = zeros(ndat,3,seq_args.nprj,seq_args.nint);
    k_out = zeros(ndat,3,seq_args.nprj,seq_args.nint);
    for iint = 1:seq_args.nint
        for iprj = 1:seq_args.nprj
            R = rot_3dtga(iprj,iint);
            k_in(:,:,iprj,iint) = k_in0 * R';
            k_out(:,:,iprj,iint) = k_out0 * R';
        end
    end

    % repeat for repetitions
    k_in = repmat(k_in,[1,1,1,1,seq_args.nrep]); % [n x 3 x nprj x nint x nrep]
    k_out = repmat(k_out,[1,1,1,1,seq_args.nrep]);
    k_in = permute(k_in,[1,3:5,2]); % [n x nprj x nint x nrep x 3]
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
    kdata = kdata(:,arg.prjs2use,arg.ints2use,arg.reps2use,:);
    k_in = k_in(:,arg.prjs2use,arg.ints2use,arg.reps2use,:);
    k_out = k_out(:,arg.prjs2use,arg.ints2use,arg.reps2use,:);

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

    % convert trajectory to spatial frequencies
    omega_in = 2*pi*seq_args.fov/seq_args.N * k_in;
    omega_out = 2*pi*seq_args.fov/seq_args.N * k_out;

    % initialize volume-wise nufft operators
    Fv_in = cell(nvol,1);
    Fv_out = cell(nvol,1);
    
    % set up spoke-wise nufft objects and volume-wise dcf objects
    nufft_args = {seq_args.N*ones(1,3), 6*ones(1,3), 2*seq_args.N*ones(1,3), ...
        seq_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};

    % loop through volumes and spokes
    for ivol = 1:nvol
        % get trajectory for current vol
        omega_vol_in = reshape(omega_in(:,:,ivol,:),[],3);
        omega_vol_out = reshape(omega_out(:,:,ivol,:),[],3);

        % assemble nufft object for whole volume
        Fv_in{ivol} = Gnufft( ...
                true(seq_args.N*ones(1,3)), ...
                [omega_vol_in, nufft_args]);
        Fv_out{vol} = Gnufft( ...
                true(seq_args.N*ones(1,3)), ...
                [omega_vol_out, nufft_args]);
    end

    % format the kspace measurements
    b = reshape(kdata,[],nvol,nc);

end