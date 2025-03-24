function [Fs_in,Fs_out,b] = setup_nuffts(kdata,k_in,k_out,seq_args,varargin)
% sets up the volume-wise nufft objects given kspace data, trajectory, and
% desired indexing parameters
% by David Frey (djfrey@umich.edu)
%
% inputs:
% prjs2use - projection indicies to use (array of indicies, default all)
% ints2use - interleaf indicies to use (array of indicies, default all)
% reps2use - repetition indicies to use (array of indicies, default all)
% rpv - rotations per volume (number of prjs/ints/reps to use per vol)
% stride - number of rotations overlapping in each volume (not yet
% implemented)
% (also reads arguments from seq_args.mat file in current directory)
%
% outputs:
% Fs_in - cell array of nufft operators per volume (echo in)
%   if spokewise = true, Fs_in will be spokewise and framewise
% Fs_out - cell array of nufft operators per volume (echo out)
%   same applies for spokewise = true
% b - formatted kspace measurements
%
%
    
    % set default arguments
    arg.prjs2use = [];
    arg.ints2use = [];
    arg.reps2use = [];
    arg.rpv = [];
    arg.stride = 0; % still need to implement

    % parse arguments
    arg = vararg_pair(arg,varargin);
    
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
    nc = size(kdata,5);
    kdata = reshape(kdata,[],nvol,nc);
    k_in = reshape(k_in,[],nvol,3);
    k_out = reshape(k_out,[],nvol,3);

    % convert trajectory to spatial frequencies
    omega_in = 2*pi*seq_args.fov/seq_args.N * k_in;
    omega_out = 2*pi*seq_args.fov/seq_args.N * k_out;

    % format the kspace measurements volume-wise
    b = reshape(kdata,[],nvol,nc);

    % initialize volume-wise nufft operators
    Fs_in = cell(nvol,1);
    Fs_out = cell(nvol,1);
    
    % set up volume-wise dcf objects
    nufft_args = {seq_args.N*ones(1,3), 6*ones(1,3), 2*seq_args.N*ones(1,3), ...
        seq_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};

    % loop through volumes
    for ivol = 1:nvol
        % get trajectory for current vol
        omegav_in = reshape(omega_in(:,ivol,:),[],3);
        omegav_out = reshape(omega_out(:,ivol,:),[],3);

        % assemble nufft object for current vol
        Fs_in{ivol} = Gnufft( ...
            true(seq_args.N*ones(1,3)), ...
            [omegav_in, nufft_args]);
        Fs_out{ivol} = Gnufft( ...
            true(seq_args.N*ones(1,3)), ...
            [omegav_out, nufft_args]);
    end

end