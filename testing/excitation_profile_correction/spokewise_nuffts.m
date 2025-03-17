%% get data
safile = '~/data/lps_rftesting/nrf3_fa3/data.h5';
[kdata,k_in,k_out,seq_args] = lps.rec.get_data(safile);

%% set up frame-wise nuffts
t = tic;
[Fs_in_vw,Fs_out_vw,b_vw] = alt_setup_nuffts(kdata,k_in,k_out,seq_args, ...
    'prjs2use', 1:3);
t_setup_vw = toc(t);

%% set up spoke-wise nuffts
t = tic;
[Fs_in_sw,Fs_out_sw,b_sw] = alt_setup_nuffts(kdata,k_in,k_out,seq_args, ...
    'prjs2use', 1:3, 'spokewise', true);
t_setup_sw = toc(t);
nspokes = size(Fs_in_sw,1);
nvol = size(Fs_in_sw,2);

%% compute volume-wise nufft
t = tic;
x0_vw = (Fs_in_vw{1} + Fs_out_vw{1})' * b_vw;
t_comp_vw = toc(t);
x0_vw = reshape(x0_vw,seq_args.N*ones(1,3));
f_vw = figure;
im(x0_vw)
title(sprintf('volume-wise\ntime to construct: %.3fs\ntime to compute inufft: %.3fs', ...
    t_setup_vw, t_comp_vw));
saveas(f_vw,'volume-wise_inufft.png')

%% get time to compute spoke-wise nufft
t = tic;
x0_sw = 0;
for ispoke = 1:nspokes
    x0_sw = x0_sw + (Fs_in_sw{ispoke,1} + Fs_in_sw{ispoke,1})' * b_sw(:,ispoke);
end
t_comp_sw = toc(t);
x0_sw = reshape(x0_sw,seq_args.N*ones(1,3));
f_sw = figure;
im(x0_sw)
title(sprintf('spoke-wise\ntime to construct: %.3fs\ntime to compute inufft: %.3fs', ...
    t_setup_sw, t_comp_sw));
saveas(f_sw,'spoke-wise_inufft.png')

%% alternative version of setup_nuffts (as of 3/17/2025)
function [Fs_in,Fs_out,b] = alt_setup_nuffts(kdata,k_in,k_out,seq_args,varargin)
% by David Frey (djfrey@umich.edu)
%
% Inputs:
% safile - scanarchive data file name
% prjs2use - projection indicies to use (array of indicies)
% ints2use - interleaf indicies to use (array of indicies)
% reps2use - repetition indicies to use (array of indicies)
% rpv - rotations per volume (number of prjs/ints/reps to use per vol)
% stride - number of rotations overlapping in each volume
% spokewise - option to compute set up spokewise nuffts
% (also reads arguments from seq_args.mat file in current directory)
%
% Outputs:
% Fs_in - cell array of nufft operators per volume (echo in)
%   if spokewise = true, Fs_in will be spokewise and framewise
% Fs_out - cell array of nufft operators per volume (echo out)
%   same applies for spokewise = true
% b - formatted kspace measurements
%
% Note: each scale is treated as a seperate frame of the reconstruction
%
    
    % set default arguments
    arg.prjs2use = [];
    arg.ints2use = [];
    arg.reps2use = [];
    arg.rpv = [];
    arg.stride = 0; % still need to implement
    arg.spokewise = false;
    arg.usepar = false; % still need to implement

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
    kdata = reshape(kdata,[],seq_args.nspokes*arg.rpv,nvol,nc);
    k_in = reshape(k_in,[],seq_args.nspokes*arg.rpv,nvol,3);
    k_out = reshape(k_out,[],seq_args.nspokes*arg.rpv,nvol,3);

    % convert trajectory to spatial frequencies
    omega_in = 2*pi*seq_args.fov/seq_args.N * k_in;
    omega_out = 2*pi*seq_args.fov/seq_args.N * k_out;

    if arg.spokewise
        % format the kspace measurements spoke-wise
        b = reshape(kdata,seq_args.nseg,seq_args.nspokes*arg.rpv,nvol,nc);

        % initialize spoke-wise nufft operators
        Fs_in = cell(seq_args.nspokes*arg.rpv,nvol);
        Fs_out = cell(seq_args.nspokes*arg.rpv,nvol);
    else
        % format the kspace measurements volume-wise
        b = reshape(kdata,seq_args.nseg*seq_args.nspokes*arg.rpv,nvol,nc);

        % initialize volume-wise nufft operators
        Fs_in = cell(nvol,1);
        Fs_out = cell(nvol,1);
    end
    
    % set up spoke-wise nufft objects and volume-wise dcf objects
    nufft_args = {seq_args.N*ones(1,3), 6*ones(1,3), 2*seq_args.N*ones(1,3), ...
        seq_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};

    % loop through volumes and spokes
    for ivol = 1:nvol
        if arg.spokewise
            for ispoke = 1:seq_args.nspokes*arg.rpv
                % get trajectory for current vol/spoke
                omegasp_in = reshape(omega_in(:,ispoke,ivol,:),[],3);
                omegasp_out = reshape(omega_out(:,ispoke,ivol,:),[],3);

                % assemble nufft object for vol/spoke
                Fs_in{ispoke,ivol} = Gnufft( ...
                    true(seq_args.N*ones(1,3)), ...
                    [omegasp_in, nufft_args]);
                Fs_out{ispoke,ivol} = Gnufft( ...
                    true(seq_args.N*ones(1,3)), ...
                    [omegasp_out, nufft_args]);
            end
        else
            % get trajectory for current vol
            omegav_in = reshape(omega_in(:,:,ivol,:),[],3);
            omegav_out = reshape(omega_out(:,:,ivol,:),[],3);

            % assemble nufft object for whole volume
            Fs_in{ivol} = Gnufft( ...
                true(seq_args.N*ones(1,3)), ...
                [omegav_in, nufft_args]);
            Fs_out{ivol} = Gnufft( ...
                true(seq_args.N*ones(1,3)), ...
                [omegav_out, nufft_args]);
        end
    end

end