function [kdata,k_in,k_out,seq_args] = lps_convert_data(safile, h5file)
% gets data from scanarchive file, and formats the kspace trajectories and
% sequence arguments based on seq_args.mat, then writes formatted data to
% h5 file for external recon if specified
% by David Frey (djfrey@umich.edu)
%
% inputs:
% safile - name of scanarchive file to read in
% h5file - name of output h5 file to write formatted data to (leave empty
% to not write to file)
%
% outputs:
% kdata - kspace data (ndat x nprj x nint x nrep x ncoil)
% k_in - spoke-in kspace trajectory (ndat x nprj x nint x nrep x 3)
% k_out - spoke-out kspace trajectory (ndat x nprj x nint x nrep x 3)
% seq_args - struct containing pulse sequence arguments
%

    % get directory and file names
    d = dir(safile);
    sadir = d(1).folder;
    safile = d(1).name;

    % load in sequence arguments
    seq_args = load([sadir,'/seq_args.mat']);
    
    % load first shot and get data size
    archive = GERecon('Archive.Load', [sadir,'/',safile]);
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
    kdata = permute(kdata,[1,5:-1:3,2]); % n x nprj x nint x nrep x nc

    % generate kspace trajectory
    [~,~,~,k_in0,k_out0] = lps.gen_lps_waveforms( ...
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
            R = lps.rot_3dtga(iprj,iint);
            k_in(:,:,iprj,iint) = k_in0 * R';
            k_out(:,:,iprj,iint) = k_out0 * R';
        end
    end

    % repeat for repetitions
    k_in = repmat(k_in,[1,1,1,1,seq_args.nrep]); % [n x 3 x nprj x nint x nrep]
    k_out = repmat(k_out,[1,1,1,1,seq_args.nrep]);
    k_in = permute(k_in,[1,3:5,2]); % [n x nprj x nint x nrep x 3]
    k_out = permute(k_out,[1,3:5,2]);

    % save h5 file
    if nargin > 1 && ~isempty(h5file)

        if isfile(h5file)
            system(sprintf('rm %s',h5file));
        end

        % save kspace data
        h5create(h5file, '/kdata/real', size(kdata), ...
            'Datatype', class(real(kdata)));
        h5write(h5file, '/kdata/real', real(kdata));
        h5create(h5file, '/kdata/imag', size(kdata), ...
            'Datatype', class(imag(kdata)));
        h5write(h5file, '/kdata/imag', imag(kdata));

        % save sampling locations
        h5create(h5file, '/ktraj/spoke_in', size(k_in), 'Datatype', class(k_in));
        h5write(h5file, '/ktraj/spoke_in', k_in);
        h5create(h5file, '/ktraj/spoke_out', size(k_out), 'Datatype', class(k_out));
        h5write(h5file, '/ktraj/spoke_out', k_out);

        % save sequence arguments
        seq_args_fields = fieldnames(seq_args);
        for i = 1:numel(seq_args_fields)
            field = seq_args_fields{i};
            val = seq_args.(field);
            if islogical(val)
                val = 1*val;
            end
            h5create(h5file, sprintf('/seq_args/%s',field), size(val), ...
                'Datatype', class(val));
            h5write(h5file, sprintf('/seq_args/%s',field), val)
        end

    end

end

