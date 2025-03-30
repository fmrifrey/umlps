% basic reconstruction example for looping star data in matlab
% by David Frey

%% load the data
safile_lp = './scanarc.h5'; % scan archive file name
[kdata,k_in,k_out,seq_args] = lps_convert_data(safile_lp);

%% let up the volume-wise NUFFT objects and data
[Fs_in,Fs_out,b] = lpsutl.setup_nuffts(kdata,k_in,k_out,seq_args);
nvol = size(b,2);

%% calculate density compensation
Ws_in = cell(nvol,1);
Ws_out = cell(nvol,1);
for ivol = 1:nvol
    Ws_in{ivol} = lpsutl.dcf_pipe(Fs_in{ivol});
    Ws_out{ivol} = lpsutl.dcf_pipe(Fs_out{ivol});
end

%% create the kspace echo-in/out filters
[Hs_in,Hs_out] = lpsutl.setup_filters(Fs_in,Fs_out, ...
    0.5, ... % kspace filter cutoff
    0.1 ... % kspace filter rolloff
    );

%% coil compress the data
nc_cc = 16;
[b_cc,~,Vr_cc] = ir_mri_coil_compress(b,'ncoil',nc_cc);

%% reconstruct with dcNUFFT & RMS coil combination
nvol = size(b,2);
xv0 = zeros([seq_args.N*ones(1,3),nc_cc,nvol]);
for ivol = 1:nvol

    % set up problem for current volume
    WAv = kronI(nc_cc, ...
        Ws_in{ivol}*Hs_in{ivol}*Fs_in{ivol} + ...
        Ws_out{ivol}*Hs_out{ivol}*Fs_out{ivol});
    Av = kronI(nc_cc, ...
        Hs_in{ivol}*Fs_in{ivol} + ...
        Hs_out{ivol}*Fs_out{ivol});
    bv = squeeze(b_cc(:,ivol,:));

    % compute dc-NUFFT solution for volume
    fprintf('reconstructing dc-NUFFT sol to vol %d/%d\n', ivol, nvol);
    xv0 = WAv' * bv; % adjoint solution
    xv0 = ir_wls_init_scale(Av, bv, xv0); % correct scale
    xv0(:,:,:,:,ivol) = reshape(xv0,[seq_args.N*ones(1,3),nc_cc]);
    
end

% compute rms
x_rms = sqrt(squeeze(sum(xv0.^2,4)));

%% estimate sensitivity maps from GRE data using eSPIRIT
safile_gre = './scanarc.h5'; % scan archive file name
[gre_kdata,gre_seq_args] = gre3d_convert_data(safile_gre);

% coil compress the gre data
tmp = reshape(gre_kdata,gre_seq_args.N^3,size(gre_kdata,4));
tmp = tmp * Vr_cc;
gre_kdata_cc = reshape(tmp, [gre_seq_args.N*ones(1,3), nc_cc]);

% get acs data
acs_idcs = floor((gre_seq_args.N-gre_seq_args.Nacs)/2):floor((gre_seq_args.N+gre_seq_args.Nacs)/2);
acs_data = gre_kdata_cc(acs_idcs,acs_idcs,acs_idcs,:);

% estimate with eSPIRIT
smaps_lr = bart('ecalib -b0 -m1', acs_data);

% upsample the sensitivity maps for gre recon
smaps = lpsutl.resample3D(smaps_lr,gre_seq_args.N*ones(1,3));

% create fft-SENSE operator and recon gre image
F = fatrix2('idim', gre_seq_args.N*ones(1,3), ...
    'odim', gre_seq_args.N*ones(1,3), ...
    'forw', @(~,x) lpsutl.fftc(x,[],1:3), ...
    'back', @(~,x) lpsutl.ifftc(x,[],1:3));
FS = Asense(F,smaps);
xv0 = zeros(gre_seq_args.N*ones(1,3));
qp = Reg1(true(gre_seq_args.N*ones(1,3)),'beta',2^-3);
x_star = qpwls_pcg1(xv0, FS, 1, gre_kdata_cc(:), qp.C, 'niter', 20);
img_gre = reshape(x_star,gre_seq_args.N*ones(1,3));

%% reconstruct LpS data with SENSE recon
niter = 10;

% upsample the sensitivity maps for gre recon
smaps = lpsutl.resample3D(smaps_lr,seq_args.N*ones(1,3));
img_lps = zeros([Fs_in{1}.idim,nvol]);
for ivol = 1%:nvol

    % set up problem for current volume
    WAv = Asense(Ws_in{ivol}*Hs_in{ivol}*Fs_in{ivol} + ...
        Ws_out{ivol}*Hs_out{ivol}*Fs_out{ivol}, smaps);
    Av = Asense(Hs_in{ivol}*Fs_in{ivol} + ...
        Hs_out{ivol}*Fs_out{ivol}, smaps);
    bv = squeeze(b_cc(:,ivol,:));

    % get initial dc-NUFFT solution for volume
    fprintf('reconstructing initial sol to vol %d/%d\n', ivol, nvol);
    xv0 = WAv' * bv; % adjoint solution
    xv0 = ir_wls_init_scale(Av, bv, xv0); % correct scale
    xv0 = reshape(xv0, seq_args.N*ones(1,3));

    % solve with CG-SENSE
    fprintf('reconstructing CG-SENSE sol for vol %d/%d\n', ivol, nvol);
    qp = Reg1(true(seq_args.N*ones(1,3)),'beta',2^-3);
    x_star = qpwls_pcg1(xv0, Av, 1, bv(:), qp.C, 'niter', niter);
    img_lps(:,:,:,ivol) = reshape(x_star,seq_args.N*ones(1,3));

end