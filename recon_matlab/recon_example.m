% basic reconstruction example for looping star data in matlab
% by David Frey

%% load the data
safile_lps = 'scanarc.h5'; % scan archive file name
[kdata,k_in,k_out,seq_args] = lps_convert_data(safile_lps);

%% let up the volume-wise NUFFT objects and data
[Fs_in,Fs_out,b] = recutl.setup_nuffts(kdata,k_in,k_out,seq_args,'prjs2use',1:200);
nvol = size(b,3);
nc = size(b,2);

%% calculate density compensation
Ws_in = cell(nvol,1);
Ws_out = cell(nvol,1);
for ivol = 1:nvol
    Ws_in{ivol} = recutl.dcf_pipe(Fs_in{ivol});
    Ws_out{ivol} = recutl.dcf_pipe(Fs_out{ivol});
end

%% create the kspace echo-in/out filters
[Hs_in,Hs_out] = recutl.setup_filters(Fs_in,Fs_out, ...
    0.8, ... % kspace filter cutoff
    0.1 ... % kspace filter rolloff
    );

%% estimate sensitivity maps from GRE data using eSPIRIT
safile_gre = './scanarc.h5'; % scan archive file name
[gre_kdata,msk,gre_seq_args] = gre3d_convert_data(safile_gre);

% coil compress the gre data
nc_cc = 16;
[gre_kdata_cc,~,Vr_cc] = ir_mri_coil_compress(gre_kdata,'ncoil',nc_cc);

% get acs data
acs_idcs = floor((gre_seq_args.N-gre_seq_args.Nacs)/2):floor((gre_seq_args.N+gre_seq_args.Nacs)/2-1);
acs_data = gre_kdata_cc(acs_idcs,acs_idcs,acs_idcs,:);

% estimate with eSPIRIT
smaps_lr = bart('ecalib -b0 -m1', acs_data);

%% recon the GRE data
% upsample sensitivity maps
smaps = recutl.resample3D(smaps_lr,gre_seq_args.N*ones(1,3));

% create fft-SENSE operator
F = fatrix2('idim', size(msk), ...
    'odim', size(msk), ...
    'omask', msk==1, ...
    'forw', @(~,x) recutl.fftc(x,[],1:3), ...
    'back', @(~,x) recutl.ifftc(x,[],1:3));
FS = Asense(F,smaps);

% create the regularizer
beta = 2^-4;
qp = Reg1(true(gre_seq_args.N*ones(1,3)), 'beta', beta);
% qpwls_psf(FS, qp.C, beta, true(gre_seq_args.N*ones(1,3)), 1, ...
%     'loop', 1, 'dx', gre_seq_args.fov/gre_seq_args.N, ...
%     'dz', gre_seq_args.fov/gre_seq_args.N);

% recon the GRE data with quadratic penalized CG-SENSE
niter = 100;
x0 = zeros(gre_seq_args.N*ones(1,3));
ytmp = reshape(gre_kdata_cc,[],nc_cc);
ytmp = ytmp(msk(:)==1,:);
x_star = qpwls_pcg1(x0, FS, 1, ytmp(:), qp.C, 'niter', niter,'isave','all');
img_gre = reshape(x_star(:,end),gre_seq_args.N*ones(1,3));

%% reconstruct LpS data with SENSE recon
niter = 30;
par_vols = false;
beta = 2^18;

% coil compress the data
b_cc = reshape(permute(b,[1,3,2]),[],nc)*Vr_cc;
b_cc = permute(reshape(b_cc,[size(b,1),nvol,nc_cc]),[1,3,2]);

% upsample the sensitivity maps for lps recon
smaps = recutl.resample3D(smaps_lr,seq_args.N*ones(1,3));

% initialize with dc-NUFFT adjoint SENSE recon
HWs_in = cell(nvol,1);
HWs_out = cell(nvol,1);
for i = 1:nvol
    HWs_in{i} = Hs_in{i}*Ws_in{i};
    HWs_out{i} = Hs_out{i}*Ws_out{i};
end
WA = recutl.A_volwise(Fs_in,Fs_out,HWs_in,HWs_out,smaps,par_vols);
x0 = WA' * b_cc;

% create the system matrix
A = recutl.A_volwise(Fs_in,Fs_out,Hs_in,Hs_out,smaps,par_vols);

% make the regularizer (quadratic differencing penalty)
qp = Reg1(true(seq_args.N*ones(1,3)),'beta',beta);
% Av1 = Asense(Hs_in{1}*Fs_in{1} + Hs_out{1}*Fs_out{1}, smaps);
% qpwls_psf(Av1, qp.C, beta, true(seq_args.N*ones(1,3)),1, ...
%     'loop', 1, 'dx', seq_args.fov/seq_args.N, 'dz', seq_args.fov/seq_args.N); % code to determine reg parameter
T = qp.C;
if nvol > 1
    T = kronI(nvol,T);
end

%% solve with CG
x0 = ir_wls_init_scale(A,b_cc,x0);
x_star = qpwls_pcg1(x0, A, 1, b_cc(:), T, 'niter', niter);
img_lps = reshape(x_star,[seq_args.N*ones(1,3),nvol]);