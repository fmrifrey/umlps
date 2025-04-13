%% load LpS data from h5 file
fname = './lps.h5';
s = recutl.loadh5struct(fname);
kdata = s.kdata.real + 1i*s.kdata.imag;
k_in = s.ktraj.spoke_in;
k_out = s.ktraj.spoke_out;
seq_args = s.seq_args;

%% set up the volume-wise NUFFT objects and data
[Fs_in,Fs_out,b] = recutl.setup_nuffts(kdata,k_in,k_out,seq_args,'rpv',32);
nvol = size(b,3);

%% calculate density compensation
Ws_in = cell(nvol,1);
Ws_out = cell(nvol,1);
parfor ivol = 1:nvol
    Ws_in{ivol} = recutl.dcf_pipe(Fs_in{ivol});
    Ws_out{ivol} = recutl.dcf_pipe(Fs_out{ivol});
end

%% create the kspace echo-in/out filters
[Hs_in,Hs_out] = recutl.setup_filters(Fs_in,Fs_out, ...
    0.8, ... % kspace filter cutoff
    0.1 ... % kspace filter rolloff
    );

%% load in the sensitivity maps and coil compress
nc = 8; % number of virtual coils to keep
fname = '../smaps.h5'; % smap file to read from

% load the smaps
s = recutl.loadh5struct(fname);
smaps = s.real + 1i*s.imag;

% compress data and get compression matrix
[tmp,~,Vr] = ir_mri_coil_compress(permute(b,[1,3,2]),'ncoil',nc);
b = permute(tmp,[1,3,2]);

% upsample smaps
smaps = recutl.resample3D(smaps,seq_args.N*ones(1,3));

% coil compress the smaps
smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[seq_args.N*ones(1,3),nc]);

%% recon the data with CG-SENSE
niter = 30;
par_vols = true;
beta = 2^20;

% initialize with dc-NUFFT adjoint SENSE recon
HWs_in = cell(nvol,1);
HWs_out = cell(nvol,1);
for i = 1:nvol
    HWs_in{i} = Hs_in{i}*Ws_in{i};
    HWs_out{i} = Hs_out{i}*Ws_out{i};
end
WA = recutl.A_volwise(Fs_in,Fs_out,HWs_in,HWs_out,smaps,par_vols);
x0 = WA' * b;

% create the system matrix
A = recutl.A_volwise(Fs_in,Fs_out,Hs_in,Hs_out,smaps,par_vols);

% make the regularizer (quadratic differencing penalty)
qp = Reg1(true(seq_args.N*ones(1,3)),'beta',beta);
% % code to determine a good regularization parameter:
% Av1 = Asense(Hs_in{1}*Fs_in{1} + Hs_out{1}*Fs_out{1}, smaps);
% qpwls_psf(Av1, qp.C, beta, true(seq_args.N*ones(1,3)),1, ...
%     'loop', 1, 'dx', seq_args.fov/seq_args.N, 'dz', seq_args.fov/seq_args.N);
T = qp.C;
if nvol > 1
    T = kronI(nvol,T);
end

% solve with CG
x0 = ir_wls_init_scale(A,b,x0);
x_star = qpwls_pcg1(x0, A, 1, b(:), T, 'niter', niter);
img_lps = reshape(x_star,[seq_args.N*ones(1,3),nvol]);