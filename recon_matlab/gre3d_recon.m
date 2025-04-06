fname = './gre3d.h5'; % gre h5 file to load

% load in the data
s = recutl.loadh5struct(fname);
kdata = s.kdata.real + 1i*s.kdata.imag;
msk = s.msk;
seq_args = s.seq_args;

%% estimate sensitivity maps with eSPIRiT
fname = './smaps.h5'; % smaps file to write to

% use bart to estimate smaps with eSPIRiT
smaps = bart(sprintf('ecalib -b0 -m1 -r%d',seq_args.Nacs), kdata);

% write out the smaps to h5 file
h5create(fname, '/real', size(smaps), ...
    'Datatype', class(real(smaps)));
h5write(fname, '/real', real(smaps));
h5create(fname, '/imag', size(smaps), ...
    'Datatype', class(imag(smaps)));
h5write(fname, '/imag', imag(smaps));

%% load in the sensitivity maps and coil compress
nc = 8; % number of virtual coils to keep
fname = './smaps.h5'; % smap file to read from

% load the smaps
s = recutl.loadh5struct(fname);
smaps = s.real + 1i*s.imag;

% compress data and get compression matrix
[kdata,~,Vr] = ir_mri_coil_compress(kdata,'ncoil',nc);

% upsample smaps
smaps = recutl.resample3D(smaps,seq_args.N*ones(1,3));

% coil compress the smaps
smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[seq_args.N*ones(1,3),nc]);

%% recon the GRE data with CG-SENSE
niter = 5; % number of CG iterations
beta = 2^-4; % regularization parameter

% create fft-SENSE operator
F = fatrix2('idim', size(msk), ...
    'odim', size(msk), ...
    'omask', msk==1, ...
    'forw', @(~,x) recutl.fftc(x,[],1:3), ...
    'back', @(~,x) recutl.ifftc(x,[],1:3));
FS = Asense(F,smaps);

% create the regularizer
qp = Reg1(true(seq_args.N*ones(1,3)), 'beta', beta);
% qpwls_psf(FS, qp.C, beta, true(seq_args.N*ones(1,3)), 1, ...
%     'loop', 1, 'dx', seq_args.fov/seq_args.N, ...
%     'dz', seq_args.fov/seq_args.N);

% recon the GRE data with quadratic penalized CG-SENSE
x0 = zeros(seq_args.N*ones(1,3));
ytmp = reshape(kdata,[],nc);
ytmp = ytmp(msk(:)==1,:);
x_star = qpwls_pcg1(x0, FS, 1, ytmp(:), qp.C, 'niter', niter,'isave','all');
img_gre = reshape(x_star(:,end),seq_args.N*ones(1,3));