%% set recon arguments
basedir = './'; % directory containing data
fname_kdata = 'lps_fmri_unfold.h5'; % name of input .h5 file (in basedir)
fname_smaps = 'smaps.h5'; % name of smaps .h5 file (in basedir)
fname_out = split(fname_kdata,'.'); fname_out = [fname_out{1},'_recon.h5']; % name of output recon .h5 file (in basedir)
ncoil_comp = 4; % number of coils to compress to
cutoff = 0.8; % kspace cutoff for echo-in/out filtering
rolloff = 0.1; % kspace rolloff for echo-in/out filtering
beta = 2^20; % regularization parameter for quadratic finite differencing penalty
niter = 30; % number of iterations for CG
M = 92; % reconstruction matrix size (leave empty for cutoff * N)
ints2use = []; % number of interleaves to use (leave empty for all)
prjs2use = []; % number of projections to use (leave empty for all)
reps2use = []; % number of repetitions to use (leave empty for all)
volwidth = 32; % number of projections per volume (leave empty for all)
par_vols = true; % option to parallelize volume-wise computations

%% load LpS data from h5 file
s = recutl.loadh5struct(fullfile(basedir,fname_kdata));
kdata = s.kdata.real + 1i*s.kdata.imag;
k_in = s.ktraj.spoke_in;
k_out = s.ktraj.spoke_out;
seq_args = s.seq_args;
if isempty(M)
    M = cutoff*seq_args.N;
end

%% set up the volume-wise NUFFT objects and data
if isempty(ints2use)
    ints2use = seq_args.nint; % use all unique in-plane rotations
end
if isempty(prjs2use)
   prjs2use = seq_args.nprj; % use all unique thru-plane rotations
end
if isempty(reps2use)
    reps2use = seq_args.nrep; % use all repetitions
end
if isempty(volwidth)
    volwidth = ints2use*prjs2use; % each rep is a vol
end
[Fs_in,Fs_out,b] = recutl.setup_nuffts(kdata,k_in,k_out,seq_args, ...
    'ints2use', ints2use, ...
    'prjs2use', prjs2use, ...
    'reps2use', reps2use, ...
    'volwidth', volwidth);
nvol = size(b,3);

%% calculate density compensation
Ws_in = cell(nvol,1);
Ws_out = cell(nvol,1);
if par_vols
    parfor ivol = 1:nvol
        Ws_in{ivol} = recutl.dcf_pipe(Fs_in{ivol});
        Ws_out{ivol} = recutl.dcf_pipe(Fs_out{ivol});
    end
else
    for ivol = 1:nvol
        Ws_in{ivol} = recutl.dcf_pipe(Fs_in{ivol});
        Ws_out{ivol} = recutl.dcf_pipe(Fs_out{ivol});
    end
end

%% create the kspace echo-in/out filters
[Hs_in,Hs_out] = recutl.setup_filters(Fs_in,Fs_out, ...
    cutoff, ... % kspace filter cutoff
    rolloff ... % kspace filter rolloff
    );

%% load in the sensitivity maps and coil compress
s = recutl.loadh5struct(fullfile(basedir,fname_smaps));
smaps = s.real + 1i*s.imag;

% compress data and get compression matrix
[tmp,~,Vr] = ir_mri_coil_compress(permute(b,[1,3,2]),'ncoil',ncoil_comp);
b = permute(tmp,[1,3,2]);

% upsample smaps
smaps = recutl.resample3D(smaps,seq_args.N*ones(1,3));

% coil compress the smaps
smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[seq_args.N*ones(1,3),ncoil_comp]);

%% recon the data with CG-SENSE
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

%% save to h5 recon file
fname = fullfile(basedir,fname_out);
if exist(fname, 'file')
    delete(fname);
end

% save solution
h5create(fname, '/sol/real', size(real(img_lps)));
h5write(fname, '/sol/real', real(img_lps));
h5create(fname, '/sol/imag', size(imag(img_lps)));
h5write(fname, '/sol/imag', imag(img_lps));

% save sequence args
h5create(fname, '/seq_args/fov', [1,1]);
h5write(fname, '/seq_args/fov', seq_args.fov);
h5create(fname, '/seq_args/tr', [1,1]);
h5write(fname, '/seq_args/tr', seq_args.tr);
h5create(fname, '/seq_args/fa', [1,1]);
h5write(fname, '/seq_args/fa', seq_args.fa);
h5create(fname, '/seq_args/dummyshots', [1,1]);
h5write(fname, '/seq_args/dummyshots', seq_args.dummyshots);
h5create(fname, '/seq_args/gmax', [1,1]);
h5write(fname, '/seq_args/gmax', seq_args.gmax);
h5create(fname, '/seq_args/smax', [1,1]);
h5write(fname, '/seq_args/smax', seq_args.smax);
h5create(fname, '/seq_args/nrf', [1,1]);
h5write(fname, '/seq_args/nrf', seq_args.nrf);
h5create(fname, '/seq_args/pislquant', [1,1]);
h5write(fname, '/seq_args/pislquant', seq_args.pislquant);
h5create(fname, '/seq_args/N', [1,1]);
h5write(fname, '/seq_args/N', seq_args.N);
h5create(fname, '/seq_args/nseg', [1,1]);
h5write(fname, '/seq_args/nseg', seq_args.nseg);
h5create(fname, '/seq_args/nspokes', [1,1]);
h5write(fname, '/seq_args/nspokes', seq_args.nspokes);
h5create(fname, '/seq_args/nint', [1,1]);
h5write(fname, '/seq_args/nint', seq_args.nint);
h5create(fname, '/seq_args/nprj', [1,1]);
h5write(fname, '/seq_args/nprj', seq_args.nprj);
h5create(fname, '/seq_args/nrep', [1,1]);
h5write(fname, '/seq_args/nrep', seq_args.nrep);

% save recon args
h5create(fname, '/recon_args/ncoil_comp', size(ncoil_comp));
h5write(fname, '/recon_args/ncoil_comp', ncoil_comp);
h5create(fname, '/recon_args/cutoff', size(cutoff));
h5write(fname, '/recon_args/cutoff', cutoff);
h5create(fname, '/recon_args/rolloff', size(rolloff));
h5write(fname, '/recon_args/rolloff', rolloff);
h5create(fname, '/recon_args/beta', size(beta));
h5write(fname, '/recon_args/beta', beta);
h5create(fname, '/recon_args/niter', size(niter));
h5write(fname, '/recon_args/niter', niter);
h5create(fname, '/recon_args/ints2use', size(ints2use));
h5write(fname, '/recon_args/ints2use', ints2use);
h5create(fname, '/recon_args/prjs2use', size(prjs2use));
h5write(fname, '/recon_args/prjs2use', prjs2use);
h5create(fname, '/recon_args/reps2use', size(reps2use));
h5write(fname, '/recon_args/reps2use', reps2use);
h5create(fname, '/recon_args/volwidth', size(volwidth));
h5write(fname, '/recon_args/volwidth', volwidth);
h5create(fname, '/recon_args/M', size(M));
h5write(fname, '/recon_args/M', M);
