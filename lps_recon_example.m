% basic reconstruction example for looping star data
% by David Frey

%% load the data
safile = './scanarc.h5'; % scan archive file name
[kdata,k_in,k_out,seq_args] = lps_convert_data(safile);

%% let up the volume-wise NUFFT objects and data
[Fs_in,Fs_out,b] = lps.setup_nuffts(kdata,k_in,k_out,seq_args);
nvol = size(b,2);

%% calculate density compensation
Ws_in = cell(nvol,1);
Ws_out = cell(nvol,1);
for ivol = 1:nvol
    Ws_in{ivol} = lps.dcf_pipe(Fs_in{ivol});
    Ws_out{ivol} = lps.dcf_pipe(Fs_out{ivol});
end

%% create the kspace echo-in/out filters
[Hs_in,Hs_out] = lps.setup_filters(Fs_in,Fs_out, ...
    0.9, ... % kspace filter cutoff
    0.1 ... % kspace filter rolloff
    );

%% coil compress the data
nc_cc = 4;
b_cc = ir_mri_coil_compress(b,'ncoil',nc_cc);

%% reconstruct with RMS coil combination
niter = 0;
nvol = size(b,2);
x0 = zeros([Fs_in{1}.idim,nc_cc,nvol]);
for ivol = 1:nvol

    % set up problem for current volume
    WAv = kronI(nc_cc, ...
        Ws_in{ivol}*Hs_in{ivol}*Fs_in{ivol} + ...
        Ws_out{ivol}*Hs_out{ivol}*Fs_out{ivol});
    Av = kronI(nc_cc, ...
        Hs_in{ivol}*Fs_in{ivol} + ...
        Hs_out{ivol}*Fs_out{ivol});
    bv = squeeze(b_cc(:,ivol,:));

    % get dc-NUFFT solution for volume
    fprintf('reconstructing initial sol to vol %d/%d\n', ivol, nvol);
    xv0 = WAv' * bv; % adjoint solution
    xv0 = ir_wls_init_scale(Av, bv, xv0); % correct scale
    x0(:,:,:,:,ivol) = reshape(xv0,[Fs_in{1}.idim,nc_cc]);

end

%% reconstruct with SENSE recon
load('smaps.mat');
niter = 10;
nvol = size(b,2);
x0 = zeros([Fs_in{1}.idim,nvol]);
for ivol = 1:nvol

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
    x0(:,:,:,ivol) = reshape(xv0,Fs_in{1}.idim);

end