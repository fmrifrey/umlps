% Very basic reconstruction example for looping star data

%% Load the data
safile = './data.h5';
[kdata,k_in,k_out,seq_args] = lps.rec.get_data(safile);

%% Set up the volume-wise NUFFT objects and data
[Fs_in,Fs_out,b] = lps.rec.setup_nuffts(kdata,k_in,k_out,seq_args);
nvol = size(b,2);

%% Calculate density compensation
Ws_in = cell(nvol,1);
Ws_out = cell(nvol,1);
for ivol = 1:nvol
    Ws_in{ivol} = lps.rec.dcf_pipe(Fs_in{ivol});
    Ws_out{ivol} = lps.rec.dcf_pipe(Fs_out{ivol});
end

%% Create the kspace echo-in/out filters
[Hs_in,Hs_out] = lps.rec.setup_filters(Fs_in,Fs_out, ...
    0.5, ... % kspace filter cutoff
    0.1 ... % kspace filter rolloff
    );

%% Coil compress the data
b_cc = ir_mri_coil_compress(b,'ncoil',1);

%% Reconstruct initial solution with dc-NUFFT
nvol = size(b,2);
x0 = zeros([Fs_in{1}.idim,nvol]);
for ivol = 1:nvol

    % get FTs and filters for current volume
    Fv_in = Fs_in{ivol};
    Fv_out = Fs_out{ivol};
    Hv_in = Hs_in{ivol};
    Hv_out = Hs_out{ivol};
    Wv_in = Ws_in{ivol};
    Wv_out = Ws_out{ivol};

    % get dc-NUFFT solution for volume
    fprintf('reconstructing initial sol to vol %d/%d\n', ivol, nvol);
    xv0 = (Wv_in*Hv_in*Fv_in + Wv_in*Hv_in*Fv_in)' * b_cc; % adjoint solution
    xv0 = ir_wls_init_scale(Hv_in*Fv_in + Hv_in*Fv_in, b_cc, xv0); % correct scale
    x0(:,:,:,ivol) = reshape(xv0,Fs_in{1}.idim);

end
