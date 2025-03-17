function [Hs_in,Hs_out] = setup_filts(Fs_in,Fs_out,cutoff,rolloff)

    nvol = length(Fs_in);

    % create kspace filter
    kfilt = @(w) (w < pi*cutoff) .* 1./(1 + exp(2*pi*(w - pi*cutoff)/(pi*rolloff)));
    
    Hs_in = cell(nvol,1);
    Hs_out = cell(nvol,1);
    for ivol = 1:nvol
        omega_in = Fs_in{ivol}.arg.arg{1};
        omega_out = Fs_out{ivol}.arg.arg{1};
        filt_in = kfilt(vecnorm(omega_in,2,2));
        filt_out = kfilt(vecnorm(omega_out,2,2));
        Hs_in{ivol} = Gdiag(filt_in);
        Hs_out{ivol} = Gdiag(filt_out);
    end

end

