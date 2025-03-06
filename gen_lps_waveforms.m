function [g,rf,k_in,k_out] = gen_lps_waveforms(varargin)

    % define default arguments
    arg.fov = 20; % fov (cm)
    arg.N = 128; % nominal matrix size
    arg.nspokes = 23; % number of lps spokes
    arg.nseg = 280; % number of samples/segment
    arg.nrf = 3; % number of samples/rf pulse
    arg.fa = 4; % rf flip angle (deg)
    arg.gmax = 4; % max gradient amplitude (G/cm)
    arg.smax = 500; % max slew rate (G/cm/s)
    arg.dt = 4e-6; % raster time (s)

    gam = 4258; % GMR of H+ (Hz/G)

    % parse inputs
    arg = vararg_pair(arg,varargin);

    % calculate loop velocity
    om = 2*pi / (arg.nspokes * arg.dt * arg.nseg);

    % calculate gradient amplitude & slew rate
    g_amp = pi*arg.N / gam / (arg.dt * arg.nseg) / ...
        (arg.nspokes * arg.fov * 2 * sin(pi/arg.nspokes));
    s_amp = om * g_amp;
    assert(g_amp <= arg.gmax, ...
        'gradient amp exceeds limit with given parameters')
    assert(s_amp <= arg.smax, ...
        'slew rate exceeds limit with given parameters')

    % construct looping star gradients
    n = 0:2*arg.nspokes*arg.nseg-1;
    gx = g_amp * cos(om * n*arg.dt);
    gy = g_amp * sin(om * n*arg.dt);

    % calculate rf amplitude
    rf_amp = arg.fa / (360 * gam * arg.nrf*arg.dt); % (G)
    
    % construct rf burst pulse
    rf = rf_amp * (mod(n,arg.nseg) < arg.nrf) .* (n < arg.nspokes*arg.nseg);

    % calculate ramp-up
    nramp = ceil(g_amp / arg.smax / arg.dt);
    ramp_up = linspace(0,1,nramp+1);
    nramp = length(ramp_up);

    % append the ramp
    gx = [ramp_up*gx(1), gx, (1-ramp_up)*gx(end)];
    gy = [ramp_up*gy(1), gy, (1-ramp_up)*gy(end)];
    g = [gx(:),gy(:)];
    rf = padarray(rf,[0,nramp],0,"both");

    % calculate the full kspace trajectory for each spoke
    k_spokes = zeros([arg.nspokes*arg.nseg,2,arg.nspokes]);
    for v = 1:arg.nspokes
        idx_spoke = nramp + (v-1)*arg.nseg + (1:arg.nspokes*arg.nseg);
        k_spokes(:,:,v) = gam*cumsum(g(idx_spoke,:),1)*arg.dt;
    end
    k_spokes = permute(k_spokes,[1,3,2]); % [nseg*nspokes x nspokes x 2]

    % isolate in/out spokes
    k_out = reshape(k_spokes(1:arg.nseg,:,:),[],2); % spoke out
    ktmp = circshift(k_spokes,arg.nseg,1); % shift along fast time to get spoke in
    ktmp = circshift(ktmp,-1,2); % shift along spokes to align spokes in time
    k_in = reshape(ktmp(1:arg.nseg,:,:),[],2); % spoke in

end