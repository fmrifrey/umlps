function write_lps_seq(varargin)

    % define default arguments
    arg.fov = 22; % fov (cm)
    arg.N = 128; % nominal matrix size
    arg.tr = 100; % repetition time (ms)
    arg.nrep = 1; % number of repetitions
    arg.nprj = 16; % number of thru-plane rotations (projections)
    arg.nint = 1; % number of in-plane rotations (interleaves)
    arg.nspokes = 23; % number of lps spokes
    arg.nseg = 280; % number of samples/segment
    arg.nrf = 3; % number of samples/rf pulse
    arg.fa = 4; % rf flip angle (deg)
    arg.gmax = 4; % max gradient amplitude (G/cm)
    arg.smax = 500; % max slew rate (G/cm/s)
    arg.dt = 4e-6; % raster time (s)
    arg.acq_delay = 1e-3; % acquisition delay (s)

    % parse inputs
    arg = vararg_pair(arg,varargin);

    % define constants/conversion factors
    gam = 4258; % GMR of H+ (Hz/G)
    gconv = 10; % gradient unit conversion factor (G/cm --> mT/m)
    sconv = 1e-2; % slew rate unit conversion factor (G/cm/s --> mT/m/ms)

    % set system limits
    sys = mr.opts('MaxGrad',arg.gmax*gconv, 'GradUnit', 'mT/m',...
        'MaxSlew', arg.smax*sconv, 'SlewUnit', 'mT/m/ms',...
        'rfDeadTime', 100e-6, ...
        'rfRingdownTime', 60e-6, ...
        'adcRasterTime', arg.dt, ...
        'gradRasterTime', arg.dt, ...
        'rfRasterTime', arg.dt, ...
        'blockDurationRaster', 4e-6, ...
        'B0', 3, ...
        'adcDeadTime', 0e-6);

    % initialize sequence
    seq = mr.Sequence(sys);
    warning('OFF', 'mr:restoreShape');
    warning('OFF', 'mr:makeArbitraryGrad');

    % create looping star waveforms
    [g,rf_wav] = gen_lps_waveforms( ...
        'fov', arg.fov, ... % fov (cm)
        'N', arg.N, ... % nominal matrix size
        'nspokes', arg.nspokes, ... % number of lps spokes
        'nseg', arg.nseg, ... % number of samples/segment
        'nrf', arg.nrf, ... % number of samples/rf pulse
        'fa', arg.fa, ... % rf flip angle (deg)
        'gmax', arg.gmax, ... % max gradient amplitude (G/cm)
        'smax', arg.smax, ... % max slew rate (G/cm/s)
        'dt', arg.dt ... % raster time (s))
        );
    % nextra = 4*ceil(size(g,1)/2 / 4) - size(g,1)/2; % make sure each half is len div by 4
    % G0 = padarray(G0,[nextra,0],0,'both');
    % rf_wav = padarray(rf_wav,[nextra,0],0,'both');

    G0 = padarray(g,[0,1],0,'post');

    % create rf (only take 1st half to play during block 1)
    rf = mr.makeArbitraryRf(rf_wav(1:end/2),arg.fa/180*pi, ...
        'delay', sys.rfDeadTime);

    % create ADC (no delay - only playing during block 2)
    acq_len = arg.nspokes*arg.nseg;
    adc = mr.makeAdc(acq_len,'Duration',arg.dt*acq_len, ...
        'Delay',arg.acq_delay,'system',sys);

    % calculate delay
    t_delay = arg.tr*1e-3 - arg.dt*size(G0,1);
    delay = mr.makeDelay(t_delay);

    % loop through reps, shots and rotations
    for repn = 1:arg.nrep
        for intn = 1:arg.nint
            for prjn = 1:arg.nprj

                % interleaf rotation: rotate by golden angle about z
                r_int = (intn-1)*pi*(3 - sqrt(5));
                R_int = eul2rotm([r_int,0,0],'ZYX');

                % projection rotation: rotate by 3D golden angles
                % ref: Generalization of three-dimensional golden-angle radial acquisition
                % to reduce eddy current artifacts in bSSFP CMR imaging (A, Fyrdahl et. al)
                phi1 = 0.4656; phi2 = 0.6823; % 3D golden ratios
                rp_prj = acos(mod((prjn-1) * phi1, 2)-1) + pi; % polar angle
                ra_prj = 2*pi*((prjn-1) * phi2); % azimuthal angle
                R_prj = eul2rotm([rp_prj,0,ra_prj],'ZYX');

                % rotate the gradients
                R = R_prj * R_int;
                iG = G0 * R.';

                % add first half (FID) gradients and rf to sequence
                gx_fid = mr.makeArbitraryGrad('x',0.99*iG(1:end/2,1),'system',sys);
                gy_fid = mr.makeArbitraryGrad('y',0.99*iG(1:end/2,2),'system',sys);
                gz_fid = mr.makeArbitraryGrad('z',0.99*iG(1:end/2,3),'system',sys);
                seq.addBlock(rf, gx_fid, gy_fid, gz_fid, ...
                    mr.makeLabel('SET', 'TRID', 1));

                % add second half (GRE) gradients and adc to sequence
                gx_gre = mr.makeArbitraryGrad('x',0.99*iG(end/2+1:end,1),'system',sys);
                gy_gre = mr.makeArbitraryGrad('y',0.99*iG(end/2+1:end,2),'system',sys);
                gz_gre = mr.makeArbitraryGrad('z',0.99*iG(end/2+1:end,3),'system',sys);
                seq.addBlock(adc, gx_gre, gy_gre, gz_gre);

                % add delay for TR
                seq.addBlock(delay);

            end
        end
    end

    % check whether the timing of the sequence is correct
    [ok, error_report] = seq.checkTiming;
    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end

    % write out the sequence
    seq.setDefinition('FOV', 1e-2*arg.fov*ones(1,3));
    seq.setDefinition('Name', 'lps');
    seq.write('lps.seq');

    % the sequence is ready, so let's see what we got
    seq.plot(); % plot sequence waveforms

end

