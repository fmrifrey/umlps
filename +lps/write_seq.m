function write_seq(varargin)
% Create the .seq file for looping star sequence
%
% by David Frey (djfrey@umich.edu)
%
% Inputs:
% tr - repetition time (ms)
% fov - field of view (cm)
% dt - raster time (s)
% N - 3D matrix size
% nrep - number of rotation sequence repetitions
% nint - number of interleaves (2D in-plane rotations)
% nprj - number of projections (3D thru-plane rotations)
% nspokes - number of looping star spokes
% nseg - number of samples/segment
% nrf - number of samples/rf hard pulse
% fa - flip angle (deg)
% gmax - max gradient amplitude (G/cm)
% smax - max slew rate (G/cm/s)
% plotseq - option to plot the sequence
% pislquant - number of TRs to use for prescan
% writepge - option to convert seq to pge file
%
% Outputs:
% lps.seq file - seq file for pulseq
% seq_args.mat - .mat file containing copy of input arguments
%

    % set default arguments
    arg.tr = 100; % repetition time (ms)
    arg.fov = 20; % fov (cm)
    arg.dt = 4e-6; % raster time (s)
    arg.N = 128; % 3D matrix size
    arg.nrep = 1; % number of rotation sequence repetitions
    arg.nint = 1; % number of interleaves (2D in-plane rotations)
    arg.nprj = 16; % number of projections (3D thru-plane rotations)
    arg.nspokes = 23; % number of lps spokes
    arg.nseg = 280; % number of samples/segment
    arg.nrf = 4; % number of samples/rf pulse
    arg.fa = 4; % rf flip angle (deg)
    arg.gmax = 4; % max gradient amplitude (G/cm)
    arg.smax = 500; % max slew rate (G/cm/s)
    arg.plotseq = false; % option to plot the sequence
    arg.pislquant = 1; % number of TRs to use for prescan
    arg.writepge = true; % option to convert seq to pge file
    
    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % set system limits
    sys = mr.opts('MaxGrad', arg.gmax*10, 'GradUnit', 'mT/m', ...
        'MaxSlew', arg.smax*1e-2, 'SlewUnit', 'mT/m/ms', ...
        'rfDeadTime', 100e-6, ...
        'rfRingdownTime', 60e-6, ...
        'adcRasterTime', arg.dt, ...
        'gradRasterTime', arg.dt, ...
        'rfRasterTime', arg.dt, ...
        'blockDurationRaster', 4e-6, ...
        'B0', 3, ...
        'adcDeadTime', 0e-6);
    
    % initialize sequence object
    seq = mr.Sequence(sys);
    warning('OFF', 'mr:restoreShape');

    % create looping star waveforms
    [g_wav,rf_wav,rf_del] = lps.seq.gen_lps_waveforms( ...
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
    g_wav = sys.gamma * 1e-2 * g_wav; % kHz/m
    G0 = padarray(g_wav.',[1,0],0,'post');

    % create rf (only take 1st half to play during block 1)
    rf = mr.makeArbitraryRf(rf_wav, arg.nspokes*arg.fa*pi/180, ...
        'delay', arg.dt*rf_del, ...
        'use', 'excitation', ...
        'system', sys);

    % create ADC
    acq_len = sys.adcSamplesDivisor*ceil(arg.nspokes*arg.nseg ...
        /sys.adcSamplesDivisor);
    adc = mr.makeAdc(acq_len, ...
        'Duration', sys.adcRasterTime*acq_len, ...
        'system', sys);

    % create TR delay
    tr_delay = mr.makeDelay(arg.tr*1e-3 - arg.dt*size(G0,2));
    
    % define sequence blocks
    for irep = 1:arg.nrep
        for iprj = 1:arg.nprj
            for iint = 1:arg.nint
        
                % rotate the gradients based on a 3DTGA rotation sequence
                R = lps.seq.rot_3dtga(iprj, iint);
                iG = R * G0;
        
                % write fID portion to sequence
                gx_fid = mr.makeArbitraryGrad('x', 0.99*iG(1,1:end/2), ...
                    'system', sys, ...
                    'first', 0.99*iG(1,1), ...
                    'last', 0.99*iG(1,end/2));
                gy_fid = mr.makeArbitraryGrad('y', 0.99*iG(2,1:end/2), ...
                    'system', sys, ...
                    'first', 0.99*iG(2,1), ...
                    'last', 0.99*iG(2,end/2));
                gz_fid = mr.makeArbitraryGrad('z', 0.99*iG(3,1:end/2), ...
                    'system', sys, ...
                    'first', 0.99*iG(3,1), ...
                    'last', 0.99*iG(3,end/2));
                seq.addBlock(rf, gx_fid, gy_fid, gz_fid, ...
                    mr.makeLabel('SET', 'TRID', 1));

                % write GRE portion to sequence
                gx_gre = mr.makeArbitraryGrad('x', 0.99*iG(1,end/2+1:end), ...
                    'system', sys, ...
                    'first', 0.99*iG(1,end/2+1), ...
                    'last', 0.99*iG(1,end));
                gy_gre = mr.makeArbitraryGrad('y', 0.99*iG(2,end/2+1:end), ...
                    'system', sys, ...
                    'first', 0.99*iG(2,end/2+1), ...
                    'last', 0.99*iG(2,end));
                gz_gre = mr.makeArbitraryGrad('z', 0.99*iG(3,end/2+1:end), ...
                    'system', sys, ...
                    'first', 0.99*iG(3,end/2+1), ...
                    'last', 0.99*iG(3,end));
                seq.addBlock(adc, gx_gre, gy_gre, gz_gre);
        
                % add tr delay
                seq.addBlock(tr_delay);
        
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
    
    % write out sequence and save args
    seq.setDefinition('FOV', arg.fov*1e-2*ones(1,3));
    seq.setDefinition('Name', 'lps');
    seq.write('lps.seq');
    if arg.writepge
        ceq = seq2ceq('lps.seq');
        writeceq(ceq, 'lps.pge', 'pislquant', arg.pislquant);
    end
    save seq_args.mat -struct arg
    
    % the sequence is ready, so let's see what we got
    if arg.plotseq
        figure
        seq.plot();
    end

end