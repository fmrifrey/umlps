% actions
createSequenceFile = 1;

if createSequenceFile

    % write lps.seq
    system('git clone git@github.com:toppeMRI/toppe.git');
    addpath toppe
    system('git clone --branch v1.5.0 git@github.com:pulseq/pulseq.git');
    addpath pulseq/matlab
    system('git clone --branch tv7 git@github.com:HarmonizedMRI/PulCeq.git');
    addpath PulCeq/matlab
    write_lps_seq;

    ceq = seq2ceq('lps.seq');
    pislquant = 1;   % number of ADC events used for receive gain calibration
    writeceq(ceq, 'lps.pge', 'pislquant', pislquant);
    
end