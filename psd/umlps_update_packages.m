% gets the current compatible packages for running code from this repo
% by David Frey (djfrey@umich.edu)

%% Pulseq
fprintf('updating pulseq... ')
system('[ -d "./pulseq" ] && rm -rf ./pulseq');
system('git clone --branch v1.5.0 git@github.com:pulseq/pulseq.git 2> /dev/null');
addpath pulseq/matlab
fprintf('done.\n')

%% PulCeq
fprintf('updating PulCeq... ')
system('[ -d "./PulCeq" ] && rm -rf ./PulCeq');
system('git clone --branch tv7 git@github.com:HarmonizedMRI/PulCeq.git 2> /dev/null');
addpath PulCeq/matlab
fprintf('done.\n')

%% toppe
fprintf('updating toppe... ')
system('[ -d "./toppe" ] && rm -rf ./toppe');
system('git clone --branch develop git@github.com:fmrifrey/toppe.git 2> /dev/null');
addpath toppe
fprintf('done.\n')

%% MIRT
fprintf('updating MIRT... ')
system('[ -d "./MIRT" ] && rm -rf ./MIRT');
system('git clone git@github.com:JeffFessler/MIRT.git 2> /dev/null');
fprintf('done.\n')
run MIRT/setup.m

%% BART
fprintf('updating bart... ')
system('[ -d "./bart" ] && rm -rf ./bart');
system('git clone git@github.com:mrirecon/bart.git 2> /dev/null');
fprintf('done.\n')
run bart/startup

%% Orchestra
addpath ~/code/packages/orchestra-sdk-2.1-1.matlab/ % replace with path to orchestra