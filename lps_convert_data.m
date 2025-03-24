% get formatted data and kspace trajectory and save it to h5 file for
% external reconstruction
% by David Frey (djfrey@umich.edu)

% read in data
safile = './scanarc.h5'; % replace with scan archive file name
[kdata,k_in,k_out,seq_args] = lps.get_data(safile);

% save h5 file
h5file = './lps_data.h5';
if isfile(h5file)
    system(sprintf('rm %s',h5file));
end

% save kspace data
h5create(h5file, '/kdata/real', size(kdata), ...
    'Datatype', class(real(kdata)));
h5write(h5file, '/kdata/real', real(kdata));
h5create(h5file, '/kdata/imag', size(kdata), ...
    'Datatype', class(imag(kdata)));
h5write(h5file, '/kdata/imag', imag(kdata));

% save sampling locations
h5create(h5file, '/ktraj/spoke_in', size(k_in), 'Datatype', class(k_in));
h5write(h5file, '/ktraj/spoke_in', k_in);
h5create(h5file, '/ktraj/spoke_out', size(k_out), 'Datatype', class(k_out));
h5write(h5file, '/ktraj/spoke_out', k_out);

% save sequence arguments
seq_args_fields = fieldnames(seq_args);
for i = 1:numel(seq_args_fields)
    field = seq_args_fields{i};
    val = seq_args.(field);
    if islogical(val)
        val = 1*val;
    end
    h5create(h5file, sprintf('/seq_args/%s',field), size(val), ...
        'Datatype', class(val));
    h5write(h5file, sprintf('/seq_args/%s',field), val)
end
