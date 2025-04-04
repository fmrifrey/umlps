# Looping Star (LpS) MRI in Pulseq
This repository contains code for generating a Looping Star sequence in Pulseq (new pge2 interpreter)

by David Frey (djfrey@umich.edu)

## File structure
- `psd_pulseq`: Contains code for pulse sequence development in Pulseq, and converting acquired data to formatted h5 files
  - `lps_write_seq.m`: Function for writing LpS pulseq files
  - `lps_convert_data.m`: Function to read LpS data from scanarchive and export to h5 for external reconstruction
  - `gre3d_write_seq.m`: Function for writing 3D GRE pulseq files
  - `gre3d_convert_data.m`: Function to read 3D GRE data from scanarchive and export to h5 for external reconstruction
  - `update_psd_packages.m`: Script to update required external packages for pulse sequence development
  - `psd_example.m`: Example script for writing sequence files and converting the data
  - `+psdutl/`: Contains core helper functions for MATLAB based pulse sequence development
- `recon_matlab`: Contains code for reconstructing the data in MATLAB
  - `lps_recon.m`: Example script for reconstructing LpS data with CG-SENSE
  - `gre3d_recon.m`: Example script for creating SENSE maps from 3D GRE data and reconstructing 3D GRE data with CG-SENSE
  - `update_recon_packages.m`: Script to update required external packages for pulse sequence development
  - `+psdutl/`: Contains core helper functions for MATLAB based image reconstruction
- `recon_pytorch`: Contains code for reconstructing the data in PyTorch
  (under development)

## PSD dependencies
This repository relies on the following external packages for MATLAB based pulse sequence development:
- [Pulseq](https://github.com/pulseq/pulseq)
- [PulCeq (tv7 branch)](https://github.com/HarmonizedMRI/PulCeq/tree/tv7)
- [toppe (develop branch)](https://github.com/toppeMRI/toppe/tree/develop)
- [MIRT](https://github.com/JeffFessler/MIRT)
- [Orchestra SDK 2.1 for MATLAB](https://weconnect.gehealthcare.com/s/feed/0D53a00008pQ1Q8CAK) (requires GE SDK access)
The packages can be automatically downloaded and updated in the current directory by running `update_psd_packages.m`

## MATLAB recon dependencies
This repository relies on the following external packages for MATLAB based image reconstruction:
- [MIRT](https://github.com/JeffFessler/MIRT)
- [BART](https://github.com/mrirecon/bart)
The packages can be automatically downloaded and updated in the current directory by running `update_recon_packages.m`

## PyTorch recon dependencies
This repository relies on the following external packages for PyTorch based image reconstruction:
- [BART](https://github.com/mrirecon/bart)
Dependencies can be handled using a package manager such as Conda with the .yml file
