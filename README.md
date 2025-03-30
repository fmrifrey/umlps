# Looping Star (LpS) MRI in Pulseq
This repository contains code for generating a Looping Star sequence in Pulseq (new pge2 interpreter)

by David Frey (djfrey@umich.edu)

## File Structure
- `psd`: Contains code for writing pulse sequence code in Pulseq, and reconstruction in MATLAB
  - `lps_write_seq.m`: Function for writing LpS pulseq files
  - `lps_convert_data.m`: Function to read LpS data from scanarchive and export to h5 for external reconstruction
  - `gre3d_write_seq.m`: Function for writing 3d GRE pulseq files
  - `gre3d_convert_data.m`: Function to read3d GRE data from scanarchive and export to h5 for external reconstruction
  - `umlps_update_packages.m`: Script to update required external packages
  - `umlps_recon_example.m`: Example script for reconstructing LpS data in MATLAB
  - `+lpsutl/`: Contains core helper functions for LpS sequence and reconstruction
- `recon`: Contains code for reconstructing the data in MIRTorch (in development)

## MATLAB Dependencies
This repository relies on the following external packages for MATLAB based development:
- [Pulseq](https://github.com/pulseq/pulseq)
- [PulCeq (tv7 branch)](https://github.com/HarmonizedMRI/PulCeq/tree/tv7)
- [toppe (develop branch)](https://github.com/toppeMRI/toppe/tree/develop)
- [MIRT](https://github.com/JeffFessler/MIRT)
- [Orchestra SDK 2.1 for MATLAB](https://weconnect.gehealthcare.com/s/feed/0D53a00008pQ1Q8CAK) (requires GE SDK access)
Run `umlps_update_packages.m` to download and set up the required packages.
