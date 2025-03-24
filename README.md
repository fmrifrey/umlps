# Looping Star MRI in Pulseq
This repository contains code for generating a Looping Star sequence in Pulseq (new pge2 interpreter)

by David Frey (djfrey@umich.edu)

## File Structure
- `lps_write_seq.m`: Function for writing LPS pulseq files
- `lps_update_packages.m`: Script to update required external packages
- `lps_convert_data.m`: Script to export data for external reconstruction
- `lps_recon_example.m`: Example script for reconstructing LPS data in MATLAB
- `+lps/`: Contains core helper functions for LPS sequence and reconstruction

## Dependencies
This repository relies on the following external packages:
- [Pulseq](https://github.com/pulseq/pulseq)
- [PulCeq (tv7 branch)](https://github.com/HarmonizedMRI/PulCeq/tree/tv7)
- [toppe (develop branch)](https://github.com/toppeMRI/toppe/tree/develop)
- [MIRT](https://github.com/JeffFessler/MIRT)
- [Orchestra SDK 2.1 for MATLAB](https://weconnect.gehealthcare.com/s/feed/0D53a00008pQ1Q8CAK) (requires GE SDK access)

Run `lps_update_packages.m` to download and set up the required packages.