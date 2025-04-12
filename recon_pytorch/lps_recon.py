# batch script for MIRTorch based Looping Star reconstruction

# import packages
import h5py
import numpy as np
import torch
from mirtorch.linear import NuSense, Diff3dgram, Diag
from mirtorch.alg.cg import CG
import os
import sys
from recutl import mri_coil_compress, resize_nd
import argparse

# parse the input arguments
parser = argparse.ArgumentParser(
    description="Looping star reconstruction code in PyTorch"
)
parser.add_argument("--basedir", required=False, type=str, default="./")
parser.add_argument("--fname_kdata", required=False, type=str, default="lps.h5")
parser.add_argument("--fname_smaps", required=False, type=str, default="smaps.h5")
parser.add_argument("--fname_out", required=False, type=str, default="recon.h5")
parser.add_argument("--device", required=False, type=str, default="cpu")
parser.add_argument("--ncoil_comp", required=False, type=int, default=8)
parser.add_argument("--cutoff", required=False, type=float, default=0.8)
parser.add_argument("--lam", required=False, type=float, default=20)
parser.add_argument("--niter", required=False, type=int, default=20)
parser.add_argument("--ints2use", required=False, type=int, default=None)
parser.add_argument("--prjs2use", required=False, type=int, default=None)
parser.add_argument("--reps2use", required=False, type=int, default=None)
parser.add_argument("--volwidth", required=False, type=int, default=None)
args = parser.parse_args()

# select device
if args.device == "cpu" or torch.cuda.is_available():
    device0 = torch.device(args.device)
else:
    raise RuntimeError("GPU not available!")
    
print(f'using device {device0}', flush=True)

# load in the data
print(f'loading data from {os.path.join(args.basedir,args.fname_kdata)}', flush=True)
with h5py.File(os.path.join(args.basedir,args.fname_kdata), 'r') as h5_file:
    kdata = h5_file['kdata/real'][:] + 1j * h5_file['kdata/imag'][:] # kspace data
    k_in = h5_file['ktraj/spoke_in'][:] # kspace spoke-in trajectory (1/cm)
    k_out = h5_file['ktraj/spoke_out'][:] # kspace spoke-out trajectory (1/cm)
    fov = h5_file['seq_args/fov'][0][0] # field of view (cm)
    tr = h5_file['seq_args/tr'][0][0] # repetition time (ms)
    fa = h5_file['seq_args/fa'][0][0] # flip angle (deg)
    dummyshots = int(h5_file['seq_args/dummyshots'][0][0]) # number of dummy shots
    gmax = h5_file['seq_args/gmax'][0][0] # max gradient amplitude (G/cm)
    smax = h5_file['seq_args/smax'][0][0] # max slew rate (G/cm/s)
    pislquant = int(h5_file['seq_args/pislquant'][0][0]) # pislquant
    nrf = int(h5_file['seq_args/nrf'][0][0]) # number of samples/rf pulse
    N = int(h5_file['seq_args/N'][0][0]) # 3D matrix size
    nseg = int(h5_file['seq_args/nseg'][0][0]) # number of points per segment
    nspokes = int(h5_file['seq_args/nspokes'][0][0]) # number of spokes
    nprj = int(h5_file['seq_args/nprj'][0][0]) # number of projections
    nint = int(h5_file['seq_args/nint'][0][0]) # number of interleaves
    nrep = int(h5_file['seq_args/nrep'][0][0]) # number of repetitions
    ncoil = int(h5_file['ncoil'][0][0]) # number of coils

# convert to tensors
kdata = torch.tensor(kdata).reshape(ncoil,nrep,nprj,nint,nseg*nspokes)
k_in = torch.tensor(k_in).reshape(3,nrep,nprj,nint,nseg*nspokes)
k_out = torch.tensor(k_out).reshape(3,nrep,nprj,nint,nseg*nspokes)

# set default values for interleaves, projections, and repetitions to use
if args.ints2use is None:
    args.ints2use = nint
if args.prjs2use is None:
    args.prjs2use = nprj
if args.reps2use is None:
    args.reps2use = nrep
if args.volwidth is None:
    args.volwidth = args.prjs2use*args.ints2use

# get number of volumes
nvol = args.reps2use*args.prjs2use*args.ints2use // args.volwidth

# arrange data into volumes
print(f'arranging {args.reps2use*args.prjs2use*args.ints2use} TRs into {nvol} sampling volumes', flush=True)
kdata2 = kdata.clone()
kdata2 = kdata2[:,np.arange(args.reps2use),:,:,:]
kdata2 = kdata2[:,:,np.arange(args.prjs2use),:,:]
kdata2 = kdata2[:,:,:,np.arange(args.ints2use),:]
kdata2 = kdata2.reshape(ncoil,args.reps2use*args.prjs2use*args.ints2use,nseg*nspokes)
kdata2 = kdata2[:,:nvol*args.volwidth,:]
kdata2 = kdata2.reshape(ncoil,nvol,args.volwidth*nseg*nspokes)
kdata2 = kdata2.permute(1,0,2)

# arrange kspace trajectory into volumes
k_in2 = k_in.clone()
k_in2 = k_in2[:,np.arange(args.reps2use),:,:,:]
k_in2 = k_in2[:,:,np.arange(args.prjs2use),:,:]
k_in2 = k_in2[:,:,:,np.arange(args.ints2use),:]
k_in2 = k_in2.reshape(3,args.reps2use*args.prjs2use*args.ints2use,nseg*nspokes)
k_in2 = k_in2[:,:nvol*args.volwidth,:]
k_in2 = k_in2.reshape(3,nvol,args.volwidth*nseg*nspokes)
k_in2 = k_in2.permute(1,0,2)

k_out2 = k_out.clone()
k_out2 = k_out2[:,np.arange(args.reps2use),:,:,:]
k_out2 = k_out2[:,:,np.arange(args.prjs2use),:,:]
k_out2 = k_out2[:,:,:,np.arange(args.ints2use),:]
k_out2 = k_out2.reshape(3,args.reps2use*args.prjs2use*args.ints2use,nseg*nspokes)
k_out2 = k_out2[:,:nvol*args.volwidth,:]
k_out2 = k_out2.reshape(3,nvol,args.volwidth*nseg*nspokes)
k_out2 = k_out2.permute(1,0,2)

# load in the sensitivity maps
print(f'loading sensitivity maps from {os.path.join(args.basedir,args.fname_smaps)}', flush=True)
with h5py.File(os.path.join(args.basedir,args.fname_smaps), 'r') as h5_file:
    smaps = torch.tensor(h5_file['/real'][:] + 1j * h5_file['/imag'][:]).unsqueeze(0).to(kdata)
smaps = resize_nd(smaps, (2,3,4), N/smaps.shape[2]) # resize to match kspace data
smaps = smaps.permute(0,1,4,3,2).repeat(nvol,1,1,1,1) # t x C x X x Y x Z

# coil compress the data
print(f'compressing to {ncoil} coils to {args.ncoil_comp} virtual coils', flush=True)
kdata_comp,Vr = mri_coil_compress(kdata2, ncoil=args.ncoil_comp)

# coil compress the sensitivity maps
smaps_comp,_ = mri_coil_compress(smaps, Vr=Vr)

# convert trajectory to spatial frequencies
om_in = 2*torch.pi * fov/N * k_in2
om_out = 2*torch.pi * fov/N * k_out2

# create filter objects
Hvec_in = 1*(torch.norm(om_in,2,dim=1,keepdim=True) <= args.cutoff*torch.pi).repeat(1,args.ncoil_comp,1)
Hvec_out = 1*(torch.norm(om_out,2,dim=1,keepdim=True) <= args.cutoff*torch.pi).repeat(1,args.ncoil_comp,1)
H_in = Diag(Hvec_in.to(device0))
H_out = Diag(Hvec_out.to(device0))

# create nufft system operators with flat sensitivity
FS_in = NuSense(smaps_comp.to(device0), om_in.to(device0))
FS_out = NuSense(smaps_comp.to(device0), om_out.to(device0))

# set up system matrices and data
print('setting up system matrices', flush=True)
A = H_in*FS_in + H_out*FS_out
AHA = A.H*A

# add L2 roughness penalty
THT = Diff3dgram(FS_in.size_in)
AHA_tikh = AHA + args.lam*THT

# set up data
y = kdata_comp.to(device0)
AHy = A.H * y

# set up the CG solver
solv = CG(AHA_tikh, max_iter=args.niter)

# solve with CG
print('solving with CG', flush=True)
x = solv.run(torch.zeros(nvol,1,N,N,N).to(kdata).to(device0), AHy)

# save the reconstruction
x = x.cpu().detach().numpy()
print(f'saving reconstruction to {os.path.join(args.basedir,args.fname_out)}', flush=True)
with h5py.File(os.path.join(args.basedir,args.fname_out), 'w') as h5_file:
    # save solution
    sol = h5_file.create_group('sol')
    sol.create_dataset('real', data=x.real.transpose(0, 1, 4, 3, 2))
    sol.create_dataset('imag', data=x.imag.transpose(0, 1, 4, 3, 2))

    # save sequence args
    seq_args = h5_file.create_group('seq_args')
    seq_args.create_dataset('fov', data=fov)
    seq_args.create_dataset('tr', data=tr)
    seq_args.create_dataset('fa', data=fa)
    seq_args.create_dataset('dummyshots', data=dummyshots)
    seq_args.create_dataset('gmax', data=gmax)
    seq_args.create_dataset('smax', data=smax)
    seq_args.create_dataset('nrf', data=nrf)
    seq_args.create_dataset('pislquant', data=pislquant)
    seq_args.create_dataset('N', data=N)
    seq_args.create_dataset('nseg', data=nseg)
    seq_args.create_dataset('nspokes', data=nspokes)
    seq_args.create_dataset('nint', data=nint)
    seq_args.create_dataset('nprj', data=nprj)
    seq_args.create_dataset('nrep', data=nrep)

    # save recon args
    recon_args = h5_file.create_group('recon_args')
    recon_args.create_dataset("basedir", data=args.basedir)
    recon_args.create_dataset("fname_kdata", data=args.fname_kdata)
    recon_args.create_dataset("fname_smaps", data=args.fname_smaps)
    recon_args.create_dataset("fname_out", data=args.fname_out)
    recon_args.create_dataset("device", data=args.device)
    recon_args.create_dataset("ncoil_comp", data=args.ncoil_comp)
    recon_args.create_dataset("cutoff", data=args.cutoff)
    recon_args.create_dataset("lam", data=args.lam)
    recon_args.create_dataset("niter", data=args.niter)
    recon_args.create_dataset("ints2use", data=args.ints2use)
    recon_args.create_dataset("prjs2use", data=args.prjs2use)
    recon_args.create_dataset("reps2use", data=args.reps2use)
    recon_args.create_dataset("volwidth", data=args.volwidth)    
