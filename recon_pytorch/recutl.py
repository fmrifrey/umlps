import torch
import numpy as np
from torch.fft import fftn, ifftn
from mirtorch.linear.util import fftshift, ifftshift

def resize_nd(x, dims, sfac=1):

    # fourier transform the data
    x = ifftshift(x, dims)
    k = fftn(x, dim=dims, norm='ortho')
    k = fftshift(k, dims)

    # zero pad or truncate the data in kspace along the specified dimensions
    if sfac > 1:
        padsize = list(np.zeros(k.ndim))
        for dim in dims:
            padsize[dim] = np.round((sfac - 1) * k.shape[dim] / 2)
        padsize = tuple(int(p) for p in padsize)
        pad = []
        for p in reversed(padsize):
            pad.extend([int(p), int(p)])
        k_resized = torch.nn.functional.pad(k, pad=pad, mode='constant', value=0)
    elif sfac < 1:
        new_size = list(k.shape)
        for dim in dims:
            new_size[dim] = np.round(k.shape[dim] * sfac)
        new_size = tuple(int(p) for p in new_size)
        k_resized = torch.zeros(new_size, dtype=k.dtype)
        k_resized = k_resized.to(k.device)
        slices = tuple(slice(0, new_size[dim]) for dim in range(k.ndim))
        k_resized[slices] = k[slices]
    else:
        k_resized = k

    # inverse fourier transform back
    k_resized = ifftshift(k_resized, dims)
    x_resized = ifftn(k_resized, dim=dims, norm='ortho')
    x_resized = fftshift(x_resized, dims)

    return x_resized

def mri_coil_compress(data, ncoil=None, Vr=None):

    # form column vectors of coil data
    data_cv = data.permute(0,*range(2,data.ndim),1).reshape(-1,data.shape[1])

    if Vr is None:
        # perform SVD to get principal component rotation
        _,_,Vt = torch.linalg.svd(data_cv, full_matrices=False)
        Vr = Vt.T

    if ncoil is None:
        # default to the number of coils in the original data
        ncoil = Vr.shape[1]
    
    # compress the data using the first ncoil columns of Vr
    data_comp_cv = data_cv @ Vr[:,:ncoil]  # (num_samples x ncoil)

    # reshape back to original dimensions with compressed coils
    data_comp = data_comp_cv.reshape(data.shape[0], -1, ncoil)
    data_comp = data_comp.permute(0, 2, 1)
    data_comp = data_comp.reshape(data.shape[0], ncoil, *data.shape[2:])  # maintain original shape except for coils

    return data_comp, Vr[:,:ncoil]

def sketched_svd(A, d):
    m, n = A.shape
    
    # generate a random Gaussian sketching matrix
    S = 1/d * torch.tensor(np.random.randn(n, d)).to(A.dtype)
    
    # compute the sketch
    B = A @ S
    
    # QR factorization
    Q, _ = torch.linalg.qr(B)
    
    # compute small SVD
    C = Q.T @ A  # (d x n)
    U_tilde, Sigma, Vt = torch.linalg.svd(C, full_matrices=True)
    
    # recover approximate left singular vectors
    U = Q @ U_tilde
    
    return U, Sigma, Vt