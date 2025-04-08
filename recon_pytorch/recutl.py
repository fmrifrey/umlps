import torch
import numpy as np

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