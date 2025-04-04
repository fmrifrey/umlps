import torch
import numpy as np

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