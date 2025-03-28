import numpy as np
import matplotlib.pyplot as plt

def im(img, viewtype='montage', rows=None, cols=None, offset=None):
    """
    Display a 3D image in different modes: montage, mid-slice, or maximum intensity projection.
    
    Parameters:
    img (numpy.ndarray): 3D image to display.
    viewtype (str): Type of view ('montage', 'mid3', 'mip3').
    rows (int): Number of rows for montage view.
    cols (int): Number of columns for montage view.
    offset (tuple): Offset for mid-slice view.
    
    Returns:
    None
    """
    
    # get image size
    Nx,Ny,Nz = img.shape

    if viewtype == 'montage':
    # 3D montage (lightbox) mode

        if rows is None:
            rows = int(np.ceil(np.sqrt(Nz)))
        if cols is None:
            cols = int(np.ceil(Nz/rows))

        # create the montage
        lightbox = np.zeros((rows*Nx,cols*Ny))
        for i in range(rows):
            for j in range(cols):
                if i*cols+j < Nz:
                    lightbox[i*Nx:(i+1)*Nx,j*Ny:(j+1)*Ny] = img[:,:,i*cols+j]
        
        # show the montage
        plt.imshow(lightbox, cmap='gray')
        plt.axis('off')
        plt.show()

    elif viewtype == 'mid3':
    # 3D mid-slice mode

        if offset is None:
            offset = (0,0,0)

        # get the slices
        xn = int(np.floor(Nx/2)) + offset[0]
        yn = int(np.floor(Ny/2)) + offset[1]
        zn = int(np.floor(Nz/2)) + offset[2]

        # check the slices
        if xn < 0 or xn >= Nx:
            raise ValueError('x slice out of bounds')
        if yn < 0 or yn >= Ny:
            raise ValueError('y slice out of bounds')
        if zn < 0 or zn >= Nz:
            raise ValueError('z slice out of bounds')
        
        # show the slices
        plt.subplot(131)
        plt.imshow(img[xn,:,:], cmap='gray')
        plt.axis('off')
        plt.title('x slice')
        plt.subplot(132)
        plt.imshow(img[:,yn,:], cmap='gray')
        plt.axis('off')
        plt.title('y slice')
        plt.subplot(133)
        plt.imshow(img[:,:,zn], cmap='gray')
        plt.axis('off')
        plt.title('z slice')
        plt.show()

    elif viewtype == 'mip3':
    # 3D maximum intensity projection mode
        
        # show the projections
        plt.subplot(131)
        plt.imshow(np.max(img,dim=0), cmap='gray')
        plt.axis('off')
        plt.title('x projection')
        plt.subplot(132)
        plt.imshow(np.max(img,dim=1), cmap='gray')
        plt.axis('off')
        plt.title('y projection')
        plt.subplot(133)
        plt.imshow(np.max(img,dim=2), cmap='gray')
        plt.axis('off')
        plt.title('z projection')
        plt.show()
