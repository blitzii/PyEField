import numpy as np


def coarsen_grid(phi, mask, epsilon):
    """
    Original MATLAB Docstring:
    Coarsens the current grid by a factor of two, symmetrically about the
    centre-line in the x-direction, and always starting from the base in the
    z-direction.  This routine is called from FieldSolverMG and GridCheck.
    """

    # Store size of (original) finer grid
    zsize_h, xsize_h = phi.shape

    # Find centre line in the x; furthest points in z and x
    """ Original MATLAB:
    cx_h = (xsize_h-1)/2+1;
    zfar_2h = 2*floor((zsize_h-1)/2)+1;
    xfar_2h = 2*floor((xsize_h-1)/4);
    """
    cx_h = (xsize_h)/2
    zfar_2h = 2*np.floor(zsize_h/2)
    xfar_2h = 2*np.floor(xsize_h/4)

    # Pre-allocate arrays
    """ Original MATLAB:
    phi_2h = zeros(size(1:2:zfar_2h,2),size(cx_h+((-xfar_2h):2:xfar_2h),2));
    mask_2h = zeros(size(1:2:zfar_2h,2),size(cx_h+((-xfar_2h):2:xfar_2h),2));

    This seems horribly inefficient as it creates an temp array in order to
    calculate the size of the new array.
    """
    phi_2h = np.zeros((np.arange(0, zfar_2h, 2).size,
                      np.arange(-xfar_2h, xfar_2h, 2).size + cx_h))
    mask_2h = np.zeros((np.arange(0, zfar_2h, 2).size,
                       np.arange(-xfar_2h, xfar_2h, 2).size + cx_h))

    # Copy across the phi and mask data from the finer grid to the coarser grid
    # Start at the first line:

    # phi_2h[0, :] = phi[0, cx_h + np.arange(-xfar_2h, xfar_2h, 2)]
    # mask_2h[0, :] = mask[0, cx_h + np.arange(-xfar_2h, xfar_2h, 2)]

    phi_2h[0, :] = phi[0, cx_h - xfar_2h:cx_h + xfar_2h:2]
    mask_2h[0, :] = mask[0, cx_h - xfar_2h:cx_h + xfar_2h:2]

    j = 0
    # Now do rest of lines in turn
    for i in range(1, zfar_2h, 2):
        # phi_2h[j, :] = phi[i, cx_h + np.arange(-xfar_2h, xfar_2h, 2)]
        # mask_2h[j, :] = mask[i, cx_h + np.arange(-xfar_2h, xfar_2h, 2)]

        phi_2h[j, :] = phi[i, cx_h - xfar_2h:cx_h + xfar_2h:2]
        mask_2h[j, :] = mask[i,  cx_h - xfar_2h:cx_h + xfar_2h:2]
        # Check whether we have missed any important information
        # (i.e. an electrode), and add any required info to new grid
        # Dave:
        # I think this checks if there was an electrode in the skipped line,
        # then adds it to the current one if that's the case.
        if np.count_nonzero(mask[i-1, cx_h - xfar_2h:cx_h + xfar_2h:2]):
            phi_2h[j, :] = phi[i-1, cx_h - xfar_2h:cx_h + xfar_2h:2]
            mask_2h[j, :] = mask[i-1, cx_h - xfar_2h:cx_h + xfar_2h:2]
        j = j + 1

    # Ensure boundary conditions retained at the top and bottom
    phi_2h[-1, :] = phi[-1, cx_h - xfar_2h:cx_h + xfar_2h:2]
    mask_2h[-1, :] = mask[-1, cx_h - xfar_2h:cx_h + xfar_2h:2]

    # Pre-allocate the permittivity array
    epsil = np.ones(zsize_h, xsize_h)

    # Interpolate the input permittivities onto the potential grid at the finer
    # grid level
    epsil[1:-1, 1:-1] = (epsilon[1:, 0:-1] + epsilon[0:-1, 0:-1] +
                         epsilon[1:, 1:] + epsilon[0:-1, 1:])/4

    epsil[0, :] = epsil[1, :]
    epsil[-1, :] = epsil[-2, :]
    epsil[:, 0] = epsil[:, 1]
    epsil[:, end] = epsil[:, -2]

    # Copy onto coarser grid
    """ MATLAB:
    epsilin_2h = epsil(2:2:(zfar_2h-1),cx_h+((-xfar_2h+1):2:(xfar_2h-1)));
    """
    # epsilon_2h = epsil[np.arange(1, zfar_2h-1, 2),
    #                    cx_h + np.arange(-xfar_2h+1, xfar_2h-1, 2)]
    epsilon_2h = epsil[1:zfar_2h - 1:2,
                       cx_h - xfar_2h + 1:cx_h + xfar_2h - 1:2]

    return phi_2h, mask_2h, epsilon_2h
