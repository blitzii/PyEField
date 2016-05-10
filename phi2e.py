def phi2E(phi):
  zsize, xsize = phi.shape
  # Returned arrays appear to be different order to MATLAB
  # Check Z and X are correct
  Ez, Ex = np.gradient(phi)
  Ex = -Ex
  Ez = -Ez

# Manually calculate gradient at x-boundaries
# This boundary is mirrored
  for i in range(1, zsize-1):
    Ex[i][0] = (phi[i][xsize-1] - phi[i][1])/2
    Ex[i, xsize-1] = (phi[i][xsize-2] - phi[i][0])/2

  return Ex, Ez
