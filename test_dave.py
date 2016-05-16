import numpy as np
from coarsen_grid import coarsen_grid


test_1_phi = np.loadtxt('./test_data/test_1_phi.csv', delimiter=',')
test_1_mask = np.loadtxt('./test_data/test_1_mask.csv', delimiter=',')
test_1_perm = np.loadtxt('./test_data/test_1_perm.csv', delimiter=',')

rar = coarsen_grid(test_1_phi, test_1_mask, test_1_perm)
