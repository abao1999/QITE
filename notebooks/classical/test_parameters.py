import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm 
from binary_functions import Int2Bas,Bas2Int,Opp2Str,Str2Opp
from Pmn import PPmunu

t_hubbard = 1
U_hubbard = 2
beta = 0.5
sigma = np.zeros([2, 2, 4], dtype=complex)
sigma[0, 0, 0] = 1.
sigma[1, 1, 0] = 1.
sigma[0, 1, 1] = 1.
sigma[1, 0, 1] = 1.
sigma[0, 1, 2] = -1j
sigma[1, 0, 2] = 1j
sigma[0, 0, 3] = 1.
sigma[1, 1, 3] = -1.
N = 10 # Number of imaginary time steps

# Defining the Hamiltonian
Hamiltonian = -t_hubbard * (np.kron(sigma[:, :, 1], sigma[:, :, 0]) + np.kron(sigma[:, :, 0], sigma[:, :, 1])) \
             + U_hubbard / 2 * (np.kron(sigma[: ,:, 0], sigma[: ,:, 0]) + np.kron(sigma[:, :, 3], sigma[:, :, 3]))


# Initialize into |00> state
psi = np.array([1,0,0,0],dtype=complex)
psi = psi/np.linalg.norm(psi)

phi = psi

# Initial expectation values
# Value of <XI>
XI = np.kron(sigma[:, :, 1], sigma[:, :, 0])
tmp = np.matmul(XI, phi)
tmp = np.real(np.matmul(np.transpose(np.conj(phi)), tmp))
print('Expectation value of <XI> is', tmp)

# Value of <IX>
IX = np.kron(sigma[:, :, 0], sigma[:, :, 1])
tmp = np.matmul(IX, phi)
tmp = np.real(np.matmul(np.transpose(np.conj(phi)), tmp))
print('Expectation value of <IX> is', tmp)

# Value of <II>
II = np.kron(sigma[:, :, 0], sigma[:, :, 0])
tmp = np.matmul(II, phi)
tmp = np.real(np.matmul(np.transpose(np.conj(phi)), tmp))
print('Expectation value of <II> is', tmp)

# Value of <ZZ>
ZZ = np.kron(sigma[:, :, 3], sigma[:, :, 3])
tmp = np.matmul(ZZ, phi)
tmp = np.real(np.matmul(np.transpose(np.conj(phi)), tmp))
print('Expectation value of <ZZ> is', tmp)

## Value of <ZY>
#ZY = np.kron(sigma[:, :, 3], sigma[:, :, 2])
#tmp = np.matmul(ZY, phi)
#tmp = np.real(np.matmul(np.transpose(np.conj(phi)), tmp))
#print('Expectation value of <ZY> is', tmp)


