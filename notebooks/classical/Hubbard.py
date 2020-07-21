import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm 
from binary_functions import Int2Bas,Bas2Int,Opp2Str,Str2Opp
from Pmn import PPmunu

t_hubbard = 1
U_hubbard = 2
beta = 1.0
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
w, v = np.linalg.eigh(Hamiltonian)
print('The exact solution is ', w[0])

# Classical trajectory

# List to store energy values
energy_classical_list = []
energy_qite_list = []

# Initialize into |00> state
psi = np.array([1,0,0,0],dtype=complex)
psi = psi/np.linalg.norm(psi)

db = beta/N
expmH = expm(-db * Hamiltonian)
phi = psi

# Store the energy for initial wavefunction
e = np.matmul(Hamiltonian, phi)
e = np.real(np.matmul(np.transpose(np.conj(phi)), e))
energy_classical_list.append(e)

# Perform imaginary time evolution
for i in range(0,N):

	# Perform imaginary time evolution
	phi = np.matmul(expmH, phi)

	# Require normalization because imaginary time evolution is non-unitary
	phi = phi / np.linalg.norm(phi)

	# Store the energy values at each step
	e = np.matmul(Hamiltonian, phi)
	e = np.real(np.matmul(np.transpose(np.conj(phi)), e))
	energy_classical_list.append(e)
print('Final energy at beta', beta, 'is ', e)

# Qite approximation

# First populate the Lie algebra rules
index = np.zeros([4**2,4**2],dtype=int)
coeff = np.zeros([4**2,4**2],dtype=complex)
	
row = 0
for i in range(4**2):
	column = 0
	for j in range(4**2):
		Pnu = Opp2Str(Int2Bas(column,4,2))
		Pmu = Opp2Str(Int2Bas(row,4,2))
		A = Pmu[0] + Pnu[0]
		B = Pmu[1] + Pnu[1]
		A,intA = PPmunu(A)
		B,intB = PPmunu(B)
		index[i,j] = Bas2Int(Str2Opp(A+B),4)
		coeff[i,j] = intA*intB
		column += 1
	row += 1

#
phi = psi
# Store the energy for initial wavefunction
e = np.matmul(Hamiltonian, phi)
e = np.real(np.matmul(np.transpose(np.conj(phi)), e))
energy_qite_list.append(e)

print('We start QITE now')
for i in range(0,N):

	# First construct Pmu_expectation matrices
	Pmu_expectation = np.zeros([16], dtype=complex)
	for Pmu in range(2**4):
		ops = Int2Bas(Pmu, 4, 2)
		operator = np.kron(sigma[:, :, ops[0]], sigma[:, :, ops[1]])
		Pmu_expectation[Pmu] = np.matmul(np.transpose(np.conj(phi)), np.matmul(operator, phi))

	# Now construct S matrix
	S = np.zeros([16, 16], dtype=complex)
	for i in range(16):
		for j in range(16):
			S[i,j] = Pmu_expectation[index[i, j]]*coeff[i, j]

	# Now construct b vector
	b = np.zeros([16], dtype=complex)
	c = 1

	# We will hardcode in the QITE step
	
	c -= 2 * db * (-t_hubbard) * Pmu_expectation[1]
	c -= 2 * db * (-t_hubbard) * Pmu_expectation[4]
	c -= 2 * db * (U_hubbard / 2) * Pmu_expectation[0]
	c -= 2 * db * (U_hubbard / 2) * Pmu_expectation[15]
	c = np.sqrt(c)
	for i in range(16):
		b[i] += (Pmu_expectation[i] / c - Pmu_expectation[i]) / db
		b[i] -= (-t_hubbard) * coeff[i, 1] * Pmu_expectation[index[i, 1]] / c
		b[i] -= (-t_hubbard) * coeff[i, 4] * Pmu_expectation[index[i, 4]] / c
		b[i] -= (U_hubbard / 2) * coeff[i, 0] * Pmu_expectation[index[i, 0]] / c
		b[i] -= (U_hubbard / 2) * coeff[i, 15] * Pmu_expectation[index[i, 15]] / c
		b[i] = 1j * b[i] - 1j * np.conj(b[i])

	# Obtain x 
	dalpha = np.eye(16) * 0.01
	x = np.linalg.lstsq(S + np.transpose(S) + dalpha, -b, rcond=-1)[0]

	# Classical evolution
	U = np.eye(4)
	for i in range(len(x)):
		ops = Int2Bas(i, 4, 2)
		operator = np.kron(sigma[:, :, ops[0]], sigma[:, :, ops[1]])
		U = np.matmul(expm(1j * db * x[i] * operator), U)
	phi = np.matmul(U, phi)
	e = np.matmul(Hamiltonian, phi)
	e = np.real(np.matmul(np.transpose(np.conj(phi)), e))
	energy_qite_list.append(e)
print('Final energy from QITE is ', e)

plt.figure()
beta_list = np.asarray(range(0, N+1)) * db
plt.plot(beta_list, energy_classical_list, '-bo', label='Exact trajectory')
plt.plot(beta_list, energy_qite_list, '-ro', label='QITE approximation')
plt.plot([beta_list[0], beta_list[-1]], [w[0], w[0]], '--k', label='Ground state')
plt.legend()
plt.show()




