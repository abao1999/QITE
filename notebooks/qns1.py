from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.noise import estimate_bitstring_probs
from scipy.optimize import minimize
import numpy as np 
import matplotlib.pyplot as plt

def single_step(step_size,p):
	p += H(1)
	p += H(2)
	p += CNOT(1,2)
	p += RZ(0.25*step_size,2)
	p += CNOT(1,2)
	p += H(1)
	p += H(2)

	p += RX(-np.pi/2,1)
	p += RX(-np.pi/2,2)
	p += CNOT(1,2)
	p += RZ(0.25*step_size,2)
	p += CNOT(1,2)
	p += RX(np.pi/2,1)
	p += RX(np.pi/2,2)

	p += CNOT(1,2)
	p += RZ(0.25*step_size,2)
	p += CNOT(1,2)

	p += RZ(0.5*3*step_size,1)
	p += RZ(0.5*3*step_size,2)


def Trotter(steps,time,p):
	for i in range(steps):
		# Times 2 here because Rz is defined as Rz(theta) = exp(-i theta/2)
		single_step(time/steps*2,p)

def ansatz(t):
	p = Program()
	p.wrap_in_numshots_loop(1000)
	ro = p.declare('ro','BIT',3)
	p += H(0)
	p += X(1)
	p += X(2)
	p += CNOT(0,1)
	Trotter(10,t,p)
	p += X(0)
	p += CNOT(0,2)
	p += X(0)
	return p,ro

time_list = [0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0]
sa_list = []
si_list = []
qc = get_qc('3q-qvm')
for c,t in enumerate(time_list):

	p1,ro = ansatz(t)
	p1 += H(0)
	p1 += MEASURE(0,ro[0])
	exe = qc.compile(p1)
	result = qc.run(exe)
	sa = estimate_bitstring_probs(result)
	sa_list.append((sa[0]-sa[1])/4)

	p2,ro = ansatz(t)
	p2 += RX(np.pi/2,0)
	p2 += MEASURE(0,ro[0])
	exe = qc.compile(p2)
	result = qc.run(exe)
	si = estimate_bitstring_probs(result)
	si_list.append((si[0]-si[1])/4)


# Classical solution


from scipy.linalg import expm 
sgm = np.zeros([2,2,4],dtype=complex)
sgm[0,0,0] = 1
sgm[1,1,0] = 1
sgmi = sgm[:,:,0]
sgm[0,1,1] = 1
sgm[1,0,1] = 1
sgmx = sgm[:,:,1]
sgm[0,1,2] = -1j
sgm[1,0,2] = 1j
sgmy = sgm[:,:,2]
sgm[0,0,3] = 1
sgm[1,1,3] = -1
sgmz = sgm[:,:,3]

Hamiltonian= np.kron(sgmx/2,sgmx/2) + np.kron(sgmy/2,sgmy/2) + np.kron(sgmz/2,sgmz/2) + 3*np.kron(sgmi,sgmz/2) + 3*np.kron(sgmz/2,sgmi)
time = np.arange(0,4,0.01)
#time = [0.0]
corr_list = []
for t in time:
	TE = expm(-1j*Hamiltonian*t)
	psi = np.zeros([4],dtype=complex)
	psi[3] = 1
	s1 = np.kron(sgmx/2,sgmi)
	s1 = np.matmul(np.transpose(np.conj(TE)),np.matmul(s1,TE))
	s2 = np.kron(sgmi,sgmx/2)
	s = np.matmul(s1,s2)
	fnc = np.matmul(np.transpose(psi),np.matmul(s,psi))
	corr_list.append(fnc)
plt.figure()
plt.plot(time_list,sa_list,'ro',label='Real PyQuil')
plt.plot(time,np.real(corr_list),'r',markersize=0.8,label='Real Exact')
plt.plot(time_list,si_list,'bo',label='Imaginary PyQuil')
plt.plot(time,np.imag(corr_list),'b',markersize=0.8,label='Imaginary Exact')

plt.grid()
plt.show()