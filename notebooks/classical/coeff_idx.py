import numpy as np
from binary_functions import Int2Bas,Bas2Int,Opp2Str,Str2Opp


def PPmunu(Oppstr):
	if Oppstr == "II":
		return "I",1
	elif Oppstr == "IX":
		return "X",1
	elif Oppstr == "IY":
		return "Y",1
	elif Oppstr == "IZ":
		return "Z",1
	elif Oppstr == "XI":
		return "X",1
	elif Oppstr == "XX":
		return "I",1
	elif Oppstr == "XY":
		return "Z",1j
	elif Oppstr == "XZ":
		return "Y",-1j
	elif Oppstr == "YI":
		return "Y",1
	elif Oppstr == "YX":
		return "Z",-1j
	elif Oppstr == "YY":
		return "I",1
	elif Oppstr == "YZ":
		return "X",1j
	elif Oppstr == "ZI":
		return "Z",1
	elif Oppstr == "ZX":
		return "Y",1j
	elif Oppstr == "ZY":
		return "X",-1j
	elif Oppstr == "ZZ":
		return "I",1
	else:
		raise ValueError



def lie_algebra(mu,nu,n):
	# Return coefficients and index for sigma mu,sigma nu
	index = ''
	coeff = 1
	for i in range(n):
		tmpA,tmpB = PPmunu(mu[i]+nu[i])
		index += tmpA
		coeff *= tmpB
	return coeff,Bas2Int(Str2Opp(index),4)

if __name__ =='__main__':
	n = 2
	index = np.zeros([4**n,4**n],dtype=int)
	coeff = np.zeros([4**n,4**n],dtype=complex)
	for i in range(4**n):
		for j in range(4**n):
			coeff[i,j],index[i,j] = lie_algebra(Opp2Str(Int2Bas(i,4,n)),Opp2Str(Int2Bas(j,4,n)),n)
			print(i,j,index[i,j])