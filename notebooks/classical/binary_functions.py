import numpy as np

def Int2Bas(n,b,nbit):
 if(n==0): return [0]*nbit
 x=[]
 while(n):
  x.append(int(n%b))
  n//=b
 return [0]*(nbit-len(x))+x[::-1]

def Bas2Int(x,b,dotme=None):
 nbit=len(x)
 if(dotme is None): dotme = [b**(nbit-1-i) for i in range(nbit)]
 return np.dot(dotme,x) 

def Psi2Str(x):
 s='|'
 for i in x: s+=str(i)
 return s+'>'

def str_op(i):
 if(i==0): return 'I'
 if(i==1): return 'X'
 if(i==2): return 'Y'
 if(i==3): return 'Z'

def Opp2Str(x):
 s=''
 for i in x: s+=str_op(i)
 return s

def Int2Str(integer):
	if integer ==0:
		return 'I'
	elif integer == 1:
		return 'X'
	elif integer == 2:
		return 'Y'
	elif integer == 3:
		return 'Z'
	else:
		raise ValueError
def Str2Opp(x):
	Opp = []
	for i in x:
		if i == 'I':
			Opp.append(0)
		elif i == 'X':
			Opp.append(1)
		elif i == 'Y':
			Opp.append(2)
		elif i == 'Z':
			Opp.append(3)
	return Opp
		
if __name__ == "__main__":
 n0 = 13
 x  = Int2Bas(n0,4,10)
 n  = Bas2Int(x,4)
 print(n0,x,n)
 
 n0 = 2**18+5
 x  = Int2Bas(n0,2,20)
 n  = Bas2Int(x,2)
 print(n0,x,n)
