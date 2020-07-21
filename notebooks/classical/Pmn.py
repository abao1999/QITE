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