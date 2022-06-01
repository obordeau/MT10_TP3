def codeGRS(q, message, V, A):
	code = []
	FqX.<X> = GF(q)['X']
	f = FqX.zero()
	for i, x in enumerate(message):
		if x < 0 or x >= q:
			print("Le message donne n'est pas valide.")
			return
		f += Integer(message[-(i+1)]) * X^i
	for i, a in enumerate(A):
		code.append(V[i] * f(a))
	return code

def polynomeLagrange(q, A, i):
	FqX.<X> = GF(q)['X'] 
	polynome = FqX.one()
	for j in range(len(A)):
		if j != i:
			polynome *= (X - A[j])  
	return polynome

def decodeGRS(q, code, V, A):
	message = []
	FqX.<X>=GF(q)['X']
	polynome = FqX.zero()
	for i, a in enumerate(A):
		polynomeLge = polynomeLagrange(q, A, i)		
		polynome += code[i] * pow(V[i] * polynomeLge(a),-1) *  polynomeLge
	message = polynome.list()
	message.reverse()
	return message

def errTrans(q, y, Nb_err):
	yprime = list(y)
	indices = random.sample(range(len(yprime)), Nb_err)
	for i in indices:
		nombresPossibles = list(range(0, q))
		nombresPossibles.remove(yprime[i])
		yprime[i] = random.choice(nombresPossibles)
	return yprime

def Syndrome(q, k, code, V, A):
	FqX.<X> = GF(q)['X']
	n = len(V)
	r = n-k
	result = FqX.zero()
	for i in range(n):
		sum = FqX.one()
		for j in range(1, r):
			sum += pow(FqX(A[i]) * X,j)
		polynome = polynomeLagrange(q, A, i)
		result += FqX(code[i]) * (pow(FqX(V[i]),-1) * pow(polynome(FqX(A[i])),-1)) * sum 
	return result

def Clef(q, synd, k, n):
	FqX.<X> = GF(q)['X']
	r = n - k
	rJ_1 = X^r
	rJ = synd
	uJ_1 = FqX.one()
	uJ = FqX.zero()
	vJ_1 = FqX.zero()
	vJ = FqX.one()
	tmp = 0
	qJ = 0
	while(FqX(rJ).degree() >= r/2):
		qJ = FqX(rJ_1) // FqX(rJ) 
		tmp = rJ
		rJ = FqX(rJ_1) % FqX(rJ) 
		rJ_1 = tmp
		tmp = uJ
		uJ = uJ_1 - uJ * qJ
		uJ_1 = tmp
		tmp = vJ
		vJ = vJ_1 - vJ * qJ
		vJ_1 = tmp
	sigma = pow(vJ(0),-1) * vJ
	omega = pow(vJ(0),-1) * rJ
	return sigma, omega