def codeGRS(q, message, V, A):
	if len(V) != len(alphas):
		print("Les listes v et alpha doivent avoir la meme longueur.")
		return
	if not Integer(q).is_prime_power():
		print("L'ordre d'un corps fini doit etre une puissance premiere.")
		return
	if not Integer(q).is_prime():
		print("L'ordre q du corps fini doit etre premier.")
		return
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
	if not Integer(q).is_prime_power():
		print("L'ordre d'un corps fini doit être une puissance première.")
		return
	if not Integer(q).is_prime():
		print("L'ordre q du corps fini doit être premier.")
		return
	FqX.<X> = GF(q)['X'] 
	polynome = FqX.one()
	for j in range(len(A)):
		if j != i:
			polynome *= (X - A[j])  
	return polynome

def decodeGRS(q, code, V, A):
	if len(V) != len(alphas):
		print("Les listes v et alpha doivent avoir la même longueur.")
		return
	if not Integer(q).is_prime_power():
		print("L'ordre d'un corps fini doit être une puissance première.")
		return
	if not Integer(q).is_prime():
		print("L'ordre q du corps fini doit être premier.")
		return
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