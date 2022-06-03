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

def test1(n):
    for i in range(1,n + 1):
        if moebius(i) not in [-1, 0, 1]:
            return False
    return True

def test2(n):
    for i in range(1,n + 1):
        sum = 0
        for d in divisors(i):
            sum += moebius(d)
        if (i == 1 and sum != 1) or (i > 1 and sum != 0):
            return False
    return True

def phiMobius(n):
    p = 0
    for d in divisors(n):
        p += moebius(n /d ) * d
    return p

def irr(p, n):
    res = 0
    for d in divisors(n):
        res += moebius(n/d) * pow(p, d)
    return res/n

def polynomes10():
    number = 0
    poly = []
    if X.is_irreducible() : 
        poly.append(X)
        number += 1
    for n in range (2,2048):
        P = 0
        i = 0
        for c in (bin(n)[2:]):
            k = int(c)
            P += (X^i)*k
            i += 1
        if not(P in poly) :
            if P.is_irreducible() :
                poly.append(P)
                number += 1
    print(number)

p = 2
n = 6
F2 = GF(p)
R = F2.polynomial_ring()
prod = R.one() #utilise pour le produit
for d in divisors(n):
  for poly in R.polynomials(of_degree = d):
    if poly.is_irreducible() and poly.is_monic(): # si irreductible et unitaire
      prod *= poly