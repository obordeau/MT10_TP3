def moebiusTest(n):
	for i in range(1, n+1):
		if moebius(i) not in [-1, 0, 1]:
			return False
	return True

def eulerTest(n):
	if moebius(1) != 1:
		return False
	for i in range(2, n+1):
		sumMobeius = 0
		divisors = Integer(i).divisors() # Here, we get the divisors of i
		for divisor in divisors:
			sumMobeius += moebius(divisor)
		if sumMobeius != 0:
			return False
	return True

def moebiusInversionFormulaTest(n):
	"""Indirectly returns the value of the Euler phi function on the integer n"""
	result = 0
	divisors = Integer(n).divisors() # Here, we get the divisors of n
	for divisor in divisors :
		result += (divisor * moebius(n / divisor))
	return result

def Irr(p, n):
	irr = 0
	divisors = Integer(n).divisors()
	for divisor in divisors :
		irr += p^divisor * moebius(n / divisor)
	irr /= n
	return int(irr)

def printF2XIrreductiblePolynomials(maxDegree):
	F2X.<X> = GF(2)['X']
	numberOfIrreductiblePolynomials = 0
	print('Polynômes irréductibles de F2[X] de degré inférieur ou égal à {} :'.format(maxDegree))
	for degree in range(0, maxDegree + 1):
		polynomialsOfSpecDegree = F2X.polynomials(degree)
		for p in polynomialsOfSpecDegree:
			if p.is_irreducible():
				print(p)
				numberOfIrreductiblePolynomials += 1
	print("Il y a {} polynômes irréductibles de F2[X] de degré inférieur ou égal à {}.".format(numberOfIrreductiblePolynomials, maxDegree))


def codeGRS(q, word, vs, alphas):
	if len(vs) != len(alphas):
		raise ValueError("List v and alpha must have the same length.")
	if not Integer(q).is_prime_power():
		raise ValueError("The order of a finite field must be a prime power.")
	if not Integer(q).is_prime():
		raise ValueError("Our implementation only works when the order 'q' of " +
			"the finite field is prime.")
	code = []
	FqX.<X> = GF(q, name='a')['X'] # It represents the polynomials in Fq[X]
	f = FqX.zero() # A polynomial which will have its coefficients in Fq
	for i, letter in enumerate(word):
		if letter < 0 or letter >= q:
			raise ValueError("The given word is not valid.")
		f += Integer(word[-(i+1)]) * X^i
	for i, alpha in enumerate(alphas):
		code.append(vs[i] * f(alpha))
	return code