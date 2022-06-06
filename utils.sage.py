def codeGRS(message, V, A):
    x = []
    for  c in message:
        x.append(C[c])
    f = R.zero()
    for i in range(k):
        f += x[i] * X**i
    return [v*f(a) for v, a in zip(V, A)]

def Li(i, A):
    L = R.one()
    for j in range(n):
        if i != j:
            L *= X - A[j]
    return L 

def decodeGRS(Y, V, A):
    f_alpha = [y/v for y, v in zip(Y, V)]
    f = R.zero()
    for i in range(n):
        L = Li(i, A)
        f += (f_alpha[i]*L)/L(A[i])
    message = []
    for c in f.coefficients() :
        message.append(C.index(c))
    return message

p = 2
n = 6
F2 = GF(p)
R = F2.polynomial_ring()
prod = R.one() #utilise pour le produit
for d in divisors(n):
  for poly in R.polynomials(of_degree = d):
    if poly.is_irreducible() and poly.is_monic(): # si irreductible et unitaire
      prod *= poly

def errTrans(y, nb_Err):
    y_prime = list(y)
    indices = random.sample(range(len(y)), nb_Err)
    for i in indices:
        e = C[randint(1, len(C)-1)] 
        y_prime[i] += e
    return y_prime
