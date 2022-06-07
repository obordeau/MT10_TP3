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

def Syndrome(y_prime, A, V):
    S = R.zero()
    for i in range(n):
        L = Li(i, A)
        S += y_prime[i]/(V[i] * L(A[i])) * sum([(A[i] * X)**j for j in range(r)])
    return S

def Clef(S):
    rr = [X**r, S]
    u = [R.one(), R.zero()]
    v = [R.zero(), R.one()]
    q = [R.zero()]
    j = 1
    while rr[j].degree() >= r/2:
        q.append(rr[j-1] // rr[j])
        rr.append(rr[j-1] % rr[j])
        u.append(u[j-1] - u[j]*q[j])
        v.append(v[j-1] - v[j]*q[j])
        j += 1
    sigma = v[j]
    omega = rr[j]
    return sigma/sigma(0), omega/sigma(0)

def Erreur(sigma, omega, A, V):
    e = [0]*n
    B = [sigma(1/A[b]) == 0 for b in range(n)]

    for b, i in enumerate(B):
        if i == True:
            L = Li(b, A)
            e[b] = -A[b] * omega(1/A[b]) * V[b] * L(A[b]) * 1/sigma.derivative()(1/A[b])
    return e

p = 2
n = 6
F2 = GF(p)
R = F2.polynomial_ring()
produit = R.one()
for d in divisors(n):
  for poly in R.polynomials(of_degree = d):
    if poly.is_irreducible() and poly.is_monic():
      produit *= poly

def test1(n):
    for i in range(1, n + 1):
        if moebius(i) not in [-1, 0, 1]:
            return False
    return True

def test2(n):
    for i in range(2, n + 1):
        if sum(map(moebius, divisors(i))) != 0:
            return False
    return True

def phiMobius(n):
    p = 0
    for d in divisors(n):
        p += moebius(n /d ) * d
    return p

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

def irr(p, n):
    res = 0
    for d in divisors(n):
        res += moebius(n / d) * pow(p, d)
    return res / n

def split(message, k):
    result = []
    while (len(message) % k != 0):
        message.append(0)
    for i in range(0, len(message), k):
        new_element = []
        for j in range(k):
            new_element.append(message[i + j + 1])
        result.append(new_element)
    return result

def gather(list):
    result = []
    for sublist in list:
        for element in sublist:
            result.append(element - 1)
    while (result[-1] == 0):
        result = result[:-1]
    return result

def decoding_with_potential_errors(y, V, A):
    S = Syndrome(y, A, V)
    if S == 0 :
        return decodeGRS(y, V, A)
    else : 
        sigma, omega = Clef(S)
        e = Erreur(sigma, omega, A, V)
        yc = [i-j for i, j in zip (y, e)]
        return decodeGRS(yc, V, A)
    
F2 = GF(2)
R.<X> = F2['X']
n = 0
for p in R.polynomials(max_degree = 10) : 
    if p.is_irreducible() : 
        n += 1
print(n)