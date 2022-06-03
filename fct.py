q = 2**8
Fq = GF(q, name='a')
R.<X> = Fq['X']
C = Fq.list()

n, k = 8, 3

x = [254, 20, 11]

#Generation de v
v = []
for i in range(n):
    v.append(C[randint(1, len(C)-1)])

#Generation de a
a = []
for i in range(n):
    c = C[randint(1, len(C)-1)]
    while c in a:
        c = C[randint(1, len(C)-1)] #generer un element unique
    a.append(c)

y = codeGRS(x, v, a)
x2 = decodeGRS(y, v, a)
yp = errTrans(y, 2)

r = n - k

S = syndrome(yp, a, v)
sigma, omega = clef(S)
e = erreur(sigma, omega, a, v)

#corriger les erreurs:
yc = [i-j for i, j in zip(yp, e)]
print (y==yc)

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

def errTrans(y, nb_Err):
    y_prime = list(y)
    indices = random.sample(range(len(y)), nb_Err)
    for i in indices:
        e = C[randint(1, len(C)-1)] 
        y_prime[i] += e
    return y_prime

def syndrome(yp, a, v):
    '''Calcule le polynome syndrome'''
    S = R.zero()
    for i in range(n):
        L = Li(i, a)
        S += yp[i]/(v[i] * L(a[i])) * sum([(a[i] * X)**j for j in range(r)])
    return S

def clef(S):
    '''Trouve les polynomes localisateur et evaluateur'''
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
    sigmat = v[j]
    omegat = rr[j]

    return sigmat/sigmat(0), omegat/sigmat(0)

def erreur(sigma, omega, a, v):
    #ui = 1/(vi Li(ai))
    e = [0]*n
    B = [sigma(1/a[b]) == 0 for b in range(n)]

    for b, i in enumerate(B):
        if i == True:
            L = Li(b, a)
            e[b] = -a[b] * omega(1/a[b]) * v[b] * L(a[b]) * 1/sigma.derivative()(1/a[b])
    return e