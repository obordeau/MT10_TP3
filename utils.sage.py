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