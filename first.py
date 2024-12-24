from itertools import product
from math import cos, sin, pi
from numpy import roots, zeros
import pickle

q = 7
tau = (2.3704694055762 - 1.1751062918847874j)


# DATA MANIPULATION
def load(filename):
    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported):", ex)


def save(obj, filename):
    try:
        with open(filename, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)


class polynomial:

    def __init__(self, coefficients):
        if isinstance(coefficients, int):
            coefficients = [coefficients % q]
        if isinstance(coefficients, str):
            coeffs = []
            coefficients = coefficients.replace('x^1', 'x')
            coefficients = coefficients.replace('*1', '')
            coefficients = coefficients.replace('*', '')
            pl = coefficients.split('x')
            if len(pl[0]) == 0:
                coeffs.append(1)
            else:
                coeffs.append((int(pl[0])))
            for i in range(1, len(pl)):
                if len(pl[i]) > 0:
                    if pl[i][0] == '^':
                        a = int(pl[i][1])
                        try:
                            if len(pl[i + 1]) > 0 and pl[i + 1][0] == '^':
                                b = int(pl[i + 1][1])
                            else:
                                b = 1
                        except:
                            b = 0
                        for _ in range(a - b - 1):
                            coeffs.append(0)
                        if pl[i][-1].isdigit():
                            coeffs.append(int(pl[i][-1]))
                        else:
                            coeffs.append(1)
                    else:
                        if len(pl[-1]) > 0:
                            if pl[-1][-1].isdigit():
                                coeffs.append(int(pl[-1][-1]))
                            else:
                                coeffs.append(1)
                else:
                    coeffs.append(0)
            coefficients = coeffs
        self.coefficients = coefficients

    def __str__(self):
        out = ''
        size = len(self.coefficients)
        for i in range(size):
            if self.coefficients[i] != 0:
                out += ' + %g*x^%d' % (self.coefficients[i], size - i - 1)
            elif self.degree() == 0:
                out += '0'
        out = out.replace('+ -', '- ')
        out = out.replace('x^0', '1')
        out = out.replace(' 1*', ' ')
        out = out.replace('x^1', 'x')
        if out[0:3] == ' + ':
            out = out[3:]
        if out[0:3] == ' - ':
            out = '-' + out[3:]
        return out

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        res = []
        for i in range(1, max(len(self.coefficients), len(other.coefficients)) + 1):
            try:
                res.append((self.coefficients[-i] + other.coefficients[-i]) % q)
            except IndexError:
                if len(self.coefficients) > len(other.coefficients):
                    res.append(self.coefficients[-i])
                else:
                    res.append(other.coefficients[-i])
        res.reverse()
        return polynomial(res)

    def __sub__(self, other):
        res = []
        for i in range(1, max(len(self.coefficients), len(other.coefficients)) + 1):
            try:
                res.append((self.coefficients[-i] - other.coefficients[-i]) % q)
            except IndexError:
                if len(self.coefficients) > len(other.coefficients):
                    res.append(self.coefficients[-i])
                else:
                    res.append((-1) * other.coefficients[-i] % q)
        res.reverse()
        return polynomial(res)

    def __mul__(self, other):
        res = []
        for i in range(len(self.coefficients) + len(other.coefficients) - 1):
            res.append(0)
        for i in range(len(self.coefficients)):
            for j in range(len(other.coefficients)):
                res[i + j] = (res[i + j] + (self.coefficients[i] * other.coefficients[j])) % q
        return polynomial(res)

    def __truediv__(self, other):
        return polynomial(self.division(other)[0])

    def __mod__(self, other):
        return polynomial(self.division(other)[1])

    def __eq__(self, other):
        if isinstance(other, polynomial):
            return self.coefficients == other.coefficients
        else:
            if self.degree() == 0:
                return self.coefficients[0] == other
            else:
                return False

    def division(self, other):
        quotient = []
        remainder = self.coefficients
        i = 0
        if isinstance(other, int):
            other = polynomial([other])
        while len(remainder) >= len(other.coefficients) and \
                len(quotient) < len(self.coefficients):
            i += 1
            loc = []
            coeff = (modInverse(other.coefficients[0], q) * remainder[0]) % q
            quotient.append(coeff)
            for j in range(1, len(remainder)):
                try:
                    a = (remainder[j] - coeff * other.coefficients[j]) % q
                    if a != 0 or len(loc) > 0:
                        loc.append(a)
                except IndexError:
                    if remainder[j] != 0 or len(loc) > 0:
                        loc.append(remainder[j])
            if len(loc) == 0:
                loc = [0]
            remainder = loc
        return quotient, remainder

    def degree(self):
        return len(self.coefficients) - 1

    def derivative(self):
        df = []
        for i in range(self.degree()):
            df.append(self.coefficients[i] * (self.degree() - i) % q)
        return polynomial(df)

    def isSquareFree(self):
        if self.degree() > 1:
            return gcd(self, self.derivative()).degree() == 0
        else:
            return True

    def isMonic(self):
        return self.coefficients[0] == 1

    def rewrite(self):
        s = self.degree() % 2
        f0_coeffs = [self.coefficients[2 * i + s] for i in range(int(self.degree() / 2) + 1)]
        f1_coeffs = [self.coefficients[2 * i + 1 - s] for i in range(int(self.degree() / 2) + s)]
        return polynomial(f0_coeffs), polynomial(f1_coeffs)


# a^n mod m; a, m are polynomials
def fastExponentiation(a, n, m):
    binn = bin(n)
    currPow = a
    if binn[-1] == '1':
        res = currPow
    else:
        res = polynomial([1])
    for i in range(2, len(binn)):
        currPow = currPow * currPow % m
        if binn[-i] == '1':
            res = (res * currPow) % m
    if res.degree() >= 1:
        print(res.degree(), len(res.coefficients), res.coefficients, res, a)
    return res.coefficients[0]


# inverse of a mod m; a, m are integers
def modInverse(A, M):
    for X in range(1, M):
        if ((A % M) * (X % M)) % M == 1:
            return X
    return -1


# cubic characters; (f/m)_3; f, m are polynomials
# 'r' in the name indicates the use of the law of reciprocity
def dirichletChar(f, m):
    f = f % m
    if f.degree() == 0:
        if f == polynomial([0]):
            return 0
        else:
            res = f.coefficients[-1] ** int(m.degree() * ((q - 1) / 3)) % q
    else:
        p = int((q ** m.degree() - 1) / 3)
        res = fastExponentiation(f, p, m)
    if res != 1:
        return complex(cos((res * pi) / 3), sin((res * pi) / 3))
    else:
        return 1


def rDirichletChar(f, m):
    if f.degree() == 0:
        res = f.coefficients[0] ** (2 * m.degree()) % q
        if res != 1 and res != 0:
            return complex(cos((res * pi) / 3), sin((res * pi) / 3))
        # else:
        #     return res
    if f.degree() > 0 and f.degree() >= m.degree():
        f = f % m
    power = (q - 1) / 3
    coeff = 1
    if f == polynomial([0]):
        return 0
    while f.degree() > 0:
        loc = m % f
        sgn = (-1) ** (power * f.degree() * m.degree())
        sgnf = f.coefficients[0] ** (power * m.degree()) % q
        sgnm = modInverse(m.coefficients[0], q) ** (power * f.degree()) % q
        coeff = coeff * sgnf * sgnm * sgn % q
        m = f
        f = loc
    res = coeff * f.coefficients[-1] ** int(m.degree() * ((q - 1) / 3)) % q
    if res != 1 and res != 0:
        return complex(cos((res * pi) / 3), sin((res * pi) / 3))
    else:
        return res


# Dirichlet character; nested function instead of a loop
def nestedChi(f, p):
    if f.degree() >= p.degree():
        f = f % p
    pow = (q - 1) / 3
    if f.degree() == 0:
        res = f.coefficients[0] ** (pow * p.degree()) % q
        return res
    else:
        sgn = (-1) ** (pow * f.degree() * p.degree())
        sgnf = f.coefficients[0] ** (pow * p.degree()) % q
        sgnm = modInverse(p.coefficients[0], q) ** (pow * f.degree()) % q
        return (sgnf * sgnm * sgn) * nestedChi(p, f) % q


def chi(f, p):
    res = nestedChi(f, p)
    if res != 1 and res != 0:
        return complex(cos((res * pi) / 3), sin((res * pi) / 3))
    else:
        return res


def dirichletChar2(f, m):
    f = f % m
    if f.degree() == 0:
        if f == polynomial([0]):
            return 0
        else:
            res = f.coefficients[-1] ** int(m.degree() * ((q - 1) / 3)) % q
    else:
        p = int((q ** m.degree() - 1) / 3)
        res = fastExponentiation(f, p, m) ** 2 % q
    if res != 1:
        return complex(cos((res * pi) / 3), sin((res * pi) / 3))
    else:
        return 1


def rDirichletChar2(f, m):
    f = f % m
    power = (q - 1) / 3
    coeff = 1
    if f == polynomial([0]):
        return 0
    while f.degree() > 0:
        loc = m % f
        sgn = (-1) ** (power * f.degree() * m.degree())
        sgnf = f.coefficients[0] ** (power * m.degree()) % q
        sgnm = modInverse(m.coefficients[0], q) ** (power * f.degree()) % q
        coeff = coeff * sgnf * sgnm * sgn % q
        m = f
        f = loc
    res = coeff * f.coefficients[-1] ** int(m.degree() * ((q - 1) / 3)) % q
    res = res ** 2 % q
    if res != 1 and res != 0:
        return complex(cos((res * pi) / 3), sin((res * pi) / 3))
    else:
        return res


# generates a list of polynomials of degree n
def polynomialsOfDegreeN(n):
    res = []
    if n != 0:
        for p in list(product(range(7), repeat=n + 1))[q ** (n - 1):][q ** (n) - q ** (n - 1):]:
            res.append(polynomial(list(p)))
    else:
        for i in range(1, q):
            res.append(polynomial([i]))
    return res


# generating&finding roots of an L-function
def LFunctionCoeffs(ch, p):
    res = []
    for d in range(p.degree() + 1):
        sum = 0
        for f in polynomialsOfDegreeN(d):
            if f.coefficients[0] == 1:
                sum = sum + ch(f, p)
        res.append(sum)
    return res


def findLRoots(ch, m):
    coeffs = LFunctionCoeffs(ch, m)[:-1]
    coeffs.reverse()
    rs = roots(coeffs)
    return rs


# normalizing a polynomial (not over a finite field)
def normPlnm(a):
    acopy = []
    for i in range(len(a)):
        acopy.append(a[i] / a[-1])
    return acopy


# gcd of two polynomials
def gcd(f, g):
    if g != 0:
        return gcd(g, f % g)
    else:
        return f


# exponential function
def eq(num, den):
    num = num % den
    if num.degree() - den.degree() != -1:
        return 1
    else:
        a = (modInverse(den.coefficients[0], q) * num.coefficients[0]) % q
        return complex(cos(2 * pi * a / q), sin(2 * pi * a / q))


def shiftedGaussSum(s, f):
    if f.degree() == 0:
        return 1
    sum = 0
    for i in range(f.degree()):
        for g in polynomialsOfDegreeN(i):
            sum += rDirichletChar(g, f) * eq(s * g, f)
    return sum


def altShGaussSum(s, f):
    if f.degree() == 0:
        return 1
    sum = 0
    for i in range(f.degree()):
        for g in polynomialsOfDegreeN(i):
            if g.isMonic():
                m = (s * g) % f
                if m.degree() + 1 == f.degree():
                    sum += rDirichletChar(g, f)
    return sum * tau


def nShGS(s, f):
    if f.degree() == 0:
        return 1
    sum = 0
    for g in polynomialsOfDegreeN(f.degree() - 1):
        sum += rDirichletChar(g, f) * eq(s * g, f)
    return sum


# sum of Gauss sums; deg sh = n
def G(sh, n):
    if n == 0:
        return 1
    s = 0
    for g in polynomialsOfDegreeN(n):
        if g.isMonic():
            s = s + shiftedGaussSum(sh, g)
    return s


def grs(sh):
    plnm = zeros(6, complex)
    for i in range(3):
        d = int((1 + sh.degree() - i) / 3)
        for j in range(d + 1):
            if j == 0:
                cf = 0
            else:
                cf = plnm[j * 3 + (i - 1)] * (q ** 3 - q ** 4)
            plnm[j * 3 + i] = G(sh, j * 3 + i) + cf
    plnm = list(reversed(plnm))
    return plnm, roots(plnm), abs(roots(plnm))


def rho(i):
    if i == 0:
        return 1
    if i == 1:
        return tau * q
    if i == 2:
        return 0


def generateCoefficients(f):
    res = []
    for a in range(0, q):
        res.append(acoeff(a, f))
    return res


def monGaussSumNew(v, f):
    s = 0
    coeffs = generateCoefficients(f)
    remX = remModX(v, f)
    for i in range(f.degree()):
        for h in polynomialsOfDegreeN(i)[:7 ** i]:
            r = 0
            for j in range(h.degree() + 1):
                r += (h.coefficients[-j - 1] * remX[-j - 1])
            r = r % q
            s += chi(h, f) * coeffs[r]
    return s




def acoeff(r, f):
    s = 0
    for a in range(1, 7):
        s += rDirichletChar(polynomial([a]), f) * complex(cos(2 * pi * a * r / q), sin(2 * pi * a * r / q))
    return s


def a1(num, den):
    num = num % den
    if num.degree() - den.degree() != -1:
        return 0
    else:
        return (modInverse(den.coefficients[0], q) * num.coefficients[0]) % q


def remModX(v, f):
    res = []
    loc = [0 for _ in range(f.degree())]
    for n in range(f.degree()):
        loc[n] = 1
        xn = polynomial(loc)
        res.append(a1(xn * v, f))
        loc[n] = 0
    return res
