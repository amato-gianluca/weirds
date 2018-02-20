"""
Software used to search for primitive abundant and weird numbers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2017, 2018 Gianluca Amato <gianluca.amato@unich.it>
"""

from sage.all import *
from subprocess import check_output

try:
    prime_pi_kimwalisch(2)
except NameError:
    print "Warning: no primecount library found"
    primecount_lib = False
else:
    primecount_lib = True

try:
    check_output(["primecount", "2"])
except OSError:
    primecount_external = False
else:
    primecount_external = True

def is_sumof_subset(n, setnum, maxnum):
    """
    Returns True if n can be obtained as the sum of the first maxnum elements of setnum.

    INPUT:

    * "n" - integer >= 0.

    * "setnum" - a list of positive (>=0) numbers, ordered from the smallest to the largest.

    * "maxnum" - integer >= 0 and <= len(setnum).
    """
    while maxnum >= 1 and setnum[maxnum-1] > n:
        maxnum -= 1
    if maxnum == 0:
        return n == 0
    if setnum[maxnum-1] == n:
        return True
    # The following is faster then the sum function
    t = 0
    for j in xrange(maxnum): t += setnum[j]
    if n >= t:
        return n == t
    if is_sumof_subset(n - setnum[maxnum-1], setnum, maxnum-1):
        return True
    return is_sumof_subset(n, setnum, maxnum-1)
 
def sigma(n, k = 1, f = None):
    """
    Computes the sum of (the k-th powers of) the divisors of n.

    INPUT:

    * "n" - integer >= 1.

    * "k" - integer >= 0 (default 1).

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    n = ZZ(n)
    k = ZZ(k)
    one = ZZ.one()
    if f == None: f = factor(n)

    if (k == ZZ.zero()):
        return prod(expt+one for p, expt in f)
    elif (k == one):
        return prod((p**(expt+one) - one).divide_knowing_divisible_by(p - one)
                        for p, expt in f)
    else:
        return prod((p**((expt+one)*k)-one).divide_knowing_divisible_by(p**k-one)
                        for p,expt in f)

def divisors(n, f = None):
    """
    Returns a list of all positive integer divisors of the nonzero integer n.

    INPUT:

    * "n" - integer != 0.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    if f == None or n < 10**21:
        return ZZ(n).divisors()
    output = [ZZ.one()]
    for p, e in f:
        prev = output[:]
        pn = ZZ.one()
        for i in range(e):
            pn *= p
            output.extend(a*pn for a in prev)
    output.sort()
    return output

def divisors_linear(n, f = None, divs = [ZZ.one()]):
    """
    If n=km is a nonzero integer and divs is the list of all positive integer divisors of k, it returns
    the list of all positive integer divisors of n.

    INPUT:

    * "n" - integer != 0.

    * "f" - optional linear factorization for "n" which may be used to speed-up computation (default None).

    * "divs" - initial list of divisors (default [1]).
    """
    if f == None or n < 10**21:
        return ZZ(n).divisors()
    output = divs[:]
    oldp = ZZ.one()
    for p in f:
        if p != oldp:
            prev = output[:]
            oldp = p
            powerp = p
        else:
            powerp *= p
        output.extend(a*powerp for a in prev)
    output.sort()
    return output

def abundance(n, f = None):
    """
    Returns the abundance of n.

    INPUT:

    * "n" - integer >= 1.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    return sigma(n, f = f) - 2*n

def deficiency(n, f = None):
    """
    Returns the deficiency of n.

    INPUT:

    * "n" - integer >= 1.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    return 2*n - sigma(n,f = f)

def abundancy(n, f = None):
    """
    Returns the abundancy of n.

    INPUT:

    * "n" - integer >= 1.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    return sigma(n, f = f)/(2*n)

def center(n, f = None):
    """
    Returns the center of n (i.e. sigma(n)/deficiency(n))

    INPUT:

    * "n" - a deficient number.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    return 2*n/deficiency(n, f = f) - 1

def center_square(n, f = None):
    """
    Let n be deficient, p a prime such that (m,p)=1 and r = center_square(n). Then:

    * if p < r, then n*p*p is abundant;

    * if p == r, then n*p*p is perfect;

    * if p > r, then n*p*p is deficient.

    INPUT:

    * "n" - a deficient number.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    sigman = sigma(n, f = f)
    defn = 2*n - sigman
    return sigman * (1 + sqrt(1 + 4*defn/sigman)) / (2*defn)

def is_primitive_abundant(n, f = None):
    """
    Returns whether n is primitive abundant (i.e., it is abundant and all its divisors are deficient).

    INPUT:

    * "n" - integer >= 1.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    if f == None: f = factor(n)
    s = sigma(n, f = f)
    if s <= 2*n: return False
    for p, e in f:
        news = s * (p ** e - 1) // (p ** (e+1) - 1)
        newn = n // p
        if news >= 2*newn: return False
    return True

def is_weird(n, f = None):
    """
    Returns whether n is weird.

    INPUT:

    * "n" - integer >= 1.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).
    """
    if f == None: f = factor(n)
    abundancen = abundance(n,f)
    if abundancen <= 0: return False
    divisorsn = divisors(n, f)
    return not is_sumof_subset(abundancen, divisorsn, len(divisorsn)-1)

def prec_prime(n, proof = None):
    """
    The largest prime strictly smaller than n. It returns 0 if there is no such a prime.

    INPUT:

    * "n" - integer.

    * "proof" - bool or None (default: None). If this is False, this function internally uses pseudo-primes.
      If this is True, real primes are used. If it is None, uses the global default (see "sage.structure.proof.proof").
    """
    if (proof is None and sage.structure.proof.all.arithmetic()) or proof:
        return previous_prime(n) if n > 2 else 0
    else:
        return Integer(pari(n-1).precprime())

def prime_pi_kimwalisch_external(n):
    """
    Counts the  number of primes less or equal than n, using the external software primecount.

    INPUT:

    * "n" - integer <= 10^31.
    """
    return Integer(check_output(["primecount", str(n)]))

def primecount(n, base=1):
    """
    Counts the number of primes from base (excluded) to n (included).

    INPUT:

    * "p" - integer

    * "base" - integer (default: 1).
    """
    if n <= base:
        return 0
    if primecount_lib:
        if n < 2**29:
            return prime_pi(n) - prime_pi(base)
        else:
            return prime_pi_kimwalisch(n) - prime_pi_kimwalisch(base)
    if primecount_external:
        if n < 10**10:
            return prime_pi(n) - prime_pi(base)
        else:
            return prime_pi_kimwalisch_external(n) - prime_pi_kimwalisch_external(base)
    return prime_pi(n) - prime_pi(base)

def primecount_bilateral(n, base=1):
    """
    Counts the number of primes from base (excluded) to n (included). If n < base,
    it returns the number of primes from n to base as a negative number.

    INPUT:

    * "p" - integer.

    * "base" - integer (default: 1).
    """
    if n > base:
        return primecount(n, base)
    elif n < base:
        tmp = -primecount(base, n)
        if base in ZZ and is_prime(int(base)):
            tmp += 1
        if n in ZZ and is_prime(int(n)):
            tmp -= 1
        return tmp
    else:
        return 0

def index_notation(n, f = None, format = "plain"):
    """
    Given a natural number n, returns its factorization and index notation in the chosen format.

    INPUT:

    * "n" - integer >= 1.

    * "f" - optional factorization for "n" which may be used to speed-up computation (default None).

    * "format" - it may be "plain", "compact" or "linear".
    """
    m = ZZ.one()
    sigmam = ZZ.one()
    res_idx = []
    res_primes = []
    if f == None: f = factor(n)
    for p, alpha in f:
        c = sigmam/(2*m - sigmam)
        idx = primecount_bilateral(p, c)
        if format == "plain":
            res_idx.append((idx, alpha))
            res_primes.append((p, alpha))
        elif format == "compact":
            res_idx.append( (idx, alpha) if alpha > 1 else idx )
            res_primes.append( (p, alpha) if alpha > 1 else p )
        elif format == "linear":
            res_idx.append(idx)
            res_idx.extend([None] * (alpha -     1))
            res_primes.extend([p] * alpha)
        palpha = p ** alpha
        m *= palpha
        sigmam *= (p * palpha - 1 ) // (p - 1)
    return (res_primes, res_idx)

def index_notation_to_factorization(psidx, proof=None):
    """
    Given an index notation in one of the supported formats, returns a factorization in
    the same format.
    """
    with sage.structure.proof.all.WithProof('arithmetic',proof):
        ps = []
        p = 1
        exp = 1
        sigma = 1
        n = 1
        for x in psidx:
            if x == None:
                ps.append(p)
                n *= p
                exp += 1
                sigma *= (p**(exp+1) - 1)
                sigma //= (p**(exp) - 1)
            else:
                if type(x) == tuple:
                    idx, exp = x
                    tupled = True
                else:
                    idx = x
                    exp = 1
                    tupled = False
                p = sigma/(2*n - sigma)
                if idx == 0:
                    p = ZZ(p)
                elif idx > 0:
                    for i in xrange(idx): p = next_prime(floor(p))
                else:
                    for i in xrange(-idx): p = prec_prime(ceil(p))
                n *= p**exp
                sigma *= (p**(exp+1) - 1) // (p-1)
                if tupled:
                    ps.append ( (p,exp) )
                else:
                    ps.append (p)
        return ps

def sfpan(omega, firstm = 1, firstp = 2, proof = None, a = None, method = "pan", verbose = -1):
    """
    The behaviour of this function depends on the parameter "method". If "method" is "pan", then a list of
    primitive abundant numbers is returned. Numbers are returned in both compact factorized form and
    compact index form. All numbers have a factorization which starts with the same factorization as
    "firstm", and have a total of "omega" prime factors (counted without multiplicity). All factors not
    included in "firstm" have exponent one. Therefore, if "firstm" is square free, all returned numbers
    are square-free, too. If the parameter "a" is None (which is the default), we return all primitive
    abundant numbers of the above form, otherwise we only do a partial check.

    When "method" is "pancount", the behaviour is similar to the case when method is "pan", but
    only the total count of primitive abunant numbers is returned.

    When "method" is "weird", the behaviour is similar to the case when method is "pan", but only
    primitive weird numbers are reported.

    INPUT:

    * "omega" - integer >= 1, the number of prime factors of the primitive abundant or weird numbers (counted
      without multiplicity).

    * "firstm" - integer (default: 1), the base number from which we start our search. It must be deficient.
       Only particular combinations of omega and firstm are allowed, those which ensure that the algorithm
       is correct and complete. Among others, we allow the cases when firstm is square-free or when firstm
       is a power of two and omega > 2.

    * "firstp" - integer (default: 2), the smallest prime number we are allowed to use.

    * "proof" - bool or None (default: None). If this is False, this function internally uses pseudo-primes.
      If this is True, real primes are used. If it is None, uses the global default (see "sage.structure.proof.proof").

    * "a" - either None OR an integer >= 0 OR a list of None and integers >= 0 (default: None). In case it
      is not a list, it is converted into a constant list of length omega. If it is a list, each element of the
      list is the maximum value of the corresponding position of the index notation for the numbers enumerated
      by the search procedure.

    * "method" - string (default: "pan"), can be one of "pan", "pancount", "weird".

    * "verbose" - integer (default: -1). When greater or equal than 0, the function prints some debugging information. If it
      is -2 or less, the function is completely silent, otherwise when method is "weird" it prints weird numbers
      as soon as they are found.
    """
    def sfpan_aux(j, prevp, m, sigmam):
        count = 0
        deficiencym = 2*m - sigmam
        amplitude = a[j]
        if j == omega - 1:
            p = prec_prime(ceil(sigmam/deficiencym))
            if method == "pancount":
                if verbose >= j+1:
                    print ps[0:j], ": counting primes from", prevp+1, "to", p
                count = primecount(p, prevp)
                if amplitude != None: count = min(count, amplitude)
            else:
                idx = -1
                while p > prevp and (amplitude == None or -idx <= amplitude):
                    ps[j] = p
                    psidx[j] = idx
                    if verbose >= j+1:
                        print ps, 1
                    if method == "weird":
                        n = m*p
                        abundancen = sigmam*(p+1) - 2*n
                        divisorsn = divisors_linear(n, ps[len(f):], divs)
                        # count numbers with small abunance
                        if abundancen < p:
                            small_abundants[0] +=1
                        if not is_sumof_subset(abundancen, divisorsn, len(divisorsn)-1):
                            res.append((list(ps), list(psidx)))
                            # weird numbers are rare: print them as soon as they are found
                            if verbose >= -1:
                                print list(ps), list(psidx), " *" if abundancen < p else ""
                    else:
                        res.append((list(ps), list(psidx)))
                    count += 1
                    idx -= 1
                    p = prec_prime(p)
        else:
            nextprime_center = next_prime(sigmam // deficiencym)
            p = max(next_prime(prevp), nextprime_center, firstp)
            if method != "pancount" or verbose >= j+1 or amplitude != None:
                idx = 1 + primecount(p, nextprime_center)
            while amplitude == None or idx <= amplitude:
                if method != "pancount" or verbose >= j+1:
                    ps[j] = p
                    psidx[j] = idx
                if method != "pancount" or verbose >= j+1 or amplitude != None:
                    idx += 1
                inner_count = sfpan_aux(j+1, p, m*p, sigmam*(p+1))
                if inner_count == 0 and amplitude == None:
                    break
                count += inner_count
                p = next_prime(p)
        if verbose >= j:
            print ps[0:j], count
        return count

    with sage.structure.proof.all.WithProof('arithmetic',proof):
        res = []
        psidx = [None] * omega
        ps = [None] * omega
        f = factor(firstm)
        (res_primes, res_idx) = index_notation(firstm, f, "compact")
        ps[0 : len(res_primes)] = res_primes
        psidx[0: len(res_idx)] = res_idx
        divs = divisors(firstm, f)
        small_abundants = [0]
        firstp = ZZ(firstp)
        if method != "pan" and method != "pancount" and method != "weird":
            raise ValueError, 'method should be one of pan, pancount or weird'
        if len(f) > 0:
            lastp = f[-1][0]
            factors_sigmaalpha =  map ( lambda x:  (x[0]**(x[1]+1) - 1) // (x[0] - 1), f )
            sigmam = prod (factors_sigmaalpha)
            defm = 2*firstm - sigmam
            if defm <= 0:
                raise ValueError, "firstm should be deficient"
            if omega <= len(f):
                raise ValueError, "omega should be longer than the factorization of firstm"
            if  not (
                  # Ensure that the algorithms is correct. This test always succeeds when firstm is square-free
                  # and whem omega > 1 and firstm is a power of two.
                  (omega == len(f)+1 and all ( next_prime(lastp) > sigmam/(defm + 2*firstm/(x-1)) for x in factors_sigmaalpha)) or
                  (omega > len(f)+1 and all ( max(next_prime(sigmam // defm), next_prime(lastp)) >= x - 1 - 2*(omega-1) for x in factors_sigmaalpha ))
               ):
                raise ValueError, 'omega = ' + str(omega) + ' and firstm = ' + str(firstm) + ' is not a valid input combination'
        else:
            lastp = ZZ.one()
            sigmam = ZZ.one()
        if not type(a) == list:
            a = omega * [a]
        count = sfpan_aux(len(res_primes), lastp, firstm, sigmam)
        if method == "weird" and verbose >= -1:
            print "Found", count, "primitive abundant numbers", small_abundants[0], "with small abundance (marked with  *) and", len(res), "weirds"
        if method == 'pancount':
            return count
        else:
            return res

def sfweird(omega, firstm = 1, firstp = 2, proof = None, a = None, verbose = -1):
    """
    Syntactic sugar for callig sfpan with parameter "method" equals to "weird".
    """
    return sfpan(omega, firstm = firstm, firstp = firstp, proof = proof, a = a, verbose = verbose, method="weird")

def pan_bigomega(bigomega, firstm = 1, firstp = 2, a = None, alsoperfect = False, omega = None, onlynsfo = False,
                 report = "all", method = "pan", proof = None, verbose = -1):
    """
    The behaviour of this function depends on the parameter "method". If "method" is "pan", then a list of
    primitive abundant numbers is returned. Numbers are returned in linear factorized form and
    linear index form. All numbers have a factorization which starts with the same factorization as
    "firstm", and have a total of "omega" prime factors, COUNTED WITH MULTIPLICTY.
    If the parameter "a" is None (which is the default), we return all primitive
    abundant numbers of the above form, otherwise we only do a partial check.

    When "method" is "weird", the behaviour is similar to the case when method is "pan", but only
    primitive weird numbers are reported.

    INPUT:

    * "bigomega" - integer >= 1, the number of prime factors of the primitive abundant or weird numbers (counted
      with multiplicity).

    * "firstm" - integer (default: 1), the base number from which we start our search. It must be deficient.

    * "firstp" - integer (default: 2), the smallest prime number we are allowed to use.

    * "a" - either None OR an integer >= 0 OR a list of None and integers >= 0 (default: None). In case it
      is not a list, it is converted into a constant list of length omega. If it is a list, each element of the
      list is the maximum value of the corresponding position of the index notation for the numbers enumerated
      by the search procedure.

    * "alsoperfect" - bool (default: False). If true, instead of primitive abundant numbers the function searches
      for primitive non-deficient numbers. This has no effect when method = "weird", if not a slight decrease
      in performance.

    * "omega" - integer > 0, restrict the search to those numbers with given value of omega.

    * "onlynsfo" - bool (default: True). If true, only number with at least a non square-free odd factor are
      considered.

    * "report" - string (default: "all"). If "all", it reports all pan or weird numbers found, accordinf to the
      method parameter. For the other choices, the value of method is is ignored. With "count" the program
      returns the number of PANs and the minimum and maximum PANs in the search space. With "simple" only the minimum
      and maximum PAN is returned. With "maxmin" minimum and maximum PANs are reported together with their factorization
      and index sequence. The behaviour of this parameter WILL be modified in the future.

    * "method" - string (default: "pan"), can be one of "pan", "pancount", "weird".

    * "verbose" - integer (default: -1). When greater or equal than 0, the function prints some debugging information. If it
      is -2 or less, the function is completely silent, otherwise when method is "weird" it prints weird numbers
      as soon as they are found.
    """
    def pan_aux(j, w, prevp, m, sigmam, sigmaprevpalpha, maxsigmapalpha, found_nosfo):
        count = 0
        minval = Infinity
        maxval = 0
        found_abundant = False
        deficiencym = 2*m - sigmam
        reporting_step = report == "all" or report == "maxmin" or verbose >= j+1
        amplitude = a[j]
        if j == bigomega - 1:
            gamma = sigmam // max(sigmaprevpalpha, maxsigmapalpha)
            plowerlimit = max(prevp, (sigmam - gamma) // (deficiencym + gamma))
            center = sigmam / deficiencym
            p = prec_prime(floor(center)+1) if alsoperfect else prec_prime(ceil(center))
            if p > prevp:
                # case when the last prime is a new prime
                found_abundant = True
                if (not onlynsfo or found_nosfo) and (omega == None or omega - w == 1) and p > plowerlimit:
                    # m * p is primitive non deficient
                    maxval = max(maxval, m*p)
                    idx = 0 if p == center else -1
                    if report != "all":
                        minval = min(minval, m*next_prime(plowerlimit))
                    if report == "count":
                        if verbose >= j+1:
                            print ps[0:j], ": counting primes from", plowerlimit+1, "to", p
                        count += primecount(p, plowerlimit)
                        if amplitude != None:
                            count = min(count, amplitude+1 if alsoperfect and center in ZZ and is_prime(floor(center)) else amplitude)
                    elif report == "maxmin":
                        ps[j] = p
                        psidx[j] = idx
                        if maxval > resmax[0]:
                            resmax[0] = maxval
                            resmax[1] = (list(ps), list(psidx))
                        if minval < resmin[0]:
                            resmin[0] = minval
                            resmin[1] = (list(ps), list(psidx))
                    elif report == "all":
                        while p > plowerlimit and (amplitude == None or -idx <= amplitude):
                            ps[j] = p
                            psidx[j] = idx
                            if verbose >= j+1:
                                print ps, 1, True
                            if method == "weird":
                                n = m*p
                                abundancen = sigmam*(p+1) - 2*n
                                divisorsn = divisors_linear(n, ps)
                                # count numbers with small abunance
                                if abundancen < p:
                                    small_abundants[0] += 1
                                if not is_sumof_subset(abundancen, divisorsn, len(divisorsn)-1):
                                    res.append((list(ps), list(psidx)))
                                    if verbose >= -1:
                                        print ps, psidx , " * " if abundancen < p else ""
                            else:
                                res.append((list(ps), list(psidx)))
                            count += 1
                            idx -= 1
                            p = prec_prime(p)
            # case when the last prime is repeated
            if prevp > 1 and (omega == None or omega - w == 0):
                tmp = prevp * sigmaprevpalpha
                if tmp <= center:
                    # m * prevp is non deficient
                    found = True
                    gamma = sigmam // maxsigmapalpha
                    if tmp > (sigmam - gamma) // (deficiencym + gamma):
                        # m * prevep is primitive non-deficient
                        minval = min(minval, m*prevp)
                        maxval = max(maxval, m*prevp)
                        if report == "count":
                            count += 1
                        if reporting_step:
                            ps[j] = prevp
                            psidx[j] = None
                        if verbose >= j+1:
                            print ps, 1, True
                        if report == "maxmin":
                            if maxval > resmax[0]:
                                resmax[0] = maxval
                                resmax[1] = (list(ps), list(psidx))
                            if minval < resmin[0]:
                                resmin[0] = minval
                                resmin[1] = (list(ps), list(psidx))
                        elif report == "all":
                            if method == "weird":
                                n = m * prevp
                                newsigmaprevpalpha = (sigmaprevpalpha * prevp) + 1
                                newsigmam = sigmam // sigmaprevpalpha * newsigmaprevpalpha
                                abundancen = newsigmam - 2*n
                                divisorsn = divisors_linear(n, ps)
                                # count numbers with small abunance
                                if abundancen < p:
                                    small_abundants[0] += 1
                                if not is_sumof_subset(abundancen, divisorsn, len(divisorsn)-1):
                                    res.append((list(ps),list(psidx)))
                                    if verbose >= -1:
                                        print ps, psidx , " * " if abundancen < p else ""
                            else:
                                res.append((list(ps), list(psidx)))
        else:
            if m > 1 and (omega == None or bigomega - j - 1 >=  omega - w):
                # case when we repeat the last prime
                newsigmaprevpalpha = (sigmaprevpalpha * prevp) + 1
                newsigmam = sigmam // sigmaprevpalpha * newsigmaprevpalpha
                newm = m*prevp
                if 2*newm - newsigmam > 0:
                    if reporting_step:
                        ps[j] = prevp
                        psidx[j] = None
                    (inner_count, inner_min, inner_max, inner_found) = pan_aux(j+1, w, prevp, newm , newsigmam, newsigmaprevpalpha, maxsigmapalpha, True if prevp > 2 else found_nosfo)
                    found_abundant = found_abundant or inner_found
                    count += inner_count
                    maxval = max(maxval, inner_max)
                    minval = min(minval, inner_min)
            if (omega == None or w < omega):
                nextprime_center = next_prime(sigmam // deficiencym)
                p = max(nextprime_center, next_prime(prevp), firstp)
                if reporting_step or amplitude != None:
                    idx = 1 + primecount(p,nextprime_center)
                while amplitude == None or idx <= amplitude:
                    if reporting_step:
                        ps[j] = p
                        psidx[j] = idx
                    if reporting_step or amplitude != None:
                        idx += 1
                    (inner_count, inner_min, inner_max, inner_found) = pan_aux(j+1, w+1, p, m*p , sigmam*(p+1), p+1, max(maxsigmapalpha, sigmaprevpalpha), found_nosfo)
                    found_abundant = found_abundant or inner_found
                    if not inner_found and amplitude == None:
                        break
                    count += inner_count
                    maxval = max(maxval, inner_max)
                    minval = min(minval, inner_min)
                    p = next_prime(p)
        if j <= verbose:
            print ps[0:j], count, minval, maxval, found_abundant
        return (count, minval, maxval, found_abundant)

    with sage.structure.proof.all.WithProof('arithmetic',proof):
        res = []
        ps = [None] * bigomega
        psidx = [None] * bigomega
        f = factor(firstm)
        (res_primes, res_idx) = index_notation(firstm, f, "linear")
        ps[0 : len(res_primes)] = res_primes
        psidx[0: len(res_idx)]= res_idx
        divs = divisors(firstm, f)
        firstp = ZZ(firstp)
        small_abundants = [0]
        found_nosfo = any ( p > 2 and e > 1 for p, e in f )
        resmax = [0, None]
        resmin = [+Infinity, None]
        if not method in ["pan", "weird"]:
            raise ValueError, "method should be one of pan or weird"
        if not report in ["all", "count", "simple", "maxmin"]:
            raise ValueError, "report should be one of all, count or simple"
        if len(f) > 0:
            lastp = f[-1][0]
            factors_sigmaalpha =  map ( lambda x:  (x[0]**(x[1]+1) - 1) // (x[0] - 1), f )
            sigmalastpalpha = factors_sigmaalpha[-1]
            firstmaxsigmapalpha = max (  factors_sigmaalpha[0:-1]  ) if len(factors_sigmaalpha) > 1 else 1
            sigmam = prod(factors_sigmaalpha)
            if 2*firstm - sigmam <= 0:
                raise ValueError, "firstm should be deficient"
            if bigomega <= len(f):
                raise ValueError, "omega should be longer than the factorization of firstm"
        else:
            lastp = ZZ.one()
            sigmalastpalpha = ZZ.one()
            firstmaxsigmapalpha = ZZ.one()
            sigmam = ZZ.one()
        if not type(a) == list:
            a = bigomega * [a]
        (count,  minval, maxval, found_abundant) = pan_aux(len(res_primes), len(f), lastp, firstm, sigmam, sigmalastpalpha, firstmaxsigmapalpha, found_nosfo)
        if method == "weird" and verbose >= -1:
            print "Found", count, "primitive abundants numbers", small_abundants[0], "with small abundance (marked with  *) and", len(res), "weirds"
        if report == "all":
            return res
        elif report == "maxmin":
            return (resmax, resmin)
        else:
            return (count, minval, maxval)

def count_opan(omega, base = None):
    """
    An hakcy procedure to count the number of opan (odd PAN) with given value of omega.

    INPUT:

    * "omega", integer > 0: the number of prime factors of PANs we want to count.

    * "base", integer > 0: the initia value of bigomega to use (default = omega).
    """
    s = 0
    totzero = 0
    i = omega if base == None else base
    while totzero < 5:
        c = pan_bigomega(i, firstp = 3, omega = omega, report="count")
        if c[0] == 0:
            totzero += 1
        else:
            totzero = 0
        s += c[0]
        print "Bigomega: " , i, " count: ", c, " total: ", s
        i += 1

def find_primitive_abundant(m):
    """
    Given m such that c(m) > largest prime factor of m, tries to determine p, q such that m*p*q
    is primitive abundant. This should always succeed (at the first attempt) when m is square-free.
    This is essentially a debugging tool used to prove that primitive abundance of m in the completion
    theorem cannot be ensured if m is not square-free.
    """
    lastp = factor(m)[-1][0]
    print "lastp: ",lastp
    c = center(m)
    print "center: ",c*1.0
    if c >= lastp:
        p = next_prime(floor(c))
        while True:
            q = prec_prime(ceil(center(m*p)))
            if q > p:
                print "p:",p, "q:", q
                if is_primitive_abundant(m*p*q):
                    print "FOUND PRIMITIVE ABUNDANT"
                p = next_prime(p)
            else:
                break

def runtest():
    """
    Test the sfpan and pan_bigomega functions. If everything works well no excpetion is thrown.
    """
    print "Testing sfpan function"
    sfpan_params = [ ((4, 1, 1), [([2, 5, 11, 53], [1, 1, 1, -1]),
                            ([2, 5, 11, 47], [1, 1, 1, -2]),
                            ([2, 5, 11, 43], [1, 1, 1, -3]),
                            ([2, 5, 11, 41], [1, 1, 1, -4]),
                            ([2, 5, 11, 37], [1, 1, 1, -5]),
                            ([2, 5, 11, 31], [1, 1, 1, -6]),
                            ([2, 5, 11, 29], [1, 1, 1, -7]),
                            ([2, 5, 11, 23], [1, 1, 1, -8]),
                            ([2, 5, 11, 19], [1, 1, 1, -9]),
                            ([2, 5, 11, 17], [1, 1, 1, -10]),
                            ([2, 5, 11, 13], [1, 1, 1, -11]),
                            ([2, 5, 13, 31], [1, 1, 2, -1]),
                            ([2, 5, 13, 29], [1, 1, 2, -2]),
                            ([2, 5, 13, 23], [1, 1, 2, -3]),
                            ([2, 5, 13, 19], [1, 1, 2, -4]),
                            ([2, 5, 13, 17], [1, 1, 2, -5]),
                            ([2, 5, 17, 19], [1, 1, 3, -1]),
                            ([2, 7, 11, 13], [1, 2, 2, -1]) ]),
                ((4,110,1), [([2, 5, 11, 53], [1, 1, 1, -1]),
                            ([2, 5, 11, 47], [1, 1, 1, -2]),
                            ([2, 5, 11, 43], [1, 1, 1, -3]),
                            ([2, 5, 11, 41], [1, 1, 1, -4]),
                            ([2, 5, 11, 37], [1, 1, 1, -5]),
                            ([2, 5, 11, 31], [1, 1, 1, -6]),
                            ([2, 5, 11, 29], [1, 1, 1, -7]),
                            ([2, 5, 11, 23], [1, 1, 1, -8]),
                            ([2, 5, 11, 19], [1, 1, 1, -9]),
                            ([2, 5, 11, 17], [1, 1, 1, -10]),
                            ([2, 5, 11, 13], [1, 1, 1, -11])]),
                ((4,10,13), [([2, 5, 13, 31], [1, 1, 2, -1]),
                            ([2, 5, 13, 29], [1, 1, 2, -2]),
                            ([2, 5, 13, 23], [1, 1, 2, -3]),
                            ([2, 5, 13, 19], [1, 1, 2, -4]),
                            ([2, 5, 13, 17], [1, 1, 2, -5]),
                            ([2, 5, 17, 19], [1, 1, 3, -1])]),
                ((3, 4, 1), [([(2, 2), 11, 19], [(1, 2), 1, -1]),
                             ([(2, 2), 11, 17], [(1, 2), 1, -2]),
                             ([(2, 2), 11, 13], [(1, 2), 1, -3])]),
                ((9, 1, 5), [([5, 7, 11, 13, 17, 19, 23, 29, 31],[3, 4, 4, 5, 5, 5, 5, 4, -1])])
    ]

    def filtera(firstm, a, res):
        if a == None: return res
        mps, mpsidx = index_notation(firstm, format="compact")
        return filter(lambda (ps, psidx): all(i == None or abs(i) <= a for i in psidx[len(mps):]), res)

    def filterweird(res):
        return filter(lambda (ps, psidx): is_weird(prod(map(lambda x: x[0]**x[1] if type(x)==tuple else x, ps))), res)

    def filternsfo(res):
        return filter(lambda (ps, psidx): any(psidx[i] == None and ps[i] > 2 for i in xrange(0, len(psidx))), res)

    def filterpan(res):
        return filter(lambda (ps, psidx): abundance(prod(ps)) > 0, res)

    for (omega, firstm, firstp), res in sfpan_params:
        for a in [None, 1,2,5]:
            resa = filtera(firstm, a, res)
            for proof in [False, True]:
                assert sfpan(omega, firstm, firstp, proof=proof, a=a) == resa, str((omega, firstm, firstp, a, resa))
                assert sfpan(omega, firstm, firstp, proof=proof, a=a, method="pancount") == len(resa), str((omega, firstm, firstp, a, resa))
                assert sfpan(omega, firstm, firstp, proof=proof, a=a, method="weird", verbose=-2) == filterweird(resa), str((omega, firstm, firstp, a, resa))

    print "Testing pan_bigomega function"
    pan_params = [ ((4, 1, 1), [([2, 2, 2, 13], [1, None, None, -1]),
                                 ([2, 2, 2, 11], [1, None, None, -2]),
                                 ([2, 2, 11, 19], [1, None, 1, -1]),
                                 ([2, 2, 11, 17], [1, None, 1, -2]),
                                 ([2, 2, 11, 13], [1, None, 1, -3]),
                                 ([2, 5, 5, 13], [1, 1, None, -1]),
                                 ([2, 5, 5, 11], [1, 1, None, -2]),
                                 ([2, 5, 11, 53], [1, 1, 1, -1]),
                                 ([2, 5, 11, 47], [1, 1, 1, -2]),
                                 ([2, 5, 11, 43], [1, 1, 1, -3]),
                                 ([2, 5, 11, 41], [1, 1, 1, -4]),
                                 ([2, 5, 11, 37], [1, 1, 1, -5]),
                                 ([2, 5, 11, 31], [1, 1, 1, -6]),
                                 ([2, 5, 11, 29], [1, 1, 1, -7]),
                                 ([2, 5, 11, 23], [1, 1, 1, -8]),
                                 ([2, 5, 11, 19], [1, 1, 1, -9]),
                                 ([2, 5, 11, 17], [1, 1, 1, -10]),
                                 ([2, 5, 11, 13], [1, 1, 1, -11]),
                                 ([2, 5, 13, 31], [1, 1, 2, -1]),
                                 ([2, 5, 13, 29], [1, 1, 2, -2]),
                                 ([2, 5, 13, 23], [1, 1, 2, -3]),
                                 ([2, 5, 13, 19], [1, 1, 2, -4]),
                                 ([2, 5, 13, 17], [1, 1, 2, -5]),
                                 ([2, 5, 17, 19], [1, 1, 3, -1]),
                                 ([2, 7, 11, 13], [1, 2, 2, -1])]),
                   ((4, 4, 1),  [([2, 2, 2, 13], [1, None, None, -1]),
                                 ([2, 2, 2, 11], [1, None, None, -2]),
                                 ([2, 2, 11, 19], [1, None, 1, -1]),
                                 ([2, 2, 11, 17], [1, None, 1, -2]),
                                 ([2, 2, 11, 13], [1, None, 1, -3])]),
                   ((4, 14, 1), [([2, 7, 11, 13], [1, 2, 2, -1])]),
                   ((3, 1, 1),  [([2, 2, 7], [1, None, 0]),
                                 ([2, 2, 5], [1, None, -1]),
                                 ([2, 5, 7], [1, 1, -1])]),
    ]

    for (bigomega, firstm, firstp), res in pan_params:
        for a in [None, 1,2,5]:
            resa = filtera(firstm, a, res)
            for proof in [False, True]:
                for alsoperfect in [False, True]:
                    respan = filterpan(resa) if alsoperfect == False else resa
                    for onlynsfo in [False, True]:
                        finalres = filternsfo(respan) if onlynsfo else respan
                        assert pan_bigomega(bigomega, firstm, firstp, alsoperfect=alsoperfect, onlynsfo=onlynsfo, a=a) == finalres, str((bigomega, firstm, firstp, a, onlynsfo, finalres))
                        assert pan_bigomega(bigomega, firstm, firstp, alsoperfect=alsoperfect, onlynsfo=onlynsfo, a=a, method="weird", verbose=-2) == filterweird(finalres), str((bigomega, firstm, firstp, a, onlynsfo, finalres))
                        assert pan_bigomega(bigomega, firstm, firstp, alsoperfect=alsoperfect, onlynsfo=onlynsfo, a=a, report="count")[0] == len(finalres), str((bigomega, firstm, firstp, a, onlynsfo, finalres))


    assert pan_bigomega(8,  firstm=2, a=3, method="weird", onlynsfo=True, verbose=-2) == [([2, 2, 2, 23, 23, 53, 691, 32587], [1, None, None, 3, None, 1, 3, -1]),
                                              ([2, 2, 13, 17, 449, 24809, 24809, 351659387],[1, None, 2, 1, 2, 1, None, -3])]

    assert pan_bigomega(5, firstm = 4, onlynsfo=True) == [([2, 2, 11, 11, 23], [1, None, 1, None, -1]),
                                              ([2, 2, 13, 13, 17], [1, None, 2, None, -1]),  ([2, 2, 13, 17, 17], [1, None, 2, 1, None])]

    print "Test index functions"
    for (omega, firstm, firstp), res in sfpan_params:
        for ps, psidx in res:
            assert index_notation_to_factorization(psidx) == ps, str((ps, psidx))
    for (omega, firstm, firstp), res in pan_params:
        for ps, psidx in res:
            assert index_notation_to_factorization(psidx) == ps, str((ps, psidx))

def write_list(filename, l):
    """
    Output the list "l" to the file "filename" one element at a time.
    """
    with open(filename,"w") as o:
        for e in l:
            o.write(str(e))
            o.write("\n")
