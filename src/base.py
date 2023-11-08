import numpy as np

# Obtain minimal gap between adyacent values for sequence X
def get_min_gap(X):
    gap = np.inf
    q = len(X)

    for i in range(q):

        d = abs(X[(i+1) % q] - X[i])

        if d < gap:
            gap = d

    return gap


# Transform number n to base b
def numberToBase(n, b):
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits[::-1]


# greatest common division
def gcd(a, b):
    if(b == 0):
        return abs(a)
    else:
        return gcd(b, a % b)


# hamming correlation with shift 0 for sequences
# u and v with the same length (assumed)
def hamming_correlation(u,v):
    u_eq_v = u == v
    return u_eq_v.sum()


# maximal hamming correlation for sequences
# u and v with the same length (assumed)
#
# as implemented in equation (3) and (4) from [1]
# [1] Lempel, A., & Greenberger, H. (1974). Families of sequences with optimal
# Hamming-correlation properties. IEEE Transactions on Information Theory, 20(1), 90-94.
def maxHC(u,v):

    # for crosscorrelation shift stars in 0
    # for autocorrelation shift stars in 1
    start = 0
    if np.array_equal(u,v): start = 1

    current_maxHC = 0
    for shift in range(start, len(u)):

        hc = hamming_correlation(u, np.roll(v, shift))
        
        if hc > current_maxHC:
            current_maxHC = hc
    
    return current_maxHC


# average hamming auto correlation for a family of
# sequences with the same length (assumed)
#
# as implemented in equation (4) from [2]
# [2] Peng, D. Y., Niu, X. H., & Tang, X. H. (2010). Average Hamming correlation 
# for the cubic polynomial hopping sequences. IET communications, 4(15), 1775-1786.
def avg_autoHC(fam):
    _avgHC = 0
    M, L = fam.shape
    for i in range(M):
        for shift in range(1, L):
            _avgHC += hamming_correlation(fam[i], np.roll(fam[i], shift))

    return _avgHC / (M * (L-1))


# average hamming cross correlation for a family of
# sequences with the same length (assumed)
#
# as implemented in equation (5) from [2]
# [2] Peng, D. Y., Niu, X. H., & Tang, X. H. (2010). Average Hamming correlation 
# for the cubic polynomial hopping sequences. IET communications, 4(15), 1775-1786.
def avg_crossHC(fam):
    _avgHC = 0
    M, L = fam.shape
    for i in range(M):
        for j in range(i):
            for shift in range(L):
                _avgHC += hamming_correlation(fam[i], np.roll(fam[j], shift))

    n = M * (M-1) / 2
    return _avgHC / (L*n)


# average maximal hamming correlation for all paris of
# sequences from a single family
def avg_maxHC(fam):
    mean = 0
    s = len(fam)
    for i in range(s):
        for j in range(i+1):
            mean += maxHC(fam[i], fam[j])

    n = s * (s+1) / 2
    return mean / n


# average maximal hamming correlation for all paris of
# sequences from two families
def avg_maxHC_2fam(fam1, fam2):
    mean = 0
    s = len(fam1)
    for i in range(s):
        for j in range(s):
            mean += maxHC(fam1[i], fam2[j])
        
    n = s**2
    return mean / n


# delete frequences above the maximum, raises error if it
# disrupts the minimum gap property
def filter_freq(seq, maxfreq, mingap):
    newseq = np.delete(seq, np.where(seq >= maxfreq)[0])
    assert get_min_gap(newseq) == mingap, "couldn't filter sequences while preserving minimum gap"
    return newseq


# split the given sequence into a family of sequences with size q
def split_seq(seq, q):
    family = []
    i=0
    j=q
    while j < len(seq):
        family.append(seq[i:j])
        i+=q
        j+=q

    return np.array(family)
    

# calculate the number of frequency hops for the payload
# as a function of its size and the coding rate
def numHops(payload_length, CR):

    assert CR==1 or CR==2, "Only CR 1/3 and CR 2/3 supported"

    length_bits = ( payload_length + 2 ) * 8 + 6
    length_bits *= (3/CR)

    nb_hops_out = ( length_bits + 47 ) // 48

    return nb_hops_out
