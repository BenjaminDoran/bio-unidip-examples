import string
import random
from itertools import chain
import numpy as np

def ecdf(x):
    """ Return CDF of array """
    return np.arange(len(x))/float(len(x))

def getInfoCntnt(cnts, is_probs=False, minval=None):
    """ Get Information Content"""
    if not is_probs:
        dat = np.array(cnts) / np.sum(cnts)
    else:
        dat = np.array(cnts)
    dat[dat == 0] = minval or (1 / (10 * len(dat)))
    return np.apply_along_axis(lambda p: -np.log2(p), 0, dat)

def shannonEntropy(cnts, is_probs=False, minval=None):
    """ Calculate Entropy across columns """
    # prepare data
    if not is_probs:
        dat = np.array(cnts) / np.sum(cnts)
    else:
        dat = np.array(cnts) 
    dat[dat == 0] = minval or (1 / (10 * len(dat)))
    # entropy calc
    return -np.sum(np.apply_along_axis(lambda p: p * np.log2(p), 0, dat))

def genStrs(pattern, times, length, number): 
    """ Generate Interspersed Motif and Noise 
        INPUT: 
            pattern (str), motif to insert
            times (int), number of times to insert motif
            length (int), length of noise regions
            number (int), number of samples to generate
            
        RETURNS: list of strings
    """
    ileave = lambda *iters: list(chain(*zip(*iters)))
    genNoise = lambda n: "".join(random.choices("ACGT", k=n))
    strGen = lambda p, t, nl: "".join([genNoise(nl)] + ileave([p for i in range(t)], [genNoise(nl) for i in range(t)]))
    
    return [strGen(pattern, times, length) for i in range(number)]

def genMutStrs(pattern="AAAAAAAAAAAAAAA", mut_num=0, insert_num=1, noise_len=50, noise_var=0, number=20):
    ileave = lambda *iters: list(chain(*zip(*iters)))
    genNoise = lambda n: "".join(random.choices("ACGT", k=n))
    mutate1 = lambda s, i: s[:i] + genNoise(1) + s[i+1:]
    def mutateN(s, n): 
        for i in [random.randint(0, len(s)-1) for x in range(n)]:
            s = mutate1(s, i) 
        return s
    def strGen(p, mn, inum, nl, nv): 
        return "".join([genNoise(nl + random.randint(0, nv))] 
                       + ileave([mutateN(p, mn) for i in range(inum)],
                                [genNoise(nl + random.randint(0, nv)) for i in range(inum)]))
    
    return [strGen(pattern, mut_num, insert_num, noise_len, noise_var) for i in range(number)]