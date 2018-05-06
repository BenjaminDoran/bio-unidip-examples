import string
import random
from itertools import chain
import numpy as np


##
## Calc information content, entropy, and ecdfs ##
##

def ecdf(x):
    """ Return CDF of array """
    return np.arange(len(x))/float(len(x))

def plot_ecdf(a):
    """ Return plot of CDF by sorted xvals """
    import matplotlib.pyplot as plt 
    sorted_ = np.sort(a)
    yvals = np.arange(len(sorted_))/float(len(sorted_))
    return plt.plot(sorted_, yvals)

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


##
## Mutate and insert string motifs ##
##

gen_noise = lambda n: "".join(random.choices("ACGT", k=n))

mutate1 = lambda s, i: s[:i] + gen_noise(1) + s[i+1:]

def mutate_n(s, n): 
    for i in [random.randint(0, len(s)-1) for x in range(n)]:
        s = mutate1(s, i) 
    return s

def insert_s(s, i, ins):
    return s[:i] + ins + s[i+1:]

def insert_n(s, n, ins, **kwargs):
    ls = len(s)
    Lins = len(ins)
    min_offset = ((ls//n)//2) if kwargs.get("min_offset") == None else kwargs.get("min_offset")
    max_offset = min_offset if kwargs.get("max_offset") == None else  kwargs.get("max_offset")
    min_gap = ls//n if kwargs.get("min_gap") == None else kwargs.get("min_gap") 
    max_gap = min_gap if kwargs.get("max_gap") == None else kwargs.get("max_gap") 
    assert min_gap <= max_gap, f"min_gap must be less than max_gap {min_gap} !<= {max_gap}"
    mutate_num = 0 if kwargs.get("mutate_num") == None else kwargs.get("mutate_num")
    
    idxs = []
    for i in range(n):
        idxs.append(random.randint(min_offset, max_offset) + i * (Lins + random.randint(min_gap, max_gap)))
        s = insert_s(s, idxs[-1], mutate_n(ins, mutate_num))

    return idxs, s

def gen_mut_strs(motif="AAAAAAAAAAAAAAAA", insert_num=1, n_strs=50, str_len=100, **kwargs):
    noise_len = str_len - ((len(motif) -1) * insert_num)
    strs = [insert_n(gen_noise(noise_len), insert_num, motif, **kwargs) for i in range(Nstrs)]
    idxs, strs = zip(*strs)
    return idxs, strs

def score(res, nent):
    scores = {}
    for i in res:
        scores[f"cluster{len(scores)+1}"] = sum(nent.iloc[i[0]:i[1]]) / (i[1] - i[0])
        print(f"Cluster {len(scores)} mean conservation", scores[f"cluster{len(scores)}"])
    return scores
