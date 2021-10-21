import numpy as np
import pymbc

np.random.seed(12345)

D = np.array([
    [1,1,1,0],
    [0,1,1,1],
    [0,1,1,1],
    [0,1,1,0],
    [0,0,1,1],
    [1,1,1,1],
    [1,1,1,1],
    [0,1,1,0],
], dtype=np.bool)

U = np.arange(D.shape[0]).astype(np.uint64)
V = np.arange(D.shape[1]).astype(np.uint64)

(U,V) = pymbc.maximum_biclique_search(U,    # U
                                      V,    # V
                                      D,    # D
                                      2,    # tau_u
                                      2,    # tau_v
                                      3,    # init_type (STAR)
                                      1,    # init_iter 1 (not used)
                                      np.random.randint(0, 32767), # random seed
                                      True, # use_star
                                      2,    # star_max_iter
                                      3)    # optimize

print("U,V",U,V)
print("size",len(U)*len(V))
assert np.all(D[np.ix_(U,V)])

