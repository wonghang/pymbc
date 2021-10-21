import numpy as np
import pymbc
import time

np.random.seed(12345)

M = 200
N = 40
d_count = int(0.01*M*N)
D = np.full((M,N),True,dtype=np.bool)
i = np.random.randint(0,M,size=d_count)
j = np.random.randint(0,N,size=d_count)

D[i,j] = False

print(M*N - np.sum(D),M*N)

st = time.time()
(U1,V1) = pymbc.maximum_biclique_search(None, # U
                                        None, # V
                                        D,    # D
                                        10,   # tau_u
                                        10,   # tau_v
                                        0,    # init_type (BEST)
                                        3,    # init_iter 3
                                        np.random.randint(0, 32767), # random seed
                                        True, # use_star
                                        2,    # star_max_iter
                                        3)    # optimize
et = time.time()
assert np.all(D[np.ix_(U1,V1)])
print("MBC* spent %.6fs" % (et-st))

st = time.time()
(U2,V2) = pymbc.maximum_biclique_search(None, # U
                                        None, # V
                                        D,    # D
                                        10,   # tau_u
                                        10,   # tau_v
                                        0,    # init_type (BEST)
                                        3,    # init_iter 3
                                        np.random.randint(0, 32767), # random seed
                                        False, # use_star
                                        2,    # star_max_iter
                                        3)    # optimize
et = time.time()
assert np.all(D[np.ix_(U2,V2)])
print("MBC spent %.6fs" % (et-st))

size1 = len(U1)*len(V1)
size2 = len(U2)*len(V2)
assert size1 == size2, (size1,size2)

print("U1,V1",U1,V1)
print("size",size1)

