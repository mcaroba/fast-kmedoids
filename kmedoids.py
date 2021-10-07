import numpy as np
import random
from src.fast_module import kmedoids_module as kmod

def kMedoids(D, k, tmax=100, init_Ms="random", n_iso=None, n_inits=1, incoherence="rel", verbosity=1):
    # determine dimensions of distance matrix D
    m, n = D.shape
    # Assert square distance matrix
    assert n == m

    if k > n:
        raise Exception('too many medoids')
    # randomly initialize an array of k medoid indices
    if str(init_Ms) == "random":
        M = np.arange(n)
        np.random.shuffle(M)
        M = np.sort(M[:k])
        n_iso = 0
    elif str(init_Ms) == "isolated":
        if n_iso is None:
            n_iso = k
            n_inits = 1
        elif n_iso > k:
            n_iso = k
            n_inits = 1
        else:
            assert type(n_iso) == int
            if n_iso == 0:
                raise Exception('You must choose a non-zero number of \
                                 isolated medoids or switch to random init.')
        if n_inits == 1: # else flag for fast Fortran execution
            M = kmod.find_isolated(D, k, n_iso)
    elif len(init_Ms) == k:
        M = init_Ms
        n_inits = 1
    elif len(init_Ms) < k:
        M = init_Ms
        if n_inits == 1: # else flag for fast Fortran execution
            for i in range(len(init_Ms), k):
                rand = np.random.random_integers(0,n-1)
                while rand in M:
                    rand = np.random.random_integers(0,n-1)
                M = np.append(M, rand)
    else:
        raise Exception("intial medoids vector has the wrong length")

#   Use Fortran module to run kmedoids algorithm fast (number of clusters is inferred
#   from size of M)
    if n_inits == 1:
        M, C = kmod.fast_kmedoids(D, tmax, M)
    else:
        if str(init_Ms) == "random" or str(init_Ms) == "isolated":
            M_init_mode = np.array(str(init_Ms), dtype='c')
            M, C = kmod.fast_kmedoids_iterative(D, k, tmax, n_inits, n_iso, incoherence, [0], M_init_mode, verbosity)
        else:
            n_iso = 0
            M_init_mode = np.array("random+", dtype='c')
            M, C = kmod.fast_kmedoids_iterative(D, k, tmax, n_inits, n_iso, incoherence, init_Ms, M_init_mode, verbosity)

#   Return results like CB's code for compatibility
    C_dict = {}
    for i in range(0, k):
        C_dict[i] = np.where( C == i )[0]

    # return results
    return M, C_dict


def cluster_incoherence(D, M, C, incoherence="rel"):
    if isinstance(C, dict):
        C_array = np.empty(len(D), dtype=int)
        for i in C:
            C_array[C[i]] = i
        C = C_array
    I = kmod.get_incoherence(D, M, C, incoherence)
    return I

