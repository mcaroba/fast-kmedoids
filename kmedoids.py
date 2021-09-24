import numpy as np
import random
from src.fast_module import kmedoids_module as kmod

def kMedoids(D, k, tmax=100, init_Ms="random", n_iso=None):
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
    elif str(init_Ms) == "isolated":
        if n_iso is None:
            n_iso = k
        elif n_iso > k:
            n_iso = k
        else:
            assert type(n_iso) == int
            if n_iso == 0:
                raise Exception('You must choose a non-zero number of \
                                 isolated medoids or switch to random init.')
        summed_distances = np.sum(D, axis = 0)
        av_distance = np.sum(D) / float(n)**2
        M = np.empty(0, dtype=int)
        # Find most isolated sample to use as first medoid
        i = np.argmax(summed_distances)
        M = np.append(M, i)
        # Find new medoids which are furthest away from previous medoids
        for i in range(1, n_iso):
            summed_distances = np.zeros(n)
            for j in range(0, len(M)):
                summed_distances += D[:, M[j]]
            # avoid repeating medoids
            summed_distances[M] = 0 
            j = np.argmax(summed_distances)
            M = np.append(M, j)
        # Randomize the rest
        for i in range(n_iso, k):
            rand = np.random.random_integers(0,n-1)
            # Make sure that the random medoid is not already chosen and
            # is not too close to another medoid by imposing the inter-
            # medoid distance to be larger than the average inter-sample
            # distance
            tries = 0
            while rand in M or any(D[rand,j] < av_distance for j in M):
                tries += 1
                rand = np.random.random_integers(0,n-1)
                if tries == 100:
                    tries = 0
                    av_distance *= 0.9
            M = np.append(M, rand)
    elif len(init_Ms) == k:
        M = init_Ms
    elif len(init_Ms) < k:
        M = init_Ms
        for i in range(len(init_Ms), k):
            rand = np.random.random_integers(0,n-1)
            while rand in M:
                rand = np.random.random_integers(0,n-1)
            M = np.append(M, rand)
    else:
        raise Exception("intial medoids vector has the wrong length")

#   Use Fortran module to run kmedoids algorithm fast (number of clusters is inferred
#   from size of M)
    M, C = kmod.fast_kmedoids(D, tmax, M)
#   Return results like CB's code for compatibility
    C_dict = {}
    for i in range(0, k):
        C_dict[i] = np.where( C == i )[0]

    # return results
    return M, C_dict
