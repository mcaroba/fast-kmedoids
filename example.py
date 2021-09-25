from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import kmedoids

data = np.random.sample([10000,2])
D = pairwise_distances(data, metric='euclidean')

M, C = kmedoids.kMedoids(D, 20, n_iso=10, init_Ms="isolated", n_inits=10)
