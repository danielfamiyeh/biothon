from random import randint
from cluster import Cluster
from copy import deepcopy

class KMedoid:
    """
    Class for KMedoid clustering object.
    """
    def __init__(self, k, data, metric):
        """
        KMedoid clustering object
        :param k:       Number of objects per cluster
        :param data:    Sequence data to cluster
        :param metric:  Distance metric to use for clustering
        """
        if k < len(data):
            self.k = k
            self.medoid_indices = set()
            while len(self.medoid_indices) < k:
                self.medoid_indices.add(randint(0, len(data)-1))
            self.clusters = [Cluster(data[i]) for i in self.medoid_indices]
            self.tentative_clusters = None
            self.tentative_objects = None
            self.objects = [d for i, d in enumerate(data) if i not in self.medoid_indices]
            self.metric = metric
        else:
            raise ValueError("Param k must be less than number of data points.")

    def solve(self, epsilon, max_iters=100):
        """
        Performs clustering.
        :param epsilon:     Fractional difference between iterations that should
                            be deemed as convergence.
        :param max_iters:   Maximum number of iterations to run for before halting.
        :return:
        """
        self._partition(self.clusters, self.objects)
        for i in range(max_iters):
            last_sse = sum([c.getSSE() for c in self.clusters])
            new_sse = self._new_medoid()

            if new_sse < last_sse:
                self.clusters = deepcopy(self.tentative_clusters)
                self.objects = deepcopy(self.tentative_objects)

            if last_sse - new_sse / last_sse < epsilon:
                print("Converged.")
                break

    def _new_medoid(self):
        """
        Genereates new tentative clusters based on randomly selected medoid.
        :return:    The squared sum-of-errors of based on the new set of
                        cluseters.
        """
        rand_cluster = Cluster("")
        while len(rand_cluster.objects) == 0:
            rand_cluster = self.clusters[randint(0, len(self.clusters)-1)]

        rand_medoid = rand_cluster.objects[randint(0, len(rand_cluster.objects)-1)]

        self.tentative_clusters = [deepcopy(c) for c in self.clusters if c != rand_cluster]
        self.tentative_objects = [deepcopy(o) for o in self.objects if o != rand_medoid.data]

        self.tentative_clusters.append(Cluster(rand_medoid.data))
        self.tentative_objects.append(rand_cluster.centroid.data)

        self._partition(self.tentative_clusters, self.tentative_objects)
        return sum([c.getSSE() for c in self.tentative_clusters])

    def _partition(self, cluster_list, others):
        """
        Partitions a set of clusters.
        :param cluster_list:    List of clusters with centroids.
        :param others:          Data points to partition.
        :return:                None
        """
        for o in others:
            dists = []
            for cluster in cluster_list:
                dists.append(self.metric(o, cluster.centroid.data))
            min_dist = min(dists)
            min_index = dists.index(min_dist)
            cluster_list[min_index].add_object(o, min_dist)

