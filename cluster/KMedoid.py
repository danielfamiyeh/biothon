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

            # Start with k random medoids
            while len(self.medoid_indices) < k:
                self.medoid_indices.add(randint(0, len(data)-1))

            # Add medoid clusters to list
            self.clusters = [Cluster(data[i]) for i in self.medoid_indices]

            self.tentative_clusters = None
            self.tentative_objects = None

            # Add other data points to object list
            self.objects = [d for i, d in enumerate(data) if i not in self.medoid_indices]
            self.metric = metric
        else:
            # Throw error if k value is invalid
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

            # If we have a better cluster configuration
            if new_sse < last_sse:
                # Copy new cluster information
                self.clusters = deepcopy(self.tentative_clusters)
                self.objects = deepcopy(self.tentative_objects)

            # Convergence test
            if last_sse - new_sse / last_sse < epsilon:
                print("Converged.")
                break

    def _new_medoid(self):
        """
        Genereates new tentative clusters based on randomly selected medoid.
        :return:    The squared sum-of-errors of based on the new set of
                        cluseters.
        """
        # Initialise a blank cluster
        rand_cluster = Cluster("")

        # Pick a random cluster
        while len(rand_cluster.objects) == 0:
            rand_cluster = self.clusters[randint(0, len(self.clusters)-1)]

        # Pick a new, random medoid in the cluster
        rand_medoid = rand_cluster.objects[randint(0, len(rand_cluster.objects)-1)]

        # Copy cluster information to object-level variables
        self.tentative_clusters = [deepcopy(c) for c in self.clusters if c != rand_cluster]
        self.tentative_objects = [deepcopy(o) for o in self.objects if o != rand_medoid.data]

        # Add new cluster info to object-level variables
        self.tentative_clusters.append(Cluster(rand_medoid.data))
        self.tentative_objects.append(rand_cluster.centroid.data)

        # Partition clusters and return sum-of-squared error
        self._partition(self.tentative_clusters, self.tentative_objects)
        return sum([c.getSSE() for c in self.tentative_clusters])

    def _partition(self, cluster_list, others):
        """
        Partitions a set of clusters.
        :param cluster_list:    List of clusters with centroids.
        :param others:          Data points to partition.
        :return:                None
        """
        # Iterate over other data points
        for o in others:
            # Initialise an empty distance list
            dists = []
            # Iterate over clusters
            for cluster in cluster_list:
                # Get data point's distance from cluster and add it to list
                dists.append(self.metric(o, cluster.centroid.data))
            # Find the minimum distance
            min_dist = min(dists)
            # Find the index of the cluster with min distance
            min_index = dists.index(min_dist)
            # Add data point to nearest cluster
            cluster_list[min_index].add_object(o, min_dist)

