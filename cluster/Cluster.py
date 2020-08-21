class ClusterObj:
    """
    ClusterObj class, to represent data in Cluster.
    """
    def __init__(self, data, dist):
        """
        ClusterObj constructor.
        :param data:    Data that represents data point in cluster.
        :param dist:    Distance from cetnroid in cluster.
        """
        self.data = data
        self.dist = dist

    def __str__(self):
        newline = '\n'
        return f"{str(self.data).replace(newline, ' ')}" +\
               f"{' : ' + str(self.dist) if self.dist is not None else ''}"

    def __repr__(self):
        return self.__str__()


class Cluster:
    """
    Cluster class, used for k-medoid clustering.
    """
    def __init__(self, centroid):
        """
        Cluster constructor.
        :param centroid: Data that represents centroid data point.
        """
        self.centroid = ClusterObj(centroid, None)
        self.objects = []

    def add_object(self, data, dist):
        """
        Adds new ClusterObj to Cluster.
        :param data:    Data that represents data point in cluster.
        :param dist:    Distance from cetnroid in cluster.
        :return:        None
        """
        self.objects.append(ClusterObj(data, dist))

    def getSSE(self):
        """
        Caclulates sum-of-squares error between data points and centroid.
        :return:    Sum-of-squares error.
        """
        return sum([obj.dist ** 2 for obj in self.objects])

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"Cluster:\nCentroid: {self.centroid}\nObjects: {self.objects}\n"

