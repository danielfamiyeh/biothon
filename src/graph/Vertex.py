from collections import deque
from src.graph.Edge import *


class Vertex:
    """
    Graph node structure
    """
    def __init__(self, data=None):
        """
        Defalt constructor
        :param data: (Object)
        """
        self.data = self if data is None else data
        self.edges = {}

    def add_edge(self, **kwargs):
        """
        Adds an undirected edge between two vertices.
            Edge can be made directed.

        :param kwargs:
            endpoint (Vertex):  End point of the edge
            weight (Float): Cost of edge
            direct ed(Boolean): Edge is directed if True, else not
        :return:
        """
        endpoint = kwargs.get("endpoint", None)
        weight = kwargs.get("weight", 0)

        if endpoint is not None:
            if endpoint not in self.edges:
                self.edges[endpoint] = []
        self.edges[endpoint].append(Edge(endpoint, weight))

    def delete_edge(self, dest, delete_all=False):
        """
        Delete first or all edges to Vertex dest
        dest (Vertex): Destination vertex key
        delete_all (Boolean): Deletes every edge to dest if true, else first
        :return: None
        """
        if dest in self.edges:
            if delete_all or len(self.edges[dest]) == 1:
                del self.edges[dest]
            else:
                del self.edges[dest][0]

    def adjacent(self, endpoint):
        """
        :param endpoint: Destination vertex
        :return: (Boolean): True if an edge exists from self to endpoint else False
        """
        return endpoint in self.edges

    def deg(self):
        return sum([len(n) for n in self.edges])

    def __str__(self):
        return f"{'{'}\ndata: {self.data}\n" \
               f"edges: {[f'{repr(k)}: {str(self.edges[k])}' for k in self.edges]}\n" \
               f"{'}'}"


