from src.graph.Vertex import *


class Graph:
    """Directed Graph ADT"""

    def __init__(self, vertices=None):
        if vertices is None:
            vertices = []
        self.g = {v: Vertex(v) for v in vertices}

    def add_vertex(self, v):
        """
        Adds vertex given by v to graph
        :param v: (Any): Vertex name/data
        :return:
        """
        if v not in self.g:
            self.g[v] = Vertex(v)

    def add_vertices(self, *args):
        for v in args:
            if v not in self.g:
                self.g[v] = Vertex(v)

    def exists(self, v):
        """
        Checks that vertex given by v exists in graph
        :param v: (String): Vertex name
        :return: (Boolean): True if vertex exists, else False
        """
        return v in self.g

    def add_edge(self, u, v, weight=0.0, directed=True):
        """
        Adds a undirected edge from vertex u to v.
        :param u: (String): Source vertex
        :param v: (String): Destination vertex
        :param weight: (Float): Edge weight
        :param directed: (Boolean): Adds a directed edge if True,
            else undirected edge
        :return: None
        """
        if u in self.g and v in self.g:
            self.g[u].add_edge(endpoint=v,
                               weight=weight)
            if not directed:
                self.g[v].add_edge(endpoint=u,
                                   weight=weight)

    def add_edges(self, edges):
        for e in edges:
            u = e[0]
            v = e[1]
            weight = e[2] if len(e) > 2 else 0.0
            directed = e[3] if len(e) > 3 else True
            self.add_edge(u, v, weight, directed)

    def get_vertex(self, v):
        """
        Gets Vertex object with name v.
        :param v: (String): Vertex name
        :return: Vertex object with name v
        """
        try:
            return self.g[v]
        except:
            print("Vertex does not exist")

    def adjacent(self, u, v):
        """
        Checks to see if an edge from u to v exists
        :param u: (String) Source vertex name
        :param v: (String) Destination vertex name
        :return: (Boolean)
        """
        if u in self.g and v in self.g:
            return self.g[u].adjacent(v)
        return False

    def neighbours(self, v):
        """
        Returns neighbourhood of vertex v
        :param v: (String): Vertex name
        :return: Set of vertices neighbouring vertex v
        """
        if v in self.g:
            source = self.g[v]
            return set([n for n in source.edges])
        return set()

    def euler_circuit(self):
        """
        Searches for eulerian circuit in graph
        :return: List in vertices
        """

        def forward(vertex, edges_explored, edges_max):
            path = []
            current = (vertex, -1)
            while True:
                print(current)
                added = False
                path.append(current)
                current_v = self.g[current[0]]
                for neighbour_name, edges in current_v.edges.items():
                    for i, edge in enumerate(edges):
                        if edges_explored[current[0]][neighbour_name] < edges_max[current[0]][neighbour_name]:
                            edges_explored[current[0]][neighbour_name] += 1
                            current = neighbour_name, i
                            added = True
                            break
                    if added:
                        break
                if not added:
                    return path

        vertices = list(self.g.keys())
        current_vertex = None
        new_cycle, old_cycle = [], []
        edges_explored = {u: {v: 0 for v in self.g[u].edges} for u in vertices}

        edges_max = {u: {v: len(self.g[u].edges[v])
                         for v in self.g[u].edges} for u in vertices}

        print(forward(vertices[1], edges_explored, edges_max))
        while True:
            for v in vertices:
                found = False
                vertex_v = self.g[v]
                edgedict = vertex_v.edges
                for neighbour, edges in edgedict.items():
                    for i, edge in enumerate(edges):
                        if edges_explored[v][neighbour][i] != \
                                edges_max[v][neighbour][i]:
                            current_vertex = v
                            found = True
                        if found:
                            break
                    if found:
                        break
                if v == vertices[len(vertices) - 1] or current_vertex is None:
                    return old_cycle

            break

    def __add__(self, other):
        """
        Returns two graphs 'glued' together in the de Brujin sense
        :param other: Other graph
        :return: Glued graph
        """
        new_graph = Graph(list(self.g.keys()) + list(other.g.keys()))
        for u in self.g:
            vertex_u = self.g[u]
            for v in vertex_u.edges:
                edgelist = vertex_u.edges[v]
                for e in edgelist:
                    new_graph.add_edge(u, v, e.weight, True)
        for u in other.g:
            vertex_u = other.g[u]
            for v in vertex_u.edges:
                edgelist = vertex_u.edges[v]
                for e in edgelist:
                    new_graph.add_edge(u, v, e.weight, True)
        return new_graph

    def __iadd__(self, other):
        """
        Glues this graph with graph other, in-place
        :param other: Other graph to be glued with
        :return: None
        """
        self.add_vertices(tuple(other.g.keys()))
        for u in other.g:
            vertex_u = other.g[u]
            for v in vertex_u.edges:
                edgelist = vertex_u.edges[v]
                for e in edgelist:
                    self.add_edge(u, v, e.weight, True)
        return self
