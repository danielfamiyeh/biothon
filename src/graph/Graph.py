from enum import Enum

class _SPF(Enum):
    """
    Arbitrary shortest path vertices enum.
    """
    SOURCE = 0
    TARGET = 1

class Vertex:
    """
    Vertex class.
    """
    def __init__(self, n, v):
        """
        Vertex constructor.
        :param n:   Name of vertex.
        :param v:   Vertex value.
        """
        self. name = n
        self.value = v
        self.edges = {} # The edge family is meant to be a nonempty set but we'll allow that ;)
        self.pred = None

    def __repr__(self):
        return f"V({self.name}:{self.value} E:{self.edges})"

class Graph:
    """
    Graph class for a simple, weighted graph.
    """
    def __init__(self):
        self._V = {}
        self.num_edges = 0
        self.num_vertices = 0

    def __repr__(self):
        as_string = "Graph\n"
        for i, v in enumerate(self._V):
            as_string += str(self._V[v]) + "\n"
        return as_string

    def adjacent(self, x, y):
        """
        Tests to see if an edge exists from vertex x to vertex y.
        :param x:   Name of source vertex.
        :param y:   Name of destination vertex.
        :return:    True if vertices are adjacent, else False.
        """
        return True if y in self._V[x].edges else False

    def neighbours(self, x):
        """
        Returns a list of neighbours of vertex x.
        :param x:   Name of verte
        :return:
        """
        return [neighbour for neighbour in self._V[x]]

    def add_vertex(self, x, v=0):
        """
        Adds a vertex with name x and value v, if vertex with name x does
            not exist in graph.
        :param x:   Name of vertex to add to graph.
        :param v:   Value of vertex x to add to graph.
        :return:    None
        """
        if x not in self._V:
            self._V[x] = Vertex(x, v)
        self.num_vertices += 1

    def remove_vertex(self, x):
        """
        Removes vertex with name x from graph.
        :param x:   Name of vertex to remove from graph.
        :return:    None
        """
        if x in self._V:
            del self._V[x]

        for v in self._V:
            for neighbour in self._V[v]:
                if neighbour[0].value == x:
                    self._V[v].edges.remove(neighbour)

    def add_edge(self, x, y, w=0, directed=True):
        """
        Adds an edge from vertex with name x to vertex with name y
            with a weight of 0 to graph. Edge can either be directed
            or undirected.
        :param x:           Name of source vertex
        :param y:           Name of destination vertex
        :param w:           Weight of edge between the two vertices
        :param directed:    Adds a directed edge if True, else an undirected edge.
        :return:            None
        """
        if x in self._V and y in self._V:
            self._V[x].edges[y] = w
            if not directed:
                self._V[y].edges[x] = w

    def remove_edge(self, x, y):
        """
        Removes an edge from vertex with name x to vertex with name y,
            should such an edge exist.
        :param x:   Name of source vertex
        :param y:   Name of destination vertex
        :return:    None
        """
        if x in self._V and y in self._V:
            source = self._V[x]
            if y in source.edges:
                del source.edges[y]

    def get_vertex_value(self, x):
        """
        Returns vertex value of vertex with name x.
        :param x:   Vertex name
        :return:    Value of vertex with name x
        """
        return self._V[x].value if x in self._V else None

    def set_vertex_value(self, x, v):
        """
        Sets vertex value of vertex with name x to value v
        :param x:   Vertex name
        :param v:   New value of vertex
        :return:    None
        """
        if x in self._V:
            self._V[x].value = v

    def get_edge_value(self, x, y):
        """
        Returns the weight of an edge from a vertex with name x to vertex
            with name y.
        :param x:   Source vertex name
        :param y:   Destination vertex name
        :return:    Edge weight
        """
        if x in self._V and y in self._V:
            if y in self._V[x].edges:
                return self._V[x].edges[y]

    def set_edge_value(self, x, y, w):
        """
        Sets the value of edge from vertex with name x to vertex with name y
            to the value given by w, should such an edge exist.
        :param x:   Name of source vertex
        :param y:   Name of destination vertex
        :param w:   Value of new weight of edge between the two/
        :return:    None
        """
        if x in self._V and y in self._V:
            if y in self._V[x].edges:
                self._V[x].edges[y] = w
            else:
                print(f"Edge does not exist between {x} and {y}.")

    def _top_sort(self, s, visited, stack):
        """
        Internal topological sorting function.
        :param s:       Vertex currently being processed.
        :param visited: List containing visited vertices.
        :param stack:   Stack containing vertices to visit.
        :return:        None
        """
        visited[s] = True
        for v in self._V[s].edges:
            if not visited[v]:
                self._top_sort(v, visited, stack)
        stack.append(s)

    def top_sort(self):
        """
        External topological sorting function.
        :return:    List contining names of veritces sorte topologically.
        """
        visited = {v: False for v in self._V}
        stack = []

        for v in self._V:
            if not visited[v]:
                self._top_sort(v, visited, stack)
        return stack

    def chaining(self):
        """
        The chaining algorithm required ClustalW MSA.
        :return:
        """
        # Initialise arbitrary source and target vertices
        source = Vertex(_SPF.SOURCE, 0)
        target = Vertex(_SPF.TARGET, 0)
        # Final chain
        chain = []

        # Iterate over vertex map
        for name, v in self._V.items():
            # Assign zero from source
            source.edges[name] = 0
            # Assign zero cost to target
            v.edges[_SPF.TARGET] = 0

        # Insert source vertex in to vertex map
        self._V[_SPF.SOURCE] = source
        # Insert target vertex in vertex map
        self._V[_SPF.TARGET] = target

        # Iterate over every vertex in vertex map
        for name, v in self._V.items():
            # Iterate over all edges of vertex
            for n, w in v.edges.items():
                # If neighbour's distance is less than
                # current vertex's distance + cost to neighbour
                if self._V[n].value < (v.value + w):
                    # Set neighbour's distance
                    self._V[n].value = v.value + w
                    # Set neighbour's predecessor
                    self._V[n].pred = v.name

        # Initialise current to target's predecessor
        current = self._V[_SPF.TARGET].pred

        # iterate backwards through predecessor chain
        while current is not None:
            # Add predeessor name vertex to chain
            chain.append(current)
            # Set current to predecessor
            current = self._V[current].pred

        return chain






