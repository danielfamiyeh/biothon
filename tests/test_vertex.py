import unittest

from src.graph.Vertex import *


class TestVertex(unittest.TestCase):
    def setUp(self):
        self.vertex_a = Vertex()
        self.vertex_b = Vertex("B")
        self.vertex_c = Vertex("C")

    def test_vertex_init(self):
        # V1.1: Vertex with no data data initialises to reference
        self.assertEqual(self.vertex_a.data, self.vertex_a)
        # V2.2: Vertex with data
        self.assertEqual(Vertex(3).data, 3)

    def test_vertex_add_edge(self):
        # V2.1: Unweighted, undirected edge (ab)<->(ba)
        self.vertex_a.add_edge(endpoint=self.vertex_b)
        # V2.2: Weighted, undirected edge (ca)<->(ac)
        self.vertex_c.add_edge(endpoint=self.vertex_a, weight=10)
        # V2.3: Weighted, directed edge (bc)
        self.vertex_b.add_edge(endpoint=self.vertex_c, weight=-5, directed=True)

        # V2.1: Test
        self.assertTrue(self.vertex_b in self.vertex_a.edges)
        self.assertTrue(self.vertex_a in self.vertex_b.edges)

        # V2.2: Test
        self.assertTrue(self.vertex_a in self.vertex_c.edges)
        self.assertTrue(self.vertex_c in self.vertex_a.edges)
        self.assertEqual(self.vertex_c.edges[self.vertex_a][0].weight, 10)
        self.assertEqual(self.vertex_a.edges[self.vertex_c][0].weight, 10)

        # V2.3: Test
        self.assertTrue(self.vertex_c in self.vertex_b.edges)
        self.assertFalse(self.vertex_b in self.vertex_c.edges)


if __name__ == '__main__':
    unittest.main()
