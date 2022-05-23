import unittest

from src.graph.Graph import *


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.graph = Graph()

    def test_add_vertex(self):
        self.graph.add_vertex("A")
        self.graph.add_vertex(1)

        self.assertTrue("A" in self.graph.g)
        self.assertTrue(1 in self.graph.g)
        self.assertFalse(-3 in self.graph.g)

    def test_add_vertices(self):
        self.graph.add_vertices("A", 0, "loool", -3)
        self.assertIn("A", self.graph.g)
        self.assertIn(0, self.graph.g)
        self.assertIn("loool", self.graph.g)
        self.assertIn(-3, self.graph.g)

    def test_exists(self):
        self.graph.add_vertex("node")
        self.assertTrue(self.graph.exists("node"))
        self.assertFalse(self.graph.exists(1))

    def test_add_edge(self):
        self.graph.add_vertex("A")
        self.graph.add_vertex("B")
        self.graph.add_vertex("C")

        self.graph.add_edge("A", "B")
        self.graph.add_edge("C", "A", -1, True)

        self.assertIn("B",
                      self.graph.g["A"].edges)
        self.assertIn("A",
                      self.graph.g["B"].edges)
        self.assertIn("A",
                      self.graph.g["C"].edges)
        self.assertNotIn("C",
                         self.graph.g["A"].edges)

    def test_adjacent(self):
        self.graph.add_vertex("A")
        self.graph.add_vertex("B")
        self.graph.add_vertex("C")

        self.graph.add_edge("A", "B")
        self.graph.add_edge("C", "A", 0, True)

        self.assertTrue(self.graph.adjacent("A", "B"))
        self.assertTrue(self.graph.adjacent("B", "A"))
        self.assertTrue(self.graph.adjacent("C", "A"))
        self.assertFalse(self.graph.adjacent("A", "C"))

    def test_neighbors(self):
        self.graph.add_vertex("A")
        self.graph.add_vertex("B")
        self.graph.add_vertex("C")

        self.graph.add_edge("A", "B")
        self.graph.add_edge("B", "C")
        self.graph.add_edge("C", "A", 0, True)

        self.assertIn("C", self.graph.neighbours("B"))
        self.assertIn("A", self.graph.neighbours("B"))

    def test_glue(self):
        other_graph = Graph(["A", "C", "D"])
        self.graph.add_vertices("A", "B", "C")

        other_graph.add_edge("A", "C", 3)

        self.graph.add_edge("A", "B", 2)
        self.graph.add_edge("A", "C", -1)

        glued_graph = self.graph + other_graph

        self.assertIn("e:3", str(glued_graph.g["A"].edges["C"]))
        self.assertIn("e:-1", str(glued_graph.g["A"].edges["C"]))

        # In-place
        self.graph += other_graph
        self.assertIn("e:3", str(self.graph.g["A"].edges["C"]))
        self.assertIn("e:-1", str(self.graph.g["A"].edges["C"]))

    def test_euler_circuit(self):
        self.graph.add_vertices("A", "B", "C", "D", "E", "F", "G")
        self.graph.add_edges([("G", "A"), ("A", "B"), ("B", "C"),
                              ("B", "G"), ("G", "E"), ("C", "G"),
                              ("E", "F"), ("E", "D"), ("D", "B"),
                              ("F", "D")])

        self.graph.euler_circuit()

if __name__ == '__main__':
    unittest.main()
