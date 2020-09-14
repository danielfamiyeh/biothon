import unittest

from src.graph.Graph import *


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.graph = Graph()
        self.vertex_names = ["A", "B", "C"]
        self.e1 = 2
        self.e2 = -3
        self.e3 = 3.25
        self.vertex_count = 0

    def testConstructor(self):
        self.assertEqual(self.vertex_count,
                         self.graph.num_vertices)

    def testAddVertex(self):
        for name in self.vertex_names:
            self.graph.add_vertex(name)
            self.vertex_count += 1

            self.assertEqual(self.vertex_count,
                             self.graph.num_vertices)
            self.assertTrue(name in self.graph.get_vertices())
            self.assertTrue(self.graph.vertex_exists(name))

    def testVertexExists(self):
        self.graph.add_vertex(self.vertex_names[0])
        self.assertTrue(self.graph.vertex_exists(
            self.vertex_names[0]
        ))

    def testAddEdge(self):
        for name in self.vertex_names:
            self.graph.add_vertex(name)

        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[1], self.e1, True)
        self.graph.add_edge(self.vertex_names[1],
                            self.vertex_names[2], self.e2, True)
        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[2], self.e3, False)

        vertices = self.graph.get_vertices()
        self.assertTrue(self.vertex_names[1] in
                        vertices[self.vertex_names[0]].edges)
        self.assertTrue(self.vertex_names[2] in
                        vertices[self.vertex_names[1]].edges)
        self.assertTrue(self.vertex_names[2] in
                        vertices[self.vertex_names[0]].edges and
                        self.vertex_names[0] in
                        vertices[self.vertex_names[2]].edges)

    def testAdjacent(self):
        for name in self.vertex_names:
            self.graph.add_vertex(name)

        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[1], self.e1, True)
        self.graph.add_edge(self.vertex_names[1],
                            self.vertex_names[2], self.e2, True)
        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[2], self.e3, False)

        self.assertTrue(self.graph.adjacent(
            self.vertex_names[0], self.vertex_names[1]
        ))
        self.assertTrue(self.graph.adjacent(
            self.vertex_names[1], self.vertex_names[2]
        ))
        self.assertFalse(self.graph.adjacent(
            self.vertex_names[1], self.vertex_names[0]
        ))
        self.assertFalse(self.graph.adjacent(
            self.vertex_names[2], self.vertex_names[1]
        ))
        self.assertTrue(self.graph.adjacent(
            self.vertex_names[0], self.vertex_names[2]
        ))
        self.assertTrue(self.graph.adjacent(
            self.vertex_names[2], self.vertex_names[0]
        ))

    def testNeighbours(self):
        for name in self.vertex_names:
            self.graph.add_vertex(name)

        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[1], self.e1, True)
        self.graph.add_edge(self.vertex_names[1],
                            self.vertex_names[2], self.e2, True)
        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[2], self.e3, False)

        v1_neighbours = self.graph.neighbours(
            self.vertex_names[0]
        )

        v2_neighbours = self.graph.neighbours(
            self.vertex_names[1]
        )

        v3_neighbours = self.graph.neighbours(
            self.vertex_names[2]
        )

        self.assertTrue(
            self.vertex_names[1] in v1_neighbours and
            self.vertex_names[2] in v1_neighbours and
            not self.vertex_names[0] in v1_neighbours
        )

        self.assertTrue(
            self.vertex_names[2] in v2_neighbours and
            not (self.vertex_names[0] in v2_neighbours and
                 self.vertex_names[1] in v2_neighbours)
        )

        self.assertTrue(
            self.vertex_names[0] in v3_neighbours and
            not (self.vertex_names[2] in v3_neighbours and
                 self.vertex_names[1] in v3_neighbours)
        )

    def testRemoveVertex(self):
        for name in self.vertex_names:
            self.graph.add_vertex(name)

        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[1], self.e1, True)
        self.graph.add_edge(self.vertex_names[1],
                            self.vertex_names[2], self.e2, True)
        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[2], self.e3, False)

        self.graph.remove_vertex(self.vertex_names[0])

        self.assertFalse(self.graph.vertex_exists(
            self.vertex_names[0]
        ))

        vertices = self.graph.get_vertices()

        for vertex in vertices:
            self.assertFalse(self.vertex_names[0] in vertices[vertex].edges)

    def testRemoveEdge(self):
        self.graph.add_vertex(self.vertex_names[0])
        self.graph.add_vertex(self.vertex_names[1])
        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[1], 3)

        self.graph.remove_edge(self.vertex_names[0],
                               self.vertex_names[1])

        self.assertFalse(self.graph.adjacent(self.vertex_names[0],
                                             self.vertex_names[1]))

    def testGetVertexValue(self):
        for i, name in enumerate(self.vertex_names):
            self.graph.add_vertex(name, i * (-1 if i % 2 == 0 else 1))
        self.assertEqual(0, self.graph.get_vertex_value(
            self.vertex_names[0]
        ))

        self.assertEqual(1, self.graph.get_vertex_value(
            self.vertex_names[1]
        ))

        self.assertEqual(-2, self.graph.get_vertex_value(
            self.vertex_names[2]
        ))

    def testSetVertexValue(self):
        for i, name in enumerate(self.vertex_names):
            self.graph.add_vertex(name, i * (-1 if i % 2 == 0 else 1))

        for i, name in enumerate(self.vertex_names):
            self.graph.set_vertex_value(name, 3-i)
            self.assertEqual(3-i, self.graph.get_vertex_value(name))

    def testGetEdgeValue(self):
        self.graph.add_vertex(self.vertex_names[0])
        self.graph.add_vertex(self.vertex_names[1])
        self.graph.add_edge(self.vertex_names[0],
                       self.vertex_names[1], 3)
        self.assertEqual(3, self.graph.get_edge_value(
            self.vertex_names[0], self.vertex_names[1]
        ))

    def testSetEdgeValue(self):
        '''
        Tests that an edge that exists between two vertices
            in the graph can have its' value changed.
        '''
        self.graph.add_vertex(self.vertex_names[0])
        self.graph.add_vertex(self.vertex_names[1])
        self.graph.add_edge(self.vertex_names[0],
                            self.vertex_names[1], -7)

        self.graph.set_edge_value(self.vertex_names[0],
                                  self.vertex_names[1], 4)

        self.assertEqual(4, self.graph.get_edge_value(self.vertex_names[0],
                                                      self.vertex_names[1]))

    def testTopologicalSort(self):
        '''
        Topological sort test using example from
            geeksforgeeks.org
        '''

    def testChainingSolver(self):
        '''
        ClustalW chaining problem solver test
        '''
        self.assertEqual(True, False)

if __name__ == '__main__':
    unittest.main()
