from typing import Sequence

from algorithms.error_corrections import CorrectedReads
from io_utils import parse_input


class DeBruijnGraph:
    """ A de Bruijn multigraph built from a collection of strings.
        User supplies strings and k-mer length k.  Nodes of the de
        Bruijn graph are k-1-mers and edges correspond to the k-mer
        that joins a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(st, k):
        """ Chop a string up into k mers of given length """
        for i in range(0, len(st) - (k - 1)):
            yield (st[i:i + k], st[i:i + k - 1], st[i + 1:i + k])

    def __init__(self, strIter, k):
        """ Build de Bruijn multigraph given string iterator and k-mer
            length k """
        self._super_string = None
        self.graph = {}  # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        for st in strIter:
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = Node(km1R)
                nodeL.number_of_outgoing += 1
                nodeR.number_of_incoming += 1
                self.graph.setdefault(nodeL, []).append(nodeR)
        # Iterate through nodes and tally how many are balanced,
        # semi-balanced, or neither
        self.nodes_semi_balanced, self.nodes_balanced, self.nodes_not_balanced = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian path (not cycle)
        self.head, self.tail = None, None
        for node in self.nodes.values():
            if node.is_balanced():
                self.nodes_balanced += 1
            elif node.is_semi_balanced():
                if node.number_of_incoming == node.number_of_outgoing + 1:
                    self.tail = node
                if node.number_of_incoming == node.number_of_outgoing - 1:
                    self.head = node
                self.nodes_semi_balanced += 1
            else:
                self.nodes_not_balanced += 1

    def has_eulerian_path(self):
        """ Return true iff graph has Eulerian path. """
        return self.nodes_not_balanced == 0 and self.nodes_semi_balanced == 2

    def has_eulerian_cycle(self):
        """ Return true iff graph has Eulerian cycle. """
        return self.nodes_not_balanced == 0 and self.nodes_semi_balanced == 0

    def is_eulerian(self):
        """ Return true iff graph has Eulerian path or cycle """
        return self.has_eulerian_path() or self.has_eulerian_cycle()

    def eulerian_path(self):
        """ Find and return Eulerian path or cycle (as appropriate) """
        assert self.is_eulerian()
        graph = self.graph
        if self.has_eulerian_path():
            graph = graph.copy()
            assert self.head is not None
            assert self.tail is not None
            graph.setdefault(self.tail, []).append(self.head)
        # graph graph has an Eulerian cycle
        tour = []
        src = next(iter(graph.keys()))  # pick arbitrary starting node

        def __visit(node):
            while len(graph[node]) > 0:
                dst = graph[node].pop()
                __visit(dst)
            tour.append(node)

        __visit(src)
        tour = tour[::-1][:-1]

        if self.has_eulerian_path():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        # Return node list
        return map(str, tour)

    @property
    def super_string(self):
        if not self._super_string:
            for point in self.eulerian_path():
                if not self._super_string:
                    self._super_string = point
                else:
                    self._super_string += point[-1]
        return self._super_string


class Node:
    """ Node in a de Bruijn graph, representing a k-1 mer.  We keep
        track of # of incoming/outgoing edges so it's easy to check
        for balanced, semi-balanced. """

    def __init__(self, km1mer):
        self.km1mer = km1mer
        self.number_of_incoming = 0
        self.number_of_outgoing = 0

    def is_semi_balanced(self):
        return abs(self.number_of_incoming - self.number_of_outgoing) == 1

    def is_balanced(self):
        return self.number_of_incoming == self.number_of_outgoing

    def __hash__(self):
        return hash(self.km1mer)

    def __repr__(self):
        return self.km1mer


def do_assembly(data: Sequence[str]):
    graph = DeBruijnGraph(data, 10)
    return graph.super_string


if __name__ == '__main__':
    # DEBUG

    list_graph = DeBruijnGraph(["konstantynopoli", "politanczykowianeczka"], 5)
    print(list_graph.super_string)

    data = parse_input('./sample_data/reads_1_percent_bad.fasta')
    for i in range(10, 20):
        corrected_reads = CorrectedReads(data, k=i, threshold=5)
        try:
            print(f"Trying {i}...")
            graph = DeBruijnGraph(corrected_reads, i)
            print(graph.super_string)
            print(f"Success!!!!!!!!!! {i}")
        except Exception:
            pass
