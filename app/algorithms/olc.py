import logging
import random
import re
from typing import Sequence

from cached_property import cached_property


class Graph:
    def __init__(self):
        self.nodes = {}

    def add_node(self, node):
        self.nodes[node.value] = node

    def __iter__(self):
        return iter(self.nodes.values())

    def __len__(self):
        return len(self.nodes)

    def get_random_node_but_not(self, not_nodes: {'Node'}) -> 'Node':
        return random.choice(list(self.node_values_set - not_nodes))

    @cached_property
    def node_values_set(self):
        return set(self.nodes.values())

    def get_node_greatest_number_of_reachable_out(self, unreachable: {'Node'}):
        return max(self.node_values_set - unreachable, key=lambda node: len(node.out_nodes_set - unreachable))


class Node:
    def __init__(self, value):
        self.value = value
        self.out = {}

    def add_edge_with_value(self, node, value):
        self.out[node] = value

    def get_next_node_and_overlap(self, visited: {'Node'}):
        for node in self.out_nodes_sorted_by_value:
            if node not in visited:
                return node, self.out[node]
        raise Exception("Integrity exception!")

    def has_non_visited_out(self, visited: {'Node'}) -> bool:
        return bool(self.out_nodes_set - visited)

    @cached_property
    def out_nodes_set(self) -> set('Node'):
        return set(self.out.keys())

    @cached_property
    def out_nodes_sorted_by_value(self):
        return [key for key, item in sorted(self.out.items(), key=lambda x: x[1], reverse=True)]

    def __eq__(self, other):
        return self.value == other.value

    def __hash__(self):
        return hash(self.value)


def olc(data: Sequence[str]):
    overlap_graph = overlap_naive(data)
    sequence = naive_graph_path(overlap_graph)
    return max(sequence, key=lambda x: len(max(x, key=len)))


def olc_suffix(data: Sequence[str]):
    overlap_graph = overlap_suffix(data)
    contigs = layout(overlap_graph)
    sequence = consensus(contigs)
    return sequence


def olc_dynamic(data: Sequence[str]):
    overlap_graph = overlap_dynamic(data)
    contigs = layout(overlap_graph)
    sequence = consensus(contigs)
    return sequence


def overlap_naive(data: Sequence[str]):
    logging.info("Building graph.")
    minimum_overlap_size = 6
    graph = Graph()
    for read in data:
        graph.add_node(Node(read))
    for read_a in graph:
        sufix = read_a.value[-minimum_overlap_size:]
        sufix_pattern = re.compile(f'.*{sufix}')
        for read_b in graph:
            if read_a == read_b:
                continue
            for prefix in re.findall(sufix_pattern, read_b.value):
                if read_a.value.endswith(prefix):
                    read_a.add_edge_with_value(read_b, len(prefix))
                    break
    logging.info("Graph has been built!")
    return graph


def naive_graph_path(graph):
    for node in graph:
        super_strings = naive_graph_path_staring_from_node(graph, node)
        if super_strings:
            yield super_strings


def naive_graph_path_staring_from_node(graph: Graph, start_node: Node) -> [str]:
    read_length = len(start_node.value)
    minimal_super_string_length = len(graph) * read_length * 0.01
    super_strings = []
    visited = set()
    node = start_node
    super_strings.append(node.value)
    visited.add(node)
    graph_len = len(graph)
    while len(visited) != graph_len:
        while node.has_non_visited_out(visited):
            node, overlap = node.get_next_node_and_overlap(visited)
            super_strings[-1] += node.value[overlap:]
            visited.add(node)
        if len(visited) != graph_len:
            node = graph.get_node_greatest_number_of_reachable_out(visited)
            super_strings.append(node.value)
            visited.add(node)
    super_strings = [super_string for super_string in super_strings if len(super_string) > minimal_super_string_length]
    logging.info(f"Generated {len(super_strings)} contigs!")
    return super_strings


def overlap_suffix(data: Sequence[str]):
    raise NotImplementedError


def overlap_dynamic(data: Sequence[str]):
    raise NotImplementedError


def layout(overlap_graph):
    raise NotImplementedError


def consensus(contigs):
    raise NotImplementedError


if __name__ == '__main__':
    # DEBUG
    logging.basicConfig(level=logging.DEBUG)
    from io_utils import parse_input, dump_output
    from algorithms.error_corrections import CorrectedReads

    data = parse_input('./sample_data/reads_1_percent_bad.fasta')
    super_string = olc(CorrectedReads(data))
    print(len(max(super_string, key=len)))
    dump_output('./sample_data/con.fasta', super_string)
