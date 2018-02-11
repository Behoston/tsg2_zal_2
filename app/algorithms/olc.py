import logging
import pickle
import random
import re
from itertools import count
from typing import Sequence

from cached_property import cached_property

from utils import timing


class Graph:
    def __init__(self):
        self.nodes = {}

    def add_node(self, node):
        self.nodes[node.id] = node

    def __iter__(self):
        return iter(self.nodes.values())

    def __len__(self):
        return len(self.nodes)

    def __bool__(self):
        return bool(self.nodes)

    def __getitem__(self, item: int or 'Node') -> 'Node':
        if isinstance(item, int):
            return self.nodes[item]
        else:
            return self.nodes[item.id]

    def remove_node(self, node: int or 'Node'):
        node = self[node]
        self.nodes.pop(node.id)
        for out_id in node.out.keys():
            self[out_id].entries.pop(node.id)
        for entry_id in node.entries.keys():
            self[entry_id].out.pop(node.id)

    def get_random_node(self) -> 'Node':
        return random.choice(list(self.nodes.values()))

    def get_node_greatest_number_of_out(self):
        return max(self.nodes.values(), key=lambda node: len(node.out))

    def get_node_with_smallest_number_of_entries(self):
        return min(self.nodes.values(), key=lambda node: len(node.entries))

    def remove_edges_can_be_inferred_1(self):
        for node in self.nodes.values():
            for maybe_inferrable_node_id in reversed(node.out_nodes_sorted_by_value):
                for following_node_id in node.out_nodes_sorted_by_value:
                    following_node = self[following_node_id]
                    if following_node_id in node.out and maybe_inferrable_node_id in following_node.out:
                        node.out.pop(maybe_inferrable_node_id)
                        self[maybe_inferrable_node_id].entries.pop(node.id)
                        break

    @property
    def average_node_value_length(self):
        return sum([len(node.value) for node in self]) / len(self)


class Node:
    _id = count()

    def __init__(self, value):
        self.id = next(self._id)
        self.value = value
        self.out = {}
        self.entries = {}

    def add_edge_with_value(self, node, value):
        self.out[node.id] = value
        node.entries[self.id] = value

    def get_next_node_id_and_overlap(self):
        return max(self.out.items(), key=lambda x: x[1])

    @property
    def has_out(self) -> bool:
        return bool(self.out)

    @cached_property
    def out_nodes_sorted_by_value(self):
        return [key for key, item in sorted(self.out.items(), key=lambda x: x[1], reverse=True)]

    def __eq__(self, other):
        return self.value == other.value

    def __hash__(self):
        return hash(self.value)


def olc_naive(data: Sequence[str]):
    overlap_graph = overlap_naive(data)
    sequence = naive_graph_path(overlap_graph)
    return max(sequence, key=lambda x: len(max(x, key=len)))


def olc(data: Sequence[str]):
    overlap_graph = overlap_naive(data)
    sequences = layout(overlap_graph)
    return sequences


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
    pickled_graph = pickle.dumps(graph)
    for node in graph:
        super_strings = naive_graph_path_staring_from_node(pickle.loads(pickled_graph), node)
        if super_strings:
            yield super_strings


def naive_graph_path_staring_from_node(graph: Graph, start_node_id: int) -> [str]:
    node = graph[start_node_id]
    read_length = len(node.value)
    minimal_super_string_length = len(graph) * read_length * 0.01
    super_strings = []

    graph.remove_node(node)
    super_strings.append(node.value)
    while len(graph):
        while node.has_out:
            node_id, overlap = node.get_next_node_id_and_overlap()
            node = graph[node_id]
            graph.remove_node(node)
            super_strings[-1] += node.value[overlap:]
        if len(graph):
            node = graph.get_node_greatest_number_of_out()
            super_strings.append(node.value)
            graph.remove_node(node)
    super_strings = [super_string for super_string in super_strings if len(super_string) > minimal_super_string_length]
    logging.info(f"Generated {len(super_strings)} contigs!")
    return super_strings


def overlap_suffix(data: Sequence[str]):
    raise NotImplementedError


def overlap_dynamic(data: Sequence[str]):
    raise NotImplementedError


def layout(overlap_graph: Graph):
    minimal_super_string_length = len(overlap_graph) * overlap_graph.average_node_value_length * 0.01
    overlap_graph.remove_edges_can_be_inferred_1()
    node = overlap_graph.get_node_with_smallest_number_of_entries()
    overlap_graph.remove_node(node)
    super_strings = [node.value]
    while overlap_graph:
        while node.has_out:
            node_id, overlap = node.get_next_node_id_and_overlap()
            node = overlap_graph[node_id]
            overlap_graph.remove_node(node)
            super_strings[-1] += node.value[overlap:]
        if overlap_graph:
            node = overlap_graph.get_node_with_smallest_number_of_entries()
            overlap_graph.remove_node(node)
            super_strings.append(node.value)
    super_strings = [super_string for super_string in super_strings if len(super_string) > minimal_super_string_length]
    return super_strings


def consensus(contigs):
    raise NotImplementedError


if __name__ == '__main__':
    # DEBUG
    logging.basicConfig(level=logging.DEBUG)
    from io_utils import parse_input, dump_output
    from algorithms.error_corrections import CorrectedReads

    data = parse_input('./sample_data/reads_1_percent_bad.fasta')
    with timing():
        super_string = olc(CorrectedReads(data))
    print(len(max(super_string, key=len)))
    dump_output('./sample_data/con.fasta', super_string)
