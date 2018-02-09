import re
from typing import Sequence


class Graph:
    def __init__(self):
        self.nodes = {}

    def add_node(self, node):
        self.nodes[node.value] = node

    def __getitem__(self, item):
        return self.nodes[item]

    def __iter__(self):
        return iter(self.nodes.values())


class Node:
    def __init__(self, value):
        self.value = value
        self.out = {}

    def add_edge_with_value(self, node, value):
        self.out[node] = value

    def get_next_node_and_overlap(self):
        return max(self.out.items(), key=lambda x: x[1])

    def __eq__(self, other):
        return self.value == other.value

    def __hash__(self):
        return hash(self.value)


def olc(data: Sequence[str]):
    overlap_graph = overlap_naive(data)
    sequence = naive_graph_path(overlap_graph)
    return max(sequence, key=len)


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
    minimum_overlap_size = 6
    graph = Graph()
    for read in data:
        graph.add_node(Node(read))
    for read_a in graph:
        sufix = read_a.value[-minimum_overlap_size:]
        sufix_pattern = f'.*{sufix}'
        for read_b in graph:
            if read_a == read_b:
                continue
            for prefix in re.findall(sufix_pattern, read_b.value):
                if read_a.value.endswith(prefix):
                    read_a.add_edge_with_value(read_b, len(prefix))
                    break
    return graph


def naive_graph_path(graph):
    for node in graph.nodes.values():
        visited = set()
        super_string = node.value
        visited.add(node)
        while node.out:
            next_node, overlap = node.get_next_node_and_overlap()
            if next_node not in visited:
                node = next_node
                super_string += node.value[overlap:]
                visited.add(node)
            else:
                node.out.pop(next_node)
        yield super_string


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
    from io_utils import parse_input, dump_output

    data = parse_input('./sample_data/reads_1_percent_bad.fasta')
    super_string = olc(data)
    print(len(super_string))
    dump_output('./sample_data/con.fasta', super_string)
