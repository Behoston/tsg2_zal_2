def olc_suffix(data: [str]):
    overlap_graph = overlap_suffix(data)
    contigs = layout(overlap_graph)
    sequence = consensus(contigs)
    return sequence


def olc_dynamic(data: [str]):
    overlap_graph = overlap_dynamic(data)
    contigs = layout(overlap_graph)
    sequence = consensus(contigs)
    return sequence


def overlap_suffix(data: [str]):
    pass


def overlap_dynamic(data: [str]):
    pass


def layout(overlap_graph):
    pass


def consensus(contigs):
    pass
