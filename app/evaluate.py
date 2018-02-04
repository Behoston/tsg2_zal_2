"""
Created on Mon Dec 18 11:58:16 2017

@author: norbert
"""

from collections import namedtuple

import pysam
from math import log

MINLEN = 300

EvaluationResult = namedtuple('EvaluationResult', [
    'reference_coverage',
    'reads_coverage',
    'alignment_errors',
    'identity_score',
    'number_of_alignments',
    'fragmentation_score',
    'overall_score',
])


def evaluate(sam_data):  # noqa C901
    print(sam_data)
    samfile = pysam.AlignmentFile(sam_data)
    reftotlen = sum(samfile.lengths)
    rdstotlen = 0

    reads = []
    fragments = []
    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            rdstotlen += read.query_length
        else:
            if not read.is_secondary:
                rdstotlen += read.infer_read_length()
            if read.reference_length >= MINLEN:
                reads.append(read)
                fragments.append((read.reference_start, read.reference_end))
    samfile.close()

    contig_overlaps = []
    for i, frag1 in enumerate(fragments):
        for frag2 in fragments[i + 1:]:
            s, e = max(frag1[0], frag2[0]), min(frag1[1], frag2[1])
            if s < e:
                contig_overlaps.append((s, e))
    contig_overlaps.sort()

    redundant = []
    ls, le = -1, -1
    for cs, ce in contig_overlaps:
        if cs <= le:
            le = max(le, ce)
        else:
            if le >= 0:
                redundant.append((ls, le))
            ls, le = cs, ce
    if le >= 0:
        redundant.append((ls, le))

    almtotlen = 0
    almmmcount = 0
    almcount = 0
    for read in reads:
        s, e = read.reference_start, read.reference_end
        alms = s
        alpairs = read.get_aligned_pairs(True, True)
        for rs, re in redundant + [(e, e)]:
            if alms < rs:
                alme = min(e, rs)
                if alme - alms >= MINLEN:
                    almtotlen += alme - alms
                    almcount += 1
                    i = alms - s
                    while alpairs[i][1] is None or alpairs[i][1] < alms:
                        i += 1
                    while i < len(alpairs) and (alpairs[i][1] is None or alpairs[i][1] < alme):
                        if alpairs[i][1] is None or alpairs[i][2] is None or alpairs[i][2].islower():
                            almmmcount += 1
                        i += 1
            if re < e:
                alms = max(alms, re)
            else:
                break

    almtotlen = almtotlen
    refcoverage = almtotlen / reftotlen if almtotlen else 0
    rdscoverage = almtotlen / rdstotlen if almtotlen else 0
    ident_score = max(0.5, 1 - 10 * almmmcount / almtotlen) if almtotlen else 0
    count_score = 1 / log(4 + almcount, 5) if almtotlen else 0

    score = refcoverage * rdscoverage * ident_score * count_score
    return EvaluationResult(
        reference_coverage=refcoverage,
        reads_coverage=rdscoverage,
        alignment_errors=almmmcount / almtotlen if almtotlen else 0,
        identity_score=ident_score,
        number_of_alignments=almcount,
        fragmentation_score=count_score,
        overall_score=score,
    )


if __name__ == '__main__':
    import sys

    result = evaluate(sys.stdin)
    print("Pokrycie referencji:", result.reference_coverage)
    print("Pokrycie odczytów:", result.reads_coverage)
    print("Błędy uliniowień:", result.alignment_errors)
    print("Ocena identyczności:", result.identity_score)
    print("Liczba uliniowień:", result.number_of_alignments)
    print("Ocena rozdrobnienia:", result.fragmentation_score)
    print("Łączna ocena:", result.overall_score)
