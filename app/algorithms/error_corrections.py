import collections
import itertools
from typing import Sequence, Iterator

from cached_property import cached_property


class CorrectedReads:
    def __init__(self, reads: Sequence[str], k=10, threshold=2):
        """
        :param reads: sequence of reads
        :param k: k-mer length
        :param threshold: if k-mer occur more often than threshold value,
               k-mer is considered correct. \n
               **Use higher values to replace k-mer more likely**
        """
        self.reads = reads
        self.k = k
        self.threshold = threshold

    def __iter__(self) -> Iterator[Sequence[str]]:
        return iter(self.corrected_reads)

    @cached_property
    def corrected_reads(self):
        result = []
        for read in self.reads:
            result.append(self.correct1mm(read))
        return result

    @cached_property
    def alphabet(self):
        return set(itertools.chain.from_iterable(self.reads))

    @cached_property
    def histogram(self):
        """ Build k-mer histogram and average # k-mer occurrences """
        histogram = {}
        for read in self.reads:
            for kmer in [read[i:i + self.k] for i in range(len(read) - (self.k - 1))]:
                histogram[kmer] = histogram.get(kmer, 0) + 1
        return histogram

    def plot_histogram(self):
        """**Require matplotlib!**"""
        try:
            from matplotlib import pyplot
        except ImportError:
            print("Please install 'matplotlib' to plot histogram!")
            return
        counter = collections.Counter(self.histogram.values())
        x, y = zip(*sorted(counter.items()))
        pyplot.plot(x, y)
        pyplot.show()

    def correct1mm(self, read):
        """ Return an error-corrected version of read. """
        # Iterate over k-mers in read
        for i in range(len(read) - (self.k - 1)):
            kmer = read[i:i + self.k]
            # If k-mer is infrequent...
            if self.histogram.get(kmer, 0) <= self.threshold:
                # Look for a frequent neighbor
                for new_kmer in self.neighbors1mm(kmer):
                    if self.histogram.get(new_kmer, 0) > self.threshold:
                        # Found a frequent neighbor; replace old kmer
                        # with neighbor
                        read = read[:i] + new_kmer + read[i + self.k:]
                        break
        # Return possibly-corrected read
        return read

    def neighbors1mm(self, kmer):
        """ Generate all neighbors at Hamming distance 1 from kmer """
        neighbors = []
        for j in range(len(kmer) - 1, -1, -1):
            oldc = kmer[j]
            for c in self.alphabet:
                if c == oldc:
                    continue
                neighbors.append(kmer[:j] + c + kmer[j + 1:])
        return neighbors
