import itertools


class CorrectedReads:
    def __init__(self, reads: [str], k=3, threshold=10):
        """
        :param reads: list of reads
        :param k: k-mer length
        :param threshold: if k-mer occur more often than threshold value,
               k-mer is considered correct. \n
               **Use higher values to replace k-mer more likely**
        """
        self.reads = reads
        self.k = k
        self.histogram = {}
        self.threshold = threshold
        self.alphabet = set()

    def __iter__(self):
        self.build_alphabet()
        self.build_kmer_histogram()
        self._iterator = iter(self.reads)
        return self

    def __next__(self):
        read = next(self._iterator)
        return self.correct1mm(read)

    def build_alphabet(self):
        self.alphabet = set(itertools.chain.from_iterable(self.reads))

    def build_kmer_histogram(self):
        """ Build k-mer histogram and average # k-mer occurrences """
        for read in self.reads:
            for kmer in [read[i:i + self.k] for i in range(len(read) - (self.k - 1))]:
                self.histogram[kmer] = self.histogram.get(kmer, 0) + 1

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
