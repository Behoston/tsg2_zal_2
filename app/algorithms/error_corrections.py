class CorrectedReads:
    def __init__(self, reads: [str], k=3, threshold=2):
        self.reads = reads
        self.k = k
        self.histogram = None
        self.threshold = threshold

    def kmerHist(self):
        """ Return k-mer histogram and average # k-mer occurrences """
        kmerhist = {}
        for read in self.reads:
            for kmer in [read[i:i + self.k] for i in range(len(read) - (self.k - 1))]:
                kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
        self.histogram = kmerhist

    def __iter__(self):
        self.kmerHist()
        self._iterator = iter(self.reads)
        return self

    def __next__(self):
        read = next(self._iterator)
        return self.correct1mm(read)

    def correct1mm(self, read):
        """ Return an error-corrected version of read.  k = k-mer length.
            kmerhist is kmer count map.  alpha is alphabet.  thresh is
            count threshold above which k-mer is considered correct. """
        # Iterate over k-mers in read
        for i in range(len(read) - (self.k - 1)):
            kmer = read[i:i + self.k]
            # If k-mer is infrequent...
            if self.histogram.get(kmer, 0) <= self.threshold:
                # Look for a frequent neighbor
                for newkmer in self.neighbors1mm(kmer, read):
                    if self.histogram.get(newkmer, 0) > self.threshold:
                        # Found a frequent neighbor; replace old kmer
                        # with neighbor
                        read = read[:i] + newkmer + read[i + self.k:]
                        break
        # Return possibly-corrected read
        return read

    def neighbors1mm(self, kmer, alpha):
        """ Generate all neighbors at Hamming distance 1 from kmer """
        neighbors = []
        for j in range(len(kmer) - 1, -1, -1):
            oldc = kmer[j]
            for c in alpha:
                if c == oldc:
                    continue
                neighbors.append(kmer[:j] + c + kmer[j + 1:])
        return neighbors
