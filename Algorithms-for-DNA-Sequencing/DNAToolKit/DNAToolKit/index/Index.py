import bisect

class Index:
    """Holds a substring index for a text"""

    def __init__(self, t: str, k: int) -> None:
        """Creates an index from all substrings of length k in text t

        Parameters
        ----------
        t : str
            text to create a substring index from
        k : int
            length of substrings
        """
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        self.index.sort()

    def query(self, p: str) -> list:
        """Returns index hits for the first k-mer of p

        Parameters
        ----------
        p : str
            pattern to match agains a text

        Returns
        -------
        list
            list of index hits
        """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
        