import bisect

class SubseqIndex:
    """Holds a subsequence index for a text t
    
    k if a subsequence of t if the characters of k appear in t in the
    same order, but not necessarily consecutively:

    Example
    -------
    t = "AACCGGTT"
    k = "AAGT"

    t[0] + t[1] + t[5] + t[7] = k -> k is a subsequence of t

    Using subsequences tends to increase the specificity of the index hits
    """
    def __init__(self, t: str, k: int, interval: int):
        """ Create index from all subsequences consisting of k characters
            spaced interval positions apart.

        Parameters
        ----------
        t : str
            text to create index for
        k : int
            number of characters to match
        interval : int
            interval between characters. 1=adjacent, 2=every other char, etc.
        """
        self.k = k
        self.interval = interval
        self.index = []
        self.span = 1 + interval * (k-1) # length of the pattern to match
        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i+self.span], i))
        
        self.index.sort()

    def query(self, p: str) -> list:
        """Return index hits for the first subsequence of pattern p

        Parameters
        ----------
        p : str
            pattern to match agains a text

        Returns
        -------
        list
            index hits
        """
        subseq = p[:self.span:self.interval]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
