import bisect

class Index(object):
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))  # tuple
        self.index.sort()
    
    def query(self, p):
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1)) # binary search
        # print(i)
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i+k:i+len(p)]:  # verify whether the rest of p match the extension of t
            offsets.append(i)
    return offsets


if __name__ == "__main__":
    t = 'GCTACGATCTAGAATCTA'
    p = 'TCTA'

    index = Index(t, 2)
    print(queryIndex(p, t, index))