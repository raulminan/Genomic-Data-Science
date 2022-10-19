#!usr/bin/env python
from DNAToolKit import DNAToolKit as dna

x = dna.assembly.overlap("CGG", "GGT", 2)

if __name__ == "__main__":
    print(x)