#!usr/bin/env python
import os
import DNAToolKit.DNAToolKit as dna

"""
Questions:

Question 1
----------
In a practical, we saw the scs function for findind the shortest common superstring
of a set of strings.

It's possible for there to be multiple different shortest common superstrings for 
the same set of input strings. Consider the input strings ABC, BCA, CAB. 
One shortest common superstring is ABCAB but another is BCABC and another is CABCA.

What is the length of the shortest common superstring of the following strings?

set = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]


Question 2
----------
How many different shortest common superstrings are there for the input strings 
given in the previous question?

Hint 1: You can modify the SCS function to keep track of this

Question 3
----------
Download this FASTQ file containing synthetic sequencing reads from a mystery virus:

All the reads are the same length (100 bases) and are exact copies of substrings from 
the forward strand of the virus genome.  You don't have to worry about sequencing 
errors, ploidy, or reads coming from the reverse strand.

Assemble these reads using one of the approaches discussed, such as greedy shortest 
common superstring.  Since there are many reads, you might consider ways to make 
the algorithm faster, such as the one discussed in the programming assignment in 
the previous module.

How many As are there in the full, assembled genome?

Hint: the virus genome you are assembling is exactly 15,894 bases long

Question 4
----------
How many Ts are there in the full, assembled genome from the previous question?

"""

if __name__ == "__main__":

    # # Question 1: -> CCTTGGATTGC

    # set = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]

    # scs, scs_list = dna.assembly.scs(set)
    # print(scs)

    # # Question 2: -> ['CCTTGCGATTGG', 'CCTTGCTGGATT', 'CCTTGGATTGC', 'GATTGCCTTGG', 'TGCCTTGGATT', 'TGGATTGCCTT']

    # print(scs_list)

    # Question 3:

    filepath = os.path.join(
        os.path.dirname(__file__),
        "data",
        "ads1_week4_reads.fq"
    )

    reads, _ = dna.read_files.read_fastq(filepath)

    # can't be bothered to implement the fast version now,
    # TODO LATER
    genome = dna.assembly.greedy_scs(reads, 30) 
    
    print(len(genome)) # 15894
    print(genome.count("A")) # 4633
    print(genome.count("T")) # 3723