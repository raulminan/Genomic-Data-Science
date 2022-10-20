#!usr/bin/env python
import os
from DNAToolKit import DNAToolKit as dna

"""
# Question 1
------------
What is the edit distance of the best match between pattern GCTGATCGATCGTACG and
and the excerpt of human chromosome 1? 

# Question 2
------------
What is the edit distance of the best match between pattern GATTTACCAGATTGAG and 
the excerpt of human chromosome 1?

# Question 3
------------
[...]
Picture the overlap graph corresponding to the overlaps just calculated.  
How many edges are in the graph?  In other words, how many distinct pairs of 
reads overlap?

# Question 4
------------

Picture the overlap graph corresponding to the overlaps computed for the previous 
question. How many nodes in this graph have at least one outgoing edge?  
(In other words, how many reads have a suffix involved in an overlap?)


"""


if __name__ == "__main__":
    
    # Question 1 -> 3

    genome = dna.read_files.read_fasta(
        os.path.join(
            os.path.dirname(__file__), 
            "data", 
            "chr1.GRCh38.excerpt.fasta"
            )
        )
    
    # pattern = "GCTGATCGATCGTACG"
    # distance = dna.utils.edit_distance_fewest_edits(pattern, genome)
    # print(f"The edit distance between the best match between the pattern and the excerpt of human chromosome 1 is {distance}")


    # Question 2 -> 2

    pattern = "GATTTACCAGATTGAG"
    # distance = dna.utils.edit_distance_fewest_edits(pattern, genome)
    # print(f"The edit distance between the best match between the pattern and the excerpt of human chromosome 1 is {distance}")


    # Question 3 -> 1729056

    reads, _ = dna.read_files.read_fastq(
        os.path.join(
            os.path.dirname(__file__),
            "data",
            "ERR266411_1.for_asm.fastq"
        )
    )

    map, graph = dna.assembly.overlap_map(reads, min_length=30)

    print(len(map))

    # Question 4 -> 
    
    print(len(graph))