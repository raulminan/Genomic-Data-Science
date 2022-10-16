#!/usr/bin/env python
"""This script contains the homework for week 1 of the course"""
import os
import matplotlib.pyplot as plt
from HWWeek1 import HWWeek1


def main():
    filename = os.path.join(os.path.dirname(__file__), "data", "lambda_virus.fa")
    week1 = HWWeek1(filename)
    genome = week1.read_genome()

    # Q1
    # How many times does AGGT or its revers complement occur in the lambda virus
    # genome?
    pattern_1 = 'AGGT'
    occurences_1 = week1.naive_with_rc(pattern_1, genome)
    
    print(len(occurences_1))
    
    # Q2
    # How many times does AGGT or its revers complement occur in the lambda virus
    # genome?

    pattern_2 = 'TTAA'
    occurences_2 = week1.naive_with_rc(pattern_2, genome)
    
    print(len(occurences_2))
    
    # Q3
    # What is the leftmost occurrence of ACTAAGT or its reverse complement in the 
    # Lambda virus genome
    
    pattern_3 = 'ACTAAGT'
    occurences_3 = week1.naive_with_rc(pattern_3, genome)
    print(occurences_3[0])

    
    # Q4
    # What is the leftmost occurrence of AGTCGA or its reverse complement in the 
    # Lambda virus genome
    
    pattern_4 = 'AGTCGA'
    occurences_4 = week1.naive_with_rc(pattern_4, genome)
    print(occurences_4[0])
    
    
    # Q5
    # How many times does TTCAAGCC occur in the Lambda virus genome when 
    # allowing up to 2 mismatches? 
    
    pattern_5 = 'TTCAAGCC'
    occurences_5 = week1.naive_2mm(pattern_5, genome)
    print(len(occurences_5))
    
    # Q6
    # What is the offset of the leftmost occurrence of AGGAGGTT
    # in the Lambda virus genome when allowing up to 2 mismatches?
    
    pattern_6 = 'AGGAGGTT'
    occurences_6 = week1.naive_2mm(pattern_6, genome)
    print(occurences_6[0])
    
    # Q7
    # This dataset has something wrong with it; one of the sequencing cycles is 
    # poor quality. Report which sequencing cycle has the problem.  Remember that a 
    # sequencing cycle corresponds to a particular offset in all the reads. For 
    # example, if the leftmost read position seems to have a problem consistently 
    # across reads, report 0. If the fourth position from the left has the problem, 
    # report 3. Do whatever analysis you think is needed to identify the bad cycle. 
    # It might help to review the "Analyzing reads by position" video.
    
    filename2 = os.path.join(os.path.dirname(__file__),"data", "ERR037900_1.first1000.fastq")
    q7 = HWWeek1(filename2)

    genome, _ = q7.read_fastq()
    gc = q7.find_GC_by_pos(genome)

    plt.plot(range(len(gc)), gc)
    plt.show()
    
    error = [i for i in range(len(gc)) if gc[i] < 0.3]
    print(error)

if __name__ == "__main__":
    main()
