#!usr/bin/env
import os
import regex as re
from DNAToolKit import DNAToolKit

filepath = os.path.join(os.path.dirname(__file__), "data/dna.fasta")
toolkit = DNAToolKit(filepath)

def main():

    # Q1
    # How many records are in the multi FASTA file?
    n = toolkit.record_count()
    print(f"There are {n} records")
    print("==================================================================\n")
    
    # Q2
    # What is the length of the longest sequence in the file?
    # Q3
    # What is the length of the shortest sequence in the file?
    shortest_sequence, longest_sequence = toolkit.seq_lenght()
    print(f"{shortest_sequence}\n")
    print(f"{longest_sequence}\n")
    print("==================================================================\n")
    # Q4
    # What is the longest ORF appearing in reading frame 2 of any of the sequences
    sequences = list(toolkit.dict.values())
    max_len = 0
    
    for sequence in sequences:
        orfs, len, _ = toolkit.orf_finder(sequence=sequence, frame=2)
        try:
            if len > max_len:
                max_len = len
        except TypeError: # raises when orf_finder() can't find an ORF
            continue
        
    print(f"The longest ORF in any sequence is {max_len}")
    print("==================================================================\n")
    # Q5
    # What is the starting position of the longest ORF in reading frame 3 in any 
    # of the sequences?
    
    sequences = list(toolkit.dict.values())
    max_len = 0
    start = 0
    
    for sequence in sequences:
        orfs, len, start_pos = toolkit.orf_finder(sequence=sequence, frame=3)
        try:
            if len > max_len:
                max_len = len
                start = start_pos
        except TypeError: # raises when orf_finder() can't find an ORF
            continue
    
    print(f"The starting position of the longest ORF in reading frame 3 is {start}")
    print("==================================================================\n") 
    
    # Q6
    # What is the length of the longest ORF appearing in any sequence and in any 
    # reading frame?
    
    sequences = sequences = list(toolkit.dict.values())
    frames = [1,2,3]
    max_len = 0
    
    for frame in frames:
        for sequence in sequences:
            orfs, len, _ = toolkit.orf_finder(sequence=sequence, frame=frame)
            try:
                if len > max_len:
                    max_len = len
            except TypeError: # raises when orf_finder() can't find an ORF
                continue
            
    print(f"The longest ORF in any sequence, in any reading frame is {max_len}")
    print("==================================================================\n") 
    
    # Q7
    # What is the length of the longest forward ORF that appears in the seq with 
    # identifier "gi|142022655|gb|EQ086233.1|16"

    
 
    # find the key with that identifier
    id = "gi|142022655|gb|EQ086233.1|16"
    sequences = [seq for key, seq in toolkit.dict.items() if id in key]
    
    # find longest ORF
    frames = [1,2,3]
    max_len = 0
    
    for frame in frames:
        for sequence in sequences:
            orfs, len, _ = toolkit.orf_finder(sequence=sequence, frame=frame)
            try:
                if len > max_len:
                    max_len = len
            except TypeError:
                continue
        
    print(f"The longest ORF in that sequence, in any reading frame is {max_len}")
    print("==================================================================\n") 
    
    # Q8
    # Find the most frequently occuring repeat of length 6 in all sequences. How
    # many times does it occur in all?

    max_repeat_global = 0
    most_frequent_repeat = ""

    # find the repeat that appears the most times in a sequence
    for id, seq in toolkit.dict.items():
        repeats = toolkit.repeats(seq, 6) # dict with all repeats of length 6 in sequence seq
        max_repeat = max(repeats.values()) # get number of times the most frequent repeat appears in that seq
        max_repeat_seq = max(repeats, key = repeats.get)
        if max_repeat > max_repeat_global:
            max_repeat_global = max_repeat
            most_frequent_repeat = max_repeat_seq
    
    # find how many times that repeat appears in all other sequences
    appearences = []
    for id, seq in toolkit.dict.items():
        appearences.append(re.findall(most_frequent_repeat, seq, overlapped=True))
    
    flat = [item for sublist in appearences for item in sublist]
    
    print(f"The most frequent repeat of length 6 is {most_frequent_repeat} repeats itself {len(flat)} times")
    print("==================================================================\n") 
    
    # Q9
    # Find all repeats of length 12 in the input file. Let's use Max to specify
    # the number of copies of the most frequent repeat of length 12. How many
    # different 12-base sequences occur Max times?
    
    Max = 0
    all_repeats = {}
    
    # find all repeats of length 12:
    for id, seq in toolkit.dict.items():
        repeats = toolkit.repeats(seq, 12)
        for repeat, n in repeats.items():
            try:
                all_repeats[repeat] += n
            except KeyError:
                all_repeats[repeat] = n # create key if it's the 1st time the repeat appears
    Max = max(all_repeats.values())
    
    # how many different repeats occur Max times?
    counter = 0
    repeats = []
    for repeat, n in all_repeats.items():
        if n == Max:
            counter += 1
            repeats.append(repeat)
        
    print(
        f"There are {counter} different repeats of length 12 that repeat themselves\
        {Max} times and they are: {repeats}"
    )
    print("==================================================================\n")         
    
    # Q10
    # Which one of the following repeats of length 7 has a maximum number of occurences
    
    candidates = {
        "AATGGCA": 0,
        "TGCGCGC": 0,
        "CATCGCC": 0,
        "CGCGCCG": 0
    }
    
    for id, seq in toolkit.dict.items():
        repeats = toolkit.repeats(seq, 7)
        for repeat, n in repeats.items():
            if repeat in candidates:
                candidates[repeat] += n
                
    print(f"The repeat that appears the most times is {max(candidates, key=candidates.get)}")
    print("==================================================================\n")

if __name__ == "__main__":
   main()