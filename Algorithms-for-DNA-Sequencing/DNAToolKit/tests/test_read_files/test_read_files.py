import pytest 
import os
from DNAToolKit.DNAToolKit import read_files

class Tests:
    """Tests for read_files module"""
    @pytest.fixture(autouse=True)
    def setup(self):
        self.fasta_file = os.path.join(os.path.dirname(__file__), "data_examples", "lambda_virus.fa")
        self.fastq_file = os.path.join(os.path.dirname(__file__), "data_examples", "ERR037900_1.first1000.fastq")
    
    def test_read_fasta_reads_file(self):
        """Test if read_fasta reads the file at all"""
        genome = read_files.read_fasta(self.fasta_file)
        assert genome is not None
    
    def test_read_fasta_reads_correctly(self):
        """Test if read_fasta reads the file correctly"""
        genome = read_files.read_fasta(self.fasta_file)
        assert genome.startswith("GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG")

    def test_read_fastq_reads_file(self):
        """Test if read_fastq reads the file at all"""
        sequences, qualities = read_files.read_fastq(self.fastq_file)
        assert sequences is not None
        assert qualities is not None

    def test_read_fastq_reads_correctly(self):
        """Test if read_fastq reads the file correctly"""
        sequences, qualities = read_files.read_fastq(self.fastq_file)
        assert sequences[0] == "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCNAACCCTAACCCTAACCCTAACCCTAACCCTAAC"
        assert qualities[0] == "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGFHHHFHFFHHHHHGHHFHEH@4#55554455HGFBF<@C>7EEF@FBEDDD<=C<E"
