"""Reference Genomes."""
import re
import numpy as np
import logging

from Bio import SeqIO

MAX_CHR = 24 #number of chromosomes. 0-21 autosomes. 22 - chrX. 23 - chrY
CHR_8 = 8-1
CHR_9 = 9-1
CHR_X = MAX_CHR-2
CHR_Y = MAX_CHR-1

def chr_name_to_no(chr_name):
    """Translate chromosome name e.g. "chr1", "chrX" to 0-indexed chromosome number. "chrX" and "chrY" are translated to 22 and 23 respectively."""
    if chr_name == "chrX":
        return CHR_X
    if chr_name == "chrY":
        return CHR_Y
    match = re.match(r"chr(\d+)", chr_name)
    if match:
        return int(match.groups()[0]) - 1
    raise ValueError(f"'{chr_name}' is not a valid chromosome name")

def chr_nos_sorted_by_name():
    """Iterate through chromosome numbers in the order of their names e.g. "chr1", "chr10", "chr11",..."""
    chr_names = [f"chr{no+1}" for no in range(CHR_X)] + ["chrX", "chrY"]
    for chr_name in sorted(chr_names):
        yield chr_name_to_no(chr_name)

class Genome:
    """A class to represent a reference genome."""

    # a matrix to quickly translate the base name into a one-hot vector
    __base2mat = np.zeros([90, 4], np.uint8)
    __base2mat[ord('A')] = np.array([1, 0, 0, 0])
    __base2mat[ord('G')] = np.array([0, 1, 0, 0])
    __base2mat[ord('C')] = np.array([0, 0, 1, 0])
    __base2mat[ord('T')] = np.array([0, 0, 0, 1])
    # reverse strand
    __base2mat_rev = np.zeros([90, 4], np.uint8)
    __base2mat_rev[ord('A')] = np.array([0, 0, 0, 1])
    __base2mat_rev[ord('G')] = np.array([0, 0, 1, 0])
    __base2mat_rev[ord('C')] = np.array([0, 1, 0, 0])
    __base2mat_rev[ord('T')] = np.array([1, 0, 0, 0])

    def __init__(self, genome_filename):
        """Create an empty object."""
        self.__genome_filename = genome_filename
        self.__seq = ["" for _ in range(MAX_CHR)]

    def read_genome(self):
        """Read the chromosome sequences into the __seq list."""
        logging.debug(f"Reading genome from {self.__genome_filename}")
        for record in SeqIO.parse(self.__genome_filename, "fasta"):
            chr_no = self.id_to_chr_no(record.id)
            if chr_no >=0:
                logging.debug(f"Reading chromosome# {chr_no+1} from {self.__genome_filename}")
                self.__seq[chr_no] = record.seq.upper()

    def id_to_chr_no(self, id):
        """Translate the chromosome id into its corresponding number. Must be overridden."""
        raise NotImplementedError("id_to_chr_no() should be overloaded")

    def fill_window(self, data, chr_no, start, window_size, complementary=False):
        """Fill one-hot vector corresponding to the genome sequence.

        Parameters
        ----------
        data : ndarray
            The array to fill
        chr_no : int
            The chromosome number
        start : int
            The position in the sequence
        window_size : int
            The size of the window to fill
        complementary : bool (default=False)
            Whether to fill complementary strand (reverse order)
        """
        for i in range(window_size):
            start = min(start, len(self.__seq[chr_no])-window_size)
            start = max(start, 0)
            if  not complementary:
                data[i]=self.__base2mat[ord(self.__seq[chr_no][start+i])]
            else:
                data[i]=self.__base2mat_rev[ord(self.__seq[chr_no][start+window_size-i-1])]
        

class Hg19(Genome):
    """The Hg19 assembly reference genome."""

    def __init__(self, genome_filename):
        """Create Hg19 empty genome object."""
        super().__init__(genome_filename)

    def id_to_chr_no(self, id):
        """Translate the record id into a corresponding chromosome number. Considers fully assembled chromosome only."""
        if id in ["chrX", "chrY"]:
            return MAX_CHR-1
        z = re.match(r"chr(\d+)$", id)
        if not z:
            return -1
        return int(z.groups()[0]) - 1
