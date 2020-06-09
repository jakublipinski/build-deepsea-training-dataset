"""Chromatin features."""
import csv
import os
import re
import logging
from collections import namedtuple

from genomes import chr_name_to_no, MAX_CHR

ChrRegion = namedtuple("ChrRegion", ["start", "end"])

class Feature:
    """A class representing a chromatin feature (e.g. CTCF, RBBP5) and its binding location in a particular reference genome."""

    def __init__(self, bin_size, accession_id):
        """Create an empty Feature object.

        Parameters
        ----------
        bin_size : int
            The size of the bin
        accession_id : string
            Feature Accession id
        """
        self.accession_id = accession_id
        self.__chr_regions = [list() for chr in range(MAX_CHR)]
        self.__bin_size = bin_size
        self.__bins = None

    def __add_region(self, chr_no, start, end):
        """Add genome position to the list of genome regions."""
        self.__chr_regions[chr_no].append(ChrRegion(start, end))

    def read_bed(self, filename, signal_threshold = None):
        """Read the narrow peaks from the bed file and add them to the list of genome regions.
        
        Parameters
        ----------
        Filename : str
            Bed filename
        signal_threshold : float (default=None)
            The signal threshold to add a region from the bed file. Add all if not provided
        """
        with open(filename) as bedfile:
            logging.debug(f"Reading regions from {filename}")
            reader = csv.reader(bedfile, delimiter='\t')
            for row in reader:
                try:
                    chr_no = chr_name_to_no(row[0])
                except ValueError:
                    continue
                if not signal_threshold or float(row[6]) >= signal_threshold:
                    self.__add_region(chr_no, int(row[1]), int(row[2]))
        self.__bins = None

    @property
    def bins(self):
        """Translate feature genome regions to the list of bins according to the following algorithm from the paper:.

        To prepare the input for the deep convolutional network model, we split the genome into 200-bp bins. For each bin 
        we computed the label for all 919 chromatin features; a chromatin feature was labeled 1 if more than half of the 
        200-bp bin is in the peak region and 0 otherwise.
        """
        if self.__bins:
            return self.__bins
        self.__bins = [set() for _ in range(MAX_CHR)]
        for chr_no in range(MAX_CHR):
            for region in self.__chr_regions[chr_no]:
                idx = region.start//self.__bin_size
                while True:
                    bin_start = idx * self.__bin_size
                    bin_end = bin_start + self.__bin_size
                    if max(0, min(bin_end, region.end) - max(bin_start, region.start)) > self.__bin_size //2:
                        self.__bins[chr_no].add(bin_start)
                    if bin_end >= region.end:
                        break
                    idx += 1
        return self.__bins

    def is_feature_in_bin(self, chr_no, bin_start):
        """Determine if the chromatin feature exists in the particular bin."""
        return bin_start in self.bins[chr_no]

    def no_of_samples(self):
        """Return number of samples analyzed."""
        return sum([len(regions) for regions in self.__chr_regions])


class Features:
    """The list of chromatin features."""

    def __init__(self, metadata_filename, bin_size):
        """Create an empty chromatin feature list based on the metadata file."""
        self.__metadata = list()
        self.__features = list()
        self.__bin_size = bin_size
        self.__bins = None
        # open the file as a tsv file in the format described at https://www.encodeproject.org/help/batch-download/
        with open(metadata_filename) as tsvfile:
            reader = csv.DictReader(tsvfile, dialect='excel-tab')
            for row in reader:
                self.__metadata.append(row)

    def read_beds(self, beds_folder, filter, signal_threshold=None):
        """Read the narrow peaks from the bed files.

        Parameters
        ----------
        beds_folder : str
            The path to the folder where all the bed files are located
        filter : set
            The set of filters to apply to the metadata file e.g. {"Lab":"Michael Snyder, Stanford", "Assembly":"hg19"}
        signal_threshold : float (default=None)
            The signal threshold to add a region from the bed file. Add all if not provided
        """ 
        for row in self.__metadata:
            if not filter or filter.items() <= row.items(): # checks if row contains filter
                accession_id = row["File accession"]
                filename = os.path.join(beds_folder, f"{accession_id}")
                feature = Feature(self.__bin_size, accession_id)
                try:
                    feature.read_bed(filename, signal_threshold=signal_threshold)
                except FileNotFoundError:
                    filename += ".bed"
                    feature.read_bed(filename, signal_threshold=signal_threshold)
                self.__features.append(feature)
        self.__bins = None

    def read_bins(self, bins_file):
        """Read list of bins from the file."""
        self.__bins = [set() for _ in range(MAX_CHR)]
        with open(bins_file) as binsfile:
            reader = csv.reader(binsfile, delimiter='\t')
            for row in reader:                
                chr_no = chr_name_to_no(row[0])
                bin_start = int(row[1])
                self.__bins[chr_no].add(bin_start)

    @property
    def bins(self):
        """Return the list of all the bins for all the chromatin features."""
        if self.__bins:
            return self.__bins
        self.__bins = [set() for _ in range(MAX_CHR)]
        for feature in self.__features:
            feature_bins = feature.bins
            for chr_no in range(MAX_CHR):
                self.__bins[chr_no] |= feature_bins[chr_no]
        return self.__bins

    def is_feature_in_bin(self, feature_no, chr_no, bin_start):
        """Determine if the particular chromatin feature exists in the particular bin."""
        return self.__features[feature_no].is_feature_in_bin(chr_no, bin_start)

    def fill_labels(self, label, chr_no, bin_start, debug_row=None):
        """Fill the label vector with the information on which chromatin features are located in the particular bin."""
        for i in range(len(self.__features)):
            label[i] = 1 if self.is_feature_in_bin(i, chr_no, bin_start) else 0
            if debug_row is not None:
                debug_row[self.__features[i].accession_id] = label[i]

    def no_of_samples(self):
        """Return number of samples analyzed."""
        return sum([feature.no_of_samples() for feature in self.__features])

    def no_of_bins(self, chromosomes=None):
        """Return number of distinct bins."""
        if not chromosomes:
            return sum([len(bin) for bin in self.__bins])
        else:
            return sum([len(self.__bins[i]) for i in chromosomes])

    def no_of_labels(self):
        """Return number of labels."""
        return len(self.__features)

    def accession_ids(self):
        """Return list of accession ids."""
        return [feature.accession_id for feature in self.__features]