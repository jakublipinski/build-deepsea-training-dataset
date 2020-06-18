"""Build the DeepSEA dataset."""

import csv

import argparse
import logging
import numpy as np

from features import Features
from genomes import Hg19, chr_nos_sorted_by_name, CHR_7, CHR_8, CHR_9
from files import save_to_mat

def fill(features, hg19, bin_start, chr_no, bin_size, window_size, complementary_sequence, \
         create_data, create_labels, complementary_shift, data, labels, idx, debug_writer=False):
    """Fill data and labels with DNA sequence and features vector.""" 
    debug_row = dict() if debug_writer else None
    if create_data:
        hg19.fill_window(data[idx], chr_no, bin_start - (window_size-bin_size)//2, window_size, debug_row=debug_row)
    if create_labels:
        features.fill_labels(labels[idx], chr_no, bin_start, debug_row=debug_row)
    if debug_writer:
        debug_writer.writerow(debug_row)

    if complementary_sequence:
        debug_row = dict() if debug_writer else None
        if create_data:
            features.fill_labels(labels[idx+complementary_shift], chr_no, bin_start, debug_row=debug_row)
        if create_labels:
            hg19.fill_window(data[idx+complementary_shift], chr_no, bin_start - (window_size-bin_size)//2, window_size, complementary=True, debug_row=debug_row)
        if debug_writer:
            debug_writer.writerow(debug_row)

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument("--metadata_file", help="Metadata (.tsv) file in the format described at https://www.encodeproject.org/help/batch-download/", required=True)
parser.add_argument("--beds_folder", help="Path to the folder containing all the .bed files", required=True)
parser.add_argument("--bin_size", type=int, default=200, help="Bin size (default: 200)")
parser.add_argument("--window_size", type=int, default=1000, help="Window size (default: 1000)")
parser.add_argument("--train_size", type=int, help="Size of the training set. Will be multiplied by 2 if complementary sequence is True (default: number of all bins minus size of valid set)")
parser.add_argument("--valid_size", type=int, help="Size of the validation set. Will be multiplied by 2 if complementary sequence is True. (default: valid_ratio * number of all bins). Validation set is located at chromosome 7 starting from 30,508,751")
parser.add_argument("--valid_ratio", type=float, default=.05, help="Ration of validation / train set. Ignored if valid_size is provided")
parser.add_argument("--test_size", type=int, help="Size of the test set. Will be multiplied by 2 if complementary sequence is True (default: number of all bins in chromosome 8 and 9)")
parser.add_argument("--pos", help="If provided loads bin positions from an external file")
parser.add_argument("--hg19", help="hg19 genome file")
parser.add_argument("--train_data_filename", help="Output train data filename (.npy format)")
parser.add_argument("--train_labels_filename", help="Output train labels filename (.npy format)")
parser.add_argument("--valid_data_filename", help="Output valid data filename (.npy format)")
parser.add_argument("--valid_labels_filename", help="Output valid labels filename (.npy format)")
parser.add_argument("--test_data_filename", help="Output test data filename (.npy format)")
parser.add_argument("--test_labels_filename", help="Output test labels filename (.npy format)" )
parser.add_argument("--train_filename", help="Output train dataset filename (.mat format)")
parser.add_argument("--valid_filename", help="Output valid dataset filename (.mat format)" )
parser.add_argument("--test_filename", help="Output test dataset filename (.mat format)")
parser.add_argument("--filter", action="append", help="Add condition to filter the metadata file eg. -f\"Lab=Michael Snyder, Stanford\" -f\"Assembly=hg19\"")
parser.add_argument("--signal_threshold", type=float, help="Signal threshold from which to add positions from the bed file")
parser.add_argument("--complementary_sequence", type=bool, default=True, help = "Add complementary sequence (default=True)")
parser.add_argument("--save_debug_info", type=bool, default=False, help="Save debug info to debug_train.tsv, debug_valid.tsv and debug_test.tsv files")
args = parser.parse_args()

# change filter from command line arguments -f\"Lab=Michael Snyder, Stanford\" -f\"Assembly=hg19\" to {"Lab":"Michael Snyder, Stanford", "Assembly":"hg19"}
filter = {f.split("=")[0]:f.split("=")[1] for f in args.filter} if args.filter else None

create_train_data = args.train_data_filename or args.train_filename
create_train_labels = args.train_labels_filename or args.train_filename
create_valid_data = args.valid_data_filename or args.valid_filename
create_valid_labels = args.valid_labels_filename or args.valid_filename
create_test_data = args.test_data_filename or args.valid_filename
create_test_labels = args.test_labels_filename or args.valid_filename
create_data = create_train_data or create_valid_data or create_test_data
create_labels = create_train_labels or create_valid_labels or create_test_labels

features = Features(args.metadata_file, args.bin_size)
features.read_beds(args.beds_folder, filter, args.signal_threshold)

hg19 = None
if create_data:
    hg19 = Hg19(args.hg19)
    hg19.read_genome()

if args.pos:
    features.read_bins(args.pos)
bins = features.bins

train_size = features.no_of_bins()
if args.train_size:
    train_size = min(train_size, args.train_size)

test_size = args.test_size or features.no_of_bins([CHR_8, CHR_9])

if args.valid_size:
    valid_size = args.valid_size
    train_size = min(train_size, features.no_of_bins() - valid_size - test_size, features.no_of_bins() - valid_size - features.no_of_bins([CHR_8, CHR_9]))
else:
    train_size -= test_size
    valid_size = int(args.valid_ratio * train_size)
    train_size -= valid_size

if train_size + valid_size + test_size > features.no_of_bins():
    print(f"Sum of training ({train_size:,}), validation ({valid_size:,}) and test ({test_size:,}) sets exceeds total number of bins ({features.no_of_bins():,})")
    exit()

revx = 2 if args.complementary_sequence else 1

labels_size = features.no_of_labels()

logging.info(f"train size: {train_size:,} valid size: {valid_size:,} test size: {test_size:,} bins: {features.no_of_bins():,} labels: {labels_size:,} samples: {features.no_of_samples():,}")

debug_train, debug_valid, debug_test = None, None, None
if args.save_debug_info:
    fieldnames=["Chr_No", "Chr", "Start", "End", "Seq"]+features.accession_ids()
    debug_train = csv.DictWriter(open("debug_train.tsv", "w"), fieldnames=fieldnames)
    debug_train.writeheader()
    debug_valid = csv.DictWriter(open("debug_valid.tsv", "w"), fieldnames=fieldnames)
    debug_valid.writeheader()
    debug_test = csv.DictWriter(open("debug_test.tsv", "w"), fieldnames=fieldnames)
    debug_test.writeheader()

train_data = np.ndarray(shape=(train_size * revx, args.window_size, 4), dtype=np.uint8)
train_labels = np.zeros(shape=(train_size * revx, labels_size), dtype=np.uint8)
valid_data = np.ndarray(shape=(valid_size * revx, args.window_size, 4), dtype=np.uint8)
valid_labels = np.zeros(shape=(valid_size * revx, labels_size), dtype=np.uint8)
test_data = np.ndarray(shape=(test_size * revx, args.window_size, 4), dtype=np.uint8)
test_labels = np.zeros(shape=(test_size * revx, labels_size), dtype=np.uint8)

train_idx = 0
valid_idx = 0
test_idx = 0

if create_labels or create_data:
    i=0
    for chr_no in chr_nos_sorted_by_name():
        for bin_start in sorted(bins[chr_no]):
            if chr_no in [CHR_8, CHR_9]:
                if test_idx < test_size:
                    fill(features, hg19, bin_start, chr_no, args.bin_size, args.window_size, args.complementary_sequence, \
                        create_test_data, create_test_labels, test_size, test_data, test_labels, test_idx, debug_test)
                    test_idx += 1
                    i+=1
            else:
                if chr_no == CHR_7 and bin_start >= 30508751 and valid_idx<valid_size:
                    fill(features, hg19, bin_start, chr_no, args.bin_size, args.window_size, args.complementary_sequence, \
                        create_valid_data, create_valid_labels, valid_size, valid_data, valid_labels, valid_idx, debug_valid)
                    valid_idx += 1
                    i+=1
                elif train_idx<train_size:
                    fill(features, hg19, bin_start, chr_no, args.bin_size , args.window_size, args.complementary_sequence, \
                        create_train_data, create_train_labels, train_size, train_data, train_labels, train_idx, debug_train)
                    train_idx += 1
                    i+=1
                
            if i%1000 == 0:
                logging.debug(f"{100*i/(train_size+valid_size+test_size):.2f}%")

    logging.debug("Saving files...")

    if args.train_filename:
        save_to_mat(args.train_filename, train_data, train_labels, "train")
    if args.valid_filename:
        save_to_mat(args.valid_filename, valid_data, valid_labels, "valid")
    if args.test_filename:
        save_to_mat(args.test_filename, test_data, test_labels, "test")

    if args.train_data_filename:
        np.save(args.train_data_filename, train_data)
    if args.train_labels_filename:
        np.save(args.train_labels_filename, train_labels)
    if args.valid_data_filename: 
        np.save(args.valid_data_filename, valid_data)
    if args.valid_labels_filename:
        np.save(args.valid_labels_filename, valid_labels)
    if args.test_data_filename:
        np.save(args.test_data_filename, test_data)
    if args.test_labels_filename:
        np.save(args.test_labels_filename, test_labels)

logging.debug("Done")
