"""Compare original DeepSEA DNA sequence with sequence built using build.py."""

import argparse
import logging
import numpy as np
from genomes import Hg19

from files import load_data_from_mat_or_npy

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument("--deepsea_data", default="", help="DeepSEA data file (.mat or .npy)", required=True)
parser.add_argument("--built_data", help="Data file built with build.py (.mat or .npy)", required=True)
args = parser.parse_args()

deepsea_data = load_data_from_mat_or_npy(args.deepsea_data)
if deepsea_data.shape[0]==1000:
    deepsea_data = np.moveaxis(deepsea_data, -1, 0)
if deepsea_data.shape[2]==1000:
    deepsea_data = np.moveaxis(deepsea_data, 1, -1)
logging.debug(f"DeepSEA data: {deepsea_data.shape}")

built_data = load_data_from_mat_or_npy(args.built_data)
logging.debug(f"Built data: {built_data.shape}")

if deepsea_data.shape != built_data.shape:
    print(f"Error. Shapes of both matrices should be same. {args.deepsea_data} shape is {deepsea_data.shape}. {args.built_data} shape is {built_data.shape}")
    exit()

n, w = built_data.shape[0], built_data.shape[1]

total_diff = np.sum(deepsea_data!=built_data)

total_diff = np.sum(deepsea_data!=built_data)
print(f"Total differences#:{total_diff:,} Total differences%: {100*total_diff/(n*w):.4f}%")




