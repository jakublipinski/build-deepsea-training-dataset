"""Compare original DeepSEA labels with labels built using build.py."""

import argparse
import logging
import numpy as np

from files import load_labels_from_mat_or_npy

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument("--deepsea_labels", default="", help="DeepSEA labels file (.mat or .npy)", required=True)
parser.add_argument("--built_labels", help="Labels file built with build.py (.mat or .npy)", required=True)
args = parser.parse_args()

deepsea_labels = load_labels_from_mat_or_npy(args.deepsea_labels)
if deepsea_labels.shape[0]==919:
    deepsea_labels = np.moveaxis(deepsea_labels, 0, -1)

logging.debug(f"DeepSEA labels: {deepsea_labels.shape}")

built_labels = load_labels_from_mat_or_npy(args.built_labels)
logging.debug(f"Built labels: {built_labels.shape}")

if deepsea_labels.shape != built_labels.shape:
    print(f"Error. Shapes of both matrices should be same. {args.deepsea_labels} shape is {deepsea_labels.shape}. {args.built_labels} shape is {built_labels.shape}")
    exit()

n, l = built_labels.shape[0], built_labels.shape[1]

total_diff = 0
for i in range(l):
    diff = np.sum(deepsea_labels[...,i]!=built_labels[...,i])
    print(f"Label idx: {i} Differences#: {diff:,} Differences%: {100*diff/n:.2f}%")
    total_diff += diff
print(f"Total differences#: {total_diff:,} Total differences%: {100*total_diff/(n*l):.2f}%")




