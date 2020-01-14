"""Converts original dataset to npy file format."""

from scipy.io import loadmat
import h5py
import numpy as np

def load_mats(filename, key):
    """Load mat file."""
    try:
        f = h5py.File(filename, 'r')
    except OSError:
        # probably an old MATLAB file. try opening with scipy
        from scipy.io import loadmat
        f = loadmat(filename)
    return f[f"{key}xdata"], f[f"{key}data"]

train_data, train_labels = load_mats("../data/train.mat", "train")
valid_data, valid_labels = load_mats("../data/valid.mat", "valid")
test_data, test_labels = load_mats("../data/test.mat", "test")

# fix data structure. based on: https://github.com/zj-zhang/deepsea-keras/blob/master/read_data.py
train_data = np.moveaxis(train_data, -1, 0)
train_labels = np.moveaxis(train_labels, -1, 0)    
valid_data = np.moveaxis(valid_data, 1, -1)
test_data = np.moveaxis(test_data, 1, -1)

np.save("train_data.npy", train_data)
np.save("train_labels.npy", train_labels)
np.save("valid_data.npy", valid_data)
np.save("valid_labels.npy", valid_labels)
np.save("test_data.npy", test_data)
np.save("test_labels.npy", test_labels)