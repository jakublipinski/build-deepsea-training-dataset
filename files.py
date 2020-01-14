"""Saves and loads data and labels from .npy and .mat files."""
import h5py
import numpy as np

def guess_key(f):
    """Guess key name from the file.keys()."""
    for key in f.keys():
        if key.endswith("xdata"):
            return key[:-5]
    return ""

def load_from_mat_or_npy(filename, key=None, load_data=True, load_labels = True):
    """Load labels and/or data from the .mat or .npy file."""
    data = None
    labels = None
    if filename.endswith(".mat"):
        try:
            f = h5py.File(filename, 'r')
        except OSError:
            # probably an old MATLAB file. try opening with scipy
            from scipy.io import loadmat
            f = loadmat(filename)
        key = key or guess_key(f)
        if load_data:
            data = f[f"{key}xdata"]
        if load_labels:
            labels = f[f"{key}data"]
            print(labels.shape)
    else:
        if load_data:
            data = np.load(filename)
        if load_labels:
            labels = np.load(filename)
    return data, labels

def load_data_from_mat_or_npy(filename, key=None):
    """Load data from the .mat or .npy file."""
    data, _ = load_from_mat_or_npy(filename, key, load_labels=False)
    return data

def load_labels_from_mat_or_npy(filename, key=None):
    """Load labels from the .mat or .npy file."""
    _, labels = load_from_mat_or_npy(filename, key, load_data=False)
    return labels

def save_to_mat(filename, data, labels, key):
    """Save data and labels to .mat file."""
    with h5py.File(filename, 'w') as f:
        f.create_dataset(f"{key}xdata",  data=data)
        f.create_dataset(f"{key}data",  data=labels)