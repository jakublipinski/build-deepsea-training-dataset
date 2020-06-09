# Build DeepSEA training dataset

[DeepSEA](http://deepsea.princeton.edu/) is a deep learning-based algorithmic framework for predicting the chromatin effects of sequence alterations with single nucleotide sensitivity. The DeepSEA authors proposed a format of the training data which become a standard for training similar algorithms. They [published the training dataset](http://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz) used to train their model. 

This project allows to rebuild the published dataset from the raw data as well as build custom datasets that can be used with DeepSea and related algorithms.

(The project grew from my curiosity as to why I <ins>can</ins> reproduce the high ROC AUC while training the DeepSEA model with the dataset provided by the authors but I <ins>cannot</ins> match it when using my own custom dataset built according to the specification from the paper.)

## The format of the training data

The training data format is described in the [Predicting effects of noncoding variants with deep learning–based sequence model](https://www.nature.com/articles/nmeth.3547) paper as follows:

> **Data for training DeepSEA**. Training labels were computed from uniformly processed ENCODE and Roadmap Epigenomics data releases. The full list of all chromatin profile files we used are provided in [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM645_ESM.xlsx).

> To prepare the input for the deep convolutional network model, we split the genome into 200-bp bins. For each bin we computed the label for all 919 chromatin features; a chromatin feature was labeled 1 if more than half of the 200-bp bin is in the peak region and 0 otherwise.

> We focused on the set of 200-bp bins with at least one TF bind- ing event, resulting in 521,636,200 bp of sequences (17% of whole genome), which was used for training and evaluating chromatin feature prediction performance. (Variant analyses were not restricted to this region.)

> Each training sample consists of a 1,000-bp sequence from the human GRCh37 reference genome centered on each 200-bp bin and is paired with a label vector for 919 chromatin features. The 1,000-bp DNA sequence is represented by a 1,000 × 4 binary matrix, with columns corresponding to A, G, C and T. The 400-bp flanking regions at the two sides provide extra contextual information to the model.

> Training and testing sets were split by chromosomes and strictly nonoverlapping. Chromosome 8 and 9 were excluded from train- ing to test chromatin feature prediction performances, and the rest of the autosomes were used for training and validation. 4,000 samples on chromosome 7 spanning the genomic coordinates 30,508,751–35,296,850 were used as the validation set. All hyper- parameters were selected on the basis of log likelihood of the validation set data. The validation set data was not used for training or testing.

The [genome coordinates](http://deepsea.princeton.edu/media/code/allTFs.pos.bed.tar.gz) of training data and the [Supplementary Table 2. DeepSEA prediction performance for each transcription factor, DNase I hypersensitive site, and histone mark profile.](https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM646_ESM.xlsx) (converted to txt) were downloaded to the [`data`](https://github.com/jakublipinski/build-deepsea-training-dataset/tree/master/data) folder as [`allTFs.pos.bed`](https://github.com/jakublipinski/build-deepsea-training-dataset/tree/master/data/allTFs.pos.bed) and [`deapsea_metadata.tsv`](https://github.com/jakublipinski/build-deepsea-training-dataset/tree/master/data/deapsea_metadata.tsv) respectively. [The Supplementary Table 1.
List of all publicly available chromatin feature profile files used for training DeepSEA](https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM645_ESM.xlsx) has been downloaded, updated (some files were unavailable under their original location) and saved as [`deepsea_data.urls`](https://github.com/jakublipinski/build-deepsea-training-dataset/tree/master/data/deapsea_data.urls).

## Downloading the original data

In order to build the original dataset you need to download the narrow peak bed files. Clone the repository, enter the `data` folder, download and uncompress the files:
```
git clone git@github.com:jakublipinski/build-deepsea-training-dataset.git
cd build-deepsea-training-dataset/data
xargs -L 1 curl -C - -O -L < deepsea_data.urls
find ./ -name \*.gz -exec gunzip {} \;
cd ..
```

## Building the dataset

After you downloaded the original data run the following script. Modify `--hg` to point to the downloaded and uncompressed [hg19 reference genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz):

```
mkdir out

python build.py \
--metadata_file data/deepsea_metadata.tsv \
--pos data/allTFs.pos.bed \
--beds_folder data/ \
--hg19 ~/BioData/hg19.fa \
--train_size 2200000 \
--valid_size 4000 \
--train_filename out/train.mat \
--valid_filename out/valid.mat \
--test_filename out/test.mat \
--train_data_filename out/train_data.npy \
--train_labels_filename out/train_labels.npy \
--valid_data_filename out/valid_data.npy \
--valid_labels_filename out/valid_labels.npy \
--test_data_filename out/test_data.npy \
--test_labels_filename out/test_labels.npy
```

The script will generate 9 files in the `out` folder.

`(train|valid|test).mat` files contain both data (DNA sequences) and labels for training, validation and testing respectively. The format of the files is same as [the ones provided](http://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz) by the paper authors.

`(train|valid|test)_(data|labels).npy` files contain data and labels for training, validation and testing saved in the `.npy` format.

## Comparing generated dataset with the original one

We can now compare the generated datasets with the ones provided by the authors. Download and uncompress the original dataset to the `data` folder first.

```
cd data
wget http://deepsea.princeton.edu/media/code/deepsea_train_bundle.v0.9.tar.gz
tar -xzf deepsea_train_bundle.v0.9.tar.gz
mv deepsea_train/train.mat deepsea_train/valid.mat deepsea_train/test.mat .
rm -f deepsea_train
cd ..
```

Compare the labels for the training set:
```
python compare_labels.py \
    --deepsea_labels data/train.mat \
    --built_labels out/train.mat
```

The results is:
```
...
Total differences#: 4,243,416 Total differences%: 0.10%
```
Comparing labels for the validation and test set shows similar level of difference. See the comments below for the discussion on the differences.

Now compare the one-hot sequence vectors in the training set:
```
python compare_sequence.py \
    --deepsea_data data/train.mat \
    --built_data out/train.mat
```

We get:
```
Total differences#:5,256 Total differences%: 0.0001%
```
We get zero differences for the valid and test sets.

## Training the DeepSEA model with custom datasets on Google Colab

You can train the DeepSEA model using the generated datasets on Google Colab [with Keras](https://github.com/zj-zhang/deepsea-keras). The following notebook allows to train the model on 5 different datasets:

* The original dataset
* The regenerated dataset with the bins from `allTFs.pos.bed` file
* The regenerated dataset with all the bins (strictly following the spec)<sup>1</sup>
* The regenerated dataset with the the bins for which chromatin feature signal enrichment ([`signalValue`](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)) is at least `7`
* The custom dataset based on ENCODE project with all the bins (strictly following the spec)
* The custom dataset based on ENCODE project with the the bins for which chromatin feature signal enrichment ([`signalValue`](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)) is at least `50`

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1Yqmhn_2FkboPXdF_dgx040m7Yjk2WEwz)

<sup>1</sup> You cannot train the model on Google Colab with this dataset due to its large size

We get the following ROC AUC after training the model on these datasets (epochs:3, batch size: 1000):

| Dataset                                   | ROC AUC |
| ----------------------------------------- | ------:|
| The original dataset                      | 0.9176 |
| The regenerated dataset (`allTFs.pod.bed`)| 0.9191 |
| The regenerated dataset (all the bins)    | 0.8548<sup>2</sup>|
| The regenerated dataset (Signal >= 7)     | 0.9026 |
| The custom ENCODE dataset                 | 0.8441 |
| The custom ENCODE dataset (Signal >= 50)  | 0.8796 |

<sup>2</sup> Run only for 1 epoch because of the large dataset size.

## Differences from the original dataset and their reasons

The 0.10% difference between the original dataset and the one generated using the tool must most likely come from some off-by-one difference in  the implementation of the algorithm which assigns chromatin features to particular bins.

The difference in the genome data is negligible.

The interesting difference appears when you regenerate the dataset strictly according to the description from the paper and without limiting the bins to the ones listed in the `allTFs.pos.bed`. In such case this tool identifies 6,575,911 training bins (instead of 4,400,000) and 641,600 testing bins (instead of 455,024):

```
mkdir out2
python build.py \
--metadata_file data/deepsea_metadata.tsv \
--beds_folder data/ \
--hg19 ~/BioData/hg19.fa \
--valid_size 4000 \
--train_filename out2/train.mat \
--valid_filename out2/valid.mat \
--test_filename out2/test.mat \
--train_data_filename out2/train_data.npy \
--train_labels_filename out2/train_labels.npy \
--valid_data_filename out2/valid_data.npy \
--valid_labels_filename out2/valid_labels.npy \
--test_data_filename out2/test_data.npy \
--test_labels_filename out2/test_labels.npy
```

It seems the selection of the bins in the original dataset was based on the signal value of the corresponding chromatin features. If you exclude the features below some specified threshold using the `--signal_threshold` option you get the dataset with a ROC AUC as in the original dataset:

```
mkdir out3
python build.py \
--metadata_file data/deepsea_metadata.tsv \
--beds_folder data/ \
--hg19 ~/BioData/hg19.fa \
--valid_size 4000 \
--train_filename out3/train.mat \
--valid_filename out3/valid.mat \
--test_filename out3/test.mat \
--train_data_filename out3/train_data.npy \
--train_labels_filename out3/train_labels.npy \
--valid_data_filename out3/valid_data.npy \
--valid_labels_filename out3/valid_labels.npy \
--test_data_filename out3/test_data.npy \
--test_labels_filename out3/test_labels.npy \
--signal_threshold 7
```

In such case the dataset contains 4,419,304 bins (original: 4,400,000) and ROC AUC (after 3 epochs) at xxx (original: 0.9176).

(You can compare the average chromatin features signal of the bins included in the `allTFs.pos.bed` with all the bins by checking out the `calculate-signal-value` branch)

## Building custom dataset with data from ENCODE

You can build your own dataset for training DeepSEA using the narrow peaks bed files from the ENCODE project.

Let's assume we want to train the model using results from the narrow peak ChIP-seq experiments targeted to transcription factors in GM12878 cell line assembled with hg19 genome. Go to [the ENCODE browser](https://www.encodeproject.org/search/?type=Experiment&status=released&searchTerm=chip&biosample_ontology.term_name=GM12878&target.investigated_as=transcription+factor&files.file_type=bed+narrowPeak&assembly=hg19) and click the `Download` button. Click `Download` button again and save the `files.txt` to the `encode_data` folder.

Download and uncompress the files:

```
cd encode_data
xargs -L 1 curl -O -L < files.txt
find ./ -name \*.gz -exec gunzip {} \;
cd ..
```

Rename file containing metadata to something more readable:

```
mv \?type=Experiment\&status=released\&biosample_ontology.term_name=GM12878\&target.investigated_as=transcription+factor\&files.file_type=bed+narrowPeak\&assembly=hg19\&searchTerm=chip metadata.tsv
```

Create datasets:
```
python build.py \
--metadata_file encode_data/metadata.tsv \
--filter "Output type=conservative IDR thresholded peaks" \
--beds_folder encode_data/ \
--hg19 ~/BioData/hg19.fa \
--valid_ratio .0.5 \
--train_filename encode_out/train.mat \
--valid_filename encode_out/valid.mat \
--test_filename encode_out/test.mat \
--train_data_filename encode_out/train_data.npy \
--train_labels_filename encode_out/train_labels.npy \
--valid_data_filename encode_out/valid_data.npy \
--valid_labels_filename encode_out/valid_labels.npy \
--test_data_filename encode_out/test_data.npy \
--test_labels_filename encode_out/test_labels.npy
```

Please note that you can filter the metadata by using the `--filter` option.

Consider excluding bins with low signal level using the `--signal_threshold` parameter:
```
python build.py \
--metadata_file encode_data/metadata.tsv \
--filter "Output type=conservative IDR thresholded peaks" \
--beds_folder encode_data/ \
--hg19 ~/BioData/hg19.fa \
--valid_ratio .0.5 \
--train_filename encode_out/train.mat \
--valid_filename encode_out/valid.mat \
--test_filename encode_out/test.mat \
--train_data_filename encode_out/train_data.npy \
--train_labels_filename encode_out/train_labels.npy \
--valid_data_filename encode_out/valid_data.npy \
--valid_labels_filename encode_out/valid_labels.npy \
--test_data_filename encode_out/test_data.npy \
--test_labels_filename encode_out/test_labels.npy
--signal_threshold 50
```

## Debug info

You can add `--save_debug_info True` command line argument to save all the data into human readable `.tsv` files: `debug_train.tsv`, `debug_valid.tsv` and `debug_test.tsv` to debug and verify the results.

## Further work

Currently only the hg19 genome assembly is supported by it should be reasonably easy to support others. I welcome your PRs.

The whole tool can be significantly sped up by using the Python `multiprocessing` module.