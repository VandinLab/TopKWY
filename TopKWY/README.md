# TopKWY: Efficient Mining of the Most Significant Patterns with Permutation Testing
### Authors: Leonardo Pellegrina (pellegri@dei.unipd.it), Fabio Vandin (fabio.vandin@unipd.it)

This folder contains the C++ source files for the BFS version of TopKWY for mining significant itemsets.
TopKWY is based on the top-k frequent closed itemset mining algorithm TopKMiner by A. Pietracaprina and F. Vandin. 

## PACKAGE DESCRIPTION

The package contains three folders:

1. `/src/` : contains the source code of TopKWY.
2. `/datasets/` : contains the datasets used for our experiments. Each folder /datasets/dat_name/ contains the files for the dataset "dat_name": the transactions (dat_name.dat), the binary class labels (dat_name.labels) and the specification file (dat_name_new.spec) to be given in input to TopKWY.
3. `/scripts/` : contains the scripts which replicates our experiments.


## REPRODUCIBILITY OF THE EXPERIMENTS DESCRIBED IN THE ARTICLE

To reproduce the experiments described in our paper, you can follow these steps:
1. TopKWY compilation: use the `make` command from inside the /src/ folder.
2. Files organisation: copy the *topkwy* executable and *runexperiments_all.py* in the main folder.
3. Run the script for experiments with `python runexperiments_all.py`. This script launches TopKWY with jp=10^4, alpha=0.05, k = [10, 10e2, 10e3, 10e4, 10e5, 10e6] for all 19 datasets, performing for each of those combination 10 runs. These parameters can be changed inside the *runexperiments_all.py* file.
4. Results: the script creates the file *all_results.csv* which contains average and variance of all the runs for every dataset and every value of k.


The file contains rows in this format:

```
k; jp; alpha; dataset_name; running_time; variance_of_running_time; peak_memory; variance_of_peak_memory;
```

All the values of all the single runs are stored in the files *final_k_jp.txt* with various values of k.

In these files the format is:
```
k; jp; alpha; dataset_name; running_time; peak_memory; number_of_tested_patterns; run_id;
```

### Output File

The significant itemsets of the "dat_name" dataset can be found in the file *dat_name_k_alpha_jp.txt* inside its folder. The itemsets are ordered by decreasing p-value, and each line is in the format:
```
items : support : as : p-value : log-p-value
```
For example the line
```
8 6 : 30 : 29 : 8.13406e-07 : -6.08969
```
indicates that the itemset with items 8, 6 appears in 30 transactions, 29 of which belongs to the minority class. The p-value p of the distribution on the classes of the itemset is 8.13406e-07. The log (in base 10) of this value is also provided, which is useful in cases where p is reported as = 0 due to double underflow.


## USAGE OF TOPKWY

TopKWY can manually be launched with the command `./topkwy [list of parameters]`


### Parameters:

1. `-s file.spec`: the name of the specification file for the input (see for example the .spec files included in /datasets/ subfolders)
2. `-k number_of_significant_results`: [optional] the desired number of significant itemsets. Default is 10e6.
3. `-j number_of_permutations`: [optional] the number of permutations to use for the Westfall-Young multiple hypothesis testing method. Default is 10e4.
4. `-a alpha`: [optional] the target upper bound to the FWER of the retrieved significant patterns. Default is 0.05.
5. `-g gfwer`: [optional] controls the *g*-FWER instead of the FWER, for any integer value of *gfwer* >= 1. The *g*-FWER is defined as the probability of reporting at least *g* false positives in the output (the FWER corresponds to *g*=1). Default is 1.
6. `-r max_ram`: [optional] the max size (in GB) of the ram allowed for the patterns' exploration (without considering the space needed for the Patricia Trie). (100 GB is the default value). This parameter allows the user to bound the memory usage of TopKWY when exploring particularly challenging datasets.

Example:
```
./topkwy -s mushroom_new.spec -k 10 -j 10000 -a 0.05
```

The specification file for the corresponding dataset contains the following information:
1. file name of the dataset transactions (the path must be relative to TopKWY executable);
2. total number of distinct items that are in the dataset;
3. maximum length of one transaction (i.e., maximum number of items in each transaction);
4. total number of transactions;
5. file name of the dataset labels (the path must be relative to TopKWY executable);

Example for mushroom dataset:
```
mushroom.dat
118
22
8124
mushroom.labels
```
See below on how to generate such a file automatically for any dataset.

### Dataset File

Each line in the dataset file represents one transaction, where every item is represented by its id, and it is space separated from other items.
As an example, the transaction `1 100 1500` contains the three items "1", "100", and "1500".

One transaction does not need to contains items in sorted order, but should not contains duplicate items.


#### DATASETS PREPARATION

The script *compute_spec_file.py* computes the .spec file of any already labelled dataset *dat_name* by using the command `python compute_spec_file.py dat_name`
In case a given dataset is unlabelled, the script *convert_labeled_datasets.py* converts it into a labeled one, by removing the item with frequency closest to 0.5 and using it as target. It also prints the correct .spec file to be used with TopKWY.


### Output File

The significant itemsets of the "dat_name" dataset can be found in the file *dat_name_k_alpha_jp.txt* inside its folder. The itemsets are ordered by decreasing p-value, and each line is in the format:
```
items : support : as : p-value : log-p-value
```
For example the line
```
8 6 : 30 : 29 : 8.13406e-07 : -6.08969
```
indicates that the itemset with items 8, 6 appears in 30 transactions, 29 of which belongs to the minority class. The p-value p of the distribution on the classes of the itemset is 8.13406e-07. The log (in base 10) of this value is also provided, which is useful in cases where p is reported as = 0 due to double underflow.
