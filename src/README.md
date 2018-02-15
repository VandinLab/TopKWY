# TopKWY: Efficient Mining of the most Significant Patterns
### Authors: Leonardo Pellegrina (pellegri@dei.unipd.it), Fabio Vandin (fabio.vandin@unipd.it)

This package contains the C++ source files for TopKWY, tested under Linux, the datasets and the scripts to reproduce experiments described in "Efficient Mining of the most Significant Patterns".


If you find any bugs with this package, please do not hesitate to contact us at pellegri@dei.unipd.it and fabio.vandin@unipd.it.


## PACKAGE DESCRIPTION

The package contains three folders:

1. /src/ : contains the source code of TopKWY.
2. /datasets/ : contains the datasets used for our experiments.
   - /datasets/dat_name/ : contains the files for the dataset "dat_name": the transactions (dat_name.dat), the binary class labels (dat_name.labels) and the specification file (dat_name_new.spec).
3. /scripts/ : contains the scripts which replicates our experiments.


## REPRODUCIBILITY OF THE EXPERIMENTS DESCRIBED IN THE ARTICLE

To reproduce the experiments described in "Efficient Mining of the most Significant Patterns", follows these steps:
1. TopKWY compilation: use the `make` command from inside the /src/ folder.
2. Files organisation: copy the *topkwy* executable and *runexperiments_all.py* in the main folder.
3. Run the script for experiments with `python runexperiments_all.py`. This script launches TopKWY with jp=10^4, alpha=0.05, k = [10,100,1000,100000,1000000] for all 19 datasets, performing for each of those combination 10 runs. These parameters can be changed inside the *runexperiments_all.py* file.
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


## USAGE OF TOPKWY

TopKWY can manually be launched with the command `./topkwy file.spec max_ram k jp alpha`


### Parameters:

1. (file.spec): the name of the specification file for the input (see for example the .spec files included in /datasets/ subfolders)
2. (max_ram): the max size (in MB) of the ram allowed for the patterns' exploration (without considering the space needed for the Patricia Trie). (100 GB is the default value) This parameter allows the user to bound the memory usage of TopKWY when exploring particularly challenging datasets.
3. (k): the maximum number of significant itemsets required.
5. (alpha): the upper bound to the FWER of the retrieved significant patterns.
4. (jp): the number of permutations to use for the Westfall-Young multiple hypothesis testing method.

Example:
```
./topkwy mushroom_new.spec 100000 10 10000 0.05
```

Among which, the specification file for the corresponding dataset contains the following information:
1. file name of the dataset transactions (the path must be relative to TopKWY executable);
2. total number of distinct items that are in the dataset;
3. maximum length of one transaction (i.e., maximum number of items in each transaction);
4. total number of transactions;
5. file name of the dataset labels (the path must be relative to TopKWY executable);
6. a default value for jp (used if it is not specified in the arguments);
7. a default value for alpha (used if it is not specified in the arguments).

Example for mushroom dataset:
```
mushroom.dat
118
22
8124
mushroom.labels
0.05
10000
```

### Dataset File

Each line in the dataset file represents one transaction, where every item is represented by its id, and it is space separated from other items.
As an example, the transaction `1 100 1500` contains the three items "1", "100", "1500".


One transaction should not contains duplicate items.


#### DATASETS CONVERSIONS

The script *convert_labeled_datasets.py* converts an unlabelled datasets into a labeled one, by removing the item with frequency closest to 0.5 and using it as target. It also print the correct .spec file to be used with TopKWY.
The script *compute_spec_file.py* computes the .spec file of any already labelled dataset *dat_name* by using the command `python compute_spec_file.py dat_name`


### Output File

The significant itemsets of the "dat_name" dataset can be found in the file *dat_name_k_alpha_jp.txt* inside its folder. The itemsets are ordered by decreasing p-value, and each line is in the format:
```
items : support : as : p-value : log-p-value
```
For example the line
```
8 6 : 30 : 29 : 8.13406e-07 : -6.08969
```
indicates that the itemset with items 8, 6 appears in 30 transactions, 29 of which belongs to the minority class. The p-value p of the distribution on the classes of the itemset is 8.13406e-07. The log of this value is also provided, which is useful in cases where p = 0.
