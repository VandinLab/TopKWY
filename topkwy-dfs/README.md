# TopKWY: Efficient Mining of the Most Significant Patterns with Permutation Testing
### Authors: Leonardo Pellegrina (pellegri@dei.unipd.it), Fabio Vandin (fabio.vandin@unipd.it)

This folder contains the C++ source files for the DFS version of TopKWY (TopKWY-dfs) for mining significant itemsets. Note that TopKWY-dfs is generally (much) slower than TopKWY. TopKWY-dfs is based on the frequent itemset mining algorithm LCM (see LCM_readme.txt).


## REPRODUCIBILITY OF THE EXPERIMENTS DESCRIBED IN THE ARTICLE

To reproduce the experiments described in our paper, first compile TopKWY-dfs using the `make` command from this folder. Then, the script `run_all_datasets.py` can be used to execute all experiments.


## USAGE OF TOPKWY-DFS

TopKWY-dfs can manually be launched with the command `./topkwy-dfs [list of parameters]`


### Parameters:

The parameters to give are:
`./topkwy-dfs output_basefilename n_perm target_fwer input_class_labels_file input_transactions_file seed [k]`

`output_basefilename` is the prefix of the path for the output files
`n_perm` the number of permutations to use for the Westfall-Young multiple hypothesis testing method. Default is 10e4.
`target_fwer` is the target upper bound to the FWER of the retrieved significant patterns. Default is 0.05.
`input_class_labels_file` is the input file containing the labels of the permutations (each line should contain either 0 or 1).
`input_transactions_file` is the input file containing the transactions (see below for the input format for the transactions).
`seed` is the random seed for generating the permuted labels.
`k` is the [optional] desired number of significant itemsets.

Example:
```
./topkwy mushroom 10000 0.05 mushroom.labels mushroom.trans 20202021 10
```


### Dataset File

Each line in the dataset file represents one transaction, where every item is represented by its id, and it is space separated from other items.
As an example, the transaction `1 100 1500` contains the three items "1", "100", and "1500".

One transaction does not need to contains items in sorted order, but should not contains duplicate items.
