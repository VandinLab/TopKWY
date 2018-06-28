# TopKWY: Efficient Mining of the Most Significant Patterns with Permutation Testing
### Authors: Leonardo Pellegrina (pellegri@dei.unipd.it), Fabio Vandin (fabio.vandin@unipd.it)

This package contains the C++ source files for TopKWY, tested under Linux, the datasets and the scripts to reproduce experiments described in "Efficient Mining of the Most Significant Patterns with Permutation Testing".


If you find any bugs with this package, please do not hesitate to contact us at pellegri@dei.unipd.it and fabio.vandin@unipd.it.


## PACKAGE DESCRIPTION

The package contains three folders:

1. /src/ : contains the source code of TopKWY.
2. /datasets/ : contains the (zipped) datasets used for our experiments. When extracted, the folder /datasets/dat_name/ contains the files for the dataset "dat_name": the transactions (dat_name.dat), the binary class labels (dat_name.labels) and the specification file (dat_name_new.spec) to be given in input to TopKWY.
3. /scripts/ : contains the scripts which replicates our experiments.


## REPRODUCIBILITY OF THE EXPERIMENTS DESCRIBED IN THE ARTICLE

To reproduce the experiments described in "Efficient Mining of the Most Significant Patterns with Permutation Testing", you can follow these steps:
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



More details about TopKWY are provided in */src/README.md*.
