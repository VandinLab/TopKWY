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


1. REPRODUCIBILITY OF THE EXPERIMENTS DESCRIBED IN THE ARTICLE

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

More details about TopKWY are provided in */src/README.txt*.
