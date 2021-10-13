# TopKWY: Efficient Mining of the Most Significant Patterns with Permutation Testing
### Authors: Leonardo Pellegrina (pellegri@dei.unipd.it), Fabio Vandin (fabio.vandin@unipd.it)

This package contains the source files for TopKWY, the first algorithm to mine the *k* most significant patterns with bounded Family-Wise Error Rate (*FWER*) using WY permutation testing.
TopKWY is currently implemented to mine significant itemsets and sugraphs.
The repository also contains all the datasets and the scripts to reproduce all experiments described in our papers "Efficient Mining of the Most Significant Patterns with Permutation Testing" (KDD 2018 https://doi.org/10.1145/3219819.3219997 and DAMI 2020 https://doi.org/10.1007/s10618-020-00687-8).

If you find any bugs with this package, please do not hesitate to contact us at pellegri@dei.unipd.it and fabio.vandin@unipd.it.


## REPOSITORY DESCRIPTION

The package contains four folders:

1. `/TopKWY/` : contains the code for the BFS version of TopKWY to mine significant itemsets
2. `/topkwy-dfs/` : contains the code for the DFS version of TopKWY to mine significant itemsets
3. `/topkwy-subgraph/` : contains the code for TopKWY to mine significant subgraphs
4. `/wylight/` : contains the code of WYlight algorithm (for itemset mining, from https://github.com/fllinares/wylight ; we made minimal changes to the Makefile for compilation)

More details on each algorithm can be found in the README files within each subfolder.
