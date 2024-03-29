# TopKWY: Efficient Mining of the Most Significant Patterns with Permutation Testing
### Authors: Leonardo Pellegrina (pellegri@dei.unipd.it), Fabio Vandin (fabio.vandin@unipd.it)

This folder contains the C++ source files for the DFS version of TopKWY for mining significant subgraphs (TopKWY-dfs).
[GASTON](http://www.liacs.nl/~snijssen/gaston/iccs.html) is used as a frequent subgraph mining algorithm:
* Nijssen, S. and Kok, J.: **A Quickstart in Frequent Structure Mining Can
  Make a Difference,** *Proceedings of the 10th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (KDD2004), 647-652, 2004.

We based our implentation on the work of Mahito Sugiyama (https://github.com/BorgwardtLab/significant-subgraph-mining) and Felipe Llinares-López (https://github.com/BorgwardtLab/Significant-Pattern-Mining-WYlight), supporting the papers
* Sugiyama, M., Llinares-López, F., Kasenburg, N., Borgwardt, K. M.: **Significant Subgraph Mining with Multiple Testing Correction,** *Proceedings of the 2015 SIAM International Conference on Data Mining* (SDM2015), 199-207, 2015.
[[PDF]](http://epubs.siam.org/doi/pdf/10.1137/1.9781611974010.5)
* Llinares-López, F., Sugiyama, M., Papaxanthos, L., Borgwardt, K. M.:
**Fast and Memory-Efficient Significant Pattern Mining via Permutation Testing,** *Proceedings of the 21st ACM SIGKDD Conference on Knowledge Discovery and Data Mining* (KDD2015), 725-734, 2015.
[[PDF]](http://dl.acm.org/ft_gateway.cfm?id=2783363)

## REPRODUCIBILITY OF THE EXPERIMENTS DESCRIBED IN THE ARTICLE

To reproduce the experiments described in our paper, first compile TopKWY-dfs (see below on detailed instructions). Then, the script `run_all_datasets.py` can be used to execute all experiments.


## USAGE OF TOPKWY-DFS

### Data format

The input is composed of two files: a graph database (a list of graphs) and a list of class labels.
In a graph database, each graph is described as follows:

```
# <comment>
t # <transaction id>
v <node id> <label id>
e <source node id> <target node id> <label id>
```

This is the same as the standard format in frequent subgraph mining softwares, such as [gSpan](https://www.cs.ucsb.edu/~xyan/software/gSpan.htm) and [GASTON](http://www.liacs.nl/~snijssen/gaston/iccs.html).
For example:

```
# start  
t # 1        // transaction id is 1
v 0 1        // node 0 with label 1
v 1 2        // node 1 with label 2
e 0 1 0      // edge from node 0 to 1 with label 0
```

In the class file, each line contains the class label of the corresponding graph.
The class label should be either "0" or "1".

### Compilation

To compile the program, go to the "src" directory and type

```
make
```
The "Boost" library is needed to compile it.

### Run

To run the algorithm, type

```
./sgmine -a <alpha> -i <input_file> -c <input_class_file> -o <output_file> -k <max_results>
```

`<alpha>` is the target FWER. It is set to be 0.05 if skipped.

List of arguments:

`-w` : Perform Westfall-Young permutation correction if specified  
`-m <maxsize>` : Maximum size of each subgraph (default: unlimited)  
`-a <alpha>` : Target FWER (default: 0.05)  
`-j <perm>` : Number of permutations in `-w` mode (default: 1000)  
`-r <seed>` : Seed for permutations in `-w` mode (default: 0)  
`-i <input_file>` : Input file of graph database  
`-c <input_class_file>` : Input file of class labels  
`-o <output_file>` : Output file of significant subgraphs  
`-k <max_results>` : Maximum number of results in output  

### Output format

Resulting significant subgraphs are written to the file `<output_file>`. The output file has the same format as the input file except for two lines: supports (start from "\# s") in minor and major classes and the *p*-value (start from "\# p") of the graph. For example:

```
# s 22 7       // support of this subgraph in minor (left) and major (right) classes
# p 1.5841e-07 // p-value of this subgraph
t 1
v 0 0
v 1 9
v 2 9
e 0 1 0
e 1 2 0
```
