#ifndef _wy_c_
#define _wy_c_

/* LIBRARY INCLUDES */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<fstream>
//#include <set>
#include <queue>

/* CODE DEPENDENCIES */
#include"var_declare.h"
#include"permutation.h"
// #include"transaction_keeping.c"
// #include"lcm_var.c"
#include"wy.h"
#include"misc.h"

//#define bounds_debug 1

/* CONSTANT DEFINES */


/* GLOBAL VARIABLES */
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;
// Target FWER
double alpha;
// Current FWER
double FWER;

// Minimum P-value for each permutation
double *min_pval;
// Region thresholds: Sigma_k = [sl1,sl2] U [N-sl2,N-sl1]
int sl1, sl2;
// Current P-value threshold
double delta;
// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
int flag;

// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute
double log_inv_binom_N_n;
// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi;
// Array for storing values of the CDF of the hypergeometric distribution for fast p-valaue computation
double *hypergeom_cdf;
double *hypergeom_cdf_inverse;

// Cell-count counter
int *a_cnt;

int last_precomputation;
int last_processed_support;
int last_processed_as;
int last_processed_max_as;
int last_processed_min_as;
int *last_processed_transactions;

// Profiling variables
long long n_pvalues_computed;
long long n_cellcounts_computed;
long long effective_total_dataset_frq;
long long tested;
long long explored;
long long explored_testable;
long long explored_closed;
//std::set<double> topkqueue;
std::priority_queue<double> topkqueue;
double last_queue_element;
int K;

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void wy_init(double target_fwer){
	int j; //Loop variable

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	alpha = target_fwer;
	// And initialise some others
	sl1 = 1; sl2 = N_over_2;
	flag = 1;
	FWER = 0;
	delta = ((double) n)/N; //$\psi(1)=\frac{n}{N}$

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	// Allocate memory for minimum p-values, raising error if it fails
	min_pval = (double *)malloc(J*sizeof(double));
	if(!min_pval){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array min_pval\n");
		exit(1);
	}
	// Initialise all p-values to 1
	for(j=0; j<J; j++) min_pval[j] = 1;

	// Allocate memory for the CDF of the hypergeometric distribution
	// (worst case memory requirement n+1), raising an error if it fails
	hypergeom_cdf = (double *)malloc((n+1)*sizeof(double));
	if(!hypergeom_cdf){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_cdf\n");
		exit(1);
	}
	hypergeom_cdf_inverse = (double *)malloc((n+1)*sizeof(double));
	if(!hypergeom_cdf_inverse){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_cdf_inverse\n");
		exit(1);
	}

	// Allocate memory for cell counts, raising an error if it fails
	a_cnt = (int *)malloc(J*sizeof(int));
	if(!a_cnt){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array a_cnt\n");
		exit(1);
	}
	for(j=0; j<J; j++) a_cnt[j] = 0;

	last_processed_transactions = (int *)malloc(N*sizeof(int));

	n_pvalues_computed = 0; n_cellcounts_computed = 0; effective_total_dataset_frq = 0; //Init profiling variables
	last_precomputation = 0; tested = 0; explored = 0; explored_testable = 0; last_queue_element = 1.0; explored_closed = 0;
	delta = 0.05; while (psi[minfreq] > delta) ++minfreq;
	cout << "<<< WY: Update minfreq to " << minfreq << " >>>" << endl;
}

/* Precompute values of log(x!) storing them in the array loggamma */
void loggamma_init(){
	int x;
	// Allocate memory for log-gamma cache, raising error if it fails
	loggamma = (double *)malloc((N+1)*sizeof(double));
	if(!loggamma){
		fprintf(stderr,"Error in function loggamma_init: couldn't allocate memory for array loggamma\n");
		exit(1);
	}
	// Initialise cache with appropriate values
	for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
	// Initialise log_inv_binom_N_n
	log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}

/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double xi1;
	int x, x_init;
	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
	if(!psi){
		fprintf(stderr,"Error in function psi_and_xi1_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
	//In this range, the recursion $\psi(x)$=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can be seen to be correct
	for(x=1; x<=n; x++) psi[x] = (((double)(n-(x-1)))/(N-(x-1)))*psi[x-1];

	// Now, start computing xi1 in the range [N-N_over_2,N] using another recursion, this time
	// starting in N
	// Note that we don't need to store all values, since this will be used only to initialise
	// psi[N_over_2]
	x_init = N-N_over_2;
	xi1 = 1;
	//In this range, the recursion $\xi_{1}(x)$=$\xi_{1}(x+1)$*[((x-1)-n)/(x+1)] can be seen to be correct
	for(x=(N-1); x>=x_init; x--) xi1 = (((double)((x+1)-n))/(x+1))*xi1;

	// Now, use the value of $\xi_{1}(N-N_over_2)$=xi1[0] to get $\psi(N_over_2)$=psi[N_over_2] using the
	// same recursion if N is odd, or simply copy the value of xi1[0] since $\xi_{1}(N-N_over_2)=\psi(N_over_2)$
	// if N is even
	if (N % 2) psi[N_over_2] = (((double)(x_init-n))/x_init)*xi1;
	else psi[N_over_2] = xi1;

	// Now, using psi[N_over_2] compute the right side of "the W", i.e. the range [n+1,N_over_2]
	// using the same recursion as for $\xi_{1}$
	for(x=(N_over_2-1); x > n; x--) psi[x] = (((double)((x+1)-n))/(x+1))*psi[x+1];

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=x_init; x<=N; x++) psi[x] = psi[N-x];
}

void compute_corrected_delta(){
	int j, idx_max;
	double delta_corrected;
	// Sort p-values
	qsort(min_pval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; delta_corrected = min_pval[idx_max];
	// Check and correct (if necessary) boundary cases
	if(delta_corrected==min_pval[idx_max+1]){
		while(min_pval[--idx_max]==delta_corrected);
		delta_corrected = min_pval[idx_max];
	}
	//DELTA = delta_corrected;
	//for(int i=0; i <20; i++) cout << "psi["<<i<<"] " << psi[i] << endl;
	/*cout << delta << endl;
	cout << delta_corrected << endl;*/
	DELTA = min(delta_corrected , last_queue_element);

}

/* Free all allocated memory and give some output for debugging purposes */
void wy_end(){
	/*int j, idx_max;
	double delta_corrected;
	// Sort p-values
	qsort(min_pval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; delta_corrected = min_pval[idx_max];
	// Check and correct (if necessary) boundary cases
	if(delta_corrected==min_pval[idx_max+1]){
		while(min_pval[--idx_max]==delta_corrected);
		delta_corrected = min_pval[idx_max];
	}
	//DELTA = delta_corrected;
	for(int i=0; i <20; i++) cout << "psi["<<i<<"] " << psi[i] << endl;
	//cout << delta << endl;
	//cout << delta_corrected << endl;
	DELTA = min(delta_corrected , last_queue_element);*/

	cout << "explored " << explored << endl;
	cout << "explored_testable " << explored_testable << endl;
	cout << "tested " << tested << endl;
	cout << "explored closed patterns " << explored_closed << endl;

	#ifdef debug_queue
	/*for(std::set<double>::iterator it = topkqueue.begin(); it != topkqueue.end(); ++it)
	{
	    cout << *it << endl;
	}*/

	cout << "content of thq queue (in reversed order)" << endl;
	std::priority_queue<double> temp = topkqueue;
  while (!temp.empty()) {
      cout << temp.top() << endl;
			temp.pop();
  }
	#endif

	// Print results
	/*
	cerr << endl << "Results:" << endl
	     << TAB << "Corrected significance threshold:         " << delta_corrected << endl
	     << TAB << "FWER at corrected significance threshold: " << floor(idx_max+1)/J << endl
	     << TAB << "Final minfreq:                            " << minfreq << endl
	     << TAB << "Final P-value lower bound:                " << delta << endl
	     << TAB << "FWER at final P-value lower bound:        " << FWER << endl;
	*/

	/*
	for(j=0;j<(J-1);j++) ofs_p << min_pval[j] << " ";
	ofs_p << endl;
	*/

	// Free allocated memory
	free(loggamma);
	free(psi);
	free(min_pval);
	free(hypergeom_cdf);
	free(hypergeom_cdf_inverse);
	free(a_cnt);
	free(last_processed_transactions);
	cout << "wy_end done " << endl;
}

/* -------------------------------FUNCTIONS TO COMPUTE FISHER EXACT TEST P-VALUES----------------------------------- */

/* This function precomputes all Fisher exact test P-values for a contingency table with margins x,n,N that is,
 * all p-values p(a,x,n,N) for a in the range [max(0,n+x-N),min(x,n)]. The results will be stored in the array
 * hypergeom_cdf such that p(a,x,n,N)=hypergeom_cdf[a]. Note that values hypergeom_cdf[a] for a outside
 * [max(0,n+x-N),min(x,n)] are undefined and could contain garbage of previous hypotheses.
 * */
void precompute_pvals(int x){
	if (last_precomputation == x){
		return;
	}
	last_precomputation = x;
	double pre_comp_xterms;
	int a, a_min, a_max;
	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)
	hypergeom_cdf[a_min] = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]));
	hypergeom_cdf_inverse[a_max] = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]));
	hypergeom_cdf[a_max] = 1; hypergeom_cdf_inverse[a_min] = 1;
	for(a=(a_min+1);a<a_max;a++) hypergeom_cdf[a] = hypergeom_cdf[a-1] + exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));
	for(a=(a_max-1);a>a_min;a--) hypergeom_cdf_inverse[a] = hypergeom_cdf_inverse[a+1] + exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));
	for(a=a_max;a>a_min;a--) hypergeom_cdf[a] = (hypergeom_cdf[a] < hypergeom_cdf_inverse[a]) ? hypergeom_cdf[a] : hypergeom_cdf_inverse[a];
}

double getPValue(int index){

	return hypergeom_cdf[index];

}

/* --------------------------------CORE FAST WESTFALL-YOUNG PERMUTATION FUNCTIONS------------------------------------ */

/* Decrease the minimum p-value threshold one level
 * Main operations that need to be performed are:
 * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
 *    is Sigma_{k} = [sl1,sl2] U [N-sl2,N-sl1], we need to figure out if Sigma_{k+1} is of the form
 *    Sigma_{k+1} = [sl1+1,sl2] U [N-sl2,N-sl1-1] (shrink left side) or Sigma_{k+1} = [sl1,sl2-1] U [N-sl2+1,N-sl1-1]
 *    (shrink right side). This is done with help of a binary flag that remembers which of the two types of region
 *    change happened the last time the threshold was decreased.
 * 2) Update variables sl1, sl2 and delta accordingly
 * 3) If sl1 has been modified, then the support of LCM has to be modified
 * 4) Since the temptative corrected significance threshold delta has changed, the FWER needs to be recomputed
 * */
void decrease_threshold(){
	int j; //Loop iterator
	int false_positives; //Number of false positives (a false positive occurs if min_pval[j] <= delta)
	// Flag==1 means the last call to decrease_threshold() shrunk "the W" on the left side
	if(flag){
		sl1++; // Shrink Sigma_k on extremes of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]) delta = min(psi[sl1] , delta);
		else{ delta = min(psi[sl2] , delta); flag = 0; }
		//Update LCM minimum support
		// LCM_th = sl1;
		minfreq = max(sl1 , (int)minfreq);
		delta = min(psi[minfreq] , delta);
		// printf("\n\n\nTHRESHOLD CHANGE!!! NEW THRESHOLD=%d\n\n\n",LCM_th);
		// cerr << minfreq << " ";
		cout << "<<< WY: Update minfreq to " << minfreq << " >>>" << endl;
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		sl2--; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]){ delta = min(psi[sl1] , delta); flag = 1; }
		else delta = min(psi[sl2] , delta);
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
	// Recompute FWER from scratch
	false_positives = 0;
	for(j=0; j<J; j++) false_positives += (min_pval[j]<=delta) ? 1 : 0;
	FWER = ((double)false_positives)/J;
}

void setK(int k_){
	K = k_;
	//cout << K << endl;
}


void addTopKElement(double new_pvalue){

	//topkqueue.insert(new_pvalue);
	topkqueue.push(new_pvalue);

	#ifdef topk_debug
	cout << "call to topk strategy " << endl;
	#endif

	if (topkqueue.size() > K){

		#ifdef topk_debug
		cout << "______________________________________ " << endl;
		cout << "increasing support due to topk strategy " << endl;
		#endif

		// remove highest pvalue
		//topkqueue.erase(--topkqueue.end());
		topkqueue.pop();
		// take the last still in the queue
		//last_queue_element = *(--topkqueue.end());
		last_queue_element = topkqueue.top();

		#ifdef topk_debug
		//cerr << "first_element " << *(topkqueue.begin()) << endl;
		cout << "last_queue_element " << last_queue_element << endl;
		#endif

		// increase the minimum support accordingly
		while (psi[minfreq] > last_queue_element){
			++minfreq;

			#ifdef topk_debug
			cout << "minfreq " << minfreq << endl;
			cout << "psi[minfreq] " << psi[minfreq] << endl;
			#endif
			cout << "<<< Top-k: Update minfreq to " << minfreq << " >>>" << endl;

		}

		delta = min(delta , last_queue_element);

		// update FWER estimate
		int false_positives = 0;
		for(int i=0; i<J; i++) false_positives += (min_pval[i]<=delta) ? 1 : 0;
		FWER = ((double)false_positives)/J;

		#ifdef topk_debug
		cout << "psi[minfreq] " << psi[minfreq] << endl;
		cout << "delta " << delta << endl;
		cout << "FWER " << FWER << endl;
		#endif

		#ifdef topk_debug
		cout << "______________________________________ " << endl;
		#endif


	}

}


/* -------------------FUNCTIONS TO PROCESS A NEWLY FOUND TESTABLE HYPOTHESIS-------------------------------------- */

/* This code contains 3 difference functions to process newly found hypotheses. All of them are virtually identical
 * and the only thing which differs is the way the function receives the list of observations (transactions) for
 * which the hypothesis has X=1.
 * LCM has a complex structure, with frequent itemsets being found at 4 different locations in the source code
 * and under 3 different circumstances. Therefore it was needed to introduce differentiation in the way the transaction
 * list is fed to the "solution processing functions" in order to keep the "transaction keeping" overhead minimal.
 *
 * To reuse this code for problems other than frequent itemset mining, the only thing that needs to be modified
 * is the line which computes the cell counts, for example, the following line in bm_process_solution:
 * 		for(i=0; i<current_trans.siz; i++) a += labels_perm[j][current_trans.list[i]];
 * 	There, current_trans.siz is the number of transactions for which the hypothesis has X=1, i.e. the margin x
 * 	of the 2x2 contingency table (note in this case it is redundant with the input argument x of the function)
 * 	Similarly, current_trans.list[i] with i ranging from 0 to (x-1) has the list of indices of the observations
 * 	for which X=1.
 * 	Simply changing those two parameters accordingly allows code reuse.
 * */

/* Process a solution involving the bitmap represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// void bm_process_solution(int x, int item, int *mask){
void bm_process_solution(int x, bounds_info& parent_info, bounds_info& current_info){
	int i,j;//Loop iterators
	double pval; //Variable to hold p-values
	char *labels_perm_aux;
	/* First, process the new hypothesis */
	++explored;

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;
	++explored_testable;

	// this is not a closed pattern and top-k strategy is not enabled, therefore nothing needs to be done
	if(x == parent_info.frequency && K <= 0) return;

	// check if monotonicity holds
	if(x > parent_info.frequency) cout << "Problem with parent_info! x = " << x << " parent_info.frequency = " << parent_info.frequency << endl;

	if(parent_info.max_as < 0){
		parent_info.max_as = n;
		parent_info.min_as = n;
		parent_info.as = n;
		parent_info.frequency = N;
	}

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	if (last_precomputation != x){
		precompute_pvals(x);
		n_pvalues_computed++; //Update profiling variable
		last_precomputation = x;
	}

	int a_s = 0;
	current_info.as = -1;
	int current_as_max_bound = min(x , n);
	int current_as_min_bound = max(0 , n - (N - x));

	// Top-K strategy
	if(K > 0){

		// compute some bounds on the pattern distribution non the non permuted class label
		current_as_max_bound = min(min(x , n) , min(x , parent_info.as) );
		current_as_min_bound = max(0 , parent_info.as - (parent_info.frequency - x));

		if(min( hypergeom_cdf[current_as_min_bound] , hypergeom_cdf[current_as_max_bound] ) <= last_queue_element){

			// compute the support on minority class of current pattern
			for(i=0; i<x; i++)
				a_s += labels[OCC_VEC[i]];

			current_info.as = a_s;

			// check if the p-value contribute to the top-k strategy
			/*cout << x  << endl;
			cout << common_transactions  << endl;
			cout << last_processed_support  << endl;
			cout << last_queue_element  << endl;
			cout << delta << endl;
			cout << hypergeom_cdf[a_s] << endl;*/
			if(hypergeom_cdf[a_s] < last_queue_element){
				// add the p-value to the queue
				addTopKElement(hypergeom_cdf[a_s]);

				/*cout << "added p-value " << hypergeom_cdf[a_s] << " at support " << x << " and as " << a_s << endl;
				cout << "x " << x << endl;
				cout << "last_processed_support " << last_processed_support << endl;
				cout << "common_transactions " << common_transactions << endl;*/
			}

		}
	}

	// permutations testing

	// this is not a closed pattern, therefore nothing needs to be done
	if(x == parent_info.frequency) return;

	// count how many closed patterns are visited
	++explored_closed;

	// compute a tighter bound on the minimum attainable p-value using parent pattern
	current_as_max_bound = min(min(x , n) , min(x , parent_info.max_as) );
	current_as_min_bound = max(0 , parent_info.min_as - (parent_info.frequency - x));
	current_info.max_as = current_as_max_bound;
	current_info.min_as = current_as_min_bound;

	if(min( hypergeom_cdf[current_as_min_bound] , hypergeom_cdf[current_as_max_bound] ) > delta)
		return; // this pattern will not improve our estimate of the lower alpha quantile of min p-values

	// compute how many transactions are in common with last tested pattern
	int common_transactions = 0;
	i=0; j=0;
	while(i<x && j<last_processed_support){
		if(OCC_VEC[i] == last_processed_transactions[j]){
			++i;
			++j;
			++common_transactions;
		}
		else{
			if(OCC_VEC[i] > last_processed_transactions[j]){
				++j;
			}
			else{
				++i;
			}
		}
	}

	if(common_transactions == x && common_transactions == last_processed_support){
		// this pattern has been already tested, so nothing needs to be done
		return;
	}

	// compute a tighter bound on the minimum attainable p-value using last processed pattern
	int current_as_max_bound_ = min(min(x , n) , min(common_transactions , last_processed_max_as) + (x - common_transactions) );
	int current_as_min_bound_ = max(0 , last_processed_min_as - (last_processed_support - common_transactions));
	current_as_max_bound = min(current_as_max_bound , current_as_max_bound_);
	current_as_min_bound = max(current_as_min_bound , current_as_min_bound_);
	current_info.max_as = current_as_max_bound;
	current_info.min_as = current_as_min_bound;
	// update the bounds of the father
	parent_info.max_as = min(n , min(parent_info.max_as , current_as_max_bound + (parent_info.frequency - x)));
	parent_info.min_as = max(parent_info.min_as , current_as_min_bound);

	int minimum_as_using_parent = current_as_min_bound;
	int maximum_as_using_parent = current_as_max_bound;

	if(min( hypergeom_cdf[current_as_min_bound_] , hypergeom_cdf[current_as_max_bound_] ) > delta)
		return; // this pattern will not improve our estimate of the lower alpha quantile of min p-values

	// if execution reaches this point, the pattern needs to be tested on the permutations

	++tested;

	// the current will be the new last processed pattern
	for(i=0; i<x; i++) {last_processed_transactions[i] = OCC_VEC[i]; }
	last_processed_support = x;
	last_processed_as = a_s;

	// Compute cell-counts for all J-permutations
	// for(i=0; i<current_trans.siz; i++){
	for(i=0; i<x; i++){
		labels_perm_aux = labels_perm[OCC_VEC[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
		for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
	}
	n_cellcounts_computed += J; //Update profiling variable
	effective_total_dataset_frq += x; // Update profiling variable

	last_processed_min_as = x;
	last_processed_max_as = 0;
	for(j=0; j<J; j++){
		// Fetch the precomputed p-value
		last_processed_min_as = min(last_processed_min_as , a_cnt[j]);
		last_processed_max_as = max(last_processed_max_as , a_cnt[j]);
		pval = hypergeom_cdf[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in bm_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
			if( (pval<=delta) && (min_pval[j]>delta)) FWER += ((double)1)/J;
			min_pval[j] = pval;
		}
	}
	current_info.max_as = last_processed_max_as;
	current_info.min_as = last_processed_min_as;

	if(last_processed_min_as < minimum_as_using_parent || last_processed_max_as > maximum_as_using_parent){
		cout << "found lower value than bound! " << endl;
		cout << "x " << x << endl;
		cout << "parent_info.frequency " << parent_info.frequency << endl;
		cout << "parent_info.max_as " << parent_info.max_as << endl;
		cout << "parent_info.min_as " << parent_info.min_as << endl;
		cout << "maximum_as_using_parent " << maximum_as_using_parent << endl;
		cout << "minimum_as_using_parent " << minimum_as_using_parent << endl;
		cout << "last_processed_min_as " << last_processed_min_as << endl;
		cout << "last_processed_max_as " << last_processed_max_as << endl;
	}

	// update the bounds of the father
	parent_info.max_as = min(n , min(parent_info.max_as , current_as_max_bound + (parent_info.frequency - x)));
	parent_info.min_as = max(parent_info.min_as , current_as_min_bound);

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
		//printf("Threshold change BM\n");
		decrease_threshold();
		// Correct possible corruption of LCM data structures due to unexpected change in minimum support
		/*
		for(i=0; i<item; i++){
			//printf("Item %d, Frq %d, Current_th %d\n",i,LCM_Ofrq[i],LCM_th);
			if(LCM_Ofrq[i]==(LCM_th-1)){
				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",i,item,LCM_th);
				LCM_BM_occurrence_delete(i);
				*mask &= ~BITMASK_1[i];
				//printf("Problem fixed!\n");
			}
		}
		*/
	}
}


/* AUXILIARY FUNCTIONS */
// Comparison function used by quicksort implementation in C
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}

#endif
