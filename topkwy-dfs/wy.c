#ifndef _wy_c_
#define _wy_c_

/* LIBRARY INCLUDES */
#include<math.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"permutation.c"
#include"transaction_keeping.c"
#include"lcm_var.c"
#include"priority_queue.c"
/* CONSTANT DEFINES */

#define min(a,b) (a < b) ? a : b
#define max(a,b) (a > b) ? a : b


/* GLOBAL VARIABLES */
FILE* results_file, *minpvals_file, *reports_file, *minsupp_file;
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
// Array for storing values of the PDF of the hypergeometric distribution for fast p-value computation
// and eventually the p-values themselves (the array is reused)
double *hypergeom_pvals;

// store statistics
int explored_itemsets;
int last_support;
int *supports;
double start_instant;
int last_explored_report;
int tested;

// top-k strategy
int K;
heap_t *topk_queue;
heap_t *alpha_quantile;

// bounds using last tested pattern
int last_precomputation;
int last_processed_support;
int last_processed_as;
int last_processed_max_as;
int last_processed_min_as;
int *last_processed_transactions;
int *temp_transaction_list;
double WY_time;
// Cell-count counter
int *a_cnt;


struct bounds_info{
  int frequency;
  int as;
  int max_as;
  int min_as;
};

void init_bounds(struct bounds_info *current_){
  current_->frequency = 0;
  current_->as = 0;
  current_->max_as = 0;
  current_->min_as = 0;
}

void reset_bounds(struct bounds_info *current_){
  current_->max_as = min(current_->frequency , n);
  current_->min_as = max(n - (N - current_->frequency) , 0);
}


/* FUNCTION DECLARATIONS */
void precompute_pvals(int);
void loggamma_init();
void psi_init();
int doublecomp(const void*,const void*);

// Profiling variables
long long n_pvalues_computed;
long long n_cellcounts_computed;
long long effective_total_dataset_frq;
long long total_common;
long long total_supp;

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void start_measuring_time(){
	start_instant = get_cpu_time();
}

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void wy_init(double target_fwer , int k_){
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

	// Allocate memory for the precomputed p-values of the hypergeometric distribution
	// (worst case memory requirement n+1), raising an error if it fails
	hypergeom_pvals = (double *)malloc((n+1)*sizeof(double));
	if(!hypergeom_pvals){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_pvals\n");
		exit(1);
	}

	// Allocate memory for cell counts, raising an error if it fails
	a_cnt = (int *)malloc(J*sizeof(int));
	if(!a_cnt){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array a_cnt\n");
		exit(1);
	}
	for(j=0; j<J; j++) a_cnt[j] = 0;

	n_pvalues_computed = 0; n_cellcounts_computed = 0; effective_total_dataset_frq = 0; //Init profiling variables

	printf("Here done\n");

	explored_itemsets = 0;
	last_support = 0;
	last_explored_report = 0;
	supports = (int *)malloc((N+1) * sizeof(int));
	for(j=0; j<=N; j++) supports[j] = 0;

	K = k_;

	reports_file = fopen("reports.txt","w");
	minsupp_file = fopen("minsupp.txt","w");
	last_processed_transactions = (int *)malloc((N+1) * sizeof(int));
	temp_transaction_list =  (int *)malloc((N+1) * sizeof(int));
	topk_queue = (heap_t *)calloc(1, sizeof (heap_t));
  alpha_quantile = (heap_t *)calloc(1, sizeof (heap_t));
  printf("size of queue: %d\n", topk_queue->len);
  push(topk_queue, 1.0, 0);
  pop(topk_queue);
  printf("size of queue: %d\n", topk_queue->len);
	last_precomputation = -1; last_processed_support = 0;
	last_processed_as = 0; last_processed_max_as = 0; last_processed_min_as = 0;
	tested = 0; WY_time = 0.0;
  /*push(topk_queue, 3.1, "i3");
  push(topk_queue, 4.1, "i4");
  push(topk_queue, 5.1, "i5");
  push(topk_queue, 1.1, "i1");
  push(topk_queue, 2.1, "i2");
  printf("size: %d\n", topk_queue->size);
  printf("len: %d\n", topk_queue->len);
  for (j = 0; j < 5; j++) {
      printf("%s\n", pop(topk_queue));
  }*/
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

	// Correct minimum attainable P-value in some edge-cases
	if((N % 2)==0){
		if (n == (N/2)) for(x=1; x<N; x++) psi[x] *= 2;
		else psi[N/2] *= 2;
	}
	printf ("N = %d , n = %d \n",N,n);
	int aj=0;
	for(aj=0;aj<20 && aj<N+1; aj++){
		printf ("psi[%d] = %e \n",aj,psi[aj]);
	}
}

/* Free all allocated memory and give some output for debugging purposes */
void wy_end(){
	int j, idx_max;
	double delta_corrected = delta;

  printf("size of queue: %d\n", topk_queue->len);
  printf("delta: %e\n", delta);
  printf("top(topk_queue): %e\n",top(topk_queue));
  
  //while(topk_queue->len > 0 && top(topk_queue) == 1.0) pop(topk_queue);

  /*double __delta = top(topk_queue);
  j = 0;
  while(psi[j] > __delta) ++j;
  printf("min supp: %d\n", j);

  for (j = 0; j < 10; j++) {
      printf("%e\n", pop(topk_queue));
  }*/

	// Sort p-values
	/*qsort(min_pval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; delta_corrected = min_pval[idx_max];
  printf("delta_corrected 1: %e\n",delta_corrected);
	// Check and correct (if necessary) boundary cases
	if(delta_corrected==min_pval[idx_max+1]){
		while(min_pval[--idx_max]==delta_corrected);
		delta_corrected = min_pval[idx_max];
    printf("delta_corrected 2: %e\n",delta_corrected);
	}

  printf("top(topk_queue): %e\n",top(topk_queue));
	delta_corrected = min(top(topk_queue) , delta_corrected);
  printf("delta_corrected 3: %e\n",delta_corrected);
*/
	// update FWER
    FWER = 0.0;
  	for(j=0; j<J; j++) FWER += (min_pval[j]<=delta_corrected) ? 1.0 : 0.0;
  	FWER = (FWER)/J;

	// find number of results
	while(topk_queue->len > 0 && top(topk_queue) > delta){
		pop(topk_queue);
	}
	int number_of_significant_patterns = topk_queue->len;

	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",delta_corrected);
	fprintf(results_file,"\t FWER at corrected significance threshold: %e\n",floor(idx_max+1)/J);
	fprintf(results_file,"\t Final LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Final P-value lower bound: %e\n",delta);
	fprintf(results_file,"\t FWER at final P-value lower bound: %e\n",FWER);
	fprintf(results_file,"\t Number of explored itemsets: %d\n",explored_itemsets);
	fprintf(results_file,"\t Number of tested itemsets: %d\n",tested);

	fprintf(results_file,"\t Supports: \n");

	int aj=0;
	for(aj=0; aj < N; aj++){
		fprintf(results_file,"%d : ",aj);
		fprintf(results_file,"%d\n",supports[aj]);
	}

  printf("\n\n");
	printf("Number of explored itemsets: %d\n",explored_itemsets);
	printf("Number of tested itemsets: %d\n",tested);
	printf("Corrected significance threshold: %e\n",delta_corrected);
	printf("FWER: %e\n",FWER);
	printf("Number of significant patterns: %d\n",number_of_significant_patterns);
	printf("Runtime for correction (s): %d\n",(int)(get_cpu_time() - start_instant));
	printf("Runtime for WY computations (s): %d\n",(int)WY_time);
	printf("Peak Memory (MB): %d\n",(int)measurePeakMemory());
	printf("\t total_common: %d\n",total_common);
	printf("\t total_supp: %d\n",total_supp);
	printf("\t total_ratio: %d\n", (double)total_common / (double)total_supp );

	fprintf(minpvals_file,"MINIMUM P-VALS (%d PERMUTATIONS)\n",J);
	for(j=0;j<(J-1);j++) fprintf(minpvals_file,"%e,",min_pval[j]);
	fprintf(minpvals_file,"%e\n",min_pval[J-1]);

	// Report the seed for reproducibility
	fprintf(minpvals_file,"\nRANDOM SEED USED FOR KEY GENERATION\n");
	fprintf(minpvals_file,"\t Seed = %u\n",(unsigned)seed);

	// Free allocated memory
	free(loggamma);
	free(psi);
	free(min_pval);
	free(hypergeom_pvals);
	free(a_cnt);
	free(supports);

	free(topk_queue);
  free(alpha_quantile);
	free(last_processed_transactions);
	free(temp_transaction_list);

	// Close files
	fclose(results_file); fclose(minpvals_file);
}

/* -------------------------------FUNCTIONS TO COMPUTE FISHER EXACT TEST P-VALUES----------------------------------- */

/* This function precomputes all Fisher exact test P-values for a contingency table with margins x,n,N that is,
 * all p-values p(a,x,n,N) for a in the range [max(0,n+x-N),min(x,n)]. The results will be stored in the array
 * hypergeom_pvals such that p(a,x,n,N)=hypergeom_pvals[a]. Note that values hypergeom_pvals[a] for a outside
 * [max(0,n+x-N),min(x,n)] are undefined and could contain garbage of previous hypotheses.
 * */
void precompute_pvals(int x){
  if (last_precomputation == x){
		return;
	}
	last_precomputation = x;
	double pre_comp_xterms, pval, p_left, p_right;
	int a, a_min, a_max;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)

	// Precompute the hypergeometric PDF in the range of interest
	for(a=a_min; a<=a_max; a++) hypergeom_pvals[a] = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));

	// The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
	// hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
	// that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. When a value is "accepted", we know its
	// respective p-value because due to the way we explore the hypergeometric pdf, there can be no other values larger than it. Therefore, everytime
	// a value is "accepted", we store the pvalue in hypergeom_pvals.
	// The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
	// that case is by "accepting" both values simultaneously.
	pval = 0;
	while(a_min<a_max){
		p_left = hypergeom_pvals[a_min]; p_right = hypergeom_pvals[a_max];
		if(p_left == p_right) { pval += (p_left+p_right); hypergeom_pvals[a_min++] = pval; hypergeom_pvals[a_max--] = pval; }
		else if(p_left < p_right){ pval += p_left; hypergeom_pvals[a_min++] = pval;}
		else{ pval += p_right; hypergeom_pvals[a_max--] = pval;}
	}
	// In this case a_min=a_max is the mode of the distribution and its p-value is 1 by definition
	if(a_min==a_max) hypergeom_pvals[a_max] = 1;
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
		if (psi[sl1] >= psi[sl2]) delta = psi[sl1];
		else{ delta = psi[sl2]; flag = 0; }
		//Update LCM minimum support
		LCM_th = sl1;
		//printf("\n\n\nTHRESHOLD CHANGE!!! NEW THRESHOLD=%d\n\n\n",LCM_th);
    printf("WY: sigma %d, delta %e\n",LCM_th,delta);
		/*printf("THRESHOLD CHANGE!!! NEW THRESHOLD=%d\n",LCM_th);
		printf("current explored patterns:%d\n",explored_itemsets);
		printf("last_support:%d\n",last_support);
		printf("elapsed:%d\n",(int)(get_cpu_time() - start_instant));*/
		fprintf(minsupp_file,"%d;%d;%d;\n",explored_itemsets,last_support,(int)(get_cpu_time() - start_instant));
		fflush(minsupp_file);
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		sl2--; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]){ delta = psi[sl1]; flag = 1; }
		else delta = psi[sl2];
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
	// Recompute FWER from scratch
	false_positives = 0;
	for(j=0; j<J; j++) false_positives += (min_pval[j]<=delta) ? 1 : 0;
	FWER = ((double)false_positives)/J;

}


void update_queue_values(){
	while(top(alpha_quantile) > min_pval[ top_value(alpha_quantile) ]){
		/*cout << " removed outdated p-value at index " << top_value(alpha_quantile) << endl;
		cout << " top(alpha_quantile) " << top(alpha_quantile) << endl;
		cout << " min_pval[ top_value(alpha_quantile) ] " << min_pval[ top_value(alpha_quantile) ] << endl;*/
		int temp_alpha_item = top_value(alpha_quantile);
		pop(alpha_quantile);
		push(alpha_quantile , min_pval[ temp_alpha_item ] , temp_alpha_item);
	}
}

void update_delta(){
	/*cout << "------------------------" << endl;
	cout << "Increasing delta to control FWER " << endl;
	cout << "alpha_quantile->len (0) " << alpha_quantile->len << endl;
	cout << "below_delta " << below_delta << endl;*/

	/* delta needs to be lowered since FWER > 0.05 */
	int queue_not_empty = alpha_quantile->len > 0;
	int not_updated_delta = top(alpha_quantile) == delta;
  int below_delta = 0; int j;
  for(j=0; j<J; j++) below_delta += (min_pval[j]<=delta) ? 1 : 0;
	int FWER_not_controlled = (double)below_delta / (double)J > alpha;

	while(queue_not_empty && (FWER_not_controlled || not_updated_delta)){

		/*cout << "alpha_quantile->len (1) " << alpha_quantile->len << endl;
		cout << "below_delta " << below_delta << endl;*/

		// remove all elements of alpha_quantile with outdated p-value

		update_queue_values();

		delta = min(delta , top(alpha_quantile));

		pop(alpha_quantile);
		--below_delta;
		FWER_not_controlled = (double)below_delta / (double)J > alpha;

		queue_not_empty = alpha_quantile->len > 0;
		if(queue_not_empty){
			update_queue_values();
			not_updated_delta = top(alpha_quantile) == delta;
		}
		else{
			printf(" alpha-quantile queue got empty! ");
      printf(" delta was %e ",delta);
			delta = 0.999 * delta;
      printf(" delta set to %e ",delta);
			//abort();
		}

	}
	/*cout << "out of main loop with below_delta reduced to " << below_delta << endl;
	cout << "not_updated_delta " << not_updated_delta << endl;
	cout << "FWER_not_controlled " << FWER_not_controlled << endl;
	cout << "queue_not_empty " << queue_not_empty << endl;*/

	if(queue_not_empty){
		update_queue_values();
		delta = min(delta , top(alpha_quantile));
	}
	else{
    printf(" alpha-quantile queue got empty! ");
    printf(" delta was %e ",delta);
    delta = 0.999 * delta;
    printf(" delta set to %e ",delta);
		//abort();
	}

		//cout << "delta (after delta lift)" << delta << endl;
		/* update minsupp if needed */
		while(psi[LCM_th] > delta){
			++LCM_th;
      printf("WY: sigma %d, delta %e\n",LCM_th,delta);
		}

		below_delta = 0;
  	for(j=0; j<J; j++) below_delta += (min_pval[j]<=delta) ? 1 : 0;
  	FWER = ((double)below_delta)/J;


    while(topk_queue->len > 0 && top(topk_queue) > delta){
      pop(topk_queue);
    }

		/*cout << "below_delta is now " << below_delta << endl;
		cout << "FWER " << ((double)below_delta / (double)jp) << endl;
		cout << "increased sigma to " << sigma << endl;
		cout << "delta is now " << delta << endl;
		cout << "psi[sigma] " << psi[sigma] << endl;//}
		cout << "alpha_quantile->len " << alpha_quantile->len << endl;*/


	//cout << "------------------------" << endl;
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

void insert_in_topk_queue(double pvalue){
   //printf("inserted %e\n",pvalue);
   push(topk_queue, pvalue, 0);
}

int count_common_transactions(int* current_transaction_list , int current_support){
  int common_transactions = 0;
	int i=0; int j=0;
  // count how many transactions are in common between last and current
	while(i<current_support && j<last_processed_support){
		if(current_transaction_list[i] == last_processed_transactions[j]){
			++i;
			++j;
			++common_transactions;
		}
		else{
			if(current_transaction_list[i] > last_processed_transactions[j]){
				++j;
			}
			else{
				++i;
			}
		}
	}
  total_common += common_transactions;
  total_supp += current_support;
  return common_transactions;
}

double get_parent_bound(struct bounds_info *parent_info , struct bounds_info *current_info){

  int x = current_info->frequency;
  //int common_transactions = x;
  int a_ = min(x , n);
  int b_ = min(x , parent_info->max_as);
  int c_ = max(0 , n - (N - x));
  int d_ = parent_info->min_as - (parent_info->frequency - x);
  int current_as_max_bound_ = min(a_ , b_);
  int current_as_min_bound_ = max(c_ , d_);
  if(current_as_min_bound_ < c_){
    printf("problem with min bound; set to %d\n", c_ );
    current_as_min_bound_ = c_;
  }
  if(current_as_max_bound_ > a_){
    printf("problem with max bound; set to %d\n", a_ );
    current_as_min_bound_ = a_;
  }
  double bound = min( hypergeom_pvals[current_as_min_bound_] , hypergeom_pvals[current_as_max_bound_] );
  /*printf("Last: %d, curr: %d, maxlast: %d, minlast=%d\n",last_processed_support,x,last_processed_max_as,last_processed_min_as);
  printf("  current_as_max_bound_: %d, current_as_min_bound_: %d\n",current_as_max_bound_,current_as_min_bound_);
  printf("  Common: %d, bound: %e, delta: %e\n",common_transactions,bound,delta);
  printf("   A: %d\n",min(common_transactions , last_processed_max_as) + (x - common_transactions));
  printf("   B: %d\n",last_processed_min_as - (last_processed_support - common_transactions));
  printf("   C: %e\n",hypergeom_pvals[current_as_min_bound_]);
  printf("   D: %e\n",hypergeom_pvals[current_as_max_bound_]);*/
  current_info->max_as = current_as_max_bound_;
  current_info->min_as = current_as_min_bound_;
  return bound;

}

double get_last_processed_bound(int current_support , int common_transactions , struct bounds_info *current_info){

  /*if(common_transactions == 0) */return psi[current_support];

  int x = current_support;
  int a_ = min(x , n);
  int b_ = min(common_transactions , last_processed_max_as + (x - common_transactions));
  int c_ = max(0 , n - (N - current_support));
  int d_ = last_processed_min_as - (last_processed_support - common_transactions);
  int current_as_max_bound_ = min(a_ , b_);
  int current_as_min_bound_ = max(c_ , d_);
  if(current_as_min_bound_ < c_){
    printf("problem with min bound; set to %d\n", c_ );
    current_as_min_bound_ = c_;
  }
  if(current_as_max_bound_ > a_){
    printf("problem with max bound; set to %d\n", a_ );
    current_as_min_bound_ = a_;
  }
  double bound = min( hypergeom_pvals[current_as_min_bound_] , hypergeom_pvals[current_as_max_bound_] );
  /*printf("Last: %d, curr: %d, maxlast: %d, minlast=%d\n",last_processed_support,x,last_processed_max_as,last_processed_min_as);
  printf("  current_as_max_bound_: %d, current_as_min_bound_: %d\n",current_as_max_bound_,current_as_min_bound_);
  printf("  Common: %d, bound: %e, delta: %e\n",common_transactions,bound,delta);
  printf("   A: %d\n",min(common_transactions , last_processed_max_as) + (x - common_transactions));
  printf("   B: %d\n",last_processed_min_as - (last_processed_support - common_transactions));
  printf("   C: %e\n",hypergeom_pvals[current_as_min_bound_]);
  printf("   D: %e\n",hypergeom_pvals[current_as_max_bound_]);*/
  current_info->max_as = min(current_info->max_as , current_as_max_bound_);
  current_info->min_as = max(current_info->min_as , current_as_min_bound_);
  return bound;

}

/* Process a solution involving the bitmap represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
void bm_process_solution(int x, int item, int *mask , struct bounds_info* parent_info , struct bounds_info* current_info){
	int i,j;//Loop iterators
	double pval; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	double temp_time_ = get_cpu_time();

  if(current_info->frequency <= 0){
    current_info->frequency = x;
    current_info->as = 0;
    reset_bounds(current_info);
  }
  if(x != current_info->frequency){
    printf("problem in bm_process_solution: supp %d, current supp %d\n",x,current_info->frequency);
  }
  if(parent_info->frequency <= 0){
    parent_info->frequency = N;
    parent_info->as = n;
    parent_info->max_as = n;
    parent_info->min_as = n;
  }
  if(x < N && x >= parent_info->frequency){
    printf("problem in bm_process_solution: supp %d, parent supp %d\n",x,parent_info->frequency);
  }

	++explored_itemsets;
	++supports[x];

	last_support = x;

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	precompute_pvals(x);
	n_pvalues_computed++; //Update profiling variable

	// Compute as
  int a_s = 0;
	for(i=0; i<current_trans.siz; i++){
		a_s += labels[current_trans.list[i]];
	}
  current_info->as = a_s;

  //printf("processing pattern with xs=%d, as=%d, pval=%e.\n",x,a_s,hypergeom_pvals[ a_s ]);
  double pval_ = hypergeom_pvals[ a_s ];
  if( pval_ <= delta || (K > 0 && topk_queue->len < K) ) insert_in_topk_queue(pval_);

  /* top-k strategy */
  while(K > 0 && topk_queue->len > K){
    pop(topk_queue);
    if(top(topk_queue) < delta){
      delta = top(topk_queue);
      // update FWER
      FWER = 0.0;
    	for(j=0; j<J; j++) FWER += (min_pval[j]<=delta) ? 1.0 : 0.0;
    	FWER = (FWER)/J;
      while(psi[LCM_th] > delta){
        ++LCM_th;
        printf("Top-k strategy: sigma %d, delta %e\n",LCM_th,delta);
        // Correct possible corruption of LCM data structures due to unexpected change in minimum support
        for(i=0; i<item; i++){
    			//printf("Item %d, Frq %d, Current_th %d\n",i,LCM_Ofrq[i],LCM_th);
    			if(LCM_Ofrq[i]==(LCM_th-1)){
    				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",i,item,LCM_th);
    				LCM_BM_occurrence_delete(i);
    				*mask &= ~BITMASK_1[i];
    				//printf("Problem fixed!\n");
    			}
    		}
    	}
    }
  }
  /* top-k strategy */

  WY_time += get_cpu_time() - temp_time_;
  temp_time_ = get_cpu_time();


  // compute a bound on the current pattern using the parent pattern
  if (get_parent_bound(parent_info , current_info) > delta )
    return;

  // compute a bound on the current pattern using the last processed pattern
  if (get_last_processed_bound(x , count_common_transactions(current_trans.list , x) , current_info) > delta )
    return;

  ++tested;

  // testing the pattern on the permutations
  // Compute cell-counts for all J-permutations
	for(i=0; i<current_trans.siz; i++){
		labels_perm_aux = labels_perm[current_trans.list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
		for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
	}
	n_cellcounts_computed += J; //Update profiling variable
	effective_total_dataset_frq += x; // Update profiling variable

  // the current will be the new last processed pattern
	for(i=0; i<x; i++) {last_processed_transactions[i] = current_trans.list[i]; }
	last_processed_support = x;
	last_processed_as = a_s;
  last_processed_min_as = x;
	last_processed_max_as = 0;

	// If not, compute permuted P-values for each permutation
	for(j=0; j<J; j++){
		last_processed_min_as = min(last_processed_min_as , a_cnt[j]);
		last_processed_max_as = max(last_processed_max_as , a_cnt[j]);
		// Fetch the precomputed p-value
		pval = hypergeom_pvals[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in bm_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
			if( (pval<=delta) && (min_pval[j]>delta)) {
        FWER += ((double)1)/J;
        push(alpha_quantile , pval , j);
      }
			min_pval[j] = pval;
		}
	}
  current_info->max_as = last_processed_max_as;
  current_info->min_as = last_processed_min_as;

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
		//printf("Threshold change BM\n");
		//decrease_threshold();
    update_delta();
		// Correct possible corruption of LCM data structures due to unexpected change in minimum support
		for(i=0; i<item; i++){
			//printf("Item %d, Frq %d, Current_th %d\n",i,LCM_Ofrq[i],LCM_th);
			if(LCM_Ofrq[i]==(LCM_th-1)){
				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",i,item,LCM_th);
				LCM_BM_occurrence_delete(i);
				*mask &= ~BITMASK_1[i];
				//printf("Problem fixed!\n");
			}
		}
	}
	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif

	if(tested - last_explored_report > 99999){
		printf("NEW REPORT:\n");
		printf("  explored patterns:%d\n",explored_itemsets);
		printf("  current tested patterns:%d\n",tested);
		printf("  last_support:%d\n",last_support);
		printf("  delta:%e\n",delta);
		printf("  elapsed:%d\n",(int)(get_cpu_time() - start_instant));
		last_explored_report = tested;
	}

	WY_time += get_cpu_time() - temp_time_;

}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x , struct bounds_info* parent_info , struct bounds_info* current_info){
	int i,j;//Loop iterators
	double pval; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	double temp_time_ = get_cpu_time();

  if(current_info->frequency <= 0){
    current_info->frequency = x;
    current_info->as = 0;
    reset_bounds(current_info);
  }
  if(x != current_info->frequency){
    printf("problem in process_solution0: supp %d, current supp %d\n",x,current_info->frequency);
  }
  if(parent_info->frequency <= 0){
    parent_info->frequency = N;
    parent_info->as = n;
    parent_info->max_as = n;
    parent_info->min_as = n;
  }
  if(x < N && x >= parent_info->frequency){
    printf("problem in process_solution0: supp %d, parent supp %d\n",x,parent_info->frequency);
  }

	++explored_itemsets;
	++supports[x];

	last_support = x;

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	precompute_pvals(x);
	n_pvalues_computed++; //Update profiling variable

	// Compute cell-counts for all permutations
  int a_s = 0;
	for(i=0; i<bm_trans_list[1].siz; i++){
		a_s += labels[bm_trans_list[1].list[i]];
	}
  current_info->as = a_s;


  //printf("processing pattern with xs=%d, as=%d, pval=%e.\n",x,a_s,hypergeom_pvals[ a_s ]);
  double pval_ = hypergeom_pvals[ a_s ];
  if( pval_ <= delta || (K > 0 && topk_queue->len < K) ) insert_in_topk_queue(pval_);

  /* top-k strategy */
  while(K > 0 && topk_queue->len > K){
    pop(topk_queue);
    if(top(topk_queue) < delta){
      delta = top(topk_queue);
      // update FWER
      FWER = 0.0;
    	for(j=0; j<J; j++) FWER += (min_pval[j]<=delta) ? 1.0 : 0.0;
    	FWER = (FWER)/J;
      while(psi[LCM_th] > delta){
        ++LCM_th;
        printf("Top-k strategy: sigma %d, delta %e\n",LCM_th,delta);
    	}
    }
  }

  WY_time += get_cpu_time() - temp_time_;
  temp_time_ = get_cpu_time();

  // compute a bound on the current pattern using the parent pattern
  if (get_parent_bound(parent_info , current_info) > delta )
    return;

  // compute a bound on the current pattern using the last processed pattern
  if (get_last_processed_bound(x , count_common_transactions(bm_trans_list[1].list , x) , current_info) > delta )
    return;

  ++tested;

  // testing the pattern on the permutations
  // Compute cell-counts for all J-permutations
  for(i=0; i<bm_trans_list[1].siz; i++){
		labels_perm_aux = labels_perm[bm_trans_list[1].list[i]]; //Avoid recomputing labels_perm[bm_trans_list[1].list[i]] all the time
		for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
	}
	n_cellcounts_computed += J; //Update profiling variable
	effective_total_dataset_frq += x; // Update profiling variable

  // the current will be the new last processed pattern
	for(i=0; i<x; i++) {last_processed_transactions[i] = bm_trans_list[1].list[i]; }
	last_processed_support = x;
	last_processed_as = a_s;
	last_processed_min_as = x;
	last_processed_max_as = 0;

	// If not, compute permuted P-values for each permutation
	for(j=0; j<J; j++){
		last_processed_min_as = min(last_processed_min_as , a_cnt[j]);
		last_processed_max_as = max(last_processed_max_as , a_cnt[j]);
		// Fetch the precomputed p-value
		pval = hypergeom_pvals[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in process_solution0!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
      if( (pval<=delta) && (min_pval[j]>delta)) {
        FWER += ((double)1)/J;
        push(alpha_quantile , pval , j);
      }
			min_pval[j] = pval;
		}
	}
  current_info->max_as = last_processed_max_as;
  current_info->min_as = last_processed_min_as;

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
		//printf("threshold change 0\n");
		//decrease_threshold();
    update_delta();
	}

	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif

	if(tested - last_explored_report > 99999){
		printf("NEW REPORT:\n");
		printf("  explored patterns:%d\n",explored_itemsets);
		printf("  current tested patterns:%d\n",tested);
		printf("  last_support:%d\n",last_support);
		printf("  delta:%e\n",delta);
		printf("  elapsed:%d\n",(int)(get_cpu_time() - start_instant));
		last_explored_report = tested;
	}

	WY_time += get_cpu_time() - temp_time_;
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask , struct bounds_info* parent_info , struct bounds_info* current_info){
	int i,j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	double pval; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	double temp_time_ = get_cpu_time();


  if(current_info->frequency <= 0){
    current_info->frequency = x;
    current_info->as = 0;
    reset_bounds(current_info);
  }
  if(x != current_info->frequency){
    printf("problem in ary_process_solution: supp %d, current supp %d\n",x,current_info->frequency);
  }
  if(parent_info->frequency <= 0){
    parent_info->frequency = N;
    parent_info->as = n;
    parent_info->max_as = n;
    parent_info->min_as = n;
  }
  if(x < N && x >= parent_info->frequency){
    printf("problem in ary_process_solution: supp %d, parent supp %d\n",x,parent_info->frequency);
  }

	++explored_itemsets;
	++supports[x];

	last_support = x;

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Sanity-check (this one is more complicated due to the way the transactions are stored)
	aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++) aux++;
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	precompute_pvals(x);
	n_pvalues_computed++; //Update profiling variable

	// Compute cell-counts for all permutations
  int a_s = 0;
  j=0;
  for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
    end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
    for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
      a_s += labels[*ptr];
      temp_transaction_list[j] = *ptr;
      ++j;
    }
  }
  current_info->as = a_s;

  //printf("processing pattern with xs=%d, as=%d, pval=%e.\n",x,a_s,hypergeom_pvals[ a_s ]);
  double pval_ = hypergeom_pvals[ a_s ];
  if( pval_ <= delta || (K > 0 && topk_queue->len < K) ) insert_in_topk_queue(pval_);

  /* top-k strategy */
  while(K > 0 && topk_queue->len > K){
    pop(topk_queue);
    if(top(topk_queue) < delta){
      delta = top(topk_queue);
      // update FWER
      FWER = 0.0;
      for(j=0; j<J; j++) FWER += (min_pval[j]<=delta) ? 1.0 : 0.0;
      FWER = (FWER)/J;
      while(psi[LCM_th] > delta){
      ++LCM_th;
      printf("Top-k strategy: sigma %d, delta %e\n",LCM_th,delta);
      // Correct possible corruption of LCM data structures due to unexpected change in minimum support
  		for(j=0; j<LCM_BM_MAXITEM; j++){
  			//printf("Item %d, Frq %d, Current_th %d\n",j,LCM_Ofrq[j],LCM_th);
  			if(LCM_Ofrq[j]==(LCM_th-1)){
  				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",j,item,LCM_th);
  				LCM_BM_occurrence_delete(j);
  				*mask &= ~BITMASK_1[j];
  				//printf("Problem fixed!\n");
          }
        }
      }
    }
  }
  /* top-k strategy */


  WY_time += get_cpu_time() - temp_time_;
  temp_time_ = get_cpu_time();


  // compute a bound on the current pattern using the parent pattern
  if (get_parent_bound(parent_info , current_info) > delta )
    return;


    // compute a bound on the current pattern using the last processed pattern
    if (get_last_processed_bound(x , count_common_transactions(temp_transaction_list , x) , current_info) > delta )
      return;

    ++tested;

    // Compute cell-counts for all permutations
  	/*for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
  		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
  		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
  			labels_perm_aux = labels_perm[*ptr];
  			for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
  		}
  	}*/
    for(i=0; i<x; i++){
  		labels_perm_aux = labels_perm[temp_transaction_list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
  		for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
  	}
  	n_cellcounts_computed += J; //Update profiling variable
  	effective_total_dataset_frq += x; // Update profiling variable

    // the current will be the new last processed pattern
  	for(i=0; i<x; i++) {last_processed_transactions[i] = temp_transaction_list[i]; }
  	last_processed_support = x;
  	last_processed_as = a_s;
    last_processed_min_as = x;
  	last_processed_max_as = 0;


	// If not, compute permuted P-values for each permutation
	for(j=0; j<J; j++){
		last_processed_min_as = min(last_processed_min_as , a_cnt[j]);
		last_processed_max_as = max(last_processed_max_as , a_cnt[j]);
		// Fetch the precomputed p-value
		pval = hypergeom_pvals[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in ary_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
      if( (pval<=delta) && (min_pval[j]>delta)) {
        FWER += ((double)1)/J;
        push(alpha_quantile , pval , j);
      }
			min_pval[j] = pval;
		}
	}
  current_info->max_as = last_processed_max_as;
  current_info->min_as = last_processed_min_as;

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
		//printf("threshold change ary\n");
		//decrease_threshold();
    update_delta();
		// Correct possible corruption of LCM data structures due to unexpected change in minimum support
		for(j=0; j<LCM_BM_MAXITEM; j++){
			//printf("Item %d, Frq %d, Current_th %d\n",j,LCM_Ofrq[j],LCM_th);
			if(LCM_Ofrq[j]==(LCM_th-1)){
				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",j,item,LCM_th);
				LCM_BM_occurrence_delete(j);
				*mask &= ~BITMASK_1[j];
				//printf("Problem fixed!\n");
			}
		}
	}

	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif

	if(tested - last_explored_report > 99999){
		printf("NEW REPORT:\n");
		printf("  explored patterns:%d\n",explored_itemsets);
		printf("  current tested patterns:%d\n",tested);
		printf("  last_support:%d\n",last_support);
		printf("  delta:%e\n",delta);
		printf("  elapsed:%d\n",(int)(get_cpu_time() - start_instant));
		last_explored_report = tested;
	}

	WY_time += get_cpu_time() - temp_time_;

}

/* AUXILIARY FUNCTIONS */
// Comparison function used by quicksort implementation in C
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}

#endif
