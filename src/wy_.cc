#ifndef _wy_c_
#define _wy_c_

#include "patricia.h"
#include "math.h"
#include <queue>
#include <set>

#define p_val_caching true
#define cache_max_memory_mb 10000.0 // 10 GByte for p-values in mb
#define cache_max_memory 10000000000.0 // 10 GByte for p-values in b

int max_cached_pvalues;
// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute
double log_inv_binom_N_n;
// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi; // montone version
double *psi_exact; // exact non-monotone version
// log of psi
double *psi_log; // montone version
double *psi_log_exact; // exact non-monotone version
// Array for storing values of the PDF of the hypergeometric distribution for fast p-value computation
// and eventually the p-values themselves (the array is reused)
double *hypergeom_pvals;
double *hypergeom_pvals_log; // log10 of p-values
double *hypergeom_pvals_buffer;
double *hypergeom_pvals_log_buffer;
int last_support;
// indexes for p-values caching
int last_minsupp;
int last_ssupp;
int current_cache_size;
int current_cache_size_mb;
// Minimum P-value for each permutation
double *min_pval;
// midpoint of interval [0,N], floor(N/2)
int N_over_2;
// corrected significance threshold
double delta;
double log_delta;
// number of p-values below delta
int below_delta;
std::set<int> below_delta_indexes;
// minimum support
int sigma;
// time needed to process solutions
double p_values_time;
// array used to return values
double* to_return;

// array storing indexes with minimum p-values < delta
bool* quantile_indexes;
double pruning_ratio;

double precomputetime;

long long cache_accesses;
long long new_precomputation;

double log10_2;
double log_10;
int log_method;

double start_instant;
double time_to_update_delta;

std::vector< std::vector<double> > p_vals_cached;
std::vector< std::vector<double> > log_p_vals_cached;

std::vector< long > itemsets_LAMP_count;
std::vector< long > itemsets_LAMP_count_naive;


struct alpha_quantile_tuple{
	double p_value;
	int index;
};

class compare_alpha_quantile
{
public:
    bool operator() (alpha_quantile_tuple a, alpha_quantile_tuple b)
    {
        return a.p_value <= b.p_value;
    }
};

alpha_quantile_tuple temp_alpha_item;

std::priority_queue<alpha_quantile_tuple, std::vector<alpha_quantile_tuple>, compare_alpha_quantile> alpha_quantile;

void benchmarklogoperations();

double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
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
	log_inv_binom_N_n = loggamma[n1] + loggamma[N-n1] - loggamma[N];
}


/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double xi1;
	int x, x_init;
	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
  psi_log = (double *)malloc((N+1)*sizeof(double));
  psi_exact = (double *)malloc((N+1)*sizeof(double));
  psi_log_exact = (double *)malloc((N+1)*sizeof(double));
	if(!psi || !psi_log){
		fprintf(stderr,"Error in function psi_and_xi1_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
  psi_log[0] = 0;
	//In this range, the recursion $\psi(x)$=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can be seen to be correct
	for(x=1; x<=n1; x++) psi[x] = (((double)(n1-(x-1)))/(N-(x-1)))*psi[x-1];

  for(x=1; x<=n1; x++) psi_log[x] = log10(((double)(n1-(x-1)))/(N-(x-1))) + psi_log[x-1];

	// Now, start computing xi1 in the range [N-N_over_2,N] using another recursion, this time
	// starting in N
	// Note that we don't need to store all values, since this will be used only to initialise
	// psi[N_over_2]
	x_init = N-N_over_2;
	xi1 = 1;
  double xi1_log = 0;
	//In this range, the recursion $\xi_{1}(x)$=$\xi_{1}(x+1)$*[((x-1)-n)/(x+1)] can be seen to be correct
	for(x=(N-1); x>=x_init; x--) xi1 = (((double)((x+1)-n1))/(x+1))*xi1;

  /*cout << "------------------" << endl;
  cout << "x_init " << x_init << endl;
  cout << "(N-1) " << (N-1) << endl;
  cout << "------------------" << endl;*/

  for(x=(N-1); x>=x_init; x--) {
	/*cout << "x " << x << endl;
	cout << "xi1_log " << xi1_log << endl;
	cout << "log10(((double)((x+1)-n1))/(x+1)) " << log10(((double)((x+1)-n1))/(x+1)) << endl;*/
	xi1_log = log10(((double)((x+1)-n1))/(x+1)) + xi1_log;
  }
  /*cout << "------------------" << endl;
  cout << "xi1_log " << xi1_log << endl;
  cout << "xi1 " << xi1 << endl;*/
	// Now, use the value of $\xi_{1}(N-N_over_2)$=xi1[0] to get $\psi(N_over_2)$=psi[N_over_2] using the
	// same recursion if N is odd, or simply copy the value of xi1[0] since $\xi_{1}(N-N_over_2)=\psi(N_over_2)$
	// if N is even
	if (N % 2) psi[N_over_2] = (((double)(x_init-n1))/x_init)*xi1;
	else psi[N_over_2] = xi1;

  if (N % 2) psi_log[N_over_2] = log10(((double)(x_init-n1))/x_init) + xi1_log;
	else psi_log[N_over_2] = xi1_log;

  /*cout << "psi_log[N_over_2] " << psi_log[N_over_2] << endl;
  cout << "psi[N_over_2] " << psi[N_over_2] << endl;*/

	// Now, using psi[N_over_2] compute the right side of "the W", i.e. the range [n+1,N_over_2]
	// using the same recursion as for $\xi_{1}$
	for(x=(N_over_2-1); x > n1; x--) psi[x] = (((double)((x+1)-n1))/(x+1))*psi[x+1];

  for(x=(N_over_2-1); x > n1; x--) psi_log[x] = log10(((double)((x+1)-n1))/(x+1)) + psi_log[x+1];

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=x_init; x<=N; x++) psi[x] = psi[N-x];

  for(x=x_init; x<=N; x++) psi_log[x] = psi_log[N-x];

	// Correct minimum attainable P-value in some edge-cases
	if((N % 2)==0){
		if (n1 == (N/2)) for(x=1; x<N; x++) psi[x] *= 2;
		else psi[N/2] *= 2;
	}

  if((N % 2)==0){
		if (n1 == (N/2)) for(x=1; x<N; x++) psi_log[x] += log10(2);
		else psi_log[N/2] += log10(2);
	}

  // monotone psi for pattern pruning
  psi_exact[0] = psi[0];
  psi_log_exact[0] = psi_log[0];

	for(int aj=1;aj<N+1; aj++){
	psi_exact[aj] = psi[aj];
	if(psi_exact[aj] > psi[aj-1])
	  psi[aj] = psi[aj-1];

	psi_log_exact[aj] = psi_log[aj];
	if(psi_log_exact[aj] > psi_log[aj-1])
	  psi_log[aj] = psi_log[aj-1];
	}

  /*for(int aj=0;aj<N+1; aj++){
		cout << "psi["<< aj <<"] = " << psi[aj] << endl;
	cout << "psi_log["<< aj <<"] = " << psi_log[aj] << endl;
	cout << "psi_exact["<< aj <<"] = " << psi_exact[aj] << endl;
	cout << "psi_log_exact["<< aj <<"] = " << psi_log_exact[aj] << endl;
}*/

	cout << "psi successfully computed "<< endl;

	last_minsupp = 0;
	last_ssupp = N;


	log10_2 = log10(2.0);
	log_10 = log(10);
	log_method = 1;

}



void initClassLabels(){
	string line;
	string end = "";
	bool t_class = true;
	N = num_tr;
	n1 = 0;
	c_labels.close();
	c_labels.open(c_fileinput.c_str(),ifstream::in); ///<class labels file;
	/*long long perm_elements = ( ((long long)N) * ((long long)(jp + 1)));
	cout << "perm_elements = " << perm_elements << endl;
	permutations = (bool *)calloc(perm_elements , sizeof(bool)); /// class label permutations
	if(!permutations){
		cout << "Error: can't allocate permutation matrix. perm_elements = " << perm_elements << endl;
		abort();
	}
	cout << "permutation matrix allocated " << endl;
	long long i=0;
	long long index = 0;
	long long jp_long = (long long)jp;*/
	long long i=0;
	labels = (bool *)malloc(((long)N+1) * sizeof(bool)); /// class labels
	if(!labels){
		cout << "Error: can't allocate labels matrix." << endl;
		abort();
	}
	//cout << "labels matrix allocated " << endl;
	N=0;
	while(!c_labels.eof()){
		getline (c_labels,line); // read the class label
		//cout << "line = " << line << endl;
		if (end.compare(( char * )line.c_str())!=0){
			t_class = (( char * )line.c_str())[0] == '1'; //strcmp(( char * )line.c_str() , "1") == 0;
			//cout << "t_class = " << t_class << endl;
			/*if(t_class)
				n1++;*/
			n1+=t_class;
			N++;
			/*for(long long j=0; j < jp_long+1; j++){
				index = i * (jp_long + 1) + j;
				permutations[index] = t_class;
			}*/
			labels[i] = t_class;
			i++;
		}
	}

	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)

	double ratio = ((double)n1)/((double)N);
	if(ratio > 0.5){
		n1 = N - n1;
		for(long long i = 0; i<N; i++){
			/*index = i * (jp_long + 1);
			t_class = !permutations[index];
			for(long long j = 0; j<jp_long+1; j++){
				index = i * (jp_long + 1) + j;
				permutations[index] = t_class;
			}*/
			labels[i] = !labels[i];
		}
	}
	cout << "Number of transactions = " << N << ", in class 1 = " << n1 << ", ratio = " << ((double)n1)/((double)N) << endl;

	loggamma_init();
	psi_init();

	// Allocate memory for minimum p-values, raising error if it fails
	/*min_pval = (double *)malloc(J*sizeof(double));
	if(!min_pval){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array min_pval\n");
		exit(1);
	}
	// Initialise all p-values to 1
	for(j=0; j<J; j++) min_pval[j] = 1;*/

	// Allocate memory for the precomputed p-values of the hypergeometric distribution
	// (worst case memory requirement n+1), raising an error if it fails
	hypergeom_pvals = (double *)malloc((n1+1)*sizeof(double));
	hypergeom_pvals_log = (double *)malloc((n1+1)*sizeof(double));
	hypergeom_pvals_buffer = hypergeom_pvals;
	hypergeom_pvals_log_buffer = hypergeom_pvals_log;
	last_support = -1;
	if(!hypergeom_pvals){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_pvals\n");
		exit(1);
	}
	min_pval = (double*)malloc((jp+1)*sizeof(double));
	for(int aj=0;aj<jp+1;aj++)
		min_pval[aj]=1.0;

	// initialize minimum support and delta
	sigma = 1;
	delta = alpha;
	log_delta = log10(delta);
	while(psi[sigma] > delta){
		++sigma;
	}
	if(sigma > suppMin){
		suppMin = sigma;
		suppMinCurr = sigma;
	}

	// output file
	std::stringstream testing1output_path_result;
  testing1output_path_result << fileinput << "_" << K_significant_patterns << "_" << jp << "_testing.txt";
	string output_path_testing = testing1output_path_result.str();

	// time for p-values
	p_values_time = 0.0;

	// allocate to return array
	to_return = (double*) malloc((jp+1) * sizeof(double));

	quantile_indexes = (bool*) malloc((jp+1) * sizeof(bool));
	for(int aj=0; aj<(jp+1); aj++)
		quantile_indexes[aj]=false;

	start_instant = get_cpu_time();

	precomputetime = 0.0;
	max_cached_pvalues = floor(cache_max_memory / ((double)n1 * 4));

	//cout << "max_cached_pvalues " << max_cached_pvalues << endl;

	cout << "WY functions initialized correctly " << endl;

	cache_accesses = 0;
	new_precomputation = 0;
	current_cache_size = 0;
	current_cache_size_mb = 0;

	pruning_ratio = 1.0 + ((double)n1 / (double)N);

	//benchmarklogoperations();
}

/* -----------------FUNCTIONS TO SAMPLE RANDOM INTEGERS AND GENERATE RANDOM PERMUTATIONS--------------------------- */

/* Sample a random integer uniformly distributed the range [0,x)
 * Default random number generator samples integers uniformly in the range [0,RAND_MAX)
 * To sample in the range [0,x) with x < RAND_MAX, one can think of sampling an integer
 * rnd from [0,RAND_MAX) and then returning rnd % x. However, this may generate a non-uniform
 * distribution unless RAND_MAX is a multiple of x. To fix that, we combine that basic idea
 * with rejection sampling, rejecting any value rnd greater or equal than RAND_MAX - RAND_MAX *x.
 * This is the same as ensuring that the "effective" RAND_MAX is an exact multiple of x, so that
 * returning rnd % x leads to a uniform sampling scheme in the range [0,x)
 * The seed of the random number generator should have been initialised externally
 * */
static int rand_int(int x){
	int rnd;
	int limit = RAND_MAX - RAND_MAX % x;

	do{
		rnd = rand();
	}while(rnd >= limit);
	return rnd % x;
}
void computePermutations(){

	long long swap = 0; // element of random index to be swapped
	long long index = 0;
	long long jp_long = (long long)jp;
	bool tmp = true;

	// Fisher-Yates algorithm - "modern" O(n) version
	// scan every permutation
	for(long long j=1; j < jp_long+1; j++){
		// main O(n) loop
		for(long long i = N-1; i > 0; i--){
			// Sample a random integer in [0,i]
			swap = rand_int(i + 1);
			swap = swap * (jp_long + 1) + j;
			index = i * (jp_long + 1) + j;
			// Swap dest[j] and dest[i]
			tmp = permutations[swap];
			permutations[swap] = permutations[index];
			permutations[index] = tmp;
		}
	}
}

void printPermutations(bool toOutput , bool toFile){

	if(toOutput || toFile){
		ofstream out_file;
		out_file.open ("permutations.txt");
		int index = 0;
		for(int i=0; i < N; i++){
			for(int j=0; j<jp+1; j++){
				index = i * (jp + 1) + j;
				if(toOutput){
					cout << " " << permutations[index];
					if(j==0)
						cout << " | ";
				}
				if(toFile){
					out_file << " " << permutations[index];
					if(j==0)
					out_file << " | ";
				}
			}
			if(toOutput)
				cout << endl;
			if(toFile)
				out_file << endl;
		}
		out_file.close();
	}

}

double sumlogs(double first_log , double second_log){

	double result = 0.0;

	if(log_method == 1){
		result = first_log + log10(1 + pow(10.0,(second_log - first_log)));
	}
	if(log_method == 2){
		result = first_log + (log1p(exp((second_log - first_log) * log_10)) / log_10);
	}
	if(log_method == 3){
		double exponent = second_log - first_log;
		if(abs(exponent) > 2){
			// approximate
			double base;
			if(exponent > 0)
				base = 10;
			else
				base = 0.1;

			result = base;
			int approx_exponent = ceil(abs(exponent));
			for(int i = 0; i < approx_exponent; i++)
				result = result * base;

			result = first_log + log10(1 + result);
		}
		else{
			// compute exact
			result = first_log + (log1p(exp(exponent * log_10)) / log_10);
		}
	}

	return result;

}


/* -------------------------------FUNCTIONS TO COMPUTE FISHER EXACT TEST P-VALUES----------------------------------- */

/* This function precomputes all Fisher exact test P-values for a contingency table with margins x,n,N that is,
 * all p-values p(a,x,n,N) for a in the range [max(0,n+x-N),min(x,n)]. The results will be stored in the array
 * hypergeom_pvals such that p(a,x,n,N)=hypergeom_pvals[a]. Note that values hypergeom_pvals[a] for a outside
 * [max(0,n+x-N),min(x,n)] are undefined and could contain garbage of previous hypotheses.
 * */
void precompute_pvals(int x){

	// perform precomputation only if x is changed from last call
	if(last_support == x){
		return;
	}


	double start = get_cpu_time();

	// remove non-necessary slots
	//cout << "s_supp = " << s_supp << endl;
	//cout << "last_ssupp = " << last_ssupp << endl;
	if(last_ssupp > s_supp){
		for(int aj = s_supp + 1; aj < p_vals_cached.size() && aj <= last_ssupp; aj++){
			if(p_vals_cached[aj].capacity() > 0){
				//cout << "removed cache for index " << aj << endl;
				std::vector<double> foo(0);
				p_vals_cached[aj].swap(foo);
				std::vector<double> foo2(0);
				log_p_vals_cached[aj].swap(foo2);
				current_cache_size -= 2*(min(aj , n1)+1);
				current_cache_size_mb = (int)((double)current_cache_size*sizeof(double)/1000000.0);
			}
		}
		//p_vals_cached.resize(s_supp + 1);
		//log_p_vals_cached.resize(s_supp + 1);
		last_ssupp = s_supp;
	}
	if(suppMin > 1 && last_minsupp < suppMin){
		int index = suppMin-1;
		while(index >= last_minsupp && index < p_vals_cached.size()){
			if(p_vals_cached[index].capacity() > 0){
				//cout << "removed cache for index " << index << endl;
				std::vector<double> foo(0);
				p_vals_cached[index].swap(foo);
				std::vector<double> foo2(0);
				log_p_vals_cached[index].swap(foo2);
				current_cache_size -= 2*(min(index , n1)+1);
				current_cache_size_mb = (int)((double)current_cache_size*sizeof(double)/1000000.0);
			}
			index--;
		}
		last_minsupp = suppMin;
	}

	//cout << "done precomp " << endl;

	bool use_cache = p_val_caching && ( (x-suppMin) < max_cached_pvalues );

	if(use_cache){

		if(p_vals_cached.size() <= x){
			p_vals_cached.resize(x+1);
			log_p_vals_cached.resize(x+1);
		}

		/*cout << "--------------------" << x << endl;
		cout << "precompute for " << x << endl;
		cout << "p_vals_cached.size() = " << p_vals_cached.size() << endl;
		cout << "p_vals_cached["<<x<<"].size() = " << p_vals_cached[x].size() << endl;*/

		if(p_vals_cached[x].size() == 0){
			if(x <= n1){
				p_vals_cached[x].resize(x+1);
				log_p_vals_cached[x].resize(x+1);
				current_cache_size += 2*(x+1);
			}
			else{
				p_vals_cached[x].resize(n1+1);
				log_p_vals_cached[x].resize(n1+1);
				current_cache_size += 2*(n1+1);
			}
			//cout << "allocated cache for index " << x << endl;
			current_cache_size_mb = (int)((double)current_cache_size*sizeof(double)/1000000.0);

			hypergeom_pvals = p_vals_cached[x].data();
			hypergeom_pvals_log = log_p_vals_cached[x].data();
		}
		else{
			//cout << "no resize"<<endl;
			/*for(int aj=0; aj < p_vals_cached[x].size(); aj++){
				hypergeom_pvals[aj] = p_vals_cached[x].at(aj);
				hypergeom_pvals_log[aj] = log_p_vals_cached[x].at(aj);
			}*/
			hypergeom_pvals = p_vals_cached[x].data();
			hypergeom_pvals_log = log_p_vals_cached[x].data();
			last_support = x;
			cache_accesses++;
			precomputetime += get_cpu_time() - start;
			return;
		}

	}
	else{

		hypergeom_pvals = hypergeom_pvals_buffer;
		hypergeom_pvals_log = hypergeom_pvals_log_buffer;

	}

	//cout << "new precomputation for " << x << endl;

	new_precomputation++;

	double pre_comp_xterms, pval, p_left, p_right , p_left_log, p_right_log , pval_log;
	int a, a_min, a_max;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n1+x-N) > 0) ? (n1+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n1) ? n1 : x;//min(x,n)

	// Precompute the hypergeometric PDF in the range of interest
	for(a=a_min; a<=a_max; a++) hypergeom_pvals[a] = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n1-a] + loggamma[x-a] + loggamma[(N-n1)-(x-a)]));

  for(a=a_min; a<=a_max; a++){
	double temp_value = pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n1-a] + loggamma[x-a] + loggamma[(N-n1)-(x-a)]);
	hypergeom_pvals_log[a] = temp_value / log_10;
  }

  /*cout << "hypergeom pdf: " << endl;

	for(int aj = 0; aj < 10 && aj < min(x,n1) ; aj++)
		cout << hypergeom_pvals_log[aj] << " ";
	cout << " ... ";
	for(int aj = min(x,n1)-10; aj < min(x,n1) ; aj++)
		cout << hypergeom_pvals_log[aj] << " ";
	cout << endl;*/

	// The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
	// hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
	// that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. When a value is "accepted", we know its
	// respective p-value because due to the way we explore the hypergeometric pdf, there can be no other values larger than it. Therefore, everytime
	// a value is "accepted", we store the pvalue in hypergeom_pvals.
	// The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
	// that case is by "accepting" both values simultaneously.
	pval = 0;
	while(a_min<a_max){

		p_left = hypergeom_pvals[a_min];
		p_right = hypergeom_pvals[a_max];

		if(p_left == p_right) {
			pval += (p_left+p_right);
			hypergeom_pvals[a_min] = pval;
			hypergeom_pvals[a_max] = pval;

	 		// movement of indices
			a_min++;
			a_max--;

		}
		else
	  if(p_left < p_right){
		pval += p_left;
		hypergeom_pvals[a_min] = pval;

		// move index
		a_min++;
	  }
		  else{
		pval += p_right;
		hypergeom_pvals[a_max] = pval;

		// move index
		a_max--;
	  }

	}
	// In this case a_min=a_max is the mode of the distribution and its p-value is 1 by definition
	if(a_min==a_max){
	hypergeom_pvals[a_max] = 1;
  }

	// log values computation
	a_min = ((n1+x-N) > 0) ? (n1+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n1) ? n1 : x;//min(x,n)
	pval_log = 0;
	while(a_min<a_max){

		p_left_log = hypergeom_pvals_log[a_min];
		p_right_log = hypergeom_pvals_log[a_max];

		if(p_left_log == p_right_log) {

			// log increment
			if(pval_log == 0){
			pval_log = p_left_log + log10_2;
			}
			else{
				double temp = p_left_log + log10_2;
				//pval_log = pval_log + log10(1 + pow(10.0,(temp - pval_log)));
				pval_log = pval_log + (log1p(exp((temp - pval_log) * log_10)) / log_10);
				//pval_log = sumlogs(pval_log , temp);
			}
			hypergeom_pvals_log[a_min] = pval_log;
			hypergeom_pvals_log[a_max] = pval_log;

			// movement of indices
			a_min++;
			a_max--;

		}
		else
		if(p_left_log < p_right_log){

			// log increment
			if(pval_log == 0){
				pval_log = p_left_log;
			}
			else{
				//pval_log = pval_log + log10(1 + pow(10.0,(p_left_log - pval_log)));
				pval_log = pval_log + (log1p(exp((p_left_log - pval_log) * log_10)) / log_10);
				//pval_log = sumlogs(pval_log , p_left_log);
			}

			hypergeom_pvals_log[a_min] = pval_log;

			// move index
			a_min++;
		}
		else{

			// log increment
			if(pval_log == 0){
				pval_log = p_right_log;
			}
			else{
				//pval_log = pval_log + log10(1 + pow(10.0,(p_right_log - pval_log)));
				pval_log = pval_log + (log1p(exp((p_right_log - pval_log) * log_10)) / log_10);
				//pval_log = sumlogs(pval_log , p_right_log);
			}

			hypergeom_pvals_log[a_max] = pval_log;

			// move index
			a_max--;
		}

	}
	// In this case a_min=a_max is the mode of the distribution and its p-value is 1 by definition
	if(a_min==a_max){
	hypergeom_pvals_log[a_max] = 0;
  }

	/*if(p_val_caching && x < max_cached_pvalues){

		for(int aj=0; aj < p_vals_cached[x].size(); aj++){
			p_vals_cached[x].at(aj) = hypergeom_pvals[aj];
		}

	}*/


	/*if(p_val_caching && ( (x-suppMin) < max_cached_pvalues || 2 * current_cache_size_mb < cache_max_memory_mb )){

		for(int aj=0; aj < p_vals_cached[x].size(); aj++){
			p_vals_cached[x].at(aj) = hypergeom_pvals[aj];
			log_p_vals_cached[x].at(aj) = hypergeom_pvals_log[aj];
		}

	}*/

	last_support = x;

	precomputetime += get_cpu_time() - start;


	  /*cout << "hypergeom after cycle: " << endl;

	for(int aj = 0; aj < 10 && aj < min(x,n1) ; aj++)
		cout << hypergeom_pvals_log[aj] << " ";
	cout << " ... ";
	for(int aj = min(x,n1)-10; aj < min(x,n1) ; aj++)
		cout << hypergeom_pvals_log[aj] << " ";
	cout << endl;*/

  /*if(abs(x - n1) < 150){
  cout << "----------------------------" << endl;
  cout << "extreme pvalues for xs = " << x << endl;
  cout << "hypergeom_pvals[0] = " << hypergeom_pvals[0] << " | hypergeom_pvals_log[0] = " << hypergeom_pvals_log[0] << endl;
  cout << "hypergeom_pvals["<<min(x,n1)<<"] = " << hypergeom_pvals[min(x,n1)] << " | hypergeom_pvals_log["<<min(x,n1)<<"] = " << hypergeom_pvals_log[min(x,n1)] << endl;
  /*for(int aj=0; aj < n1+1; aj++){
	cout << "hypergeom_pvals["<<aj<<"] = " << hypergeom_pvals[aj] << " | hypergeom_pvals_log["<<aj<<"] = " << hypergeom_pvals_log[aj] << endl;
  }*//*
  }*/

	//p_values_time += get_cpu_time() - start;
}

double computePValue(int x_S , int a_S){

	precompute_pvals(x_S);

	return hypergeom_pvals[a_S];

}

double computeLogPValue(int x_S , int a_S){

	precompute_pvals(x_S);

	if(SOFT_DEBUG){
	cout << "processed itemset with support " << x_S  << " and a_S " << a_S << endl;

	cout << "its log pvalue is " << hypergeom_pvals_log[a_S] << endl;

	for(int aj = 0; aj < 5 && aj < a_S ; aj++)
		cout << hypergeom_pvals_log[aj] << " ";
	cout << " ... ";
	for(int aj = a_S-5; aj < a_S ; aj++)
		cout << hypergeom_pvals_log[aj] << " ";
	cout << endl;
	}

	return hypergeom_pvals_log[a_S];

}


double* computePValues(int x_S , int* a_S){

	//double start = get_cpu_time();

	/*bool updated_min_pvals = false;*/

	precompute_pvals(x_S);

	double s_inst = get_cpu_time();

	for(int aj=0; aj<jp+1; aj++){
		//cout << "aj " << aj<<endl;
		to_return[aj]=hypergeom_pvals[a_S[aj]];
		if(to_return[aj] < min_pval[aj]){
			min_pval[aj] = to_return[aj];
		}
		if(aj > 0 && to_return[aj] <= delta){
			if(quantile_indexes[aj] == false){
				temp_alpha_item.p_value = to_return[aj];
				temp_alpha_item.index = aj;
				alpha_quantile.push(temp_alpha_item);
				/*updated_min_pvals = true;
				cout << "improved min p-value of permutation with index " << aj << endl;
				cout << "x_S = " << x_S << " a_S[aj] = " << a_S[aj] << endl;*/
				below_delta++;
				quantile_indexes[aj] = true;
				/*below_delta_indexes.insert(aj);*/
			}
		}
	}
	/*
	if(updated_min_pvals){
		// print distribution of as
		std::vector<int> top5(0);
		std::vector<int> low5(0);
		int temp = 0;
		int temp_index = 0;
		for(int aj=0; aj<jp+1; aj++){
			// don't consider as of lower alpha quantile
			if(below_delta_indexes.find(aj) == below_delta_indexes.end()){
			// top operations
			if(top5.size() < 5){
				top5.push_back(a_S[aj]);
				if(top5.size() > 2){
				// find new min in top5 and swap it
				temp = top5[0];
				temp_index = 0;
				for(int h=1; h<top5.size(); h++){
					if(top5[h] < temp){
						temp_index = h;
						temp = top5[h];
					}
				}
				// swap new element with the new min
				top5[temp_index] = top5[0];
				top5[0] = temp;
			}
			}
			else{
				if(a_S[aj] > top5[0]){
					//swap min with new value
					top5[0] = a_S[aj];
					temp_index = 0;
				// find new min in top5 and swap it
				temp = top5[0];
				for(int h=1; h<top5.size(); h++){
					if(top5[h] < temp){
						temp_index = h;
						temp = top5[h];
					}
				}
				// swap new element with the new min
				top5[temp_index] = top5[0];
				top5[0] = temp;
			}
			}
			// low operations
			if(low5.size() < 5){
				low5.push_back(a_S[aj]);
				if(low5.size() > 2){
				// find new max in low5 and swap it
				temp = low5[0];
				temp_index = 0;
				for(int h=1; h<low5.size(); h++){
					if(low5[h] > temp){
						temp_index = h;
						temp = low5[h];
					}
				}
				// swap new element with the new max
				low5[temp_index] = low5[0];
				low5[0] = temp;
				}
			}
			else{
				if(a_S[aj] < low5[0]){
					//swap max with new value
					low5[0] = a_S[aj];
				// find new max in low5 and swap it
				temp = low5[0];
				temp_index = 0;
				for(int h=1; h<low5.size(); h++){
					if(low5[h] > temp){
						temp_index = h;
						temp = low5[h];
					}
				}
				// swap new element with the new max
				low5[temp_index] = low5[0];
				low5[0] = temp;
			}
			}
		}
		}

		// print the vectors
		for(int h=0; h<low5.size(); h++){
			cout << low5[h] << " , ";
		}
		cout << " , ... " ;
		for(int h=0; h<top5.size(); h++){
			cout << " , " << top5[h];
		}
		cout << endl;

	}*/

	if((double)below_delta / (double)jp > alpha){

		/*cout << "------------------------" << endl;
		cout << "Increasing delta to control FWER " << endl;
		cout << "alpha_quantile.size() (0) " << alpha_quantile.size() << endl;
		cout << "below_delta " << below_delta << endl;*/

		/* delta needs to be lowered since FWER > 0.05 */
		while((double)below_delta / (double)jp > alpha || alpha_quantile.top().p_value == delta){

			/*cout << "alpha_quantile.size() (1) " << alpha_quantile.size() << endl;
			cout << "below_delta " << below_delta << endl;*/

			// remove all elements of alpha_quantile with outdated p-value

			while(alpha_quantile.top().p_value > min_pval[ alpha_quantile.top().index ]){
				/*cout << " removed outdated p-value at index " << alpha_quantile.top().index << endl;
				cout << " alpha_quantile.top().p_value " << alpha_quantile.top().p_value << endl;
				cout << " min_pval[ alpha_quantile.top().index ] " << min_pval[ alpha_quantile.top().index ] << endl;*/
				temp_alpha_item = alpha_quantile.top();
				temp_alpha_item.p_value = min_pval[ temp_alpha_item.index ];
				alpha_quantile.pop();
				alpha_quantile.push(temp_alpha_item);
			}

			delta = min(delta , alpha_quantile.top().p_value);
			alpha_quantile.pop();
			--below_delta;

		}
		//cout << "below_delta reduced to " << below_delta << endl;

		while(alpha_quantile.top().p_value > min_pval[ alpha_quantile.top().index ])
			alpha_quantile.pop();

			delta = min(delta , alpha_quantile.top().p_value);
			log_delta = min( (delta > 0.0) ? log10(delta) : 0.0 , log_delta);

			//cout << "delta (after delta lift)" << delta << endl;
			/* update minsupp if needed */
			while(psi[sigma] > delta){
				++sigma;
			}

			below_delta = 0;
			/*below_delta_indexes.clear();*/

			for(int aj=1; aj<jp+1; aj++){
					quantile_indexes[aj] = min_pval[aj] <= delta;
					below_delta+=quantile_indexes[aj];
					/*below_delta_indexes.insert(aj);*/
			}



			/*cout << "below_delta is now " << below_delta << endl;
			cout << "FWER " << ((double)below_delta / (double)jp) << endl;
			cout << "increased sigma to " << sigma << endl;
			cout << "delta is now " << delta << endl;
			cout << "psi[sigma] " << psi[sigma] << endl;//}
			cout << "alpha_quantile.size() " << alpha_quantile.size() << endl;*/


		if(sigma > suppMin){
			cout << "WY: updated suppMin to " << sigma << endl;
			suppMinCurr = sigma;
			suppMin = sigma;
			//cout << "  below_delta is now " << below_delta << endl;
			cout << "  FWER " << ((double)below_delta / (double)jp) << endl;
			cout << "  delta " << delta << endl;
			cout << "  psi[sigma] " << psi[sigma] << endl;//}
			//cout << "alpha_quantile.size() " << alpha_quantile.size() << endl;
		}

		//cout << "------------------------" << endl;

	}

	time_to_update_delta += get_cpu_time() - s_inst;

	//p_values_time += get_cpu_time() - start;

	return to_return;

}

double getPsi(int x_S){

	return psi[x_S];

}

double getLogPsi(int x_S){

	return psi_log[x_S];

}

double getPsiExact(int x_S){

	return psi_exact[x_S];

}

double getLogPsiExact(int x_S){

	return psi_log_exact[x_S];

}

double getExactPsiBound(int x_S , int x_S_anchestor , int a_S_anchestor){

	if(x_S > N){
		cout << "problem in getExactPsiBound! "<< x_S << " " << x_S_anchestor << " " << a_S_anchestor <<endl;
		abort();
	}

	double to_return_ = 0.0;
	int max_a_S = x_S;
	int min_a_S = 0;

	if(a_S_anchestor < max_a_S)
		max_a_S = a_S_anchestor;

	if(x_S - (x_S_anchestor - a_S_anchestor) > min_a_S)
		min_a_S = x_S - (x_S_anchestor - a_S_anchestor);

	// perform computation of p-values only if the bound will be inferior to standard psi value
	if(max_a_S < x_S || min_a_S > 0){

		precompute_pvals(x_S);

		if(hypergeom_pvals[min_a_S] < hypergeom_pvals[max_a_S])
			to_return_ = hypergeom_pvals[min_a_S];
		else
			to_return_ = hypergeom_pvals[max_a_S];

	}
	else{
		to_return_ = psi[x_S];
	}

	//cout << " exact psi bound = " << to_return_ << endl;
	#ifdef LAMP

	if(itemsets_LAMP_count_naive.size()<=x_S)
		itemsets_LAMP_count_naive.resize(x_S+1);
	itemsets_LAMP_count_naive[x_S]++;


	int aj=1;
	while((psi[aj] >= to_return_ && to_return_ > 0.0) || (psi[aj] > to_return_)){
		aj++;
	}
	if(itemsets_LAMP_count.size() <= aj)
		itemsets_LAMP_count.resize(aj+1);
	itemsets_LAMP_count[aj]++;

	#endif
	//cout << " itemset of support " << x_S << " has hyp of psi["<< (aj-1) <<"], count " << itemsets_LAMP_count[aj-1] <<endl;

	return to_return_;

}

double getPsiBound(int x_S , int x_S_anchestor , int min_a_S_anchestor , int max_a_S_anchestor){

	double to_return_ = 0.0;
	int max_a_S = min(x_S , max_a_S_anchestor);
	int min_a_S = max(0 , min_a_S_anchestor - (x_S_anchestor - x_S));

	// perform computation of p-values only if the bound will be inferior to standard psi value
	if(max_a_S < x_S || min_a_S > 0){

		precompute_pvals(x_S);

		if(hypergeom_pvals[min_a_S] < hypergeom_pvals[max_a_S])
			to_return_ = hypergeom_pvals[min_a_S];
		else
			to_return_ = hypergeom_pvals[max_a_S];

	}
	else{
		to_return_ = psi_exact[x_S];
	}

	return to_return_;

}

double getLogPsiBound(int x_S , int x_S_anchestor , int min_a_S_anchestor , int max_a_S_anchestor){

	double to_return_ = 0.0;
	int max_a_S = min(x_S , max_a_S_anchestor);
	int min_a_S = max(0 , min_a_S_anchestor - (x_S_anchestor - x_S));

	// perform computation of p-values only if the bound will be inferior to standard psi value
	if(max_a_S < x_S || min_a_S > 0){

		precompute_pvals(x_S);

		/*if(hypergeom_pvals[min_a_S] < hypergeom_pvals[max_a_S])
			to_return_ = hypergeom_pvals[min_a_S];
		else
			to_return_ = hypergeom_pvals[max_a_S];*/
		to_return_ = min(hypergeom_pvals_log[min_a_S] , hypergeom_pvals_log[max_a_S]);


		/*cout << "getLogPsiBound: x_S_anchestor =  " << x_S_anchestor <<
		" min_a_S_anchestor = " << min_a_S_anchestor <<  " max_a_S_anchestor = " << max_a_S_anchestor << endl;
		cout << "getLogPsiBound: x_S = " << x_S << " min_a_S =  " << min_a_S <<
		" max_a_S = " << max_a_S <<  " max_a_S_anchestor = " << max_a_S_anchestor << endl;
		cout << "PsiBound: to_return_ =  " << to_return_ << endl;*/

	}
	else{
		to_return_ = psi_log_exact[x_S];
	}

	return to_return_;

}

int canHaveFPChild(int x_S , int min_a_S , int max_a_S){

	if(suppMin <= max_a_S || min((int)(pruning_ratio * suppMin) , suppMin + 100) <= x_S || p_vals_cached.size() <= x_S){
		return 1;
	}

	// perform a first test to check most extreme value
	int x_X = suppMin;
	int max_a_S_child = min(x_X , max_a_S);
	int min_a_S_child = max(0 , min_a_S - (x_S - x_X));
	double p_val = 0.0;
	double psi_curr = getLogPsi(suppMin);
	if(min_a_S_child != 0 || max_a_S_child != x_X){
		precompute_pvals(x_X);
		p_val = min(hypergeom_pvals_log[min_a_S_child] , hypergeom_pvals_log[max_a_S_child]);
		if(p_val < getLogPsi(suppMin)){
				//cout << "it can be significant with just one precomputation" << endl;
				return x_X;
		}
	}
	else {
		//cout << "it can be significant with NO precomputation" << endl;
		return x_X;
	}


	/*if(!p_val_caching)
		return 1;*/

	x_X = x_S - 1;

	//cout << " processing pattern with x_S = " << x_S << " and min_a_S = " << min_a_S << " and max_a_S = " << max_a_S << endl;

	if(x_X < suppMin){
		//cout << "halt here, no children can be a FP" << endl;
		return -1;
	}

	while(true){

		precompute_pvals(x_X);


		max_a_S_child = min(x_X , max_a_S);
		min_a_S_child = max(0 , min_a_S - (x_S - x_X));

		//cout << " processing sub with x_X = " << x_X << endl;
		//cout << " computed min_a_S_child = " << min_a_S_child << " and max_a_S_child = " << max_a_S_child << endl;

		p_val = min(hypergeom_pvals_log[min_a_S_child] , hypergeom_pvals_log[max_a_S_child]);

		//cout << " min attainable pval = " << p_val << " and current delta = " << psi_curr << endl;

		if(p_val < psi_curr){
			//cout << " min attainable pval = " << p_val << " and current delta = " << psi_curr << endl;
			//cout << "it can be FP, halting here " << endl;
			return x_X;
		}
		else{
			//cout << "it can not be FP, let's continue" << endl;
			x_X--;
			if(x_X < suppMin){
				//cout << " min attainable pval = " << p_val << " and current delta = " << psi_curr << endl;
				//cout << " pruned pattern with x_S = " << x_S << " and min_a_S = " << min_a_S << " and max_a_S = " << max_a_S << endl;
				//cout << "halt here, no children can be FP" << endl;
				return -1;
			}
		}
	}
}

int canHaveSignificantChild(int x_S , int a_S){

	/*if(!p_val_caching)
		return 1;*/

	if(suppMin <= a_S || min((int)(pruning_ratio * suppMin) , suppMin + 100) <= x_S || p_vals_cached.size() <= x_S){
		return 1;
	}

	// perform a first test to check most extreme value
	int x_X = suppMin;
	int max_a_S = min(x_X , a_S);
	int min_a_S = max(0 , a_S - (x_S - x_X));
	double p_val = 0.0;
	double psi_curr = getLogPsi(suppMin);
	if(min_a_S != 0 || max_a_S != x_X){
		precompute_pvals(x_X);
		p_val = min(hypergeom_pvals_log[min_a_S] , hypergeom_pvals_log[max_a_S]);
		if(p_val < getLogPsi(suppMin)){
				//cout << "can be significant with just one precomputation" << endl;
				return x_X;
		}
	}
	else {
		//cout << "can be significant with NO precomputation" << endl;
		return x_X;
	}

	x_X = x_S - 1;

	max_a_S = a_S;
	min_a_S = 0;
	p_val = 0.0;
	psi_curr = getPsi(suppMin);

	//cout << " processing pattern with x_S = " << x_S << " and a_S = " << a_S << endl;

	if(x_X < suppMin){
		//cout << "halt here, no children can be significant" << endl;
		return -1;
	}

	while(true){

		precompute_pvals(x_X);

		max_a_S = min(x_X , a_S);
		min_a_S = max(0 , a_S - (x_S - x_X));

		//cout << " processing sub with x_X = " << x_X << endl;
		//cout << " computed min_a_S = " << min_a_S << " and min_a_S = " << max_a_S << endl;

		p_val = min(hypergeom_pvals_log[min_a_S] , hypergeom_pvals_log[max_a_S]);

		//cout << " min attainable pval = " << p_val << " and current delta = " << psi_curr << endl;

		if(p_val < psi_curr){
			//cout << "it can be significant, halting here" << endl;
			return x_X;
		}
		else{
			//cout << "it can not be significant, let's continue" << endl;
			x_X--;
			if(x_X < suppMin){
				//cout << " pruned pattern with x_S = " << x_S << " and a_S = " << a_S << endl;
				//cout << "halt here, no children can be significant" << endl;
				return -1;
			}
		}
	}

}

void printMinPvalues(){
	for(int aj=0; aj<jp+1; aj++)
		cout << min_pval[aj] << " ";
	cout << endl;
}

void printTimes(){
	cout << "Time for p-values precomputation = " << precomputetime << endl;
	cout << "p-values cache accesses = " << cache_accesses << endl;
	cout << "p-values new precomputations = " << new_precomputation << endl;
	cout << "Time to update delta = " << time_to_update_delta << endl;
}

void benchmarklogoperations(){

	cout << "Benchmark for log operations begins" << endl;

	int max_value = min(N , 10000);

	log_method = 1;

	double start_1 = get_cpu_time();

	for(int aj = 1; aj < max_value ; aj++)
		precompute_pvals(aj);

	start_1 = get_cpu_time() - start_1;
	cout << "time for first method: " << start_1 << endl;

	log_method = 2;

	double start_2 = get_cpu_time();

	for(int aj = 1; aj < max_value ; aj++)
		precompute_pvals(aj);

	start_2 = get_cpu_time() - start_2;
	cout << "time for second method: " << start_2 << endl;

	log_method = 3;

	double start_3 = get_cpu_time();

	for(int aj = 1; aj < max_value ; aj++)
		precompute_pvals(aj);

	start_3 = get_cpu_time() - start_3;
	cout << "time for third method: " << start_3 << endl;


	log_method = 0;

	double start_0 = get_cpu_time();

	for(int aj = 1; aj < max_value ; aj++)
		precompute_pvals(aj);

	start_0 = get_cpu_time() - start_0;

	cout << "time for no method: " << start_0 << endl;


}

double computeEmpiricalFWER(){
	below_delta = 0;
	/*below_delta_indexes.clear();*/

	for(int aj=1; aj<jp+1; aj++){
			quantile_indexes[aj] = min_pval[aj] <= delta;
			below_delta+=quantile_indexes[aj];
			/*below_delta_indexes.insert(aj);*/
	}

	return (double)below_delta / (double)jp;
}



#endif
