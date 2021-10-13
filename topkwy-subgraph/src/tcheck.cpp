#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <boost/math/distributions/hypergeometric.hpp>
#include "legoccurrence.h"
#include "closeleg.h"
#include "database.h"
#include "graphstate.h"
#include "tcheck.h"
#include "wy.h"

// read class labels
void readClass(char *class_file) {
  ifstream ifs(class_file);
  if (ifs.fail()) {
    cerr << "File do not exist." << endl;
    exit(0);
  }

  int c;
  string str;
  while (getline(ifs, str)) {
    sscanf(str.data(), "%d", &c);
    CLASS_VEC.push_back(c);
    N_TOTAL++;
    if (c == 1) N_SMALL++;
    else if (c != 0) {
      cerr << "class labels should be 0 or 1" << endl;
      exit(0);
    }
  }

  if ((double)N_SMALL > (double)N_TOTAL / 2.0) {
    N_SMALL = N_TOTAL - N_SMALL;

    for (vector<int>::iterator itr = CLASS_VEC.begin(), end = CLASS_VEC.end(); itr != end; ++itr) {
      *itr == 1 ? *itr = 0 : *itr = 1;
    }
  }

  ifs.close();
}

void outputSubgraph(int count, Frequency frequency, double p_value) {
  OFS << "# s " << count << " " << frequency - count << endl;
  OFS << "# p " << p_value << endl;
  graphstate.print_ofs();
}



// measure cpu time
double get_cpu_time(){
	return (double)clock() / CLOCKS_PER_SEC;
}

// appear in "patterntree.cpp" (2 times), "path.cpp" (2 times), "patterngraph.cpp"
// data types are defined in "legoccurrence.h"
void checkTestable(vector<LegOccurrence>& elements, Frequency frequency, bounds_info& parent_info, bounds_info& current_info) {
  #ifdef deepdebug
  cout << "call to checkTestable with support " << (int)frequency << endl;
  #endif
  if (!RERUN) {
    if (WYLIGHT) {
      double temp_time = get_cpu_time();
      checkConditionWY<LegOccurrence>(elements, frequency, parent_info, current_info);
      wy_processing_time += get_cpu_time() - temp_time;
    }
    else checkCondition();
  }
  else computePvalue<LegOccurrence>(elements, frequency);
}

void checkTestableCl(vector<CloseLegOccurrence>& elements, Frequency frequency, bounds_info& parent_info, bounds_info& current_info) {
  #ifdef deepdebug
  cout << "call to checkTestableCl with support " << (int)frequency << endl;
  #endif
  if (!RERUN) {
    if (WYLIGHT) {
        double temp_time = get_cpu_time();
        checkConditionWY<CloseLegOccurrence>(elements, frequency, parent_info, current_info);
        wy_processing_time += get_cpu_time() - temp_time;
    }
    else checkCondition();
  }
  else computePvalue<CloseLegOccurrence>(elements, frequency);
}

void checkCondition() {
  Tid id_prev = NOTID;

  COUNT += 1.0;
  if (COUNT > THRESHOLD) {
    minfreq++;
    cerr << "<<< Update minfreq to " << minfreq << " >>>" << endl;
    // compute the minimum achievable P value, which is "binom(N_SMALL, minfreq) / binom(N_TOTAL, minfreq)"
    double minp = 1;
    for (int i = 0; i < minfreq; i++) {
      minp *= (double)(N_SMALL - i) / (double)(N_TOTAL - i);
    }
    // update threshold
    THRESHOLD = ALPHA / minp;
    cerr << "    Admissible number of subgraphs: " << floor(THRESHOLD) << endl;
  }
}



template<typename T> void checkConditionWY(vector<T>& elements, Frequency frequency, bounds_info& parent_info, bounds_info& current_info) {
  Tid id_prev = NOTID, i = 0;
  typename vector<T>::iterator itr, end;

  #ifdef bounds_debug
  cout << "call to checkConditionWY" << endl;
  cout << "_____________________________ " << endl;
  cout << "pattern of support " << (int)frequency << endl;
  #endif

  for (itr = elements.begin(), end = elements.end(); itr != end; ++itr) {
    if (id_prev != itr->tid) {
      // retrieving occurrences in graphs
      OCC_VEC[i] = itr->tid;
      #ifdef bounds_debug
      cout << "trans " << i << " " << OCC_VEC[i] << endl;
      #endif
      id_prev = itr->tid;
      i++;
    }
  }

  #ifdef bounds_debug
  cout << "_____________________________ " << endl;
  #endif

  bm_process_solution((int)frequency, parent_info, current_info);
  // OUTPUT(frequency);
}

template<typename T> void computePvalue(vector<T>& elements, Frequency frequency) {
  typename vector<T>::iterator itr, end;

  Tid id_prev = NOTID;
  int count = 0;
  for (itr = elements.begin(), end = elements.end(); itr != end; ++itr) {
    if (id_prev != itr->tid) {
      // retrieving occurrences in graphs
      if (CLASS_VEC[itr->tid] == 1) count++;
      id_prev = itr->tid;
    }
  }

  /* this seems to be inconsistent with p-values from WY code */

  // compute p-value from "count", "frequency", "N_SMALL", and "N_TOTAL"
  /*boost::math::hypergeometric_distribution<> phyper(frequency, N_SMALL, N_TOTAL);
  double p_value_L = 0, p_value_R = 0, p_value_D = 0;
  if (count == 0) count = 1;
  else if (count == N_SMALL) count = N_SMALL - 1;
  if ((int)frequency - count == 0) count = frequency - 1;
  else if ((int)frequency - count == N_TOTAL - N_SMALL) count = (int)frequency - N_TOTAL + N_SMALL + 1;
  if (0 < count && ((int)frequency + N_SMALL - N_TOTAL) < count && count < N_SMALL && count < (int)frequency) {
    p_value_L = cdf(phyper, count);
    p_value_R = cdf(complement(phyper, count - 1));
    p_value_D = p_value_L < p_value_R ? 2 * p_value_L : 2 * p_value_R;
  } else {
    p_value_D = 1;
  }*/

  /*if (p_value_D <= DELTA) {
    COUNT += 1.0;
    //cout << "rejected p-value: " << p_value_D << endl;
    outputSubgraph(count, frequency, p_value_D);
  }*/

  precompute_pvals(frequency);
  double pvalue = getPValue(count);

  if (pvalue <= DELTA) {
    COUNT += 1.0;
    #ifdef debug_results
    cout << "rejected p-value: " << pvalue << endl;
    cout << "frequency: " << (int)frequency << endl;
    cout << "count: " << (int)count << endl;
    #endif
    outputSubgraph(count, frequency, pvalue);
  }


}
