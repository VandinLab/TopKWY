/*
 * Copyright (C) 2017 Leonardo Pellegrina, Fabio Vandin
 * This file is part of TopKWY and it is based on TopKMiner
 * Copyright (C) 2008 Andrea Pietracaprina, Fabio Vandin
 * Copyright (C) 2008 Advanced Computing Group, University of Padova, Italy
 *
 * TopKWY is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#include "queueheap.h"
#include "wy.cc"
#include <queue>




// test: novel data structure for results ordering

void QueueMinMax_results_std::insert(ResultPattern* res){
	#ifdef trace_state
	queue_program_state = 9;
	#endif
	results_p_values.push(res);

}

void QueueMinMax_results_std::insert_observed(double log_observed_pval){
	#ifdef trace_state
	queue_program_state = 91;
	#endif
	observed_p_values.push(log_observed_pval);
}

// removal of non-significant patterns
void QueueMinMax_results_std::removeMax(double log_delta){

	#ifdef trace_state
	queue_program_state = 10;
	#endif

	if(results_p_values.size() == 0)
		return;

	ResultPattern* removed = results_p_values.top();
	free(removed->itemset);
	free(removed);
	results_p_values.pop();

	#ifdef trace_state
	queue_program_state = 0;
	#endif

}

// output a list of significant patterns
std::vector< ResultPattern* > QueueMinMax_results_std::removeMin(double log_delta){

	#ifdef trace_state
	queue_program_state = 11;
	#endif

	std::vector< ResultPattern* > toreturn;
	if(q_min_pvalue > log_delta || results_p_values.size() == 0){
		return toreturn;
	}

	std::priority_queue< ResultPattern*, std::vector<ResultPattern*>, Compare_Results_toplarge > new_results_p_values;
	q_min_pvalue = 1.0;
	while(results_p_values.size() > 0){
		ResultPattern* removed = results_p_values.top();
		results_p_values.pop();
		if(removed->log_p_value <= log_delta){
			toreturn.push_back(removed);
		}
		else{
			new_results_p_values.push(removed);
			if(removed->log_p_value < q_min_pvalue){
				q_min_pvalue = removed->log_p_value;
			}
		}
	}
	results_p_values = new_results_p_values;

	#ifdef trace_state
	queue_program_state = 0;
	#endif

	return toreturn;
}


/*double QueueMinMax_results_std::getMin(){
	return q_min_pvalue;
}*/


double QueueMinMax_results_std::getLogMin(){
	return q_min_pvalue;
}

/*double QueueMinMax_results_std::getMax(){
	return results_p_values.top();
}*/

int QueueMinMax_results_std::getLogMax(){
	return results_p_values.top()->log_p_value;
}

tipoInt QueueMinMax_results_std::getInqueue(){
	return results_p_values.size();
}

tipoInt QueueMinMax_results_std::getKElements(){
	return observed_p_values.size();
}

void QueueMinMax_results_std::printElements(){
	cout << "printElements call!" << endl;
	/* todo */
	cout << "printElements call end!" << endl;
}

QueueMinMax_results_std::QueueMinMax_results_std(){

	q_min_pvalue = 1.0;

}
