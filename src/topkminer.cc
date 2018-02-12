/*
 * Copyright (C) 2017 Leonardo Pellegrina, Fabio Vandin
 * This file is part of SignificantMiner and it is based on TopKMiner
 * Copyright (C) 2008 Andrea Pietracaprina, Fabio Vandin
 * Copyright (C) 2008 Advanced Computing Group, University of Padova, Italy
 *
 * SignificantMiner is free software; you can redistribute it and/or
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

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "patricia.cc"
#include "wy_.cc"
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

#define initial_to_del_size 128
#define debug_todel false
#define maxNLRAM max_ram//256000 // in Mbyte
#define maxNLLen max_nl_size//120
#define midNLLen 60
#define inital_temp_size 40
#define max_uint_8 256
#define max_uint_16 65536

#define trace_state 1
#define use_secondary_memory false

#define write_results_to_file true
#define enable_pruning true
#define enable_bounds true
//#define debug_memory_usage true


//double time_keeper;
double temp_time;
double as_processing_time;
std::vector<tipoInt*>* todelete;
int nextfree_a_S;
int nextfree_nl_list;
std::vector< std::vector<tipoInt> > nl_lists;
int allocated;
double memory_peak;
double saved;
int maxlivinginfos;
int maxinqueue;
PatriciaTrie* p;
bool resized;
tipoInt* temp_nl_a_S;
uint8_t* temp_nl_a_S_8;
uint16_t* temp_nl_a_S_16;
uint32_t* temp_nl_a_S_32;
tipoInt temp_size;
tipoInt temp_size_8;
tipoInt temp_size_16;
tipoInt temp_size_32;
//tipoInt s_supp;
char *file_path;
double in_use_ram;
int memory_baseline;
double queue_ram;
ofstream out_file_stat_2;
int inferior_psi_bound_happened;
int inferior_psi_bound_happened_neglect;
int program_state;
int debug_info_1;
int debug_info_2;
int max_nl_size_;

double time_pat;

bool store_results;

long list_length_inqueue;
long nl_nodes_length_inqueue;

std::vector<int> list_stat;
std::vector<int> list_op_stat;
int compressed_nl_lists;

int pruned_patterns;

bool memory_bounded;
long memory_bounded_max_inqueue;
long memory_bounded_max_nl_length;
long memory_bounded_max_nl_min_length;

std::vector<int> perm_supp;

typedef QueueMinMax_test* typeQueue;
typedef QueueMinMax_Restest* typeQueueRes;
//std::vector<double> timekeepers;

/**raise minimum support using util node heuristic
*/

// Measure peak memory usage
double measurePeakMemory(){
  struct rusage t;
  getrusage(RUSAGE_SELF, &t);
  return ((double)((size_t)t.ru_maxrss))/1000.0;
}

tipoInt utilNodeHeuristic(){
	tipoInt* forHeur=(tipoInt*)calloc(maximumSupport+1, sizeof(tipoInt));
	PatriciaNode* patNode;
	for (int i=0; i<numNodes; i++){
		patNode=addressOf[i];
		if (patNode->lastItem==patNode->lastInter){
			forHeur[patNode->count]++;
		}
	}
	tipoInt cumul=0;
	for (int j=maximumSupport+1; j>=0; j--){
		cumul=cumul+forHeur[j];
		if (cumul>=k_max){
			return j;
		}
	}
	return 1;
}

void notFreqInPat(tipoInt support, ItemList* IL){
	tipoInt notFreqTot=0; ///total number of items that are not frequent in Patricia
	tipoInt notFreqInter=0; ///items that are not frequent in intersections of Patricia
	PatriciaNode* patNode;
	for (int i=0; i<numNodes; i++){
		patNode=addressOf[i];
		tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));
		for (int j=0; j<numItems; j++){
			if (IL[(&patNode->item)[j]].count<support){
				notFreqTot++;
			}
		}
		tipoInt numInter=(tipoInt) (patNode->lastInter-patNode->lastItem);
		for (int l=0; l<numInter; l++){
			if (IL[(&patNode->item)[numItems+l]].count<support){
				notFreqTot++;
				notFreqInter++;
			}
		}
	}
}

void toDelFromVisited(HeaderTable* HT, ItemList* IL){
	tipoInt overhead=0;
	for (int i=num_item-1;i>=0; i--){
		if (IL[i].count<suppMin&&HT[i].visitedList!=NULL){
			overhead=overhead+(numNPI[i]+1);
		}
	}
}

/**raise minimum support using the trivial heuristic ; used in computation
NB: it requires the explicit update of utilNodeArray
*/

void trivialHeur(){
	return;
	///trivial heuristic on suppMinCurr
	for (tipoInt i=suppMinCurr; i<=maximumSupport; i++){
		if (utilNodeArray[i]>=k){
			suppMinCurr=i;
		}
	}
	///trivial heuristic on suppMin
	for (tipoInt i=suppMin; i<=maximumSupport; i++){
		if (utilNodeArray[i]>=k_max){
			suppMin=i;
		}
	}
}

void initListIndex(){
	nl_lists.resize(64);
}

tipoInt getNewListIndex(){

	nextfree_nl_list++;
	if(nextfree_nl_list >= nl_lists.size()){
		nl_lists.resize(nextfree_nl_list*2);
	}
	return (nextfree_nl_list-1);

}

void resetListIndex(){

	for(int i=0; i < nextfree_nl_list; i++){
		nl_lists[i].clear();
	}
	nextfree_nl_list=0;

}

void resizeToDel(){
	if(debug_todel)
	cout << "resize of todelete " << endl;
	resized = true;
	//int oldsize = todelete->size();
	todelete->resize(todelete->size()*2);
	/*for(int aj=oldsize; aj<todelete->size(); aj++){
		todelete->at(aj)=(tipoInt*)malloc((jp+1) * sizeof(tipoInt));
		if(debug_todel)
		cout << "new memory " << aj << " allocated at " << todelete->at(aj) << endl;
	}*/
}

tipoInt getNewTodeleteIndex(){

	if(debug_todel)
	cout << "getting new a_S index..." << endl;
	if(debug_todel)
	cout << "nextfree_a_S = " << nextfree_a_S << "/" << todelete->size() << endl;

	if(nextfree_a_S >= todelete->size())
		resizeToDel();
	if(nextfree_a_S >= allocated){
		todelete->at(allocated)=(tipoInt*)malloc((jp+1) * sizeof(tipoInt));
		if(debug_todel)
		cout << "new memory " << allocated << " allocated at " << todelete->at(allocated) << endl;
		allocated++;
	}

	if(debug_todel)
	cout << "assigned " << todelete->at(nextfree_a_S) << endl;
	nextfree_a_S++;
	return nextfree_a_S-1;

}

void resetToDel(){

	if(debug_todel)
	cout << "resetting todel..." << endl;
	nextfree_a_S=0;
	/*if(resized){
		for(int aj=0; aj<todelete->size()/2; aj++){
			free(todelete->at(aj));
			todelete->at(aj)=(tipoInt*)malloc((jp+1) * sizeof(tipoInt));
			cout << "new memory " << aj << " allocated at " << todelete->at(aj) << endl;
		}
		resized = false;
	}*/
	if(debug_todel)
	cout << "resetted nextfree_a_S" << endl;
}

void freeToDel(){
	for(int aj=0; aj<todelete->size(); aj++){
		if(debug_todel)
		cout << "freed " << aj << " of todel allocated at " << todelete->at(aj) << endl;
		free(todelete->at(aj));
	}
	free(todelete);
}

void allocateToDel(){
	todelete = new std::vector<tipoInt*>();
	todelete->resize(initial_to_del_size);
	for(int aj=0; aj < todelete->size(); aj++){
		todelete->at(aj)=(tipoInt*)malloc((jp+1)*sizeof(tipoInt));
		if(debug_todel)
		cout << "allocated space " << aj << " at " << todelete->at(aj) << endl;
	}
	allocated = todelete->size()-1;
	resized = false;
}

/** find the ppce of Clo(emptyset) passing the IL and the queue;
 */

//void findExtCloEmpty(ItemList* IL, QueueMinMax* q){
void findExtCloEmpty(ItemList* IL, typeQueue q , typeQueueRes q_res){

  //cout << "perm_supp.size() = " << perm_supp.size() << endl;
  if(perm_supp.size() < jp+1)
    perm_supp.resize(jp+1);

	///now we insert in the heap the ppc-e that are frequent (at the moment); we can begin from i=0 and go down the HT because we now that the items in Clo(emptyset) are eliminated from the dataset (and added at the end);
	Itemset* tmpIts;	///temporary itemset used to insert the itemset in the max queue
	Itemset* minIts;	///temporary itemset used to insert the itemset in the min queue
	ItsInfo* tmpInfo;	///temporary info used to memorize the infos for Itemset tmpIts and link them
	NodeIDElem* prev;
	for (tipoInt i=0; i<num_item; i++){
		/*cout << "item i " << i << endl;
		cout << "IL[i].count = " << IL[i].count << endl;
		cout << "here ok -2" << endl;
    cout << "suppMin " << suppMin << endl;
    cout << "HT[i].intersection[0] = " << HT[i].intersection[0] << endl;*/
		///test if the itemset is frequent and if it is prefix-preserving
		if ((IL[i].count>=suppMin)&&(HT[i].intersection[0]==0)) {

			/*cout << "here ok -1" << endl;
			cout << "HT[i].intersection[0] = " << HT[i].intersection[0] << endl;

			cout << "here ok 0" << endl;*/

			///create the Itemset for the max Queue and its ItsInfo
			tmpIts=(Itemset*)calloc(1,sizeof(Itemset));
			//cout << "created itemset (1) with adress "<< tmpIts <<endl;
			living_itemsets++;
			tmpIts->support=IL[i].count;
			tmpInfo=new ItsInfo();
			tmpInfo->a_S_=IL[i].a_S;
			living_infos++; if(living_infos > maxlivinginfos) maxlivinginfos = living_infos;
			///create the nodeIDList
			tmpInfo->IDlist=(tipoInt*)malloc((numNPI[i]+1)*sizeof(tipoInt));
			tmpInfo->sizelista=numNPI[i]+1;
			tmpInfo->IDlist[0]=numNPI[i];
			for(tipoInt s=1; s<=numNPI[i];s++){
				tmpInfo->IDlist[s]=HT[i].visitedList[s].IDnode;
			}
			tmpInfo->minItem=i;
			tmpInfo->coreIndex=i;
			tmpInfo->typeList=0;
			tmpInfo->prefixCI=(tipoInt*)malloc(sizeof(tipoInt));
			tmpInfo->prefixCI[0]=0;

			//cout << "New empty ppce generation with support " << tmpIts->support << " and as " << tmpInfo->a_S_ << endl;

			//tmpInfo->nl_nodes.resize(numNPI[i]+1);
			tmpInfo->minnodes.reserve(numNPI[i]+1);
			tmpInfo->minnodes_indexes.reserve(numNPI[i]+1);

		    if(numNPI[i]+1 > max_nl_size_){
		        max_nl_size_ = numNPI[i]+1;
		        cout << " new max_nl_size_ = " << max_nl_size_ << endl;
		    }
		    //if(tmpInfo->nl_nodes.size() < 0){
		    if(numNPI[i]+1 < 0){
		        //cout << "ERROR WITH tmpInfo->nl_nodes.size(): " << tmpInfo->nl_nodes.size() << endl;
		        cout << "ERROR WITH numNPI[i]+1: " << numNPI[i]+1 << endl;
		    }

		    #ifdef trace_state
			program_state = -13;
			#endif

			/*cout << "-13" << endl;
			cout << "tmpInfo->nl_nodes.size() " << tmpInfo->nl_nodes.size() << endl;
			cout << "numNPI[i] " << numNPI[i] << endl;*/


			for(tipoInt s=1; s<=numNPI[i];s++){
				//tmpInfo->nl_nodes[s].push_back(HT[i].visitedList[s].IDnode);
				tmpInfo->minnodes.push_back(HT[i].visitedList[s].IDnode);
				tmpInfo->minnodes_indexes.push_back(s-1);
				/*cout << "node " << s << " has nl_list with node id " << HT[i].visitedList[s].IDnode << endl;*/
			}



		    /*
			int support_debug = 0;
			for(int i_=0; i_ < tmpInfo->minnodes.size(); i_++){
			    support_debug+= addressOf[ tmpInfo->minnodes[i_] ]->count;
			}
		    if(support_debug != tmpIts->support){
		      	cout << "something wrong with min nodes list of cloextension : " << endl;
		      	cout << "support_debug = " << support_debug << endl;
		       	cout << "tmpIts->support = " << tmpIts->support << endl;cout << "minnodes: " << endl;
				for(int i_=0; i_ < tmpInfo->minnodes.size(); i_++){
					cout << tmpInfo->minnodes[i_] << " " << endl;
				}
				cout << endl;
				cout << "minnodes_indexes: " << endl;
				for(int i_=0; i_ < tmpInfo->minnodes_indexes.size(); i_++){
					cout << tmpInfo->minnodes_indexes[i_] << " " << endl;
				}
				cout << endl;


		    }*/

			/* SignificantMiner */
			//temp_time = get_cpu_time();
			#ifdef trace_state
			program_state = -14;
			#endif
			/*cout << "-14" << endl;*/


		    tmpInfo->max_a_S = 0;
		    tmpInfo->min_a_S = tmpIts->support;

		    if(IL[i].a_S_permutations){

			        //cout << "cached a_s of ppce of item i = " << i << endl;

			    #ifdef trace_state
	  			program_state = -140;
	  			#endif

	            IL[i].a_S_permutations[0] = IL[i].a_S;
              #ifdef testexpanded
      				computePValues(tmpIts->support , IL[i].a_S_permutations); // torna qui
              tmpInfo->p_value=to_return[0];
              #endif
              #ifndef testexpanded
              tmpInfo->p_value=computePValue(tmpIts->support , IL[i].a_S);
              #endif

	       			for(int aj=1; aj < (jp+1); aj++){
	       				if(IL[i].a_S_permutations[aj] > tmpInfo->max_a_S){
	       					tmpInfo->max_a_S = IL[i].a_S_permutations[aj];
	       				}
	       				else{
	       					if(IL[i].a_S_permutations[aj] < tmpInfo->min_a_S){
	       						tmpInfo->min_a_S = IL[i].a_S_permutations[aj];
	       					}
	       				}
	       			}
		    }
		    else{

		        //cout << "generating a_s of ppce of item i = " << i << endl;

		  			#ifdef trace_state
		  			program_state = -141;
		  			#endif

		        // reset perm_supp
		        for(int aj=0; aj<perm_supp.size(); aj++)
		        		perm_supp[aj]=0;

		        perm_supp[0]=IL[i].a_S;

		        // compute the support on the permutations of current ppce
		        tested++;
		        /*for(int i_=0; i_ < tmpInfo->nl_nodes.size(); i_++){
		          //cout << "Node  " << i_ << " points to: " << endl;
		          for(int j=0; j <  tmpInfo->nl_nodes[i_].size(); j++){
		            //cout << " " << j <<" - node  " << ppce->info->nl_nodes.at(i_).at(j) << endl;
		            if(addressOf[tmpInfo->nl_nodes[i_][j]]->count >= 128){
		              for(int aj=1;aj<(jp+1); aj++){
		                perm_supp[aj] += addressOf[ tmpInfo->nl_nodes[i_][j] ]->a_S_[aj];
		              }
		            }
		            else{
		              for(int aj=1;aj<(jp+1); aj++){
		                perm_supp[aj] += addressOf[ tmpInfo->nl_nodes[i_][j] ]->a_S[aj];
		              }
		            }
		          }
		        }*/
		        int support_debug = 0;
            uint8_t *temp_a_S;
            tipoInt *temp_a_S_;
		        for(int i_=0; i_ < tmpInfo->minnodes.size(); i_++){
		            if(addressOf[tmpInfo->minnodes[i_]]->count >= 128){
                  temp_a_S_ = addressOf[ tmpInfo->minnodes[i_] ]->a_S_;
		              for(int aj=1;aj<(jp+1); aj++){
		                perm_supp[aj] += temp_a_S_[aj];
		              }
		            }
		            else{
                  temp_a_S = addressOf[ tmpInfo->minnodes[i_] ]->a_S;
		              for(int aj=1;aj<(jp+1); aj++){
		                perm_supp[aj] += temp_a_S[aj];
		              }
		            }

		            support_debug+= addressOf[ tmpInfo->minnodes[i_] ]->count;
		        }

		        if(support_debug != tmpIts->support){
		        	cout << "something wrong with min nodes list of cloextension : " << endl;
		        	cout << "support_debug = " << support_debug << endl;
		        	cout << "tmpIts->support = " << tmpIts->support << endl;
		        }

            #ifdef testexpanded
		        computePValues(tmpIts->support , perm_supp.data()); // torna qui
		   			tmpInfo->p_value=to_return[0];
            #endif

            #ifndef testexpanded
            tmpInfo->p_value=computePValue(tmpIts->support , IL[i].a_S);
            #endif



		   			for(int aj=1; aj < (jp+1); aj++){
		   				tmpInfo->max_a_S = max( tmpInfo->max_a_S , perm_supp[aj] );
		   				tmpInfo->min_a_S = min( perm_supp[aj] , tmpInfo->min_a_S );
		   			}

		    }

      //cout << "tmpInfo->p_value = " << tmpInfo->p_value << endl;
      //cout << "tmpInfo->min_a_S = " << tmpInfo->min_a_S << endl;
      //cout << "tmpInfo->max_a_S = " << tmpInfo->max_a_S << endl;
			//tmpInfo->a_S=(tipoInt*)malloc((jp+1)*sizeof(tipoInt));


			/*cout << "-141" << endl;*/

			//cout << "ppce of Clo(empty) with support " << tmpIts->support << " has max_a_S " << tmpInfo->max_a_S << " and min " << tmpInfo->min_a_S << endl;

			/*if(tmpInfo->min_a_S > tmpInfo->max_a_S){
				cout << "ERROR HERE 1" << endl;
				cout << "tmpIts->support " << tmpIts->support << endl;
				cout << "tmpInfo->min_a_S " << tmpInfo->min_a_S << endl;
				cout << "tmpInfo->max_a_S " << tmpInfo->max_a_S << endl;
			}*/

			#ifdef trace_state
			program_state = -142;
			#endif

			/*cout << "-142" << endl;*/



			//cout << "(1) allocated jp+1 things at " << tmpInfo->a_S << endl;

			//cout << "here ok 1" << endl;
			//cout << "IL[i].item = " << IL[i].item << endl;
			//cout << "IL[i].a_S = " << IL[i].a_S << endl;
			//cout << "IL[i].a_S[1] = " << IL[i].a_S[1] << endl;
			//cout << "IL[i].count = " << IL[i].count << endl;
			//for(int aj=0; aj<jp+1; aj++)
			//tmpInfo->a_S[aj]=IL[i].a_S[aj];
			//time_keeper += get_cpu_time() - temp_time;
			//cout << "not set a_S for " << tmpIts->support << endl;
			/* SignificantMiner */
			tmpIts->info=tmpInfo;

			/* SignificantMiner */
			/*if(DEEP_DEBUG){
				cout << "support of item set to " << tmpIts->support << endl;
				cout << "a_S for item set to " << tmpIts->info->a_S << endl;
			}*/
			/* SignificantMiner */

			///create the Itemset for the min Queue: its ItsInfo is the same of the maxQueue
			//minIts=(Itemset*)calloc(1,sizeof(Itemset));
			//cout << "created itemset  (2) with adress "<< minIts <<endl;
			living_itemsets++;
			//minIts->support=IL[i].count;
			//minIts->info=tmpInfo;

			/*cout << "-1421" << endl;*/

			q->insert(tmpIts);	///insertion in the max queue
			//qMin->insert(minIts);	///insertion in the min queue
			if(q->getInqueue() > maxinqueue) maxinqueue = q->getInqueue();

			/*cout << "-1422" << endl;*/

      //cout << "tmpIts generated! tmpIts->support " << tmpIts->support << " tmpIts->info->a_S_ " << tmpIts->info->a_S_ << endl;

			      // if the extracted itemset can be significant
			      bool can_be_significant = getLogPsiExact(tmpIts->support) <= getLogPsi(suppMin);
			      if(can_be_significant){
			      double log_p_value = computeLogPValue(tmpIts->support , tmpIts->info->a_S_);
			      if( store_results && log_p_value <= getLogPsi(suppMin) ){
			        //cout << "tmpIts p_value: " << tmpIts->info->p_value << " tmpIts support " << tmpIts->support << " tmpIts a_S " << tmpIts->info->a_S_ << endl;
			        //cout << "RES QUEUE OBSERVED INSERTION: it can be significant since getPsi(suppMin) = " << getPsi(suppMin) << endl;
			        //double log_p_value = computeLogPValue(tmpIts->support , tmpIts->info->a_S_);

			        //cout << "   pattern p_value: " << tmpIts->info->p_value << " pattern log_p_value " << log_p_value << endl;
			        //cout << "   getPsi(suppMin): " << getPsi(suppMin) << " getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;

			        //cout << "toOutput_pattern->log_p_value = " << toOutput_pattern->log_p_value << endl;

			        q_res->insert_observed(log_p_value);

                  /* Top-K strategy */
                  while(K_significant_patterns > 0 && q_res->getKElements() >= K_significant_patterns){
                    int new_significance_level = q_res->getMaxIndex();
                    new_significance_level++;
                    if(new_significance_level > suppMin){
                      suppMin = new_significance_level;
                      suppMinCurr = new_significance_level;
                    }

                    cout << "Top-K strategy update! ";
                    /*cout << "K_significant_patterns " << K_significant_patterns << endl;
                    cout << "q_res->getKElements() returned " << q_res->getKElements() << endl;
                    cout << "q_res->getMaxIndex() returned " << q_res->getMaxIndex() << endl;*/
                    cout << " s_supp " << s_supp;
                    cout << " new suppMin " << suppMin << endl;
                    /*cout << "getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;
                    cout << "getPsi(suppMin) " << getPsi(suppMin) << endl;*/

                    delta = getPsi(suppMin);
                    q_res->removeMax();
                  }

      		}
      }



      /*cout << "-1423" << endl;*/

			free(tmpIts);

			/* cout << "-1424" << endl;*/
			//free(minIts);
			//living_itemsets--;
			living_itemsets--;

			//cout << "here ok 2" << endl;

			///heuristic on suppMinCurr

			/*if (inserted<k-1){
				utilNodeArray[tmpIts->support]++;
			}
			else{
				if (inserted==k-1){
					///cumulate over the utilNodeArray
					for(tipoInt c=maximumSupport-1;c>=suppMin;c--){
						utilNodeArray[c]+=utilNodeArray[c+1];
					}
					firstBound=true;
				}
				///raise minimum support using utilNodeArray: this is for the bound on suppMinCurr
				trivialHeurCurr(tmpIts->support);
			}*/
			inserted++;

			///heuristic on suppMin

			/*itsOfSupp[tmpIts->support]++;

			if (inserted>k_max){
				///we must check if we can have a new bound on suppMin
				if (k_max<inqueue-itsOfSupp[suppMin]+1){
					///we have a new bound on suppMin and we can delete non frequent itemsets from the queue
					Itemset* todel;
					for (tipoInt t=0; t<itsOfSupp[suppMin]; t++){
						todel=qMin->removeMin();
						inqueue--;
						free(todel->info->prefixCI);
						free(todel->info->IDlist);
						delete todel->info;
						cout << "removed infrequent itemset with support " << todel->support << endl;
					}
					suppMin=qMin->heap->support;
				}
			}
			else{
				if (inserted==k_max){
					///this is the first time we see k_max itemset in queue --> bound on suppMin
					suppMin=qMin->heap->support;
				}
			}*/


	#ifdef trace_state
	program_state = -142;
	#endif
	/*cout << "-142" << endl;*/

			/* SignificantMiner */
			// check if infrequent itemsets can be removed from queue
			if(DEEP_DEBUG)
			cout << "checking to remove infrequent itemsets..." << endl;
			//Itemset* todel;
			//todel = qMin->getMin();
			//temp_time = get_cpu_time();
			//while(todel != NULL && todel->support < suppMin){
			debug_info_1 = q->getMin();
			debug_info_2 = q->getInqueue();
			while(q->getMin() > 0 && q->getMin() < suppMin && q->getMax() >= suppMin){
				//.cout << "free itemset (11) with adress "<< todel <<endl;
				//living_itemsets--;
				//free(todel);
				/*cout << "remove infrequent itemsets..." << endl;*/

				//todel=qMin->removeMin();
				//q->printLowestHeapLevel();
				q->removeMin();
				//free(todel->info->prefixCI);
				//free(todel->info->IDlist);
				//delete todel->info;
				//todel->info = NULL;
				//if(DEEP_DEBUG)
				//cout << "removed infrequent itemset with support " << todel->support << endl;

				//.cout << "free itemset (1) with adress "<< todel <<endl;
				//free(todel);
				//living_itemsets--;
				//todel = qMin->getMin();
			}

			/*cout << "-142 bis" << endl;*/
			//if(todel!=NULL){
				//.cout << "free itemset (12) with adress "<< todel <<endl;
				//living_itemsets--;
				//free(todel);
			//}
			//time_keeper += get_cpu_time() - temp_time;
			/* SignificantMiner */

			//cout << "here ok 3" << endl;


			///in visitedList[0].IDnode we put the length of the current list
			HT[i].visitedList[0].IDnode=0;

			//cout << "-143" << endl;


		#ifdef trace_state
		program_state = -143;
		#endif

			//cout << "here ok 4" << endl;
		}
	}
}


/** fill the HT: it requires the StartNodeList, the Itemset (from which it generates the ppc-e) and the Queue
Hp: it uses min-item nodes
*/

void fillHT(StartNodeList startList, Itemset* father, typeQueue q){
//void fillHT(StartNodeList startList, Itemset* father, QueueMinMax* q){

	///set the right lastVisitedNode pointer and lastNodePointer
	for (int r=0; r<father->info->minItem; r++){
		HT[r].lastNodePnt=&(HT[r].pointerList);
	}
	///read the startList and initialize the lists of HT

	PatriciaNode* patNode;		///PatriciaNode pointed by current
	Transaction transTmp;		///maximum intersection recoverable by patNode and current intersection

	int i=1;	///pointer to elements of inter
	int j=1;	///pointer to element of interTmp
	int res=1;	///pointer to next element of inter that will contain

	for (tipoInt m=1; m<=startList[0].ID; m++){
		///read the node and insert it in all the visitedLists of items there are lower than coreIndex of Itemset; we must also update count and intersection
		patNode=addressOf[startList[m].ID];	///node pointed by current
		tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode
		tipoInt numInter=(tipoInt) (patNode->lastInter-patNode->lastItem);

		if (nextfree+numItems+numInter>lastvalid){
      cout << "call manual mem 10" << endl;
			resizeManualMem();
		}
		transTmp=nextfree;
		nextfree=nextfree+numItems+numInter;	///not "...+1;" beacuase minItem will not stay in the intersection

		///now set the intersection
		for(int n=1; n<=numInter;n++){
			transTmp[n]=(patNode->lastItem+1)[numInter-n];
		}

		transTmp[0]=numInter;
		///current item ID
		tipoInt index;
		///we must begin from the last item because we can update intersection
		for (int l=numItems-1; l>=0; l--){
			index=(&patNode->item)[l];
			if (index<father->info->minItem){
				///update count
				HT[index].count+=patNode->count;
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				/*for(int aj=0; aj<jp+1; aj++)
				HT[index].a_S[aj]+=patNode->a_S[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				//cout << "increased a_S (1) for patricia node of " << patNode->a_S << " to " << HT[index].a_S << endl;
				/* SignificantMiner */

				///we update intersection and visitedList only of item<coreindex
				if (index<father->info->coreIndex){
					///update the intersection in HT
					if (HT[index].intersection!=NULL){

						///update intersection: we intersect the transactions and memorize the intersection in the "inter" transaction
						i=1;	///pointer to elements of HT[l].intersection
						j=1;	///pointer to element of transLast
						res=1;	///pointer to next element of intersection that will contain an item in the intersection
						while(i<=HT[index].intersection[0]&&j<=transTmp[0]){
							if (HT[index].intersection[i]==transTmp[j]){
									///intersection[i] is in the intersection
								HT[index].intersection[res]=HT[index].intersection[i];
								res++;
								i++;
								j++;
							}
							else{
								if(HT[index].intersection[i]>transTmp[j]){
										///intersection[i] is not in the intersection
									i++;
								}
								else{
										///intersection[i] can still be in the intersection
									j++;
								}
							}
						}
						HT[index].intersection[0]=res-1;
					}
					else{
						///set the intersection
						if (nextfree+transTmp[0]+1>lastvalid){
              cout << "call manual mem 11" << endl;
							resizeManualMem();
						}
						HT[index].intersection=nextfree;
						nextfree=nextfree+transTmp[0]+1;
						for (int o=0;o<=transTmp[0]; o++){
							HT[index].intersection[o]=transTmp[o];
						}
					}

					///update visitedList and lastVisitedNode
					tipoInt numElem=HT[index].visitedList[0].IDnode;
					HT[index].visitedList[0].IDnode++;
					}
			}
			///update transLast if item is different from minItem
			if (index!=father->info->minItem){
				transTmp[0]++;
				transTmp[transTmp[0]]=index;
			}
		}

		///now transTmp contains all the intersection that must be passed to the father of patNode

		///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
		PatriciaNode* upper=patNode->father;	///upper is the father of patNode
		if (upper!=NULL){	///the node isn't a root's children
			if (upper->pntToNode==NULL){

				///create the pointer to the node
				if ((nextfree+sizePnt+transTmp[0])>lastvalid){
          cout << "call manual mem 12" << endl;
					resizeManualMem();
				}
				upper->pntToNode=(NodePointer*)nextfree;
				//todelete->push_back(upper->pntToNode);
				living_nodepointers++;
				nextfree=nextfree+transTmp[0]+sizePnt;	///update nextfree
				upper->pntToNode->count=patNode->count;
				//temp_time = get_cpu_time();

				//upper->pntToNode->a_S_index = getNewTodeleteIndex();

				//cout << "(2) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

				//cout << "access (1) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				todelete->at(upper->pntToNode->a_S_index)[aj]=patNode->a_S[aj];*/
			/*cout << "done 1.2 " << nextfree_a_S << endl;
			for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
				//time_keeper += get_cpu_time() - temp_time;
				/* SignificantMiner */
				/*if(upper->pntToNode->a_S!=addressOf[upper->pntToNode->node]->a_S){
					if(addressOf[upper->pntToNode->node]->perm_a_S)
					cout << "todelete->at(upper->pntToNode->a_S_index) = " << todelete->at(upper->pntToNode->a_S_index) << " , addressOf[upper->pntToNode->node]->a_S = " << addressOf[upper->pntToNode->node]->a_S << endl;
					return;
				}
				upper->pntToNode->a_S=patNode->a_S;*/
				/* SignificantMiner */
				upper->pntToNode->node=upper->ID;
				upper->pntToNode->nextNodePointer=-1;
				///we copy transTmp into the new intersection
				(&upper->pntToNode->intersection)[0]=transTmp[0];
				for (int o=1;o<=transTmp[0]; o++){
					(&upper->pntToNode->intersection)[o]=transTmp[o];
				}
				///insert the node in the pointerList of its greater item
				*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
				///update lastNodePnt
				HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
			}
			else{
				///we must update the field of the pointer to upper
				upper->pntToNode->count+=patNode->count;
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//cout << "access (2) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				todelete->at(upper->pntToNode->a_S_index)[aj]+=patNode->a_S[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				/*upper->pntToNode->a_S+=patNode->a_S;
				if(DEEP_DEBUG)
				cout << "increased a_S (2) for patricia node of " << patNode->a_S << " to " << todelete->at(upper->pntToNode->a_S_index) << endl;*/
				/* SignificantMiner */
				///as intersection we can use transTmp
				i=1;	///pointer to elements of upper->pntToNode->intersection
				j=1;	///pointer to element of transTmp
				res=1;	///pointer to next element of upper->pntToNode->intersection that will contain an item in the intersection
				while(i<=(&upper->pntToNode->intersection)[0]&&j<=transTmp[0]){
					if ((&upper->pntToNode->intersection)[i]==transTmp[j]){
									///upper->pntToNode->intersection[i] is in the intersection
						(&upper->pntToNode->intersection)[res]=(&upper->pntToNode->intersection)[i];
						res++;
						i++;
						j++;
					}
					else{
						if((&upper->pntToNode->intersection)[i]>transTmp[j]){
										///upper->pntToNode->intersection[i] is not in the intersection
							i++;
						}
						else{
										///upper->pntToNode->intersection[i] can still be in the intersection
							j++;
						}
					}
				}
				(&upper->pntToNode->intersection)[0]=res-1;
			}
		}
	}

	///now we must consider all the list of items greater than father->info->minItem and visit the element in their pointerList; we must begin from the last item (the greatest)

	NodePointer* currentpnt;	///temporary NodePointer used for the elements in the pointedList
	tipoInt off;			///offset of currentpnt in manualmem

	for (tipoInt l=father->info->minItem-1;l>=0;l--){

		off=HT[l].pointerList;
		while(off!=-1){
			currentpnt=(NodePointer*)&(manualmem[off]);
			patNode=addressOf[currentpnt->node];	///PatriciaNode pointed by currentpnt

			tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode
			if ((nextfree+numItems+(&currentpnt->intersection)[0]+1)>lastvalid){
        cout << "call manual mem 13" << endl;
				resizeManualMem();
			}
			transTmp=nextfree;
			nextfree=nextfree+numItems+(&currentpnt->intersection)[0]+1;
			transTmp[0]=(&currentpnt->intersection)[0];
			for (int o=1;o<=transTmp[0];o++){
				transTmp[o]=(&currentpnt->intersection)[o];
			}

			///now consider all the items in the node and update the corresponding HT rows
			tipoInt index;	///current item ID
			for (tipoInt s=numItems-1;s>=0; s--){
				index=(&patNode->item)[s];

				///we know that all the items in the node are all lower than father->info->minItem
				HT[index].count+=currentpnt->count;	///update the count with the count of the pointer
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//cout << "access (3) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				HT[index].a_S[aj]+=todelete->at(currentpnt->a_S_index)[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				//HT[index].a_S+=currentpnt->a_S;
				//if(DEEP_DEBUG)
				//cout << "increased a_S (3) for patricia node of " << currentpnt->a_S << " to " << HT[index].a_S << endl;
				/* SignificantMiner */

				///we update intersection and visitedList only of item<coreindex

				if (index<father->info->coreIndex){
					///update the intersection in HT
					if (HT[index].intersection!=NULL){
						///update intersection: we intersect the transactions and memorize the intersection in the "inter" transaction
						i=1;	///pointer to elements of HT[l].intersection
						j=1;	///pointer to element of transLast
						res=1;	///pointer to next element of intersection that will contain an item in the intersection
						while(i<=HT[index].intersection[0]&&j<=transTmp[0]){
							if (HT[index].intersection[i]==transTmp[j]){
									///intersection[i] is in the intersection
								HT[index].intersection[res]=HT[index].intersection[i];
								res++;
								i++;
								j++;
							}
							else{
								if(HT[index].intersection[i]>transTmp[j]){
										///intersection[i] is not in the intersection
									i++;
								}
								else{
										///intersection[i] can still be in the intersection
									j++;
								}
							}
						}
						HT[index].intersection[0]=res-1;
					}
					else{
						///set the intersection
						if (nextfree+transTmp[0]+1>lastvalid){
              cout << "call manual mem 15" << endl;
							resizeManualMem();
						}
						HT[index].intersection=nextfree;
						nextfree=nextfree+transTmp[0]+1;
						for (int o=0;o<=transTmp[0];o++){
							HT[index].intersection[o]=transTmp[o];
						}
					}

					///update visitedList and lastVisitedNode
					tipoInt numElem=HT[index].visitedList[0].IDnode;
				}
				///here there is no need to verify if the item is different from minItem because we are not considering nodes that contains minItem
				transTmp[0]++;
				transTmp[transTmp[0]]=index;
			}

			///now transTmp contains all the intersection that must be passed to the pointer to the father of patNode

			///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
			PatriciaNode* upper=patNode->father;	///upper is the father of patNode
			if (upper!=NULL){	///the node isn't a root's children
				if (upper->pntToNode==NULL){
				///create the pointer to the node
				if (nextfree+sizePnt+transTmp[0]>lastvalid){
          cout << "call manual mem 16" << endl;
					resizeManualMem();
				}
				//upper->pntToNode=(NodePointer*)nextfree;
				//todelete->push_back(upper->pntToNode);
				living_nodepointers++;
				nextfree=nextfree+transTmp[0]+sizePnt;	///update nextfree
				upper->pntToNode->count=currentpnt->count;
				/* SignificantMiner */
				//temp_time = get_cpu_time();

				//upper->pntToNode->a_S_index = getNewTodeleteIndex();

				//cout << "(3) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

				//cout << "access (4) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				//cout << "access (4) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				todelete->at(upper->pntToNode->a_S_index)[aj]=todelete->at(currentpnt->a_S_index)[aj];*/
			/*cout << "done 2.2 " << nextfree_a_S << endl;
			for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
				//time_keeper += get_cpu_time() - temp_time;
				/* SignificantMiner */
				upper->pntToNode->node=upper->ID;
				upper->pntToNode->nextNodePointer=-1;
				///we copy transTmp into the new intersection
				(&upper->pntToNode->intersection)[0]=transTmp[0];
				for (int o=1;o<=transTmp[0]; o++){
					(&upper->pntToNode->intersection)[o]=transTmp[o];
				}
				///insert the node in the pointerList of its greater item
				*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
				///update lastNodePnt
				HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
				///create the pointer to the node
				}
				else{
					///we must update the field of the pointer to upper
					upper->pntToNode->count+=currentpnt->count;
					/* SignificantMiner */
					//temp_time = get_cpu_time();
					//cout << "access (5) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (5) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
					/*for(int aj=0; aj<jp+1; aj++)
					todelete->at(upper->pntToNode->a_S_index)[aj]+=todelete->at(currentpnt->a_S_index)[aj];*/
					//time_keeper += get_cpu_time() - temp_time;
					/*if(DEEP_DEBUG)
					cout << "increased a_S (4) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;*/
					/* SignificantMiner */
					///as intersection we can use transTmp
					i=1;	///pointer to elements of upper->pntToNode->intersection
					j=1;	///pointer to element of transTmp
					res=1;	///pointer to next element of upper->pntToNode->intersection that will contain an item in the intersection
					while(i<=(&upper->pntToNode->intersection)[0]&&j<=transTmp[0]){
						if ((&upper->pntToNode->intersection)[i]==transTmp[j]){
										///upper->pntToNode->intersection[i] is in the intersection
							(&upper->pntToNode->intersection)[res]=(&upper->pntToNode->intersection)[i];
							res++;
							i++;
							j++;
						}
						else{
							if((&upper->pntToNode->intersection)[i]>transTmp[j]){
											///upper->pntToNode->intersection[i] is not in the intersection
								i++;
							}
							else{
											///upper->pntToNode->intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					(&upper->pntToNode->intersection)[0]=res-1;
				}
			}
			off=currentpnt->nextNodePointer;
			patNode->pntToNode=NULL;
		}
	}
}


///clear the first id+1 rows of the header table
void clearHT(tipoInt id){
	for (tipoInt i=0; i<=id; i++){
		HT[i].count=0;			///count of the extension
		HT[i].a_S=0;			///a_S of the extension
		/* SignificantMiner */
		//temp_time = get_cpu_time();
		//free(HT[i].a_S);
		//HT[i].a_S = (tipoInt*)malloc((jp+1) * sizeof(tipoInt));
		//for(int aj=0; aj<jp+1; aj++)
		//HT[i].a_S[aj]=0;			///a_S of the extension
		//time_keeper += get_cpu_time() - temp_time;
		/* SignificantMiner */
		HT[i].intersection=NULL;	///intersection on the prefix
		HT[i].visitedList[0].IDnode=0;		///number of vidited elements
		HT[i].pointerList=-1;		///empty list
		HT[i].lastNodePnt=&(HT[i].pointerList);
	}
}

///main function

/**find the closure of the (closed) itemset from which the HT was filled and insert its (frequent) ppc-e in the queue;
Hp: it uses min-item nodes
*/

Transaction findPpce(Itemset* father, typeQueue q, StartNodeList startListFather){
//Transaction findPpce(Itemset* father, QueueMinMax* q, StartNodeList startListFather){
	Transaction closure;	///here we insert all the items that are in Clo(father)
	if (nextfree+father->info->minItem+2>lastvalid){
    cout << "call manual mem 17" << endl;
		resizeManualMem();
	}
	closure=nextfree;
	nextfree=nextfree+father->info->minItem+2;
	closure[0]=1;
	closure[1]=father->info->minItem;
	///now we check all the rows of HT to find the items in Clo(father) and the ppce of father

	///we must begin from i=minItem-1 because so we can update intersection deleting the items in Clo(father)
	Itemset* ppce;	///used to create new Itemset to insert in queue
	int parity;			///if parity%2==1 then a StartNode must be inserted in the ppce->info->IDlist; also all the element with label!=0 must be inserted in the list
	for (int i=father->info->minItem-1; i>=0; i--){
		///if the item has count<suppMin, there is nothing to do
		if (HT[i].count>=suppMin){
			///if the item has the same count of father, it is in closure=Clo(father)
			if (HT[i].count==father->support){
				///insert i in closure
				closure[0]++;
				closure[closure[0]]=i;
			}
			else{
				///if it isn't in closure and it is < coreIndex, we must check if it can generate a ppce: conditions on coreindex and prefix-preservation
				if(i<father->info->coreIndex){
					///to verify if this extension is prefix-preserving, we check if the intersection in HT[i].intersection contains only items in closure-minItem; this is real IFF the length are the same (because HT[i].intersection surely contains the item in closure that are lower than i)
					if(HT[i].intersection[0]==closure[0]-1){
						///i is the core index of a ppc-e of father
						ppce=(Itemset*)calloc(1,sizeof(Itemset));
						//.cout << "created itemset (3) with adress "<< ppce <<endl;
						living_itemsets++;
						ppce->support=HT[i].count;
						ppce->info=new ItsInfo();
						living_infos++; if(living_infos > maxlivinginfos) maxlivinginfos = living_infos;
						ppce->info->minItem=father->info->minItem;
						ppce->info->coreIndex=i;
						/* SignificantMiner */
						//temp_time = get_cpu_time();
						//ppce->info->a_S=(tipoInt*)malloc((jp+1) * sizeof(tipoInt));

						//cout << "(4) allocated jp+1 things at " << ppce->info->a_S << endl;


						//for(int aj=0; aj<jp+1; aj++)
						//ppce->info->a_S[aj]=HT[i].a_S[aj]; // this can be used now and not keep in the queue

						#ifdef trace_state
						program_state = -15;
						#endif


						//computePValues(ppce->support , HT[i].a_S); // torna qui


						//time_keeper += get_cpu_time() - temp_time;
						/* SignificantMiner */

						///now we must find the NodeIDList for the ppc-e

						///first of all, we mark the StartNode in the father's StartNodeList
						tipoInt numNod=0;		///number of nodes that will be inserted in the nodeIDList

						///create the nodeIDList
						ppce->info->IDlist=(tipoInt*)malloc((numNod+1)*sizeof(tipoInt));
						ppce->info->sizelista=numNod+1;

						///now scan the startListFather and insert the ID in the nodeIDList
						int ne=1;	///next position for insertion in nodeIDList
						parity=0;	///summation on "label"
						for (tipoInt j=1;j<=startListFather[0].ID;j++){
							if ((startListFather[j].label!=0)||(parity%2==1)){
								ppce->info->IDlist[ne]=startListFather[j].ID;
								parity=parity+startListFather[j].label;
								startListFather[j].label=0;
								ne++;
							}
						}
						ppce->info->IDlist[0]=numNod;


						ppce->info->a_S_=HT[i].a_S;
						ppce->info->p_value=to_return[0];

						q->insert(ppce);
						if(q->getInqueue() > maxinqueue) maxinqueue = q->getInqueue();

						#ifdef trace_state
						program_state = -16;
						#endif


						/*if (inserted<k_max-1){
							utilNodeArray[ppce->support]++;
						}
						else{
							if (inserted==k_max-1){
								///cumulate over the utilNodeArray
								for(tipoInt c=maximumSupport-1;c>=suppMin;c--){
									utilNodeArray[c]+=utilNodeArray[c+1];
								}
							}
							///raise minimum support using utilNodeArray
							trivialHeur(ppce->support);
						}
						inserted++;*/

					}
				}
			}
		}
	}
	return closure;
}

/** given an Itemset, find its closure without generating its ppc-e; used in computation to produce closed itemset with the same support of the k-th.
Hp: the StartNodeList is referred to min-item nodes
*/

Transaction findClo(StartNodeList startList, Itemset* father){
	///set the right lastNodePointer (we're not interested in lastVisitedNodes because we don't generate ppc-e)
	for (int r=0; r<num_item; r++){
		HT[r].lastNodePnt=&(HT[r].pointerList);
	}
	///read the startList and initialize the lists of HT
	PatriciaNode* patNode;		///PatriciaNode pointed by current

	for (tipoInt m=1; m<=startList[0].ID; m++){
		///read the node and update count of the items it contains

		patNode=addressOf[startList[m].ID];	///node pointed by current
		tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode
		tipoInt numInter=(tipoInt) (patNode->lastInter-patNode->lastItem);
		///current item ID
		tipoInt index;
		///we can begin from the last item
		for (int l=numItems-1; l>=0; l--){
			index=(&patNode->item)[l];
			if (index<father->info->minItem){
				///update count
				HT[index].count+=patNode->count;
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				/*for(int aj=0; aj<jp+1; aj++)
				HT[index].a_S[aj]+=patNode->a_S[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				//if(DEEP_DEBUG)
				//cout << "increased a_S (5) for patricia node of " << patNode->a_S << " to " << HT[index].a_S << endl;
				/* SignificantMiner */
			}
		}

		///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
		PatriciaNode* upper=patNode->father;	///upper is the father of patNode
		if (upper!=NULL){	///the node isn't a root's children
			if (upper->pntToNode==NULL){
				///create the pointer to the node
				if ((nextfree+sizePnt)>lastvalid){
          cout << "call manual mem 18" << endl;
					resizeManualMem();
				}
				upper->pntToNode=(NodePointer*)nextfree;
				//todelete->push_back(upper->pntToNode);
				living_nodepointers++;
				nextfree=nextfree+sizePnt;	///update nextfree
				upper->pntToNode->count=patNode->count;
				/* SignificantMiner */
				//temp_time = get_cpu_time();

				//upper->pntToNode->a_S_index = getNewTodeleteIndex();

				//cout << "(5) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;


				//cout << "access (6) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				todelete->at(upper->pntToNode->a_S_index)[aj]=patNode->a_S[aj];*/
			/*cout << "done 3.2 " << nextfree_a_S << endl;
			for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
				//time_keeper += get_cpu_time() - temp_time;
				/* SignificantMiner */
				upper->pntToNode->node=upper->ID;
				upper->pntToNode->nextNodePointer=-1;
				(&upper->pntToNode->intersection)[0]=0;
				///insert the node in the pointerList of its greater item
				*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
				///update lastNodePnt
				HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
			}
			else{	///we must update the field of the pointer to upper
				upper->pntToNode->count+=patNode->count;
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//cout << "access (7) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				todelete->at(upper->pntToNode->a_S_index)[aj]+=patNode->a_S[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				/*if(DEEP_DEBUG)
				cout << "increased a_S (6) for patricia node of " << patNode->a_S << " to " << todelete->at(upper->pntToNode->a_S_index) << endl;*/
				/* SignificantMiner */
			}
		}
	}


	///now we must consider all the list of items greater than father->info->minItem and visit the element in their pointerList; we must begin from the last item (the greatest)

	NodePointer* currentpnt;	///temporary NodePointer used for the elements in the pointedList
	tipoInt off; 	///offset of the pointer in manualmem
	for (tipoInt l=father->info->minItem-1;l>=0;l--){
		off=HT[l].pointerList;
		///for every element in the pointedList we must visit the node
		while(off!=-1){
			currentpnt=(NodePointer*)&(manualmem[off]);
			patNode=addressOf[currentpnt->node];	///PatriciaNode pointed by currentpnt
			tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode

			///now consider all the items in the node and update the corresponding HT rows
			tipoInt index;	///current item ID
			for (tipoInt s=numItems-1;s>=0; s--){
				index=(&patNode->item)[s];
				///we know that all the items in the node are all lower than father->info->minItem
				HT[index].count+=currentpnt->count;	///update the count with the count of the pointer
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//cout << "access (8) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				HT[index].a_S[aj]+=todelete->at(currentpnt->a_S_index)[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				/*if(DEEP_DEBUG)
				cout << "increased a_S (7) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << HT[index].a_S << endl;
				*//* SignificantMiner */
			}

			///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
			PatriciaNode* upper=patNode->father;	///upper is the father of patNode
			if (upper!=NULL){	///the node isn't a root's children
				if (upper->pntToNode==NULL){
					///create the pointer to the node
					if (nextfree+sizePnt>lastvalid){
            cout << "call manual mem 19" << endl;
						resizeManualMem();
					}
					upper->pntToNode=(NodePointer*)nextfree;
					//todelete->push_back(upper->pntToNode);
					living_nodepointers++;
					nextfree=nextfree+sizePnt;	///update nextfree
					upper->pntToNode->count=currentpnt->count;
					/* SignificantMiner */
					//temp_time = get_cpu_time();

					//upper->pntToNode->a_S_index = getNewTodeleteIndex();

					//cout << "(6) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

					//cout << "access (9) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (9) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
					/*for(int aj=0; aj<jp+1; aj++)
					todelete->at(upper->pntToNode->a_S_index)[aj]=todelete->at(currentpnt->a_S_index)[aj];*/
				/*cout << "done 4.2 " << nextfree_a_S << endl;
				for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
					//time_keeper += get_cpu_time() - temp_time;
					/* SignificantMiner */
					upper->pntToNode->node=upper->ID;
					upper->pntToNode->nextNodePointer=-1;
					///insert the node in the pointerList of its greater item
					*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
					///update lastNodePnt
					HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
				}
				else{	///we must update the field of the pointer to upper
					upper->pntToNode->count+=currentpnt->count;
					/* SignificantMiner */
					//temp_time = get_cpu_time();
					//cout << "access (10) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (10) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
					/*for(int aj=0; aj<jp+1; aj++)
					todelete->at(upper->pntToNode->a_S_index)[aj]+=todelete->at(currentpnt->a_S_index)[aj];*/
					//time_keeper += get_cpu_time() - temp_time;
					/*if(DEEP_DEBUG)
					cout << "increased a_S (8) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;*/
					/* SignificantMiner */
				}
			}
			off=currentpnt->nextNodePointer;
			addressOf[currentpnt->node]->pntToNode=NULL;	///delete the reference to the current node
		}
	}

	///now find the items in Clo(father)
	Transaction closure;	///here we insert all the items that are in Clo(father)
	if (nextfree+father->info->minItem+2>lastvalid){
    cout << "call manual mem 20" << endl;
		resizeManualMem();
	}
	closure=nextfree;
	nextfree=nextfree+father->info->minItem+2;
	closure[0]=1;
	closure[1]=father->info->minItem;
	for (int i=father->info->minItem-1; i>=0; i--){
				///if the item has the same count of father, it is in closure=Clo(father)
		if (HT[i].count==father->support){
			///insert i in closure
			closure[0]++;
			closure[closure[0]]=i;
		}
	}
	return closure;
}

/** fill the HT starting form nodes that contains the coreIndex: it requires the NodeIDList, the Itemset (from which generate the ppc-e) and the Queue; it returns the coreIndex-prefix of the father (without minItem)
*/

//void fillHTCI(NodeIDList nodeList, Itemset* father, QueueMinMax* q){
void fillHTCI(NodeIDList nodeList, Itemset* father, typeQueue q){

	//cout << "FILLHTCI CALL" << endl;

	///read the startList and initialize the lists of HT
	PatriciaNode* patNode;		///PatriciaNode pointed by current
	Transaction transTmp;		///maximum intersection recoverable by patNode and current intersection
	tipoInt* temp_as;
	tipoInt* temp_as_;

	int i=1;	///pointer to elements of inter
	int j=1;	///pointer to element of interTmp
	int res=1;	///pointer to next element of inter that will contain

	int x,y;

	tipoInt offIDList=1;	///this pointer indicate where the next pointer of the nodeList begins in nodeList
	for (tipoInt m=1; m<=nodeList[0]; m++){

		patNode=addressOf[nodeList[offIDList]];	///node pointed by current

		tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode

		///now find the position of the coreIndex in the node
		tipoInt CI=0;
		while ((&patNode->item)[CI]!=father->info->coreIndex){
			CI++;
		}

		if(nodeList[offIDList+3]==-1){

			///the node will not have a pntToNode because it also contains the minItem, so the intersecton associated with the pointer to the node is the same memorized in the node itself

			tipoInt numInter=(tipoInt) (patNode->lastInter-patNode->lastItem);	///number of items in the intersection (that is memorized in the node)

			if (nextfree+numItems+numInter>lastvalid){
        cout << "call manual mem 1" << endl;
				resizeManualMem();
			}
			transTmp=nextfree;
			nextfree=nextfree+numItems+numInter;	///not "...+1;" becuase coreIndex will not stay in the intersection

			///now set the intersection: we now that all the items in the intersection are not in father because the node also contains the minItem of father
			for(int n=1; n<=numInter;n++){
				transTmp[n]=(patNode->lastItem)[numInter-n+1];
			}

			transTmp[0]=numInter;
			///current item ID
			tipoInt index;

			///now insert in the intersection all the items in the node that are not in P
			x=numItems-1;	///pointer to a item in the node
			y=1;		///poiter to a item in father->info->prefixCI
			while (x>CI && y <= father->info->prefixCI[0]){
				if ((&patNode->item)[x]==father->info->prefixCI[y]){
					///the item is in the closure, so we must skip it
					y++;
					x--;
				}
				else{
					if ((&patNode->item)[x]>father->info->prefixCI[y]){
						///the item is not in the closure, so we must add it in intersection of pntToNode
						transTmp[0]++;
						transTmp[transTmp[0]]=(&patNode->item)[x];
						x--;
					}
					else{
						///we are not sure that the item is not in the closure
						y++;
					}
				}
			}

			///now we must add all the item in the patNode that we have not compared with the closure in transTmp
			while(x>CI){
				transTmp[0]++;
				transTmp[transTmp[0]]=(&patNode->item)[x];
				x--;
			}

			///we skip the coreIndex because it will not stay in the intersection
			for (int l=CI-1; l>=0; l--){
				index=(&patNode->item)[l];

				///update count
				HT[index].count+=nodeList[offIDList+1];
				/* SignificantMiner */
				HT[index].a_S+=nodeList[offIDList+2];
				/* TODO */

				/*cout << "HT["<<index<<"].count="<<HT[index].count<<endl;
				cout << "HT["<<index<<"].a_S="<<HT[index].a_S<<endl;
				cout << "nodeList[offIDList+1]="<<nodeList[offIDList+1]<<endl;
				cout << "nodeList[offIDList+2]="<<nodeList[offIDList+2]<<endl;*/


				//temp_time = get_cpu_time();
				/*if(addressOf[nodeList[offIDList]]->a_S[0]!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
				//cout << "Error with nodelist!(1a)"<<endl;
				//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
				//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
				//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
				//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
				//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
					//cout << "using father nl_a_S "<<endl;
					int offset = (m-1)*(jp+1);
					if(father->info->nl_a_S_type==8){
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_8[offset+aj];
					}
					else if(father->info->nl_a_S_type==16){
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_16[offset+aj];
					}
					else{
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_32[offset+aj];
					}
				}
				else{
					//cout << "using node nl_a_S "<<endl;
					for(int aj=0; aj<jp+1; aj++)
					HT[index].a_S[aj]+=addressOf[nodeList[offIDList]]->a_S[aj];
				}*/
				//time_keeper += get_cpu_time() - temp_time;


				if(DEEP_DEBUG){
					cout << "nodelist: increased count (1) for patricia node of " << nodeList[offIDList+1] << " to " << HT[index].count << endl;
					cout << "nodelist: updated a_S of " << nodeList[offIDList+2] << " to " << HT[index].a_S << endl;
				}
				/* SignificantMiner */

				///update the intersection in HT
				if (HT[index].intersection!=NULL){

					///update intersection: we intersect the transactions and memorize the intersection in the "inter" transaction
					i=1;	///pointer to elements of HT[l].intersection
					j=1;	///pointer to element of transLast
					res=1;	///pointer to next element of intersection that will contain an item in the intersection
					while(i<=HT[index].intersection[0]&&j<=transTmp[0]){
						if (HT[index].intersection[i]==transTmp[j]){
							///intersection[i] is in the intersection
							HT[index].intersection[res]=HT[index].intersection[i];
							res++;
							i++;
							j++;
						}
						else{
							if(HT[index].intersection[i]>transTmp[j]){
								///HT.intersection[i] is not in the intersection
								i++;
							}
							else{
								///HT.intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					HT[index].intersection[0]=res-1;
				}
				else{
					///set the intersection
					if (nextfree+transTmp[0]+1>lastvalid){
            cout << "call manual mem 2" << endl;
						resizeManualMem();
					}
					HT[index].intersection=nextfree;
					nextfree=nextfree+transTmp[0]+1;
					for (int o=0;o<=transTmp[0]; o++){
						HT[index].intersection[o]=transTmp[o];
					}
				}

				///double visitedList if necessary
				if(HT[index].visitedList[0].IDnode==HT[index].visitedList[0].numInInter){
					tipoInt newsizelist=HT[index].visitedList[0].numInInter*2;
					HT[index].visitedList=(VisitedNodeElem*) realloc(HT[index].visitedList, (1+newsizelist)*sizeof(VisitedNodeElem));
					HT[index].visitedList[0].numInInter=newsizelist;
				}

				///update visitedList
				tipoInt numElem=HT[index].visitedList[0].IDnode;
				HT[index].visitedList[numElem+1].offPntNode=-1;
				HT[index].visitedList[numElem+1].IDnode=patNode->ID;
				HT[index].visitedList[numElem+1].numInInter=-1;		///so we know that the intersection that must be considered (for future fillHT) is in the node
				HT[index].visitedList[0].IDnode++;

				///update transLast
				transTmp[0]++;
				transTmp[transTmp[0]]=index;
			}

			///now transTmp contains all the intersection that must be passed to the father of patNode

			///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
			PatriciaNode* upper=patNode->father;	///upper is the father of patNode
			if (upper!=NULL){	///the node isn't a root's children
				if (upper->pntToNode==NULL){

					tipoInt numItemUpper=(tipoInt) 1+(upper->lastItem-&(upper->item));	///number of items in patNode

					///create the pointer to the node
					if ((nextfree+sizePnt+numItemUpper+transTmp[0])>lastvalid){
            cout << "call manual mem 3" << endl;
						resizeManualMem();
					}
					upper->pntToNode=(NodePointer*)nextfree;
					//todelete->push_back(upper->pntToNode);
					living_nodepointers++;
					nextfree=nextfree+numItemUpper+transTmp[0]+sizePnt;	///update nextfree
					upper->pntToNode->count=nodeList[offIDList+1];
					upper->pntToNode->count_a_S=nodeList[offIDList+2];
					/* SignificantMiner */
					//temp_time = get_cpu_time();

					//upper->pntToNode->a_S_index = getNewTodeleteIndex();
					upper->pntToNode->nl_list_index = getNewListIndex();

				/*cout << "here ok IL? ins4 " << endl;
				cout << p->IL[1].item << endl;*/

					//cout << "(7) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

          int temp_as_value = (addressOf[nodeList[offIDList]]->count >= 128) ? addressOf[nodeList[offIDList]]->a_S_[0] : addressOf[nodeList[offIDList]]->a_S[0];

					if(temp_as_value!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
					//cout << "Error with nodelist!(1)"<<endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
					//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
					//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
					//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
						//int offset = (m-1)*(jp+1);
						//cout << "access (11) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
						//temp_as = todelete->at(upper->pntToNode->a_S_index);
						//cout << "using father nl_a_S "<<endl;

						nl_lists[upper->pntToNode->nl_list_index].clear();
						int upper_index = (m < father->info->minnodes_indexes.size()) ? father->info->minnodes_indexes[m] : father->info->minnodes.size();
						nl_lists[upper->pntToNode->nl_list_index].reserve(upper_index - father->info->minnodes_indexes[m-1]);
						for(int aj = father->info->minnodes_indexes[m-1]; aj < upper_index; aj++){
							nl_lists[upper->pntToNode->nl_list_index].push_back(father->info->minnodes[aj]);
						}
						//nl_lists[upper->pntToNode->nl_list_index] = father->info->nl_nodes[m-1];

						/*if(father->info->nl_a_S_type==8){
						for(int aj=0; aj<jp+1; aj++)
							temp_as[aj]=father->info->nl_a_S.nl_a_S_8[offset+aj];
						}
						else if(father->info->nl_a_S_type==16){
							for(int aj=0; aj<jp+1; aj++)
								temp_as[aj]=father->info->nl_a_S.nl_a_S_16[offset+aj];
						}
						else{
							for(int aj=0; aj<jp+1; aj++)
								temp_as[aj]=father->info->nl_a_S.nl_a_S_32[offset+aj];
						}*/

					}
					else{
						//cout << "access (12) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
						//temp_as = todelete->at(upper->pntToNode->a_S_index);

						nl_lists[upper->pntToNode->nl_list_index].clear();
						nl_lists[upper->pntToNode->nl_list_index].push_back(nodeList[offIDList]);

						//cout << "using node nl_a_S "<<endl;
						/*for(int aj=0; aj<jp+1; aj++)
							temp_as[aj]=addressOf[nodeList[offIDList]]->a_S[aj]; */
					}
					/*cout << "done 5.2 " << nextfree_a_S << endl;

					for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/


					//time_keeper += get_cpu_time() - temp_time;

					/*if(DEEP_DEBUG){
						cout << "updated count (2) to " << nodeList[offIDList+1] << endl;
						cout << "a_S updated to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;
					}*/
					/* SignificantMiner */
					upper->pntToNode->node=upper->ID;
					upper->pntToNode->nextNodePointer=-1;
					///we copy transTmp into the new intersection
					(&upper->pntToNode->intersection)[0]=transTmp[0];
					for (int o=1;o<=transTmp[0]; o++){
						(&upper->pntToNode->intersection)[o]=transTmp[o];
					}
					///insert the node in the pointerList of its greater item
					*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
					///update lastNodePnt
					HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
				}
				else{
					///we must update the field of the pointer to upper
					upper->pntToNode->count+=nodeList[offIDList+1];
					upper->pntToNode->count_a_S+=nodeList[offIDList+2];
					/* SignificantMiner */
					//temp_time = get_cpu_time();

          int temp_as_value = (addressOf[nodeList[offIDList]]->count >= 128) ? addressOf[nodeList[offIDList]]->a_S_[0] : addressOf[nodeList[offIDList]]->a_S[0];

					if(temp_as_value!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
					//cout << "Error with nodelist!(1b)"<<endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
					//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
					//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
					//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
						//int offset = (m-1)*(jp+1);
						//temp_as = todelete->at(upper->pntToNode->a_S_index);
						//cout << "access (14) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
						//cout << "using father nl_a_S "<<endl;

						int upper_index = (m < father->info->minnodes_indexes.size()) ? father->info->minnodes_indexes[m] : father->info->minnodes.size();
						for(int aj = father->info->minnodes_indexes[m-1]; aj < upper_index; aj++){
							nl_lists[upper->pntToNode->nl_list_index].push_back(father->info->minnodes[aj]);
						}

						/*for(int aj=0; aj < father->info->nl_nodes[m-1].size(); aj++){
							nl_lists[upper->pntToNode->nl_list_index].push_back(father->info->nl_nodes[m-1][aj]);
						}*/

						/*if(father->info->nl_a_S_type==8){
						for(int aj=0; aj<jp+1; aj++)
							temp_as[aj]+=father->info->nl_a_S.nl_a_S_8[offset+aj];
							//temp_as[aj]+=father->info->nl_a_S[offset+aj];
						}
						else if(father->info->nl_a_S_type==16){
							for(int aj=0; aj<jp+1; aj++)
								temp_as[aj]+=father->info->nl_a_S.nl_a_S_16[offset+aj];
						}
						else{
							for(int aj=0; aj<jp+1; aj++)
								temp_as[aj]+=father->info->nl_a_S.nl_a_S_32[offset+aj];
						}*/

					}
					else{
						//cout << "access (15) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
						//temp_as = todelete->at(upper->pntToNode->a_S_index);

						nl_lists[upper->pntToNode->nl_list_index].push_back(nodeList[offIDList]);

						//cout << "using node nl_a_S "<<endl;
						/*for(int aj=0; aj<jp+1; aj++)
							temp_as[aj]+=addressOf[nodeList[offIDList]]->a_S[aj]; */
					}
					//time_keeper += get_cpu_time() - temp_time;


				/*cout << "here ok IL? ins6 " << endl;
				cout << p->IL[1].item << endl;*/

					/*if(DEEP_DEBUG){
						cout << "nodelist: updated count of " << nodeList[offIDList+1] << " to " << upper->pntToNode->count << endl;
						cout << "nodelist: updated a_S of " << nodeList[offIDList+2] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;
					}*/
					/* SignificantMiner */
					///as intersection we can use transTmp
					i=1;	///pointer to elements of upper->pntToNode->intersection
					j=1;	///pointer to element of transTmp
					res=1;	///pointer to next element of upper->pntToNode->intersection that will contain an item in the intersection
					while(i<=(&upper->pntToNode->intersection)[0]&&j<=transTmp[0]){
						if ((&upper->pntToNode->intersection)[i]==transTmp[j]){
										///upper->pntToNode->intersection[i] is in the intersection
							(&upper->pntToNode->intersection)[res]=(&upper->pntToNode->intersection)[i];
							res++;
							i++;
							j++;
						}
						else{
							if((&upper->pntToNode->intersection)[i]>transTmp[j]){
								///upper->pntToNode->intersection[i] is not in the intersection
								i++;
							}
							else{
								///upper->pntToNode->intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					(&upper->pntToNode->intersection)[0]=res-1;
				}
			}
			///update the offIDList
			/* SignificantMiner */
			tipoInt nf = 4; // number of nodelist fields: 3 + 1 (added a_S)
			/* SignificantMiner */
			offIDList=offIDList+nf; ///we now that the length of the intersection is 0 (as number of tipoInt that memorize the items in the intersection associated to the pointer)
		}

		else{
			///the node will have a pntToNode
			///we must create a NodePointer to this node because there can be some item that is in the node and will be the CI for a ppce -> we must remember the intersection of this node
			tipoInt numInter=nodeList[offIDList+3];	///number of items in the intersection (that is memorized in the nodeList)

			if (nextfree+numItems+numInter+sizePnt>lastvalid){
        cout << "call manual mem 4" << endl;
				resizeManualMem();
			}
			patNode->pntToNode=(NodePointer*)nextfree;
			//todelete->push_back(patNode->pntToNode);
			living_nodepointers++;
			nextfree=nextfree+numItems+numInter+sizePnt;	///not "...+1;" because coreIndex will not stay in the intersection
			patNode->pntToNode->count=nodeList[offIDList+1];
			patNode->pntToNode->count_a_S=nodeList[offIDList+2];
			/* SignificantMiner */
			//temp_time = get_cpu_time();

			//patNode->pntToNode->a_S_index = getNewTodeleteIndex();
			patNode->pntToNode->nl_list_index = getNewListIndex();

			//cout << "(8) allocated jp+1 things at " << todelete->at(patNode->pntToNode->a_S_index) << endl;

      int temp_as_value = (addressOf[nodeList[offIDList]]->count >= 128) ? addressOf[nodeList[offIDList]]->a_S_[0] : addressOf[nodeList[offIDList]]->a_S[0];

			if(temp_as_value!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
				//cout << "Error with nodelist!(1)"<<endl;
				//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
				//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
				//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
				//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
				//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
				//int offset = (m-1)*(jp+1);
				//cout << "access (16) to index patNode->pntToNode->a_S_index = " << patNode->pntToNode->a_S_index << endl;
				//temp_as = todelete->at(patNode->pntToNode->a_S_index);
				//cout << "using father nl_a_S "<<endl;

				nl_lists[patNode->pntToNode->nl_list_index].clear();
				int upper_index = (m < father->info->minnodes_indexes.size()) ? father->info->minnodes_indexes[m] : father->info->minnodes.size();
				nl_lists[patNode->pntToNode->nl_list_index].reserve(upper_index - father->info->minnodes_indexes[m-1]);
				for(int aj = father->info->minnodes_indexes[m-1]; aj < upper_index; aj++){
					nl_lists[patNode->pntToNode->nl_list_index].push_back(father->info->minnodes[aj]);
				}

				//nl_lists[patNode->pntToNode->nl_list_index] = father->info->nl_nodes[m-1];

				/*if(father->info->nl_a_S_type==8){
					for(int aj=0; aj<jp+1; aj++)
						temp_as[aj]=father->info->nl_a_S.nl_a_S_8[offset+aj];
						//temp_as[aj]=father->info->nl_a_S[offset+aj];
				}
				else if(father->info->nl_a_S_type==16){
					for(int aj=0; aj<jp+1; aj++)
						temp_as[aj]=father->info->nl_a_S.nl_a_S_16[offset+aj];
				}
				else{
					for(int aj=0; aj<jp+1; aj++)
						temp_as[aj]=father->info->nl_a_S.nl_a_S_32[offset+aj];
				}*/

			}
			else{
				//cout << "access (17) to index patNode->pntToNode->a_S_index = " << patNode->pntToNode->a_S_index << endl;
				//temp_as = todelete->at(patNode->pntToNode->a_S_index);

				nl_lists[patNode->pntToNode->nl_list_index].push_back(nodeList[offIDList]);

				//cout << "using node nl_a_S "<<endl;
				/*for(int aj=0; aj<jp+1; aj++)
					temp_as[aj]=addressOf[nodeList[offIDList]]->a_S[aj];*/
			}
			/*cout << "done 6.2 " << nextfree_a_S << endl;

			for(int aj=0; aj<jp+1; aj++)

				cout << todelete->at(patNode->pntToNode->a_S_index)[aj] << " ";
			cout << endl;*/

			//time_keeper += get_cpu_time() - temp_time;

			/*if(DEEP_DEBUG){
				cout << "updated count (1) to " << nodeList[offIDList+1] << endl;
				//cout << "a_S was " << todelete->at(patNode->pntToNode->a_S_index) << endl;
				cout << "updated a_S to " << todelete->at(patNode->pntToNode->a_S_index) << endl;
			}*/
			/* SignificantMiner */
			patNode->pntToNode->node=patNode->ID;
			patNode->pntToNode->nextNodePointer=-1;

			///now set the intersection
			for(int n=1; n<=numInter;n++){
				(&patNode->pntToNode->intersection)[n]=nodeList[offIDList+3+n];
			}
			(&patNode->pntToNode->intersection)[0]=numInter;

			///now insert in the intersection all the items in the node that are not in P
			x=numItems-1;	///pointer to a item in the node
			y=1;		///poiter to a item in father->info->prefixCI
			while (x>CI && y <= father->info->prefixCI[0]){
				if ((&patNode->item)[x]==father->info->prefixCI[y]){
					///the item is in the closure, so we must skip it
					y++;
					x--;
				}
				else{
					if ((&patNode->item)[x]>father->info->prefixCI[y]){
						///the item is not in the closure, so we must add it in intersection of pntToNode
						(&patNode->pntToNode->intersection)[0]++;
						(&patNode->pntToNode->intersection)[(&patNode->pntToNode->intersection)[0]]=(&patNode->item)[x];
						x--;
					}
					else{
						///we are not sure that the item is not in the closure
						y++;
					}
				}
			}

			///now we must add all the item in the patNode that we have not compared with the closure in transTmp
			while(x>CI){
				(&patNode->pntToNode->intersection)[0]++;
				(&patNode->pntToNode->intersection)[(&patNode->pntToNode->intersection)[0]]=(&patNode->item)[x];
				x--;
			}

			///current item ID
			tipoInt index;

			///we must begin from the CI-1 item because we can update intersection
			for (int l=CI-1; l>=0; l--){
				index=(&patNode->item)[l];

				///update count
				HT[index].count+=nodeList[offIDList+1];
				HT[index].a_S+=nodeList[offIDList+2];


				/*cout << "HT["<<index<<"].count="<<HT[index].count<<endl;
				cout << "HT["<<index<<"].a_S="<<HT[index].a_S<<endl;
				cout << "nodeList[offIDList+1]="<<nodeList[offIDList+1]<<endl;
				cout << "nodeList[offIDList+2]="<<nodeList[offIDList+2]<<endl;*/

				/* SignificantMiner */

        int temp_as_value = (addressOf[nodeList[offIDList]]->count >= 128) ? addressOf[nodeList[offIDList]]->a_S_[0] : addressOf[nodeList[offIDList]]->a_S[0];

				//temp_time = get_cpu_time();
				if(temp_as_value!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
					//cout << "Error with nodelist!(1c)"<<endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
					//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
					//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
					//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
					//int offset = (m-1)*(jp+1);
					//cout << "using father nl_a_S "<<endl;


					/*if(father->info->nl_a_S_type==8){
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_8[offset+aj];
							//HT[index].a_S[aj]+=father->info->nl_a_S[offset+aj];
					}
					else if(father->info->nl_a_S_type==16){
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_16[offset+aj];
					}
					else{
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_32[offset+aj];
					}*/

				}
				else{
					//cout << "using node nl_a_S "<<endl;
					/*for(int aj=0; aj<jp+1; aj++)
						HT[index].a_S[aj]+=addressOf[nodeList[offIDList]]->a_S[aj]; */
				}
				//time_keeper += get_cpu_time() - temp_time;

				if(DEEP_DEBUG){
					cout << "nodelist: increased count (2) for patricia node of " << nodeList[offIDList+1] << " to " << HT[index].count << endl;
					cout << "nodelist: increased a_S (2) for patricia node of " << nodeList[offIDList+2] << " to " << HT[index].a_S << endl;
				}
				/* SignificantMiner */

				///update the intersection in HT
				if (HT[index].intersection!=NULL){
					///update intersection: we intersect the transactions and memorize the intersection in the "inter" transaction
					i=1;	///pointer to elements of HT[l].intersection
					j=1;	///pointer to element of transLast
					res=1;	///pointer to next element of intersection that will contain an item in the intersection
					while(i<=HT[index].intersection[0]&&j<=(&patNode->pntToNode->intersection)[0]){
						if (HT[index].intersection[i]==(&patNode->pntToNode->intersection)[j]){
							///intersection[i] is in the intersection
							HT[index].intersection[res]=HT[index].intersection[i];
							res++;
							i++;
							j++;
						}
						else{
							if(HT[index].intersection[i]>(&patNode->pntToNode->intersection)[j]){
								///intersection[i] is not in the intersection
								i++;
							}
							else{
								///intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					HT[index].intersection[0]=res-1;
				}
				else{
					///set the intersection
					if (nextfree+(&patNode->pntToNode->intersection)[0]+1>lastvalid){
            cout << "call manual mem 4" << endl;
						resizeManualMem();
					}
					HT[index].intersection=nextfree;
					nextfree=nextfree+(&patNode->pntToNode->intersection)[0]+1;
					for (int o=0;o<=(&patNode->pntToNode->intersection)[0]; o++){
						HT[index].intersection[o]=(&patNode->pntToNode->intersection)[o];
					}
				}

				///double visitedList if necessary

				if(HT[index].visitedList[0].IDnode==HT[index].visitedList[0].numInInter){
					tipoInt newsizelist=HT[index].visitedList[0].numInInter*2;
					HT[index].visitedList=(VisitedNodeElem*) realloc(HT[index].visitedList, (1+newsizelist)*sizeof(VisitedNodeElem));
					HT[index].visitedList[0].numInInter=newsizelist;
				}

				///update visitedList and lastVisitedNode
				tipoInt numElem=HT[index].visitedList[0].IDnode;
				HT[index].visitedList[0].IDnode++;
				HT[index].visitedList[numElem+1].offPntNode=((tipoInt*)patNode->pntToNode)-manualmem;
				HT[index].visitedList[numElem+1].IDnode=patNode->ID;

				HT[index].visitedList[numElem+1].numInInter=numInter;	///only items not in the node must be memorized

				///update transTmp
				(&patNode->pntToNode->intersection)[0]++;
				(&patNode->pntToNode->intersection)[(&patNode->pntToNode->intersection)[0]]=index;
			}

			///now transTmp contains all the intersection that must be passed to the father of patNode

			///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
			PatriciaNode* upper=patNode->father;	///upper is the father of patNode
			if (upper!=NULL){	///the node isn't a root's children
				if (upper->pntToNode==NULL){
					tipoInt numItemUpper=(tipoInt) 1+(upper->lastItem-&(upper->item));	///number of items in patNode
					///create the pointer to the node
					if ((nextfree+sizePnt+numItemUpper+(&patNode->pntToNode->intersection)[0])>lastvalid){
            cout << "call manual mem 5" << endl;
						resizeManualMem();
					}
					upper->pntToNode=(NodePointer*)nextfree;
					//todelete->push_back(upper->pntToNode);
					living_nodepointers++;
					nextfree=nextfree+numItemUpper+(&patNode->pntToNode->intersection)[0]+sizePnt;	///update nextfree
					upper->pntToNode->count=patNode->pntToNode->count;
					upper->pntToNode->count_a_S=patNode->pntToNode->count_a_S;
					/* SignificantMiner */
					//temp_time = get_cpu_time();
					//upper->pntToNode->a_S_index = getNewTodeleteIndex();
					upper->pntToNode->nl_list_index = getNewListIndex();

					//cout << "(9) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

					//cout << "access (17) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (17) to index upper->pntToNode->a_S_index = " << patNode->pntToNode->a_S_index << endl;
					/*temp_as = todelete->at(upper->pntToNode->a_S_index);
					temp_as_ = todelete->at(patNode->pntToNode->a_S_index);
					for(int aj=0; aj<jp+1; aj++)
					temp_as[aj]=temp_as_[aj];*/

					nl_lists[upper->pntToNode->nl_list_index] = nl_lists[patNode->pntToNode->nl_list_index];


					/*cout << "done 7.2 " << nextfree_a_S << endl;
					for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/

					//time_keeper += get_cpu_time() - temp_time;
					/*if(DEEP_DEBUG){
						cout << "count set to " << patNode->pntToNode->count << endl;
						cout << "a_S set to " << todelete->at(patNode->pntToNode->a_S_index)[0] << endl;
					}*/
					/* SignificantMiner */
					upper->pntToNode->node=upper->ID;
					upper->pntToNode->nextNodePointer=-1;
					///we copy transTmp into the new intersection
					(&upper->pntToNode->intersection)[0]=(&patNode->pntToNode->intersection)[0];
					for (int o=1;o<=(&patNode->pntToNode->intersection)[0]; o++){
						(&upper->pntToNode->intersection)[o]=(&patNode->pntToNode->intersection)[o];
					}
					///insert the node in the pointerList of its greater item
					*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
					///update lastNodePnt
					HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
				}
				else{
					///we must update the field of the pointer to upper
					upper->pntToNode->count+=patNode->pntToNode->count;
					upper->pntToNode->count_a_S+=patNode->pntToNode->count_a_S;
					/* SignificantMiner */
					//temp_time = get_cpu_time();
					//cout << "access (18) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (18) to index upper->pntToNode->a_S_index = " << patNode->pntToNode->a_S_index << endl;
					/*temp_as = todelete->at(upper->pntToNode->a_S_index);
					temp_as_ = todelete->at(patNode->pntToNode->a_S_index);
					for(int aj=0; aj<jp+1; aj++)
					temp_as[aj]+=temp_as_[aj];*/

					for(int aj=0; aj<nl_lists[patNode->pntToNode->nl_list_index].size(); aj++){
						nl_lists[upper->pntToNode->nl_list_index].push_back(nl_lists[patNode->pntToNode->nl_list_index][aj]);
					}

					//time_keeper += get_cpu_time() - temp_time;
					//if(DEEP_DEBUG)
					//cout << "increased a_S (9) for patricia node of " << todelete->at(patNode->pntToNode->a_S_index)[0] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;
					/* SignificantMiner */
					///as intersection we can use transTmp
					i=1;	///pointer to elements of upper->pntToNode->intersection
					j=1;	///pointer to element of transTmp
					res=1;	///pointer to next element of upper->pntToNode->intersection that will contain an item in the intersection
					while(i<=(&upper->pntToNode->intersection)[0]&&j<=(&patNode->pntToNode->intersection)[0]){
						if ((&upper->pntToNode->intersection)[i]==(&patNode->pntToNode->intersection)[j]){
										///upper->pntToNode->intersection[i] is in the intersection
							(&upper->pntToNode->intersection)[res]=(&upper->pntToNode->intersection)[i];
							res++;
							i++;
							j++;
						}
						else{
							if((&upper->pntToNode->intersection)[i]>(&patNode->pntToNode->intersection)[j]){
											///upper->pntToNode->intersection[i] is not in the intersection
								i++;
							}
							else{
											///upper->pntToNode->intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					(&upper->pntToNode->intersection)[0]=res-1;
				}
			}

			///delete the reference to the pointer into patNode
			patNode->pntToNode=NULL;

			/* SignificantMiner */
			tipoInt nf = 4; // number of nodelist fields: 3 + 1 (added a_S)
			/* SignificantMiner */
			///update the offIDList
			offIDList=offIDList+nodeList[offIDList+3]+nf; ///nodeList[offIDList+2]= length of the intersection, 3 is number of field that are always present (ID, count and length of intersection)

		}

	}


	///now we must consider all the list of items greater than father->info->minItem and visit the element in their pointerList; we must begin from the last item (the greatest)

	NodePointer* currentpnt;	///temporary NodePointer used for the elements in the pointedList
	tipoInt off;			///offset of currentpnt in manualmem

	for (tipoInt l=father->info->coreIndex-1;l>=0;l--){
		off=HT[l].pointerList;
		while(off!=-1){
			currentpnt=(NodePointer*)&(manualmem[off]);
			patNode=addressOf[currentpnt->node];	///PatriciaNode pointed by currentpnt
			tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode
			tipoInt numInter=(&currentpnt->intersection)[0];

			///now consider all the items in the node and update the corresponding HT rows
			tipoInt index;	///current item ID
			for (tipoInt s=numItems-1;s>=0; s--){
				index=(&patNode->item)[s];

				///we know that all the items in the node are all lower than father->info->minItem
				HT[index].count+=currentpnt->count;	///update the count with the count of the pointer
				HT[index].a_S+=currentpnt->count_a_S;


				/*cout << "HT["<<index<<"].count="<<HT[index].count<<endl;
				cout << "HT["<<index<<"].a_S="<<HT[index].a_S<<endl;
				cout << "currentpnt->count="<<currentpnt->count<<endl;
				cout << "currentpnt->count_a_S="<<currentpnt->count_a_S<<endl;*/


				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//cout << "access (18) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				/*temp_as = todelete->at(currentpnt->a_S_index);
				for(int aj=0; aj<jp+1; aj++)
				HT[index].a_S[aj]+=temp_as[aj];	*/
				//time_keeper += get_cpu_time() - temp_time;
				/*if(DEEP_DEBUG){
					cout << "l = " << l << " off = " << off << endl;
					cout << "increased count (10) for patricia node of " << currentpnt->count << " to " << HT[index].count << endl;
					cout << "increased a_S (10) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << HT[index].a_S[0] << endl;
				}*/
				/* SignificantMiner */
				///we update intersection and visitedList only of item<coreindex

				///update the intersection in HT
				if (HT[index].intersection!=NULL){
					///update intersection: we intersect the transactions and memorize the intersection in the "inter" transaction
					i=1;	///pointer to elements of HT[l].intersection
					j=1;	///pointer to element of transLast
					res=1;	///pointer to next element of intersection that will contain an item in the intersection
					while(i<=HT[index].intersection[0]&&j<=(&currentpnt->intersection)[0]){
						if (HT[index].intersection[i]==(&currentpnt->intersection)[j]){
							///intersection[i] is in the intersection
							HT[index].intersection[res]=HT[index].intersection[i];
							res++;
							i++;
							j++;
						}
						else{
							if(HT[index].intersection[i]>(&currentpnt->intersection)[j]){
								///intersection[i] is not in the intersection
								i++;
							}
							else{
								///intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					HT[index].intersection[0]=res-1;
				}
				else{
					///set the intersection
					if (nextfree+(&currentpnt->intersection)[0]+1>lastvalid){
            cout << "call manual mem 6" << endl;
						resizeManualMem();
					}
					HT[index].intersection=nextfree;
					nextfree=nextfree+(&currentpnt->intersection)[0]+1;
					for (int o=0;o<=(&currentpnt->intersection)[0];o++){
						HT[index].intersection[o]=(&currentpnt->intersection)[o];
					}
				}

				///double visitedList if necessary

				if(HT[index].visitedList[0].IDnode==HT[index].visitedList[0].numInInter){
					tipoInt newsizelist=HT[index].visitedList[0].numInInter*2;
					HT[index].visitedList=(VisitedNodeElem*) realloc(HT[index].visitedList, (1+newsizelist)*sizeof(VisitedNodeElem));
					HT[index].visitedList[0].numInInter=newsizelist;
				}

				///update visitedList and lastVisitedNode
				tipoInt numElem=HT[index].visitedList[0].IDnode;
				HT[index].visitedList[numElem+1].offPntNode=((tipoInt*)patNode->pntToNode)-manualmem;
				HT[index].visitedList[numElem+1].IDnode=patNode->ID;
				HT[index].visitedList[numElem+1].numInInter=numInter;
				HT[index].visitedList[0].IDnode++;

				///here there is no need to verify if the item is different from coreIndex because we are not considering nodes that contains minItem
				(&currentpnt->intersection)[0]++;
				(&currentpnt->intersection)[(&currentpnt->intersection)[0]]=index;
			}

			///now transTmp contains all the intersection that must be passed to the pointer to the father of patNode

			///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
			PatriciaNode* upper=patNode->father;	///upper is the father of patNode
			if (upper!=NULL){	///the node isn't a root's children
				if (upper->pntToNode==NULL){
				///create the pointer to the node
				tipoInt numItemUpper=(tipoInt) 1+(upper->lastItem-&(upper->item));
				if (nextfree+sizePnt+numItemUpper+(&currentpnt->intersection)[0]>lastvalid){
          cout << "call manual mem 7" << endl;
					resizeManualMem();
				}
				upper->pntToNode=(NodePointer*)nextfree;
				//todelete->push_back(upper->pntToNode);
				living_nodepointers++;
				nextfree=nextfree+numItemUpper+(&currentpnt->intersection)[0]+sizePnt;	///update nextfree
				upper->pntToNode->count=currentpnt->count;
				upper->pntToNode->count_a_S=currentpnt->count_a_S;
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//upper->pntToNode->a_S_index = getNewTodeleteIndex();
				upper->pntToNode->nl_list_index = getNewListIndex();

				//cout << "(9) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

				//cout << "access (19) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				//cout << "access (19) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				/*temp_as = todelete->at(upper->pntToNode->a_S_index);
				temp_as_ = todelete->at(currentpnt->a_S_index);
				for(int aj=0; aj<jp+1; aj++)
				temp_as[aj]=temp_as_[aj];*/

				for(int aj=0; aj<nl_lists[currentpnt->nl_list_index].size(); aj++){
					nl_lists[upper->pntToNode->nl_list_index].push_back(nl_lists[currentpnt->nl_list_index][aj]);
				}

			/*cout << "done 8.2 " << nextfree_a_S << endl;
			for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
				//time_keeper += get_cpu_time() - temp_time;
				/*if(DEEP_DEBUG){
					cout << "count set to " << upper->pntToNode->count << endl;
					cout << "a_S set to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;
				}*/
				/* SignificantMiner */
				upper->pntToNode->node=upper->ID;
				upper->pntToNode->nextNodePointer=-1;
				///we copy transTmp into the new intersection
				(&upper->pntToNode->intersection)[0]=(&currentpnt->intersection)[0];
				for (int o=1;o<=(&currentpnt->intersection)[0]; o++){
					(&upper->pntToNode->intersection)[o]=(&currentpnt->intersection)[o];
				}
				///insert the node in the pointerList of its greater item
				*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
				///update lastNodePnt
				HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
				///create the pointer to the node
				}
				else{
					///we must update the field of the pointer to upper
					upper->pntToNode->count+=currentpnt->count;
					upper->pntToNode->count_a_S+=currentpnt->count_a_S;
					/* SignificantMiner */
					//temp_time = get_cpu_time();

					//cout << "access (20) to index upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (20) to index currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
					/*temp_as = todelete->at(upper->pntToNode->a_S_index);
					temp_as_ = todelete->at(currentpnt->a_S_index);
					for(int aj=0; aj<jp+1; aj++)
					temp_as[aj]+=temp_as_[aj];*/

					for(int aj=0; aj<nl_lists[currentpnt->nl_list_index].size(); aj++){
						nl_lists[upper->pntToNode->nl_list_index].push_back(nl_lists[currentpnt->nl_list_index][aj]);
					}

					//time_keeper += get_cpu_time() - temp_time;
					/*if(DEEP_DEBUG){
						cout << "increased count (11) for patricia node of " << currentpnt->count << " to " << upper->pntToNode->count << endl;
						cout << "increased a_S (11) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;
					}*/
					/* SignificantMiner */
					///as intersection we can use transTmp
					i=1;	///pointer to elements of upper->pntToNode->intersection
					j=1;	///pointer to element of transTmp
					res=1;	///pointer to next element of upper->pntToNode->intersection that will contain an item in the intersection
					while(i<=(&upper->pntToNode->intersection)[0]&&j<=(&currentpnt->intersection)[0]){
						if ((&upper->pntToNode->intersection)[i]==(&currentpnt->intersection)[j]){
										///upper->pntToNode->intersection[i] is in the intersection
							(&upper->pntToNode->intersection)[res]=(&upper->pntToNode->intersection)[i];
							res++;
							i++;
							j++;
						}
						else{
							if((&upper->pntToNode->intersection)[i]>(&currentpnt->intersection)[j]){
											///upper->pntToNode->intersection[i] is not in the intersection
								i++;
							}
							else{
											///upper->pntToNode->intersection[i] can still be in the intersection
								j++;
							}
						}
					}
					(&upper->pntToNode->intersection)[0]=res-1;
				}
			}
			off=currentpnt->nextNodePointer;
			patNode->pntToNode=NULL;
		}
	}
}

void testPattern(double psi_bound_log , Itemset* ppce, Itemset* father){

	/*cout << "testPattern call, ppce = " << ppce << " father = " << father << endl;
	cout << "testPattern call, ppce->support = " << ppce->support << " ppce->info->max_a_S = " << ppce->info->max_a_S << " ppce->info->min_a_S = " << ppce->info->min_a_S << endl;
	if(father) cout << "testPattern call, father->support = " << father->support << " father->info->max_a_S = " << father->info->max_a_S << " father->info->min_a_S = " << father->info->min_a_S << endl;*/

	/*cout << "minnodes: " << endl;
	for(int i_=0; i_ < ppce->info->minnodes.size(); i_++){
		cout << ppce->info->minnodes[i_] << " " << endl;
	}
	cout << endl;
	cout << "minnodes_indexes: " << endl;
	for(int i_=0; i_ < ppce->info->minnodes_indexes.size(); i_++){
		cout << ppce->info->minnodes_indexes[i_] << " " << endl;
	}
	cout << endl;*/

	// reset perm_supp
	/*for(int aj=0; aj<perm_supp.size(); aj++)
		perm_supp[aj]=0;*/

	#ifdef trace_state
	program_state = -121;
	#endif

    double temp_as_time = get_cpu_time();

	// compute the support on the permutations of current ppce
	int support_debug = 0;
	tested++;
	/*for(int i_=0; i_ < ppce->info->nl_nodes.size(); i_++){
		//cout << "Node  " << i_ << " points to: " << endl;
		for(int j=0; j <  ppce->info->nl_nodes[i_].size(); j++){
			//cout << " " << j <<" - node  " << ppce->info->nl_nodes.at(i_).at(j) << endl;
			if(addressOf[ ppce->info->nl_nodes[i_][j] ]->count >= 128){
				for(int aj=1;aj<(jp+1); aj++){
					perm_supp[aj] += addressOf[ ppce->info->nl_nodes[i_][j] ]->a_S_[aj];
				}
			}
			else{
				for(int aj=1;aj<(jp+1); aj++){
					perm_supp[aj] += addressOf[ ppce->info->nl_nodes[i_][j] ]->a_S[aj];
				}
			}
			support_debug += addressOf[ ppce->info->nl_nodes[i_][j] ]->count;
		}
	}*/



  tipoInt *pattern_a_S_;
  uint8_t *pattern_a_S;

  // first element
  int pattern_support = addressOf[ ppce->info->minnodes[0] ]->count;
  if(pattern_support >= 128){
    pattern_a_S_ = addressOf[ ppce->info->minnodes[0] ]->a_S_;
    for(int aj=1;aj<(jp+1); aj++){
      perm_supp[aj] = pattern_a_S_[aj];
    }
  }
  else{
    pattern_a_S = addressOf[ ppce->info->minnodes[0] ]->a_S;
    for(int aj=1;aj<(jp+1); aj++){
      perm_supp[aj] = pattern_a_S[aj];
    }
  }
  support_debug += pattern_support;

  // loop
	for(int i_=1; i_ < ppce->info->minnodes.size(); i_++){
    int pattern_support = addressOf[ ppce->info->minnodes[i_] ]->count;
		if(pattern_support >= 128){
      pattern_a_S_ = addressOf[ ppce->info->minnodes[i_] ]->a_S_;
			for(int aj=1;aj<(jp+1); aj++){
				perm_supp[aj] += pattern_a_S_[aj];
			}
		}
		else{
      pattern_a_S = addressOf[ ppce->info->minnodes[i_] ]->a_S;
			for(int aj=1;aj<(jp+1); aj++){
				perm_supp[aj] += pattern_a_S[aj];
			}
		}
		support_debug += pattern_support;
	}

	as_processing_time += get_cpu_time() - temp_as_time;

	if(support_debug != ppce->support){
	cout << "ERROR! something wrong with nllists!" << endl;
	cout << "support_debug = " << support_debug << endl;
	cout << "ppce->support = " << ppce->support << endl;

	cout << "minnodes: " << endl;
	for(int i_=0; i_ < ppce->info->minnodes.size(); i_++){
		cout << ppce->info->minnodes[i_] << " " << endl;
	}
	cout << endl;
	cout << "minnodes_indexes: " << endl;
	for(int i_=0; i_ < ppce->info->minnodes_indexes.size(); i_++){
		cout << ppce->info->minnodes_indexes[i_] << " " << endl;
	}
	cout << endl;

	abort();
	}

  						//cout << "ppce with support "<<ppce->support << " has as: " << endl;
  						//for(int aj=0;aj<(jp+1); aj++){
  							//cout << " " << perm_supp[aj];
  						//}
  						//cout << endl;


              #ifdef trace_state
    					program_state = -122;
    					#endif


  						computePValues(ppce->support , perm_supp.data()); // torna qui


              #ifdef enable_bounds
            	ppce->info->max_a_S = 0;
            	ppce->info->min_a_S = ppce->support;
              for(int aj=0; aj<(jp+1); aj++){
                ppce->info->max_a_S = max(ppce->info->max_a_S , perm_supp[aj]);
              	ppce->info->min_a_S = min(ppce->info->min_a_S , perm_supp[aj]);
              }
              #endif



              #ifdef trace_state
    					program_state = -123;
    					#endif

  						// update min and max a_S of the ppce across permutations
    						/*int max_non_alpha = 0;
    						int min_non_alpha = ppce->support; */

              #ifndef enable_bounds
              ppce->info->min_a_S = 0;
    					ppce->info->max_a_S = ppce->support;
              #endif


  						/*if(max_non_alpha < ppce->info->max_a_S || min_non_alpha > ppce->info->min_a_S){
  						cout << "New max and min computed for support : " << ppce->support << endl;
  						cout << "   max and min computed on all indexes : " << ppce->info->min_a_S << " " << ppce->info->max_a_S << endl;
  						cout << "   max and min computed on non-alpha indexes : " << min_non_alpha << " " << max_non_alpha << endl;}*/

              #ifdef trace_state
    					program_state = -124;
    					#endif

  						/*if(ppce->info->min_a_S > ppce->info->max_a_S){
  							cout << "ERROR HERE 1" << endl;
  							cout << "ppce->support " << ppce->support << endl;
  							cout << "ppce->info->min_a_S " << ppce->info->min_a_S << endl;
  							cout << "ppce->info->max_a_S " << ppce->info->max_a_S << endl;
  						}*/

  						// update the pruning condition
              bool can_have_FP_child = true;
              #ifdef enable_pruning
  						can_have_FP_child = canHaveFPChild(ppce->support , ppce->info->min_a_S , ppce->info->max_a_S) != -1;
              #endif


              #ifdef enable_bounds
              #ifdef testexpanded
              if(father){
  						// update the bounds on the father itemset
  						if(ppce->info->min_a_S > father->info->min_a_S){
  							/*if(ppce->info->min_a_S < 49){
  							cout << "updated father min_a_S: \nfather->support "
  							<< father->support << " old min " << father->info->min_a_S << " new min " <<
  							ppce->info->min_a_S << endl;}*/

  							father->info->min_a_S = ppce->info->min_a_S;
  						}
  						if(ppce->info->max_a_S + (father->support - ppce->support) < father->info->max_a_S){
  							/*if(ppce->info->min_a_S < 49){
  							cout << "updated father max_a_s: \nfather->support "
  							<< father->support << " old max " << father->info->max_a_S << " new max " <<
  							ppce->info->max_a_S << endl;}*/

  							father->info->max_a_S = ppce->info->max_a_S + (father->support - ppce->support);
  						}

  						if(father->info->min_a_S > father->info->max_a_S){
  							cout << "ERROR HERE 2" << endl;
  							cout << "father->support " << father->support << endl;
  							cout << "father->info->min_a_S " << father->info->min_a_S << endl;
  							cout << "father->info->max_a_S " << father->info->max_a_S << endl;
  							cout << "ppce->support " << ppce->support << endl;
  							cout << "ppce->info->min_a_S " << ppce->info->min_a_S << endl;
  							cout << "ppce->info->max_a_S " << ppce->info->max_a_S << endl;
  							abort();
  						}
              }
              #endif
              #endif

}

/**find the closure aftre coreIndex of the (closed) itemset from which the HT was filled and insert its (frequent) ppc-e in the queue; the father's conditional dataset is represented by coreIndex nodes; also the coreIndex is inserted in the transaction returned
*/

//Transaction findPpceCI(Itemset* father, QueueMinMax* q){
Transaction findPpceCI(Itemset* father, typeQueue q , typeQueueRes q_res){

	//cout << "FINDPPCECI CALL" << endl;

	Transaction closure;	///here we insert all the items that are in Clo(father) and that are <= CI
	if (nextfree+father->info->coreIndex+2>lastvalid){
    cout << "call manual mem 8" << endl;
		resizeManualMem();
	}
	closure=nextfree;
	nextfree=nextfree+father->info->coreIndex+2;
	closure[0]=1;
	closure[1]=father->info->coreIndex;
	///now we check all the rows of HT to find the items in Clo(father) and the ppce of father

	///we must begin from i=minItem-1 because so we can update intersection deleting the items in Clo(father)
	Itemset* ppce;	///used to create new Itemset to insert in max queue

	Itemset* minPpce; ///used to insert the ppce in the min queue

	bool can_have_significant_child, can_have_FP_child, prune_current_ppce;
  prune_current_ppce = false;

	for (int i=father->info->coreIndex-1; i>=0; i--){
		///if the item has count<suppMin, there is nothing to do
		if (HT[i].count>=suppMin){
			///if the item has the same count of father, it is in closure=Clo(father)
			if (HT[i].count==father->support){
				///insert i in closure
				closure[0]++;
				closure[closure[0]]=i;
			}
			else{
				///if it isn't in closure we must check if it can generate a ppce: conditions on coreindex and prefix-preservation

				///to verify if this extension is prefix-preserving, we check if the intersection in HT[i].intersection contains only items in closure-minItem; this is real IFF the length are the same (because HT[i].intersection surely contains the item in closure that are lower than i)
				if(HT[i].intersection[0]==closure[0]-1){


          #ifdef trace_state
					program_state = -23;
					#endif

          /*cout << "-23" << endl;
          cout << "ppce->support " << ppce->support << endl;
          cout << "ppce->info->a_S_ " << ppce->info->a_S_ << endl;
          cout << "ppce->info->p_value " << ppce->info->p_value << endl;
          cout << "getPsi(suppMin) " << getPsi(suppMin) << endl;
          cout << "getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;
          cout << "computeLogPValue(ppce->support , ppce->info->a_S_) " << computeLogPValue(ppce->support , ppce->info->a_S_) << endl;*/

          int ppce_support = HT[i].count;
          int ppce_a_S = HT[i].a_S;

          bool can_be_significant = true;

          #ifdef enable_bounds
          can_be_significant = getLogPsiExact(ppce_support) <= getLogPsi(suppMin);
          #endif

          bool observed = false;

          //cout << "ppce->support " << ppce_support << " ppce->info->a_S_ " << ppce_a_S << endl;
          // if the extracted itemset can be significant
          if(can_be_significant){
            double log_p_value = computeLogPValue(ppce_support , ppce_a_S);
          if( store_results && log_p_value <= getLogPsi(suppMin) ){
            //cout << "ppce p_value: " << log_p_value << " ppce support " << ppce_support <<" ppce a_S " << ppce_a_S << endl;
            //cout << "RES QUEUE OBSERVED INSERTION: it can be significant since getLogPsi(suppMin) = " << getLogPsi(suppMin) << endl;
            //double log_p_value = computeLogPValue(ppce_support , ppce_a_S);
            //cout << "RES QUEUE OBSERVED INSERTION: it can be significant since log_p_value = " << log_p_value << endl;

            //cout << "   pattern log_p_value " << log_p_value << endl;
            //cout << "   getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;

            //cout << "toOutput_pattern->log_p_value = " << toOutput_pattern->log_p_value << endl;

            //q_res->printLowestHeapLevel();
            // insert it in the queue
            q_res->insert_observed(log_p_value);
            observed = true;

                          /* Top-K strategy */
                          while(K_significant_patterns > 0 && q_res->getKElements() >= K_significant_patterns){
                            int new_significance_level = q_res->getMaxIndex();
                            new_significance_level++;
                            if(new_significance_level > suppMin){
                              suppMin = new_significance_level;
                              suppMinCurr = new_significance_level;
                            }

                            cout << "Top-K strategy update! ";
                            /*cout << "K_significant_patterns " << K_significant_patterns << endl;
                            cout << "q_res->getKElements() returned " << q_res->getKElements() << endl;
                            cout << "q_res->getMaxIndex() returned " << q_res->getMaxIndex() << endl;*/
                            cout << " s_supp " << s_supp;
                            cout << " new suppMin " << suppMin << endl;
                            /*cout << "getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;
                            cout << "getPsi(suppMin) " << getPsi(suppMin) << endl;*/

                            delta = getPsi(suppMin);
                            q_res->removeMax();
                          }
            //cout << "inserted!" << endl;
            //q_res->printLowestHeapLevel();

          }

      		}

          #ifdef enable_pruning

          if(!observed){

          #ifdef trace_state
          program_state = -161;
          #endif

          #ifdef LAMP
					getExactPsiBound(HT[i].count , father->support , father->info->a_S_);
					#endif

          //cout << "HT[i].count of ppce = " << HT[i].count << endl;

					// before creating everything for the ppce, test if we can prune it
					can_have_significant_child = canHaveSignificantChild(HT[i].count , HT[i].a_S) != -1;

					#ifdef trace_state
					program_state = -162;
					#endif

					int max_a_S_ = min(HT[i].count , father->info->max_a_S);
					int min_a_S_ = max(0 , HT[i].count - (father->support - father->info->min_a_S));

					/*if(father->info->max_a_S < max_a_S_)
						max_a_S_ = father->info->max_a_S;

					if(HT[i].count - (father->support - father->info->min_a_S) > min_a_S_)
						min_a_S_ = HT[i].count - (father->support - father->info->min_a_S);*/

					can_have_FP_child = canHaveFPChild(HT[i].count , min_a_S_ , max_a_S_) != -1;

          }

					prune_current_ppce = !observed && !can_have_significant_child && !can_have_FP_child;

					pruned_patterns+=prune_current_ppce;

					#ifdef trace_state
					program_state = -162;
					#endif

          #endif

					if(!prune_current_ppce){

					#ifdef trace_state
					program_state = -163;
					#endif

          //cout << "non pruned "<<endl;


					///i is the core index of a ppc-e of father
					ppce=(Itemset*)calloc(1,sizeof(Itemset));
					//.cout << "created itemset (4) with adress "<< ppce <<endl;
					living_itemsets++;
					ppce->support=HT[i].count;
					ppce->info=new ItsInfo();
					ppce->info->a_S_=HT[i].a_S;
					living_infos++; if(living_infos > maxlivinginfos) maxlivinginfos = living_infos;
					/* SignificantMiner */
					//temp_time = get_cpu_time();
					//ppce->info->a_S=(tipoInt*)malloc((jp+1)*sizeof(tipoInt));

					//cout << "(10) allocated jp+1 things at " << ppce->info->a_S << endl;

					//for(int aj=0; aj<jp+1; aj++)
					//ppce->info->a_S[aj]=HT[i].a_S[aj];


					//time_keeper += get_cpu_time() - temp_time;
					/* SignificantMiner */
					ppce->info->minItem=father->info->minItem;
					ppce->info->coreIndex=i;

					///now we must find the NodeIDList for the ppc-e

					///first of all, we scan the visitedList to find the size of the memory that must be allocated for the NodeIDList

					tipoInt sizeList=0;
					for (int j=1; j<=HT[i].visitedList[0].IDnode;j++){
						if (HT[i].visitedList[j].numInInter!=-1){
							sizeList+=HT[i].visitedList[j].numInInter;	///for the intersection of the node;
						}
					}
					/* SignificantMiner */
					tipoInt nf = 4; // number of nodelist fields: 3 + 1 (added a_S)
					/* SignificantMiner */
					//sizeList=sizeList+3*HT[i].visitedList[0].IDnode;	///for the ID of the node, the count associated and the length of every intersection
					sizeList=sizeList+nf*HT[i].visitedList[0].IDnode;	///for the ID of the node, the count associated and the length of every intersection
					sizeList++;	///for the length of the list

					//if (sizeList*sizeof(tipoInt)>maxListSize){
					//if (ppce->info->nl_a_S_size>maxNLSize){
          if(use_secondary_memory){
					queue_ram = (double)(q->getInqueue() * ( sizeof(Itemset) + sizeof(ItsInfo) )) / 1000000.0;
					//in_use_ram = (double)(nlinqueue_8 * (double)sizeof(uint8_t) + nlinqueue_16 * (double)sizeof(uint16_t) + nlinqueue_32 * (double)sizeof(uint32_t)) / 1000000.0 * (double)(jp+1);
					//double uncompressed_ram = (double)(nlinqueue_8 + nlinqueue_16 + nlinqueue_32) / 1000000.0 * (double)(jp+1) * (double)sizeof(tipoInt);
					//saved = in_use_ram / uncompressed_ram * 100.0;
					queue_ram += (list_length_inqueue * 6 * sizeof(tipoInt))/ 1000000.0;
					queue_ram += (nl_nodes_length_inqueue * sizeof(tipoInt))/ 1000000.0;
					in_use_ram = 0.0;
					in_use_ram += pat_tree_ram + queue_ram;
					/*if(in_use_ram > memory_peak)
						memory_peak = in_use_ram;*/
          memory_peak = max(memory_peak , in_use_ram);
					bool new_allocation = true;
					if ( use_secondary_memory && (in_use_ram > maxNLRAM || (in_use_ram > max_ram_l1 && ppce->info->nl_a_S_size > maxNLLen) || (in_use_ram > max_ram_l2 && ppce->info->nl_a_S_size > midNLLen)))
						new_allocation = false;
          }
					/*if ( in_use_ram > maxNLRAM || ppce->info->nl_a_S_size > maxNLLen){
						cout<<"store in file 1!"<<endl;


						cout << nlinqueue << endl;
						cout << (jp+1) << endl;
						cout << sizeof(tipoInt) << endl;
						cout << (double)sizeof(tipoInt) << endl;
						cout << in_use_ram << endl;


						ppce->info->fileName=inserted;
						ppce->info->typeList=2;
					}		*/

					///now create the NodeIDList that must be memorized
					ppce->info->IDlist=(tipoInt*)malloc((sizeList)*sizeof(tipoInt));
					ppce->info->sizelista=sizeList;
					/* SignificantMiner */
					//temp_time = get_cpu_time();
					//int are_zeros = 0;
					ppce->info->nl_a_S_size=0;//HT[i].visitedList[0].IDnode;
					ppce->info->nl_a_S.nl_a_S_t = NULL;
					ppce->info->nl_a_S_type = 0;

					#ifdef trace_state
					program_state = -164;
					#endif

					/*ppce->info->nl_nodes.resize(HT[i].visitedList[0].IDnode);
					if(ppce->info->nl_nodes.size() > max_nl_size_){
			            max_nl_size_ = ppce->info->nl_nodes.size();
			            cout << " new max_nl_size_ = " << max_nl_size_ << endl;
					}
					if(ppce->info->nl_nodes.size() < 0){
			            cout << "ERROR WITH ppce->info->nl_nodes.size(): " << ppce->info->nl_nodes.size() << endl;
			        }*/


					std::vector< std::vector<int> > temp_nl_list(HT[i].visitedList[0].IDnode);
			        ppce->info->minnodes.reserve(HT[i].visitedList[0].IDnode);
			        ppce->info->minnodes_indexes.reserve(HT[i].visitedList[0].IDnode);
					if(HT[i].visitedList[0].IDnode > max_nl_size_){
			            max_nl_size_ = HT[i].visitedList[0].IDnode;
			            cout << " new max_nl_size_ = " << max_nl_size_ << endl;
					}
					if(HT[i].visitedList[0].IDnode < 0){
			            cout << "ERROR WITH HT[i].visitedList[0].IDnode: " << HT[i].visitedList[0].IDnode << endl;
			        }


					/*if(ppce->support < max_uint_8)
						ppce->info->nl_a_S_type = 8;
					else
						if(ppce->support < max_uint_16)
							ppce->info->nl_a_S_type = 16;
						else
							ppce->info->nl_a_S_type = 32;*/

					//cout << "current ppce with support " << ppce->support << " is set to type " << (int)ppce->info->nl_a_S_type << endl;
					//cout << "father s= " << father->support << " list sz= " << father->info->nl_a_S_size << endl;
					//cout << "gen ppce s= " << ppce->support << " list sz= " << ppce->info->nl_a_S_size << endl;

					#ifdef trace_state
					program_state = -16;
					#endif

					/*if(new_allocation){
						if(ppce->info->nl_a_S_type == 8)
							ppce->info->nl_a_S.nl_a_S_8=(uint8_t*)malloc((HT[i].visitedList[0].IDnode * (jp+1)) * sizeof(uint8_t));
						else
							if(ppce->info->nl_a_S_type == 16)
								ppce->info->nl_a_S.nl_a_S_16=(uint16_t*)malloc((HT[i].visitedList[0].IDnode * (jp+1)) * sizeof(uint16_t));
							else
								ppce->info->nl_a_S.nl_a_S_32=(uint32_t*)malloc((HT[i].visitedList[0].IDnode * (jp+1)) * sizeof(uint32_t));
						nlasmalloc++;
						//cout << "allocated nl_a_S of size " << ppce->info->nl_a_S_size << " at " << ppce->info->nl_a_S << endl;
					}
					else{
						if(ppce->info->nl_a_S_type == 8){
							while(temp_size_8 < ppce->info->nl_a_S_size){
								temp_size_8 = temp_size_8*2;
								free(temp_nl_a_S_8);
								temp_nl_a_S_8=(uint8_t*)malloc(temp_size_8*(jp+1)*sizeof(uint8_t));
							}
							ppce->info->nl_a_S.nl_a_S_8=temp_nl_a_S_8;
						}
						else if(ppce->info->nl_a_S_type == 16){
							cout << " temp_size_16 " << temp_size_16 << endl;
							cout << " ppce->info->nl_a_S_size " << ppce->info->nl_a_S_size << endl;
							cout << " temp_nl_a_S_16 " << temp_nl_a_S_16 << endl;
							while(temp_size_16 < ppce->info->nl_a_S_size){
								temp_size_16 = temp_size_16*2;
								free(temp_nl_a_S_16);
								temp_nl_a_S_16=(uint16_t*)malloc(temp_size_16*(jp+1)*sizeof(uint16_t));
								cout << " temp_size_16 " << temp_size_16 << endl;
								cout << " temp_size_16 " << ppce->info->nl_a_S_size << endl;
								cout << " temp_nl_a_S_16 " << temp_nl_a_S_16 << endl;
							}
							ppce->info->nl_a_S.nl_a_S_16=temp_nl_a_S_16;
						}
						else if(ppce->info->nl_a_S_type == 32){
							while(temp_size_32 < ppce->info->nl_a_S_size){
								temp_size_32 = temp_size_32*2;
								free(temp_nl_a_S_32);
								temp_nl_a_S_32=(uint32_t*)malloc(temp_size_32*(jp+1)*sizeof(uint32_t));
							}
							ppce->info->nl_a_S.nl_a_S_32=temp_nl_a_S_32;
						}

					}*/

					#ifdef trace_state
					program_state = -17;
					#endif

					//nodelist_inqueue+=HT[i].visitedList[0].IDnode;
					//time_keeper += get_cpu_time() - temp_time;
					/*if(DEEP_DEBUG){
					cout << "created nl_a_S (1) for ppce with support: " << ppce->support << endl;
					cout << "created nl_a_S (1) with size: " << HT[i].visitedList[0].IDnode << " * " << (jp+1) << " = " << (HT[i].visitedList[0].IDnode * (jp+1)) << endl;
					}*/
					/* SignificantMiner */
					int ne=1;	///next position for insertion in nodeIDList
					NodePointer* pnttmp;
					for (int j=1; j<=HT[i].visitedList[0].IDnode;j++){
						ppce->info->IDlist[ne]=HT[i].visitedList[j].IDnode;	///ID of the node
						ne++;
						if(HT[i].visitedList[j].offPntNode==-1){
							ppce->info->IDlist[ne]=addressOf[HT[i].visitedList[j].IDnode]->count;
							/* SignificantMiner */
							ne++;
							int temp_as_value = (addressOf[HT[i].visitedList[j].IDnode]->count >= 128) ? addressOf[HT[i].visitedList[j].IDnode]->a_S_[0] : addressOf[HT[i].visitedList[j].IDnode]->a_S[0];
							ppce->info->IDlist[ne]=temp_as_value;
							if(DEEP_DEBUG){
								cout << "created a_S (1) as " << ppce->info->IDlist[ne] << endl;
							}
							tipoInt index = 0;
							tipoInt offset = (j-1)*(jp+1);
							//temp_time = get_cpu_time();

							//ppce->info->nl_nodes[j-1].push_back(HT[i].visitedList[j].IDnode);
							temp_nl_list[j-1].push_back(HT[i].visitedList[j].IDnode);

							/*if(ppce->info->nl_a_S_type==8){
								for(int aj=0; aj<jp+1; aj++){
									index = offset + aj;
									ppce->info->nl_a_S.nl_a_S_8[index] = addressOf[HT[i].visitedList[j].IDnode]->a_S[aj];
									//index = offset + aj;
									//ppce->info->nl_a_S[index] = addressOf[HT[i].visitedList[j].IDnode]->a_S[aj];
								}
							}
							else if(ppce->info->nl_a_S_type==16){
								for(int aj=0; aj<jp+1; aj++){
									index = offset + aj;
									ppce->info->nl_a_S.nl_a_S_16[index] = addressOf[HT[i].visitedList[j].IDnode]->a_S[aj];
								}
							}
							else{
								for(int aj=0; aj<jp+1; aj++){
									index = offset + aj;
									ppce->info->nl_a_S.nl_a_S_32[index] = addressOf[HT[i].visitedList[j].IDnode]->a_S[aj];
								}
							}*/
							//cout << endl<<endl;
							//time_keeper += get_cpu_time() - temp_time;
							//cout << "allocated  nl_a_S (1)" << endl;
							/* SignificantMiner */
						}
						else{
							pnttmp=(NodePointer*)&(manualmem[HT[i].visitedList[j].offPntNode]);	///nodePointer of the visited node
							ppce->info->IDlist[ne]=pnttmp->count;	///for the count of the pointer
							/* SignificantMiner */
							ne++;

							//cout << "access (21) to pnttmp->a_S_index = " << pnttmp->a_S_index << endl;
							ppce->info->IDlist[ne]=pnttmp->count_a_S;
							if(DEEP_DEBUG){
								cout << "created a_S (2) as " << ppce->info->IDlist[ne] << endl;
							}
							//tipoInt index = 0;
							//tipoInt offset = (j-1)*(jp+1);
							//temp_time = get_cpu_time();
							//cout << "todelete->at(pnttmp->a_S_index)[aj]="<<endl;

              int *temp_prt = nl_lists[pnttmp->nl_list_index].data();
							for(int aj=0; aj < nl_lists[pnttmp->nl_list_index].size(); aj++){
								//ppce->info->nl_nodes[j-1].push_back(nl_lists[pnttmp->nl_list_index][aj]);
								temp_nl_list[j-1].push_back(temp_prt[aj]);
							}

							/*if(ppce->info->nl_a_S_type==8){
								for(int aj=0; aj<jp+1; aj++){
									index = offset + aj;
									ppce->info->nl_a_S.nl_a_S_8[index] = todelete->at(pnttmp->a_S_index)[aj];
									//index = offset + aj;
									//ppce->info->nl_a_S[index] = todelete->at(pnttmp->a_S_index)[aj];
								}
							}
							else if(ppce->info->nl_a_S_type==16){
								for(int aj=0; aj<jp+1; aj++){
									index = offset + aj;
									ppce->info->nl_a_S.nl_a_S_16[index] = todelete->at(pnttmp->a_S_index)[aj];
								}
							}
							else{
								for(int aj=0; aj<jp+1; aj++){
									index = offset + aj;
									ppce->info->nl_a_S.nl_a_S_32[index] = todelete->at(pnttmp->a_S_index)[aj];
								}
							}*/
							//cout <<endl<<endl;
							//time_keeper += get_cpu_time() - temp_time;
							//cout << "allocated  nl_a_S (2)" << endl;
							/* SignificantMiner */
						}
						ne++;
						ppce->info->IDlist[ne]=HT[i].visitedList[j].numInInter;
						ne++;
						///for the length of the intersection
						///now copy the intersection: we must remove all the items in the intersection that are in the closure of the father and that are <	CI of father and < Ci of ppce (these items are surely in the intersection)
						tipoInt inseriti=0;
						if (closure[0]==1){///there are no items that must be erased
							for (int s=0;s<HT[i].visitedList[j].numInInter;s++){
								ppce->info->IDlist[ne+s]=(&pnttmp->intersection)[s+1];
							}
						}
						else{///we must erase something
							tipoInt jb=2;	///pointer in the closure of father: we begin from the first item that is not the CI of the father
							inseriti=0;
							for (int s=0;s<HT[i].visitedList[j].numInInter;s++){
								if ((&pnttmp->intersection)[s+1]!=closure[jb]){
									ppce->info->IDlist[ne+inseriti]=(&pnttmp->intersection)[s+1];
									inseriti++;
								}
								else{
									if (jb<closure[0])
										jb++;
									ppce->info->IDlist[ne-1]--;
								}
							}
						}
						if (HT[i].visitedList[j].offPntNode!=-1){
							///update ne if the numInInter is >=0
							if(closure[0]==1){
								ne=ne+ppce->info->IDlist[ne-1];
							}
							else {ne=ne+inseriti;}
						}
					}

					/*if(false){
					double save = ((double)are_zeros / (double)((jp+1)*ppce->info->nl_a_S_size));
					cout << " new list. supp " << ppce->support << " size " << ppce->info->nl_a_S_size << ". zeros " << are_zeros << " save " << save << endl;
        }*/

					#ifdef trace_state
					program_state = -18;
					#endif

					ppce->info->IDlist[0]=HT[i].visitedList[0].IDnode;
					ppce->info->typeList=1;
					//if (sizeList*sizeof(tipoInt)>maxListSize){
					//if (ppce->info->nl_a_S_size>maxNLSize){
					if ( use_secondary_memory && (in_use_ram > maxNLRAM || (in_use_ram > max_ram_l1 && ppce->info->nl_a_S_size > maxNLLen) || (in_use_ram > max_ram_l2 && ppce->info->nl_a_S_size > midNLLen))){
					//if ( in_use_ram > maxNLRAM || ppce->info->nl_a_S_size > maxNLLen){

						ppce->info->fileName=inserted;

						#ifdef trace_state
						program_state = -19;
						#endif


						//cout<<"store in file 2!"<<endl;
						cout<<"WR " << inserted ;

						//cout<<"nlinqueue " << nlinqueue <<endl;
						//cout<<"nlasmalloc " << nlasmalloc <<endl;
						cout<<" s_supp " << s_supp ;
						//cout<<"inqueue " << inqueue <<endl;
						cout<<" suppMin " << suppMin ;
						//cout<<"produced " << produced <<endl;
						cout<<" ppce->support " << ppce->support;
						cout<<" ppce->info->nl_a_S_size " << ppce->info->nl_a_S_size;
						cout<<" ppce->info->nl_a_S_type " << (int)ppce->info->nl_a_S_type;
						cout<<" ram " << in_use_ram <<endl;

						//cout<<"filename: " << ppce->info->fileName <<endl;
						ppce->info->typeList=2;
						FILE * pFile;

						for(int aj=0; aj<20; aj++) file_path[aj]=' ';
						sprintf(file_path, "%d.l", ppce->info->fileName);

  						pFile = fopen ( file_path , "wb" );
						fwrite (&sizeList, sizeof(tipoInt), 1, pFile );
  						fwrite (ppce->info->IDlist , sizeof(tipoInt), sizeList, pFile );
  						if(ppce->info->nl_a_S_type==8){
  							fwrite (ppce->info->nl_a_S.nl_a_S_8, sizeof(uint8_t), ppce->info->nl_a_S_size*(jp+1), pFile );
  							//free(ppce->info->nl_a_S.nl_a_S_8);
  						}
  						else
  							if(ppce->info->nl_a_S_type==16){
  								fwrite (ppce->info->nl_a_S.nl_a_S_16, sizeof(uint16_t), ppce->info->nl_a_S_size*(jp+1), pFile );
  								//free(ppce->info->nl_a_S.nl_a_S_16);
  							}
  							else{
  								fwrite (ppce->info->nl_a_S.nl_a_S_32, sizeof(uint32_t), ppce->info->nl_a_S_size*(jp+1), pFile );
  								//free(ppce->info->nl_a_S.nl_a_S_32);
  							}
  						fclose (pFile);

  						/*for(int aj=0; aj<ppce->info->nl_a_S_size*(jp+1); aj++)
							cout << ppce->info->nl_a_S[aj] << " ";
						cout << endl;*/

						free(ppce->info->IDlist);
						//if(new_allocation){
						//	nlasmalloc--;
							//cout << "DEallocated nl_a_S of size " << ppce->info->nl_a_S_size << " at " << ppce->info->nl_a_S << endl;
						//	free(ppce->info->nl_a_S);
						//}
						ppce->info->nl_a_S.nl_a_S_t=NULL;
						ppce->info->IDlist=NULL;

						//cout<<" done writing " << endl ;

					}
					/*else{
						if(ppce->info->nl_a_S_type==8){
							nlinqueue_8+=ppce->info->nl_a_S_size;
						}
						else if(ppce->info->nl_a_S_type==16){
							nlinqueue_16+=ppce->info->nl_a_S_size;
						}
						else{
							nlinqueue_32+=ppce->info->nl_a_S_size;
						}
					}*/

					#ifdef trace_state
					program_state = -20;
					#endif
					/*if(nlinqueue > maxnlinqueue){
						maxnlinqueue = nlinqueue;
						//cout << "maxnlinqueue reached with suppMin " << suppMin << endl;
					}*/

					///now copy the CI prefix of father into prefixCI of the ppce
					ppce->info->prefixCI=(tipoInt*)malloc((father->info->prefixCI[0]+closure[0]+1)*sizeof(tipoInt));

					#ifdef trace_state
					program_state = -26;
					#endif

					for (int w=1;w<=father->info->prefixCI[0];w++){
						ppce->info->prefixCI[w]=father->info->prefixCI[w];
					}

					#ifdef trace_state
					program_state = -27;
					#endif

					for (int w=1;w<=closure[0];w++){
						ppce->info->prefixCI[w+father->info->prefixCI[0]]=closure[w];
					}
					ppce->info->prefixCI[0]=father->info->prefixCI[0]+closure[0];

					#ifdef trace_state
					program_state = -22;
					#endif

					/* Print nl_nodes of ppce */

					//cout << "-------------------------- " << endl;
					//cout << "Nodes of new ppce with support " << ppce->support << endl;
					/*int support_with_nl_nodes = 0;
					std::vector<int> a_S_with_nl_nodes((jp+1),0);
					for(int i=0; i < ppce->info->nl_nodes.size(); i++){
						//cout << "Node  " << i << " points to: " << endl;
						for(int j=0; j <  ppce->info->nl_nodes.at(i).size(); j++){
							//cout << " " << j <<" - node  " << ppce->info->nl_nodes.at(i).at(j) << endl;
							support_with_nl_nodes += addressOf[ppce->info->nl_nodes.at(i).at(j)]->count;
							for(int aj=0;aj<(jp+1); aj++){
								a_S_with_nl_nodes[aj] += addressOf[ppce->info->nl_nodes.at(i).at(j)]->a_S[aj];
							}
						}
					}
					//cout << "computed count: " << support_with_nl_nodes << endl;
					//cout << "-------------------------- " << endl;
					if(ppce->support != support_with_nl_nodes){
						cout << "ERROR! support does not match with node list sum!" << endl;
					}

					for(int aj=0;aj<(jp+1); aj++){
						if(HT[i].a_S[aj] != a_S_with_nl_nodes[aj])
							cout << "ERROR! a_S["<< aj <<"] does not match with node list sum!" << endl;
					}*/

					// compute size of minnodes
					int counter_ = 0;
					for(int i_=0; i_ < temp_nl_list.size(); i_++){
						counter_ = counter_ + temp_nl_list[i_].size();
					}

					ppce->info->minnodes.reserve(counter_);
					ppce->info->minnodes_indexes.reserve(temp_nl_list.size());
					// copy the list of min nodes
					ppce->info->minnodes_indexes.push_back(0);

          int *temp_prt;
					for(int i_=0; i_ < temp_nl_list.size(); i_++){
						if(i_ > 0) ppce->info->minnodes_indexes.push_back(ppce->info->minnodes_indexes[i_-1] + temp_nl_list[i_-1].size());
            temp_prt = temp_nl_list[i_].data();
						for(int j=0; j <  temp_nl_list[i_].size(); j++){
							ppce->info->minnodes.push_back(temp_prt[j]);

						}
					}

					/*cout << "temp_nl_list: " << endl;
					for(int i_=0; i_ < temp_nl_list.size(); i_++){
						for(int j=0; j <  temp_nl_list[i_].size(); j++){
							cout << temp_nl_list[i_][j] << " " << endl;
						}
					}
					cout << endl;

					cout << "minnodes: " << endl;
					for(int i_=0; i_ < ppce->info->minnodes.size(); i_++){
						cout << ppce->info->minnodes[i_] << " " << endl;
					}
					cout << endl;
					cout << "minnodes_indexes: " << endl;
					for(int i_=0; i_ < ppce->info->minnodes_indexes.size(); i_++){
						cout << ppce->info->minnodes_indexes[i_] << " " << endl;
					}
					cout << endl;*/

					#ifdef trace_state
					program_state = -12;
					#endif

					int x_S = ppce->support;
					int x_S_anchestor = father->support;
					int min_a_S_anchestor = father->info->min_a_S;
					int max_a_S_anchestor = father->info->max_a_S;

					int num_sum = 0;
					//int cutoff_nl_list_length = 50;

					// compute list statistics
					for(int h=0; h < temp_nl_list.size(); h++){
						//cout << "Node  " << h << " points to: " << endl;
						/*for(int j=0; j <  ppce->info->nl_nodes.at(h).size(); j++){
							//cout << " " << j <<" - node  " << ppce->info->nl_nodes.at(h).at(j) << endl;
							perm_supp[0] += addressOf[ppce->info->nl_nodes.at(h).at(j)]->a_S[0];
						}*/

						// list statistic
						//if(ppce->info->nl_nodes.at(h).size() < cutoff_nl_list_length){

						#ifdef print_list_stats
						if(list_stat.size() <= temp_nl_list[h].size()){
							list_stat.resize(temp_nl_list[h].size() + 1);
						}
						list_stat[temp_nl_list[h].size()]++;
						num_sum+=temp_nl_list[h].size();
						#endif


						nl_nodes_length_inqueue+=temp_nl_list[h].size();

					}



					list_length_inqueue+=ppce->info->sizelista;

					perm_supp[0] = ppce->info->a_S_;

					// number of operations statistic
          #ifdef print_list_stats
					if(list_op_stat.size() <= num_sum){
						list_op_stat.resize(num_sum + 1);
					}
					list_op_stat[num_sum]++;
          #endif

					//cout << "ppce with support "<<ppce->support << " has as: " << perm_supp[0] << endl;

					#ifdef trace_state
					program_state = -120;
					#endif

					ppce->info->p_value=1.0;
					double psi_bound_log = getLogPsi(ppce->support);

					if(can_be_significant){
            			ppce->info->p_value=computePValue(ppce->support , ppce->info->a_S_);
                  #ifdef enable_bounds
                  #ifdef testexpanded
            					psi_bound_log = getLogPsiBound(x_S , x_S_anchestor , min_a_S_anchestor , max_a_S_anchestor);

            					if(psi_bound_log >= getLogPsi(ppce->support)){
            						inferior_psi_bound_happened++;
            					}
                  #endif
                  #endif
					}

          #ifdef testexpanded
          if(psi_bound_log >= getLogPsi(suppMin)){
            inferior_psi_bound_happened_neglect++;
            // we do not update the as of the current ppce, we use the one of the anchestor
            ppce->info->max_a_S = min(max_a_S_anchestor , ppce->support);
            ppce->info->min_a_S = max(0 , min_a_S_anchestor - (x_S_anchestor - ppce->support));

            /*if(ppce->info->min_a_S > ppce->info->max_a_S){
              cout << "ERROR HERE 3" << endl;
              cout << "ppce->support " << ppce->support << endl;
              cout << "ppce->info->min_a_S " << ppce->info->min_a_S << endl;
              cout << "ppce->info->max_a_S " << ppce->info->max_a_S << endl;
            }*/
            //cout << "ppce with support " << ppce->support << " has psi_bound " << psi_bound << " vs standard " << getPsi(ppce->support) << endl;
            //cout << "father has support " << father->support << " has max_a_S " << max_a_S_anchestor << " and min " << min_a_S_anchestor << endl;
            //cout << "ppce with support " << ppce->support << " has max_a_S " << ppce->info->max_a_S << " and min " << ppce->info->min_a_S << endl;
          }
          else{
            testPattern(psi_bound_log , ppce , father);
          }
          #endif

          #ifndef testexpanded
          ppce->info->max_a_S = min(max_a_S_anchestor , ppce->support);
          ppce->info->min_a_S = max(0 , min_a_S_anchestor - (x_S_anchestor - ppce->support));
          #endif


					#ifdef trace_state
					program_state = -13;
					#endif

          #ifdef enable_pruning


					if(observed || can_have_significant_child || can_have_FP_child){
						q->insert(ppce);
					}
					else{
						pruned_patterns++;
						free(ppce->info->IDlist);
						ppce->info->IDlist = NULL;
						free(ppce->info->prefixCI);
						ppce->info->prefixCI=NULL;
						/*for(int aj=0; aj<ppce->info->nl_nodes.size(); aj++){
							std::vector<int> foo(0);
							ppce->info->nl_nodes[aj].swap(foo);
						}
						std::vector< std::vector<int> > foo(0);
						ppce->info->nl_nodes.swap(foo);*/

						{
						std::vector<int> foo1(0);
						std::vector<int> foo2(0);
						ppce->info->minnodes_indexes.swap(foo1);
						ppce->info->minnodes.swap(foo2);
						}

						delete ppce->info;
            //cout << "pattern is pruned!" << endl;
            if(observed){
              cout <<  "ERROR: one observed pattern has been pruned! " << endl;
            }
					}

          #endif

          #ifndef enable_pruning
          q->insert(ppce);
          #endif

          free(ppce);



          		#ifdef trace_state
          		program_state = -54;
          		#endif











					///create the Itemset for the min Queue: its ItsInfo is the same if the maxQueue
					//minPpce=(Itemset*)calloc(1,sizeof(Itemset));

					#ifdef trace_state
					program_state = -24;
					#endif

					//.cout << "created itemset (5) with adress "<< minPpce <<endl;
					//living_itemsets++;
					//minPpce->support=HT[i].count;
					//minPpce->info=ppce->info;
					/* SignificantMiner */
					/*temp_time = get_cpu_time();
					for(int aj=0; aj<jp+1; aj++)
					minPpce->info->a_S[aj]=HT[i].a_S[aj];
					//time_keeper += get_cpu_time() - temp_time;*/
					/* SignificantMiner */

					#ifdef trace_state
					program_state = -25;
					#endif


					//qMin->insert(minPpce);
					if(q->getInqueue() > maxinqueue) maxinqueue = q->getInqueue();

					///heuristic on suppMinCurr

					/*if (inserted<k-1){
						utilNodeArray[ppce->support]++;
					}
					else{
						if (inserted==k-1&&!firstBound) {
							///cumulate over the utilNodeArray
							for(tipoInt c=maximumSupport-1;c>=suppMin;c--){
								utilNodeArray[c]+=utilNodeArray[c+1];
							}
							firstBound=true;
						}
						///raise minimum support using utilNodeArray
						trivialHeurCurr(ppce->support);
					}*/
					inserted++;

					#ifdef trace_state
					program_state = -21;
					#endif

					///heuristic on suppMin

					/*itsOfSupp[ppce->support]++;

					if (produced+inqueue>=k_max){
					///we must check if we can have a new bound on suppMin
						if (k_max<inqueue+produced-itsOfSupp[suppMin]+1){
						///we have a new bound on suppMin and we can delete non frequent itemsets from the queue
							Itemset* todel;
							for (tipoInt t=0; t<itsOfSupp[suppMin]; t++){
								todel=qMin->removeMin();
								free(todel->info->prefixCI);
								if (todel->info->typeList!=2){
									free(todel->info->IDlist);
								}
								else{
									remove(todel->info->fileName);
								}
								delete todel->info;
								inqueue--;
								cout << "removed infrequent itemset with support " << todel->support << endl;
							}
							suppMin=qMin->heap->support;
						}
					}*/


					/* SignificantMiner */
					// check if infrequent itemsets can be removed from queue
					if(DEEP_DEBUG)
					cout << "checking to remove infrequent itemsets..." << endl;
					//Itemset* todel;
					//todel = qMin->getMin();
					//temp_time = get_cpu_time();
					//while(todel != NULL && todel->support < suppMin){
					while(q->getMin() > 0 && q->getMin() < suppMin && q->getMax() >= suppMin){
						//.cout << "free itemset (11) with adress "<< todel <<endl;
						//living_itemsets--;
						//free(todel);

						//todel=qMin->removeMin();
						//q->printLowestHeapLevel();
						q->removeMin();
						//free(todel->info->prefixCI);
						//free(todel->info->IDlist);
						//delete todel->info;
						//todel->info = NULL;
						//if(DEEP_DEBUG)
						//cout << "removed infrequent itemset with support " << todel->support << endl;

						//.cout << "free itemset (1) with adress "<< todel <<endl;
						//free(todel);
						//living_itemsets--;
						//todel = qMin->getMin();
					}
					//if(todel!=NULL){
						//.cout << "free itemset (12) with adress "<< todel <<endl;
						//living_itemsets--;
						//free(todel);
					//}
					//time_keeper += get_cpu_time() - temp_time;
					/* SignificantMiner */

					#ifdef trace_state
					program_state = -22;
					#endif
					//time_keeper += get_cpu_time() - temp_time;
					/* SignificantMiner */
				}

				}
			}
		}
	}
	return closure;
}

/** given an Itemset, find its closure (with the coreIndex inside) without generating its ppc-e; used in computation to produce closed itemset with the same support of the k-th when an itemset is memorized through coreIndex nodes;
*/

Transaction findCloCI(NodeIDList nodeList, Itemset* father){

	//cout << "findCloCI CALLED" << endl;

	///set the right lastVisitedNode pointer and lastNodePointer
	PatriciaNode* patNode;		///PatriciaNode pointed by current
	Transaction transTmp;		///maximum intersection recoverable by patNode and current intersection

	tipoInt offIDList=1;	///this pointer indicate where the next pointer of the nodeList begins in nodeList
	for (tipoInt m=1; m<=nodeList[0]; m++){

		patNode=addressOf[nodeList[offIDList]];	///node pointed by current

		tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode

		///now find the position of the coreIndex in the node
		tipoInt CI=numItems-1;
		while ((&patNode->item)[CI]!=father->info->coreIndex){
			CI--;
		}

		///the node will not have a pntToNode because it also contains the minItem, so the intersecton associated with the pointer to the node is the same memorized in the node itself

		tipoInt index;

		///we skip the coreIndex because it is yet know to be in the closure
		for (int l=CI-1; l>=0; l--){
			index=(&patNode->item)[l];
			///update count
			HT[index].count+=nodeList[offIDList+1];
			/* SignificantMiner */
			HT[index].a_S+=nodeList[offIDList+2];

			//temp_time = get_cpu_time();

			/*if(addressOf[nodeList[offIDList]]->a_S[0]!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
				//cout << "Error with nodelist!(11)"<<endl;
				//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
				//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
				//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
				//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
				//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
				int offset = (m-1)*(jp+1);

				if(father->info->nl_a_S_type==8){
					for(int aj=0; aj<jp+1; aj++)
						HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_8[offset+aj];
  				}
  				else
  					if(father->info->nl_a_S_type==16){
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_16[offset+aj];
  					}
  					else{
						for(int aj=0; aj<jp+1; aj++)
							HT[index].a_S[aj]+=father->info->nl_a_S.nl_a_S_32[offset+aj];
  					}

				//for(int aj=0; aj<jp+1; aj++)
					//HT[index].a_S[aj]+=father->info->nl_a_S[offset+aj];
			}
			else{
				for(int aj=0; aj<jp+1; aj++)
					HT[index].a_S[aj]+=addressOf[nodeList[offIDList]]->a_S[aj];
			}*/
			//time_keeper += get_cpu_time() - temp_time;

			if(DEEP_DEBUG){
				cout << "nodelist: increased count (3) for patricia node of " << nodeList[offIDList+1] << " to " << HT[index].count << endl;
				cout << "nodelist: increased a_S of " << nodeList[offIDList+2] << " to " << HT[index].a_S << endl;
			}
			/* SignificantMiner */
		}

		///now transTmp contains all the intersection that must be passed to the father of patNode
		///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
		PatriciaNode* upper=patNode->father;	///upper is the father of patNode
		if (upper!=NULL){	///the node isn't a root's children
			if (upper->pntToNode==NULL){

				///create the pointer to the node
				if ((nextfree+sizePnt)>lastvalid){
          cout << "call manual mem 9" << endl;
					resizeManualMem();
				}
				upper->pntToNode=(NodePointer*)nextfree;
				//todelete->push_back(upper->pntToNode);
				living_nodepointers++;
				nextfree=nextfree+sizePnt;	///update nextfree
				upper->pntToNode->count=nodeList[offIDList+1];
				/* SignificantMiner */
				upper->pntToNode->count_a_S=nodeList[offIDList+2];
				//upper->pntToNode->a_S_index = getNewTodeleteIndex();
				upper->pntToNode->nl_list_index = getNewListIndex();

				//cout << "(11) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

				//temp_time = get_cpu_time();

        int temp_as_value = (addressOf[nodeList[offIDList]]->count >= 128) ? addressOf[nodeList[offIDList]]->a_S_[0] : addressOf[nodeList[offIDList]]->a_S[0];

				if(temp_as_value!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
					//cout << "Error with nodelist!(12)"<<endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
					//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
					//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
					//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
					//int offset = (m-1)*(jp+1);
					//cout << "access (22) to upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;

					//nl_lists[upper->pntToNode->nl_list_index] = father->info->nl_nodes[m-1];

					nl_lists[upper->pntToNode->nl_list_index].clear();
					int upper_index = (m < father->info->minnodes_indexes.size()) ? father->info->minnodes_indexes[m] : father->info->minnodes.size();
					nl_lists[upper->pntToNode->nl_list_index].reserve(upper_index - father->info->minnodes_indexes[m-1]);
					for(int aj = father->info->minnodes_indexes[m-1]; aj < upper_index; aj++){
						nl_lists[upper->pntToNode->nl_list_index].push_back(father->info->minnodes[aj]);
					}

					/*if(father->info->nl_a_S_type==8){
						for(int aj=0; aj<jp+1; aj++)
							todelete->at(upper->pntToNode->a_S_index)[aj]=father->info->nl_a_S.nl_a_S_8[offset+aj];
	  				}
	  				else
	  					if(father->info->nl_a_S_type==16){
							for(int aj=0; aj<jp+1; aj++)
								todelete->at(upper->pntToNode->a_S_index)[aj]=father->info->nl_a_S.nl_a_S_16[offset+aj];
	  					}
	  					else{
							for(int aj=0; aj<jp+1; aj++)
								todelete->at(upper->pntToNode->a_S_index)[aj]=father->info->nl_a_S.nl_a_S_32[offset+aj];
	  					}*/


					//for(int aj=0; aj<jp+1; aj++)
						//todelete->at(upper->pntToNode->a_S_index)[aj]=father->info->nl_a_S[offset+aj];
				}
				else{
					//cout << "access (23) to upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					/*for(int aj=0; aj<jp+1; aj++)
						todelete->at(upper->pntToNode->a_S_index)[aj]=addressOf[nodeList[offIDList]]->a_S[aj]; */
					nl_lists[upper->pntToNode->nl_list_index].push_back(nodeList[offIDList]);
				}
				/*cout << "done 8.2 " << nextfree_a_S << endl;
				for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
				//time_keeper += get_cpu_time() - temp_time;

				/*if(DEEP_DEBUG)
				cout << "nodelist: set a_S to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;*/
				/* SignificantMiner */
				upper->pntToNode->node=upper->ID;
				upper->pntToNode->nextNodePointer=-1;
				(&upper->pntToNode->intersection)[0]=0;
				///insert the node in the pointerList of its greater item
				*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
				///update lastNodePnt
				HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
			}
			else{
				upper->pntToNode->count+=nodeList[offIDList+1];
				/* SignificantMiner */
				upper->pntToNode->count_a_S+=nodeList[offIDList+2];

				//temp_time = get_cpu_time();

        int temp_as_value = (addressOf[nodeList[offIDList]]->count >= 128) ? addressOf[nodeList[offIDList]]->a_S_[0] : addressOf[nodeList[offIDList]]->a_S[0];

				if(temp_as_value!=nodeList[offIDList+2] || addressOf[nodeList[offIDList]]->count!=nodeList[offIDList+1]) {
					//cout << "Error with nodelist!(13)"<<endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->a_S[0]=" << addressOf[nodeList[offIDList]]->a_S[0] << endl;
					//cout << "Error with nodelist: nodeList[offIDList+2]=" << nodeList[offIDList+2] << endl;
					//cout << "Error with nodelist: father->info->nl_a_S[(m-1)*(jp+1)]=" << father->info->nl_a_S[(m-1)*(jp+1)] << endl;
					//cout << "Error with nodelist: addressOf[nodeList[offIDList]]->count=" << addressOf[nodeList[offIDList]]->count << endl;
					//cout << "Error with nodelist: nodeList[offIDList+1]=" << nodeList[offIDList+1] << endl;
					//int offset = (m-1)*(jp+1);
					//cout << "access (24) to upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;

					/*if(father->info->nl_a_S_type==8){
						for(int aj=0; aj<jp+1; aj++)
							todelete->at(upper->pntToNode->a_S_index)[aj]+=father->info->nl_a_S.nl_a_S_8[offset+aj];
	  				}
	  				else
	  					if(father->info->nl_a_S_type==16){
							for(int aj=0; aj<jp+1; aj++)
								todelete->at(upper->pntToNode->a_S_index)[aj]+=father->info->nl_a_S.nl_a_S_16[offset+aj];
	  					}
	  					else{
							for(int aj=0; aj<jp+1; aj++)
								todelete->at(upper->pntToNode->a_S_index)[aj]+=father->info->nl_a_S.nl_a_S_32[offset+aj];
	  					}*/

					int upper_index = (m < father->info->minnodes_indexes.size()) ? father->info->minnodes_indexes[m] : father->info->minnodes.size();
					for(int aj = father->info->minnodes_indexes[m-1]; aj < upper_index; aj++){
						nl_lists[upper->pntToNode->nl_list_index].push_back(father->info->minnodes[aj]);
					}

					//for(int aj=0; aj < father->info->nl_nodes[m-1].size(); aj++)
	  					//nl_lists[upper->pntToNode->nl_list_index].push_back( father->info->nl_nodes[m-1][aj] );


					//for(int aj=0; aj<jp+1; aj++)
					//todelete->at(upper->pntToNode->a_S_index)[aj]+=father->info->nl_a_S[offset+aj]; //cout << "increased to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;
				}
				else{
					//cout << "access (25) to upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					/*for(int aj=0; aj<jp+1; aj++)
						todelete->at(upper->pntToNode->a_S_index)[aj]+=addressOf[nodeList[offIDList]]->a_S[aj]; */
					nl_lists[upper->pntToNode->nl_list_index].push_back(nodeList[offIDList]);
				}
				//time_keeper += get_cpu_time() - temp_time;

				/*if(DEEP_DEBUG)
				cout << "nodelist: increased a_S of " << nodeList[offIDList+2] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;*/
				/* SignificantMiner */
				(&upper->pntToNode->intersection)[0]=0;
			}
		}

		///update the offIDList
		/* SignificantMiner */
		tipoInt nf = 4; // number of nodelist fields: 3 + 1 (added a_S)
		/* SignificantMiner */
		if (nodeList[offIDList+3]==-1){
			offIDList=offIDList+nf; ///we now that the length of the intersection is 0 (as number of tipoInt that memorize the items in the intersection associated to the pointer)
		}
		else{
			offIDList=offIDList+nodeList[offIDList+3]+nf;
		}
	}


	///now we must consider all the list of items greater than father->info->minItem and visit the element in their pointerList; we must begin from the last item (the greatest)

	NodePointer* currentpnt;	///temporary NodePointer used for the elements in the pointedList
	tipoInt off;			///offset of currentpnt in manualmem

	for (tipoInt l=father->info->coreIndex-1;l>=0;l--){
		off=HT[l].pointerList;
		while(off!=-1){
			currentpnt=(NodePointer*)&(manualmem[off]);
			patNode=addressOf[currentpnt->node];	///PatriciaNode pointed by currentpnt

			tipoInt numItems=(tipoInt) 1+(patNode->lastItem-&(patNode->item));	///number of items in patNode

			///now consider all the items in the node and update the corresponding HT rows
			tipoInt index;	///current item ID
			for (tipoInt s=numItems-1;s>=0; s--){
				index=(&patNode->item)[s];
				///we know that all the items in the node are all lower than father->info->CI
				HT[index].count+=currentpnt->count;	///update the count with the count of the pointer
				/* SignificantMiner */
				//temp_time = get_cpu_time();
				//cout << "access (25) to currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				HT[index].a_S+=currentpnt->count_a_S;
				/*for(int aj=0; aj<jp+1; aj++)
				HT[index].a_S[aj]+=todelete->at(currentpnt->a_S_index)[aj];*/
				//time_keeper += get_cpu_time() - temp_time;
				/*if(DEEP_DEBUG)
				cout << "increased a_S (12) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << HT[index].a_S << endl;*/
				/* SignificantMiner */
			}

			///now we must create the NodePointer to the father of the patNode and insert it in the pointerList of its greater item
			PatriciaNode* upper=patNode->father;	///upper is the father of patNode
			if (upper!=NULL){	///the node isn't a root's children
				if (upper->pntToNode==NULL){
				///create the pointer to the node
				if (nextfree+sizePnt>lastvalid){
          cout << "call manual mem 9" << endl;
					resizeManualMem();
				}
				upper->pntToNode=(NodePointer*)nextfree;
				//todelete->push_back(upper->pntToNode);
				living_nodepointers++;
				nextfree=nextfree+sizePnt;	///update nextfree
				upper->pntToNode->count=currentpnt->count;
				/* SignificantMiner */
				upper->pntToNode->count_a_S=currentpnt->count_a_S;
				//temp_time = get_cpu_time();
				//upper->pntToNode->a_S_index = getNewTodeleteIndex();
				upper->pntToNode->nl_list_index = getNewListIndex();

				//cout << "(12) allocated jp+1 things at " << todelete->at(upper->pntToNode->a_S_index) << endl;

				//cout << "access (26) to upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
				//cout << "access (26) to currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
				/*for(int aj=0; aj<jp+1; aj++)
				todelete->at(upper->pntToNode->a_S_index)[aj]=todelete->at(currentpnt->a_S_index)[aj];*/

				nl_lists[upper->pntToNode->nl_list_index]=nl_lists[currentpnt->nl_list_index];
			/*cout << "done 10.2 " << nextfree_a_S << endl;
			for(int aj=0; aj<jp+1; aj++)
						cout << todelete->at(upper->pntToNode->a_S_index)[aj] << " ";
					cout << endl;*/
				//time_keeper += get_cpu_time() - temp_time;
				/* SignificantMiner */
				upper->pntToNode->node=upper->ID;
				upper->pntToNode->nextNodePointer=-1;
				///insert the node in the pointerList of its greater item
				*(HT[*upper->lastItem].lastNodePnt)=((tipoInt*)upper->pntToNode)-manualmem;
				///update lastNodePnt
				HT[*upper->lastItem].lastNodePnt=&upper->pntToNode->nextNodePointer;
				///create the pointer to the node
				}
				else{
					///we must update the field of the pointer to upper
					upper->pntToNode->count+=currentpnt->count;
					/* SignificantMiner */
					upper->pntToNode->count_a_S+=currentpnt->count_a_S;
					//temp_time = get_cpu_time();
					//cout << "access (27) to upper->pntToNode->a_S_index = " << upper->pntToNode->a_S_index << endl;
					//cout << "access (27) to currentpnt->a_S_index = " << currentpnt->a_S_index << endl;
					/*for(int aj=0; aj<jp+1; aj++)
					todelete->at(upper->pntToNode->a_S_index)[aj]+=todelete->at(currentpnt->a_S_index)[aj];*/

					for(int aj=0; aj < nl_lists[currentpnt->nl_list_index].size(); aj++)
						nl_lists[upper->pntToNode->nl_list_index].push_back(nl_lists[currentpnt->nl_list_index][aj] );

					//time_keeper += get_cpu_time() - temp_time;
					/*if(DEEP_DEBUG)
					cout << "increased a_S (13) for patricia node of " << todelete->at(currentpnt->a_S_index)[0] << " to " << todelete->at(upper->pntToNode->a_S_index)[0] << endl;*/
					/* SignificantMiner */
					(&upper->pntToNode->intersection)[0]=0;
				}
			}
			off=currentpnt->nextNodePointer;
			patNode->pntToNode=NULL;
		}
	}

	///now set to null the pntToNode of nodes that contains the CI
	for (tipoInt h=1; h<=HT[father->info->coreIndex].visitedList[0].IDnode; h++){
		addressOf[HT[father->info->coreIndex].visitedList[h].IDnode]->pntToNode=NULL;
	}

	///now find items in closure(father)
	Transaction closure;	///here we insert all the items that are in Clo(father)
	if (nextfree+father->info->coreIndex+2>lastvalid){
    cout << "call manual mem 9" << endl;
		resizeManualMem();
	}
	closure=nextfree;
	nextfree=nextfree+father->info->coreIndex+2;
	closure[0]=1;
	closure[1]=father->info->coreIndex;
	for (int i=father->info->coreIndex-1; i>=0; i--){
		///if the item has the same count of father, it is in closure=Clo(father)
		if (HT[i].count==father->support){
			///insert i in closure
			closure[0]++;
			closure[closure[0]]=i;
		}
	}
	return closure;
}

///given a NodeIDList with typeList==0, return the corrisponding NodeIDList with typeList==1, that is set to -1 all the intersection length in the new NodeIDList

NodeIDList changeNIDL(NodeIDList nil){

	/* SignificantMiner */
	tipoInt nf = 4; // number of nodelist fields: 3 + 1 (added a_S)
	/* SignificantMiner */

	NodeIDList toret;	///NodeIDList to return
	///find the size for the NodeIDList to return
	tipoInt newsize=nf*nil[0]+1;
	toret=(tipoInt*)malloc(newsize*sizeof(tipoInt));
	toret[0]=nil[0];
	for(tipoInt j=0;j<nil[0]; j++){
		toret[j*nf+1]=nil[j+1];	///ID of the node
		toret[j*nf+2]=addressOf[nil[j+1]]->count;	///count of the node
		/* SignificantMiner */
		//cout << "a_S lost (4) for support " << toret[j*3+2] << endl;
		//cout << "having value " << addressOf[nil[j+1]]->a_S << endl;
		/* SignificantMiner */
    int temp_as_value = (addressOf[nil[j+1]]->count >= 128) ? addressOf[nil[j+1]]->a_S_[0] : addressOf[nil[j+1]]->a_S[0];
		toret[j*nf+3]=temp_as_value;	///a_S of the node

		if(DEEP_DEBUG){
			cout << "changeNIDL: ID of linked node = " << toret[j*nf+1] << endl;
			cout << "   changeNIDL: support inside list = " << toret[j*nf+2] << endl;
			cout << "   changeNIDL: a_S of linked node = " << toret[j*nf+3] << endl;
			/*for(int aj=0;aj<jp+1;aj++)
				cout<< addressOf[toret[j*nf+1]]->a_S[aj] <<" | ";
			cout<<endl;*/
		}

		toret[j*nf+4]=-1;	///intersection
	}
	return toret;
}

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  cout << "Error: signal handler " << sig <<  endl;
  #ifdef trace_state
  cout << "program state: " << (int)program_state << endl;
  #endif
  cout << "queue program state: " << (int)queue_program_state << endl;
  cout << "debug_info_1: " << debug_info_1 << endl;
  cout << "debug_info_2: " << debug_info_2 << endl;

  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);

  kill(getpid(), SIGKILL);
  exit(1);
}


/** main procedure
*/

int main(int argc, char **argv){

	signal(SIGTERM, handler);
	signal(SIGSEGV, handler);
	signal(SIGABRT, handler);

	cout << "\n*******************************************************************************" << endl;
	cout <<   "*                                                                             *" << endl;
	cout <<   "*   SignificantMiner: Efficient Incremental Mining of Significant Itemsets    *" << endl;
	cout <<   "*           Copyright (C) 2017 Leonardo Pellegrina, Fabio Vandin              *" << endl;
	cout <<   "*                                                                             *" << endl;
	cout <<   "*              This program comes with ABSOLUTELY NO WARRANTY;                *" << endl;
	cout <<   "*              see the COPYING file in the package for details.               *" << endl;
	cout <<   "*       This is free software, and you are welcome to redistribute it         *" << endl;
  cout <<   "*             under certain conditions; see the COPYING file in               *" << endl;
	cout <<   "*                         the package for details.                            *" << endl;
	cout <<   "*                                                                             *" << endl;
	cout <<   "*******************************************************************************" << endl;

	double elapsed = get_cpu_time();
	double total_time = elapsed;
  double time_withoutpatricia = 0;
	nlinqueue_8 = 0;
	nlinqueue_16 = 0;
	nlinqueue_32 = 0;
	inferior_psi_bound_happened = 0;
	inferior_psi_bound_happened_neglect = 0;
	list_length_inqueue = 0;
	nl_nodes_length_inqueue = 0;
  max_nl_size_ = 0;

	memory_bounded = false;

  K_significant_patterns = 100;

	in_use_ram = 0.0;

  as_processing_time = 0.0;

	pruned_patterns = 0;
	tested = 0;

	cache_accesses = 0;
	new_precomputation = 0;

	int timesnum = 3;
	double temptime = 0.0;
	//timekeepers.resize(timesnum);

	///to memorize statistics on output
	tipoInt* inOutput;

	///read input parameters
	/*PatriciaTrie* */p=new PatriciaTrie();

	///read input parameters
	if ( readSpecFile(argc,argv) == 0){
		return 0;
	}

	store_results = K_significant_patterns > 0;

	//time_keeper = 0.0;
	nodelist_inqueue = 0;
	living_itemsets=0;
	info_deleted=0;
	queue_program_state=0;

	/* SignificantMiner */
	// read all class labels and compute the permutations
	c_labels.open(c_fileinput.c_str(),ifstream::in); ///<class labels file;

  cout << "Reading labels... "<< endl;
  double reading_dataset_time = get_cpu_time();
	initClassLabels();
  reading_dataset_time = get_cpu_time() - reading_dataset_time;
  cout << "done in " << reading_dataset_time << " seconds" << endl;

  cout << "Memory peak after permutation matrix allocation: " << (int)measurePeakMemory() << " MByte" << endl;

  //cout << "Computing permutations... "<< endl;
  double permutations_time = get_cpu_time();
	//computePermutations();
  permutations_time = get_cpu_time() - permutations_time;
  //cout << "done in " << permutations_time << " seconds" << endl;


	printPermutations(DEEP_DEBUG, SOFT_DEBUG);
	cout.precision(10);
	/* SignificantMiner */

	/**first scan of the dataset: collect the support of singletons in the count_element array;
	*/

  cout << "First scan of dataset... "<< endl;
  double firstscan_dataset_time = get_cpu_time();
	firstScan_new(p);
  firstscan_dataset_time = get_cpu_time() - firstscan_dataset_time;
  cout << "done in " << firstscan_dataset_time << " seconds" << endl;

  cout << "Memory peak after count_a_S allocation: " << (int)measurePeakMemory() << " MByte" << endl;

	#ifdef trace_state
	program_state = 8;
	#endif

  cout << "IL operations... "<< endl;
  double ILoperations_time = get_cpu_time();

	///now count the number of effective items in the dataset

	num_item=0;
	cloempty=string();
	string tmpString;
	tipoInt cloemptyLenght=0;	///lenght of the closure of emptyset

	///remove the items in Clo(empty) and memorize them
	bool emptyIsNotClosed=false;	///it indicates if the emptyset is closed or not
	for (int j=0; j<=maxID; j++){
		if ((p->count_element[j])>0) num_item++;
		if (p->count_element[j]==effect_num_tr){
			///update emptyIsNotClosed
			emptyIsNotClosed=true;
			///update length of the closure of emptyset
			cloemptyLenght++;
//			cloempty=cloempty+" ";
			tmpString=itos(j);
			cloempty=cloempty+tmpString+" ";
			///now we can consider the item as not frequent
			cout << j << " is set infrequent"<<endl;
			p->count_element[j]=0;
			/*for(int aj=0; aj<jp+1; aj++)
			p->count_a_S_element[j*(jp+1)+aj]=0;*/
		}
	}

	///create, order IL and use it to raise minimum supports
	p->makeIL();

	srand(time(0));
	p->IL_sort(0,maxID);

	#ifdef trace_state
	program_state = 7;
	#endif

	///for output
  #ifdef write_results_to_file
	ofstream outfile;
  std::stringstream ssoutput_path_result;
  ssoutput_path_result << fileinput << "_" << K_significant_patterns << "_" << alpha<< "_" << jp << "_result.txt";
	string output_path_result = ssoutput_path_result.str();

	outfile.open (output_path_result.c_str());
  #endif
  std::stringstream testingoutput_path_result;
  testingoutput_path_result << fileinput << "_" << K_significant_patterns << "_" << jp << "_testing2.txt";
	string output_path_testing2 = testingoutput_path_result.str();
	out_file_stat_2.open(output_path_testing2.c_str());
	out_file_stat_2 << "elapsed;nlinqueue_8;nlinqueue_16;nlinqueue_32;s_supp;inqueue;suppMin;produced;underbound;ram;" << endl;

	///print Clo(empty) on output
	if ( emptyIsNotClosed ){
		//outfile << cloempty << " : " << num_tr << endl;
		/* SignificantMiner */
		//outfile << cloempty << " : " << num_tr << " : " << n1 << endl;
		/* SignificantMiner */
	}

	///create the array for the util node heuristic

	maximumSupport=p->IL[0].count;	///a util node will have a count lower than the maximum count of an item
	utilNodeArray=(tipoInt*)calloc((maximumSupport+1),sizeof(tipoInt));
	itsOfSupp=(tipoInt*)calloc((maximumSupport+1),sizeof(tipoInt));

	inOutput=new tipoInt [max_trans_length+1];
	for(int i=0;i<max_trans_length+1;i++){
		inOutput[i]=0;
	}
	p->ILRaise();
	if(DEEP_DEBUG)
	p->printIL();
	dataset.close();

  ILoperations_time = get_cpu_time() - ILoperations_time;
  cout << "done in " << ILoperations_time << " seconds" << endl;

	///array for the number of nodes in which an item is
	numNPI=(tipoInt*)calloc(num_item,sizeof(tipoInt));

	/**second scan: "transform" transactions and insert them in the PatriciaTrie
	*/

  cout << "Memory peak after IL operations allocation: " << (int)measurePeakMemory() << " MByte" << endl;

  cout << "Second scan of dataset and Patricia Trie creation... "<< endl;
  double secondscan_dataset_time = get_cpu_time();

	p->secondScan_new();

	#ifdef trace_state
	program_state = 6;
	#endif

  memory_baseline = (int)measurePeakMemory();

  cout << "Memory peak after Patricia Trie allocation: " << (int)measurePeakMemory() << " MByte" << endl;

  time_withoutpatricia = get_cpu_time();

	///allocate the header table that is used in every ppc-e generation

	HT=(HeaderTable*) calloc((num_item),sizeof(HeaderTable));

	///create the addressOf array
	addressOf=new PatriciaNode* [numNodes];

	///create the array of lists for the NodeIDList of Clo(emptyset)
	CloEmptyNodeIDList=new NodeIDElem*[num_item+1];
	lastCloEmpty= new NodeIDElem**[num_item+1];
	for (int i=0; i<=num_item; i++){
		lastCloEmpty[i]=&(CloEmptyNodeIDList[i]);
	}

	///memory (re)used for computation in every generation of ppce
	//tipoInt bigmem=totItem*5;
  tipoInt bigmem=totItem*100;
	manualmem=(tipoInt*)malloc(bigmem*sizeof(tipoInt));
	nextfree=&(manualmem[0]);
	lastvalid=&(manualmem[bigmem-1]);
	int prov=(tipoInt)(lastvalid-nextfree);

	///initialization of the space that will be used for the visited lists

	for (int h=0; (h<num_item)&&(p->IL[h].count>=suppMin); h++){
		///create the visited node list that will be used in every computation of fillHT: as upper bound for the number of elements of the list, wee use the number of nodes that contains the item "i"
		HT[h].visitedList=(VisitedNodeElem*)malloc((numNPI[h]+1)*sizeof(VisitedNodeElem));

		///in visitedList[0].IDnode we put the length of the current list
		HT[h].visitedList[0].IDnode=0;
	}

	///now assign ID and intersect to the nodes of the PatriciaTrie
	p->assignIDAndIntersect();
	if(DEEP_DEBUG)
	p->printAfter();

	///now the lastCloEmpty array is not util
	delete [] lastCloEmpty ;


  secondscan_dataset_time = get_cpu_time() - secondscan_dataset_time;
  cout << "done in " << secondscan_dataset_time << " seconds" << endl;

  //p->printCountDistribution();

	#ifdef trace_state
	program_state = 71;
	#endif

	// free not used anymore memory
	//free(permutations);
	//free(p->count_a_S_element);

	#ifdef trace_state
	program_state = 72;
	#endif


	///now reinitialize the utilNodeArray for the trivial heuristic
	memset(utilNodeArray,0,sizeof(tipoInt)*(maximumSupport+1));

	///if Clo(emptyset)!=emptyset, there is a closed itemset that is already in q (ideally): we subtract 1 by k and k_max so that there is no need to resize utilNodeArray
	/*if (emptyIsNotClosed){
		k--;
		k_max--;
		if (k==0){
			suppMinCurr=effect_num_tr;
		}
		if (k_max==0){
			suppMin=effect_num_tr;
		}
		///update statistics
		inOutput[cloemptyLenght]++;
	}*/

	///we have subtract 1 by k and k_max, so produced is always =0

	produced=0;
	inqueue=0;
	inserted=0;
	significant_itemsets=0;

	//suppMinCurr = 1;
	//suppMin = 1;

	s_supp = effect_num_tr;

	cout << "s_supp = " << s_supp << endl;

	///create the priority queue for the max-extraction
	//Queue* q=new Queue();

	///create the priority queue for the min-extraction
	//QueueMin* qMin=new QueueMin();

	// queue for SignificantMiner
	//typeQueue q=new QueueMinMax();
	typeQueue q=new QueueMinMax_test(effect_num_tr);
	//QueueMinMaxResult* q_res=new QueueMinMaxResult();
	typeQueueRes q_res=new QueueMinMax_Restest(n1);


	temp_size = inital_temp_size;
	temp_size_8 = inital_temp_size;
	temp_size_16 = inital_temp_size;
	temp_size_32 = inital_temp_size;
	temp_nl_a_S = (tipoInt*)malloc(temp_size*(jp+1)*sizeof(tipoInt));
	temp_nl_a_S_8 = (uint8_t*)malloc(temp_size_8*(jp+1)*sizeof(uint8_t));
	temp_nl_a_S_16 = (uint16_t*)malloc(temp_size_16*(jp+1)*sizeof(uint16_t));
	temp_nl_a_S_32 = (uint32_t*)malloc(temp_size_32*(jp+1)*sizeof(uint32_t));
	nlasmalloc=0;

	file_path=new char[20];




	/*cout << "Allocated first part of memory" << endl;
	cout << "Number of nodes " << numNodes << endl;
	cout << "sizeof(std::vector<tipoInt>) " << sizeof(std::vector<tipoInt>) << endl;
	cout << "sizeof(tipoInt) " << sizeof(tipoInt) << endl;
	cout << "sizeof(NodePointer) " << sizeof(NodePointer) << endl;*/
	//sleep(5);

	#ifdef trace_state
	program_state = 5;
	#endif

	//pat_tree_ram = ((((long long)numNodes*((long long)jp+1)*sizeof(tipoInt) + sizeof(PatriciaNode) * ((long long)numNodes) )/1000000.0));

	double closed_ext_time = get_cpu_time();

	findExtCloEmpty(p->IL, q, q_res);//, qMin);

	closed_ext_time = get_cpu_time() - closed_ext_time;

	cout << "Time to generate closed extensions: " << closed_ext_time << endl;

  cout << "Memory peak after findExtCloEmpty: " << (int)measurePeakMemory() << " MByte" << endl;

	#ifdef trace_state
	program_state = 4;
	#endif

	perm_supp.resize(jp+1);
	list_stat.resize(1);
	list_op_stat.resize(1);

	nextfree=manualmem;

	///now clear the HT for next extractions
	tipoInt firstNotF=0;
	while(firstNotF<num_item&&p->IL[firstNotF].count>=suppMin){
		HT[firstNotF].count=0;			///count of the extension
		/* SignificantMiner */
		//free(HT[firstNotF].a_S);
		//HT[firstNotF].a_S=(tipoInt*)malloc((jp+1)*sizeof(tipoInt));
		//for(int aj=0; aj<jp+1; aj++)
		HT[firstNotF].a_S=0;
		/* SignificantMiner */
		HT[firstNotF].intersection=NULL;	///intersection on the prefix
		HT[firstNotF].visitedList[0].IDnode=0;		///number of vidited elements
		HT[firstNotF].pointerList=-1;		///empty list
		HT[firstNotF].lastNodePnt=&(HT[firstNotF].pointerList);

		///NEW: DYNAMIC VISITEDLIST
		if (numNPI[firstNotF]>10){
			HT[firstNotF].visitedList=(VisitedNodeElem*) realloc(HT[firstNotF].visitedList, (1+(numNPI[firstNotF]/5))*sizeof(VisitedNodeElem));
			HT[firstNotF].visitedList[0].numInInter=numNPI[firstNotF]/5;		///number of vidited elements that can be stored in visitedList
		}
		else{
			HT[firstNotF].visitedList[0].numInInter=numNPI[firstNotF];		///number of vidited elements that can be stored in visitedList
		}
		HT[firstNotF].visitedList[0].IDnode=0;		///number of visited elements
		firstNotF++;
	}

	int fals;

	Itemset* max;
	Transaction toOutput;
	StartNodeList sttmp;
	NodeIDList tmpnidl;
	bool finish;

	//if (k_max==k) finish=true;
	//else finish=false;
	finish=false;
	bool firstOfSuppMin=false;	///it indicates if we had found the first itemset of support=suppMinCurr
	/*if (q->last!=-1){
		if (suppMinCurr!=q->heap->support)	firstOfSuppMin=false;
		else firstOfSuppMin=true;
	}*/
	///extract max and generate ppc-e

	bool created;	///it is true if tmpnidl is different from max->info->IDlist




	elapsed = get_cpu_time();

	//for(int aj=0; aj<timesnum; aj++)
	//	timekeepers[aj] = 0.0;

	allocateToDel();

	max_ram_l1 = 0.8 * (double)max_ram;
	max_ram_l2 = 0.9 * (double)max_ram;
	cout << "Max allowed ram levels: " << max_ram << " L1: " << max_ram_l1 << " L2: " << max_ram_l2 << endl;

	long last_itemset_report = 0;
	long last_itemset_memorycheck = 0;
	int last_time_report = 0;

  time_pat = get_cpu_time() - total_time;
  cout << "Total initialization time: " << time_pat << endl;
  time_pat = secondscan_dataset_time;

  cout << "Memory peak after entire initialization: " << (int)measurePeakMemory() << " MByte" << endl;

  // PATTERN EXPLORATION BEGINS

	//while((produced<k) && (q->last!=-1) && !(finish&&firstOfSuppMin)){
	while(q->getMax()>=suppMin){

    /*cout << "q->getMax() = " << q->getMax() << endl;
    cout << "q->getInqueue() = " << q->getInqueue() << endl;
    cout << "suppMin = " << suppMin << endl;*/

		if(memory_bounded == false && (produced - last_itemset_memorycheck) > 9999){

			#ifdef debug_memory_usage
			cout << " q->getQueueMemory() = " << q->getQueueMemory() << endl;
			#endif

			last_itemset_memorycheck = produced;

			if(measurePeakMemory()-memory_baseline > max_ram_l1){
				memory_bounded = true;
				memory_bounded_max_inqueue = q->getInqueue() * 0.9;
				memory_bounded_max_nl_length = list_length_inqueue * 0.9;
				memory_bounded_max_nl_min_length = nl_nodes_length_inqueue * 0.9;
				max_ram =  ((int)measurePeakMemory()-memory_baseline)+500;
				cout << "MEMORY BOUNDED MODE ACTIVATE! Memory limit = " << max_ram << " Max elements in queue = " << memory_bounded_max_inqueue << endl;
				cout << "memory_bounded_max_nl_length = " << memory_bounded_max_nl_length << " memory_bounded_max_nl_min_length = " << memory_bounded_max_nl_min_length << endl;
			}

		}

		if(memory_bounded == true && (produced - last_itemset_memorycheck) > 9999){

			#ifdef debug_memory_usage
			cout << " q->getQueueMemory() = " << q->getQueueMemory() << endl;
			#endif

			last_itemset_memorycheck = produced;

			if(measurePeakMemory()-memory_baseline > max_ram){
			//if(in_use_ram > max_ram){
				max_ram = ((int)measurePeakMemory()-memory_baseline)+500;
				memory_bounded_max_inqueue = memory_bounded_max_inqueue * 0.9;
				memory_bounded_max_nl_length = memory_bounded_max_nl_length * 0.9;
				memory_bounded_max_nl_min_length = memory_bounded_max_nl_min_length * 0.9;
				if(q->low_factor > 1.0)
					q->low_factor -= 0.25;
				cout << "MEMORY BOUNDED MODE UPDATE! Memory limit = " << max_ram << " Max elements in queue are now = " << memory_bounded_max_inqueue << endl;
				cout << "memory_bounded_max_nl_length = " << memory_bounded_max_nl_length << " memory_bounded_max_nl_min_length = " << memory_bounded_max_nl_min_length << " low_factor = " << q->low_factor << endl;
			}

		}

		if( (produced - last_itemset_report) > 99999 ){

		last_itemset_report = produced;

		tipoInt elapsed_since_last_report = (tipoInt)(get_cpu_time() - last_time_report);


		if( (elapsed_since_last_report > 99 ) || (elapsed_since_last_report > 999) ) {

			// compute list statistics
			double tot = 0.0;
			double avg = 0.0;

			#ifdef print_list_stats
			for(int aj=0; aj < list_stat.size(); aj++){
				tot += list_stat[aj];
			}
			for(int aj=0; aj < list_stat.size(); aj++){
				avg += list_stat[aj] * aj / tot;
			}
			tot = 0.0;
			double avg_op = 0.0;
			for(int aj=0; aj < list_op_stat.size(); aj++){
				tot += list_op_stat[aj];
			}
			for(int aj=0; aj < list_op_stat.size(); aj++){
				avg_op += list_op_stat[aj] * aj / tot;
			}
			#endif

			cout << "------ NEW REPORT for " << fileinput<<" ------- "  << endl;


			cout << "--- EXPLORATION INFO ---"  << endl;
			cout << "   elapsed: " << (tipoInt)(get_cpu_time() - elapsed)  << endl;
			cout << "   s_supp " << s_supp  << endl;
			cout << "   suppMin " << suppMin  << endl;
			cout << "   explored " << produced  << endl;
			cout << "   tested " << tested  << endl;
			cout << "   underbound " << inferior_psi_bound_happened_neglect  << endl;
			cout << "   pruned " << pruned_patterns << endl;
      cout << "   as processing time " << as_processing_time << endl;
			cout << "   time for p-values precomputation " << precomputetime << endl;
			cout << "   cache accesses " << cache_accesses << endl;
			cout << "   new precomputation count " << new_precomputation << endl;
      cout << "   current p-values cache memory " << current_cache_size_mb << endl;
			cout << "   peak memory: " << measurePeakMemory() << endl;
      cout << "   memory bounded mode ";
      if(memory_bounded){
        cout << "active!"<< endl;
        cout << "   low factor: "<< q->low_factor << endl;
      }
      else{
        cout << "non active" << endl;
      }
			cout << "   list_length_inqueue " << list_length_inqueue  << endl;
			cout << "   nl_nodes_length_inqueue " << nl_nodes_length_inqueue << endl;
			#ifdef print_list_stats
			cout << "   nllist_avglength " << avg  << endl;
			cout << "   nllist_avgop " << avg_op  << endl;
			#endif
		    #ifdef debug_memory_usage
		    cout << "   q->getQueueMemory() = " << q->getQueueMemory() << endl;
		    #endif

			cout << "--- ITEMSET QUEUE INFO ---"  << endl;
			cout << "   inqueue " << q->getInqueue()   << endl;
			cout << "   q->getMin() " << q->getMin()   << endl;

			cout << "--- RESULTS INFO ---"  << endl;
			cout << "   significant patterns extracted: " << significant_itemsets << endl;
			cout << "   patterns in results queue: " << q_res->getInqueue() << endl;
		    cout << "   q_res->getMax() " << q_res->getMax() << endl;
		    cout << "   q_res->getLogMin() " << q_res->getLogMin() << endl;
		    cout << "   getLogPsi(s_supp) " << getLogPsi(s_supp) << endl;
		    cout << "------------------------------- " << endl;

			//cout << " ram " << (tipoInt)in_use_ram
			//out_file_stat_2 << (tipoInt)(get_cpu_time() - elapsed) << ";" << nlinqueue_8 << ";" << nlinqueue_16 << ";" << nlinqueue_32 << ";" << s_supp << ";" << q->getInqueue() << ";" << suppMin << ";" << produced << ";" << inferior_psi_bound_happened_neglect << ";" << (tipoInt)in_use_ram << endl;
			last_itemset_report = produced;
			last_time_report = (tipoInt)(get_cpu_time());

		}

		}

		#ifdef trace_state
		program_state = 0;
		#endif


		if(memory_bounded == true &&
			(q->getInqueue() > memory_bounded_max_inqueue ||
				in_use_ram > max_ram ||
				memory_bounded_max_nl_length < list_length_inqueue  ||
				memory_bounded_max_nl_min_length < nl_nodes_length_inqueue )) {

			#ifdef trace_state
			program_state = 1;
			#endif


			max = q->removeLow();
      if(max == NULL){
          cout << "Extracted max IS NULL! " << endl;
          abort();
      }
		}
		else{

		#ifdef trace_state
		program_state = 11;
		#endif

		max=q->removeMax();	///remove the first itemset from the priority queue

    if(max == NULL){
        cout << "Extracted max IS NULL! " << endl;
        abort();
    }

		// update higher support bound
		s_supp = max->support;

		}

    if(max == NULL){
    		cout << "Extracted max IS NULL! " << endl;
        abort();
    }

		/*cout << "here ok IL? 1 "<<endl;
		cout << p->IL[1].item << endl;*/

    #ifdef trace_state
    program_state = 111;
    #endif

		bool new_allocation=true;


		#ifdef trace_state
		program_state = 1111;
		debug_info_1 = -1;
		debug_info_2 = 0;
		//debug_info_1 = max->info->nl_nodes.size();
		//if(debug_info_1 != 0)
		  //debug_info_2 = max->info->nl_nodes[0].size();
		#endif


		// update lists stats
    /*if(max->info->nl_nodes.size() > 0 && max->info->nl_nodes.size() <= max_nl_size_){
		list_length_inqueue-=max->info->sizelista;
		for(int aj=0; aj<max->info->nl_nodes.size(); aj++){
			nl_nodes_length_inqueue-=max->info->nl_nodes[aj].size();
			debug_info_2 += max->info->nl_nodes[aj].size();
		}
    }
    else{
      cout << "ERROR IN UPDATE LISTS STATS WITH max->info->nl_nodes.size(): " << max->info->nl_nodes.size() << endl;
    }*/

    if(max->info->minnodes.size() > 0 && max->info->minnodes_indexes.size() <= max_nl_size_){
		list_length_inqueue-=max->info->sizelista;
		nl_nodes_length_inqueue-=max->info->minnodes.size();
    }
    else{
      cout << "ERROR IN UPDATE LISTS STATS WITH max->info->minnodes_indexes.size(): " << max->info->minnodes_indexes.size() << endl;
    }

    #ifdef trace_state
    program_state = 112;
    #endif


		if(SOFT_DEBUG)
		cout << "Extracted max with typelist = " << max->info->typeList << endl;

		if (max->info->typeList==0){
			tmpnidl=changeNIDL(max->info->IDlist);
			created=true;
		}
		else{
			if (max->info->typeList==2){
				#ifdef trace_state
				program_state = -1;
				#endif

				//cout << "reading from file 1" << endl;
				cout<<"RD " << max->info->fileName;
				cout << " max->support " << max->support;
				cout << " max->info->nl_a_S_type " << (int)max->info->nl_a_S_type;
				cout << " max->info->nl_a_S_size " << max->info->nl_a_S_size << endl;

				FILE * pFile;
				for(int aj=0; aj<20; aj++) file_path[aj]=' ';
				sprintf(file_path, "%d.l", max->info->fileName);
				pFile = fopen ( file_path , "rb" );
				//cout << "file opened " << endl;
				tipoInt length;
				fread(&length,sizeof(tipoInt),1,pFile);
				tmpnidl=(tipoInt*) malloc(length*sizeof(tipoInt));
				fread(tmpnidl,sizeof(tipoInt),length,pFile);

				//max->info->nl_a_S = (tipoInt*)malloc(max->info->nl_a_S_size*(jp+1)*sizeof(tipoInt));
				//nlasmalloc++;

				//cout << " file opened ";

				if(max->info->nl_a_S_type == 8){
					while(temp_size_8 < max->info->nl_a_S_size){
						temp_size_8 = temp_size_8*2;
						free(temp_nl_a_S_8);
						temp_nl_a_S_8=(uint8_t*)malloc(temp_size_8*(jp+1)*sizeof(uint8_t));
					}
					max->info->nl_a_S.nl_a_S_8=temp_nl_a_S_8;
					fread(max->info->nl_a_S.nl_a_S_8,sizeof(uint8_t),max->info->nl_a_S_size*(jp+1),pFile);
				}
				else if(max->info->nl_a_S_type == 16){
					while(temp_size_16 < max->info->nl_a_S_size){
						temp_size_16 = temp_size_16*2;
						free(temp_nl_a_S_16);
						temp_nl_a_S_16=(uint16_t*)malloc(temp_size_16*(jp+1)*sizeof(uint16_t));
					}
					max->info->nl_a_S.nl_a_S_16=temp_nl_a_S_16;
					fread(max->info->nl_a_S.nl_a_S_16,sizeof(uint16_t),max->info->nl_a_S_size*(jp+1),pFile);
				}
				else if(max->info->nl_a_S_type == 32){
					while(temp_size_32 < max->info->nl_a_S_size){
						temp_size_32 = temp_size_32*2;
						free(temp_nl_a_S_32);
						temp_nl_a_S_32=(uint32_t*)malloc(temp_size_32*(jp+1)*sizeof(uint32_t));
					}
					max->info->nl_a_S.nl_a_S_32=temp_nl_a_S_32;
					fread(max->info->nl_a_S.nl_a_S_32,sizeof(uint32_t),max->info->nl_a_S_size*(jp+1),pFile);
				}

				//cout << "done reading";

				new_allocation=false;
				fclose (pFile);
				created=false;
				remove(file_path);

				//cout << " done removing " << endl;

			}
			else{
				tmpnidl=max->info->IDlist;
				created=false;
			}
		}

		/*cout << "here ok IL? 2.1 "<<endl;
		cout << p->IL[1].item << endl;*/

    #ifndef testexpanded
        // test pattern on the permutations
        double psi_bound_log = getLogPsiExact(max->support);

        #ifdef enable_bounds

            psi_bound_log = min( computeLogPValue(max->support , max->info->max_a_S ) , computeLogPValue(max->support , max->info->min_a_S) );

            if(psi_bound_log >= getLogPsi(max->support)){
                  inferior_psi_bound_happened++;
            }
        #endif

        if(psi_bound_log >= getLogPsi(suppMin)){
            inferior_psi_bound_happened_neglect++;

                /*cout << "PsiBound: max->support = " << max->support <<
                " max->info->max_a_S = " << max->info->max_a_S <<  " max->info->min_a_S = " << max->info->min_a_S << endl;
                cout << "PsiBound: psi_bound_log =  " << psi_bound_log << " getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;*/

        }
        else{
            testPattern(psi_bound_log , max , NULL);

            /*cout << "tested: max->support = " << max->support <<
            " max->info->max_a_S = " << max->info->max_a_S <<  " max->info->min_a_S = " << max->info->min_a_S << endl;
            cout << "tested: psi_bound_log =  " << psi_bound_log << " getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;*/
        }
    #endif

		if (max->info->coreIndex==0 || suppMin == max->support){	///there is no need to generate the ppce

			#ifdef trace_state
			program_state = -2;
			#endif

			toOutput=findCloCI(tmpnidl, max);	///find the closure of max
			produced++;




		}
		else{

			//temptime = get_cpu_time();

			#ifdef trace_state
			program_state = -3;
			#endif

			fillHTCI(tmpnidl, max, q);		///fill the header table to generate the ppc-e
			produced++;

		/*cout << "here ok IL? A " << endl;
		cout << p->IL[1].item << endl;*/
			#ifdef trace_state
			program_state = -4;
			#endif

			toOutput=findPpceCI(max, q, q_res);//, qMin);	///closure of max

			//timekeepers[0]+= get_cpu_time() - temptime;

		}


		//temptime = get_cpu_time();

		/*cout << "here ok IL? B " << endl;
		cout << p->IL[1].item << endl;*/


		#ifdef trace_state
		program_state = -5;
		#endif

		clearHT(max->info->minItem);		///now clear the HT for next ppc-e generation

		/*cout << "here ok IL? 2.2 " << endl;
		cout << p->IL[1].item << endl;


		cout << "here ok 7" << endl;*/


		//timekeepers[1]+= get_cpu_time() - temptime;

		#ifdef trace_state
		program_state = -51;
		#endif


		/* SignificantMiner */

		/*for(int aj=0; aj < todelete->size(); aj++){
			cout << "(2) freed memory at " << todelete->at(aj)->a_S << endl;
			free(todelete->at(aj)->a_S);
			living_nodepointers--;
			nextfree_a_S=0;
		}
		todelete->clear();*/

		resetToDel();

		#ifdef trace_state
		program_state = -52;
		#endif

		resetListIndex();

		/* SignificantMiner */

		/*cout << "here ok 7.5" << endl;*/

		#ifdef trace_state
		program_state = -53;
		#endif



		nextfree=manualmem;	///free the manualmem


		/*cout << "here ok 7.75 " << toOutput[1] << endl;
		cout << "here ok 7.75 " << p->IL << endl;
		cout << "here ok 7.75 " << p->IL[toOutput[1]].item << endl;*/

		for (int i=1; i<=toOutput[0];i++){
			toOutput[i]=p->IL[toOutput[i]].item;
		}


		//q_res->printLowestHeapLevel();

		// if the extracted itemset can be significant
    //cout << "max extracted! max->support " << max->support << " max->info->a_S_ " << max->info->a_S_ << endl;
		if(store_results && (max->info->p_value < getPsi(suppMin) || (getPsi(suppMin) == 0.0 && computeLogPValue(max->support , max->info->a_S_) < getLogPsi(suppMin) ))){
			//cout << "max p_value: " << max->info->p_value << " max a_S " << max->info->a_S_ << endl;
			//cout << "RES QUEUE INSERTION: it can be significant since getPsi(suppMin) = " << getPsi(suppMin) << endl;
			// compose the full result itemset
			ResultPattern* toOutput_pattern=(ResultPattern*)malloc(sizeof(ResultPattern));
			toOutput_pattern->itemset = (tipoInt*)calloc((toOutput[0] + max->info->prefixCI[0] + 1),sizeof(tipoInt));
			toOutput_pattern->support = max->support;
			toOutput_pattern->a_S = max->info->a_S_;
			toOutput_pattern->p_value = max->info->p_value;
			toOutput_pattern->log_p_value = computeLogPValue(max->support , max->info->a_S_);
			toOutput_pattern->itemset[0] = toOutput[0] + max->info->prefixCI[0];

      //cout << "   pattern p_value: " << toOutput_pattern->p_value << " pattern log_p_value " << toOutput_pattern->log_p_value << endl;
      //cout << "   getPsi(suppMin): " << getPsi(suppMin) << " getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;

			//cout << "toOutput_pattern->log_p_value = " << toOutput_pattern->log_p_value << endl;

			for (int i=1; i<=max->info->prefixCI[0];i++){
				toOutput_pattern->itemset[i]=max->info->prefixCI[i];
			}
			for (int i=1; i<=toOutput[0];i++){
				toOutput_pattern->itemset[i+max->info->prefixCI[0]]=toOutput[i];
			}

			//q_res->printLowestHeapLevel();
			// insert it in the queue
			q_res->insert(toOutput_pattern);
			//cout << "inserted!" << endl;
			//q_res->printLowestHeapLevel();
		}

		#ifdef trace_state
		program_state = -54;
		#endif

		// output all itemset in the queue which are surely significant
		//while(q_res->getMin() > -0.5 && q_res->getMin() <= getPsi(s_supp)){
		while(q_res->getInqueue() > 0 && q_res->getLogMin() <= getLogPsi(s_supp)){

			//cout << "getLogPsi(s_supp) " << getLogPsi(s_supp) << endl;
			//cout << "get min returned " << q_res->getLogMin() << endl;

			//q_res->printLowestHeapLevel();

			ResultPattern* toOutput_pattern = q_res->removeMin();

			//cout << " remove min done! " << endl;
			//q_res->printLowestHeapLevel();

      if(SOFT_DEBUG){
        cout << "New significant pattern found: " << endl;
        cout << "getLogPsi(s_supp): " << getLogPsi(s_supp) << endl;
        cout << cloempty ;
        for (int i=1; i<=toOutput_pattern->itemset[0];i++){
          cout << toOutput_pattern->itemset[i] << " ";
        }
        cout << ": " << toOutput_pattern->support << " : " << toOutput_pattern->a_S << " : " << toOutput_pattern->p_value << " : " << toOutput_pattern->log_p_value << endl;
      }

			//cout << ": " << toOutput_pattern->support << " : " << toOutput_pattern->a_S << " : " << toOutput_pattern->p_value << endl;

			///print the closed itemset in output file
      #ifdef write_results_to_file
			outfile << cloempty ;	///in every closed itemset there is Clo(emptyset)
			for (int i=1; i<=toOutput_pattern->itemset[0];i++){
				outfile << toOutput_pattern->itemset[i] << " ";
			}
			outfile << ": " << toOutput_pattern->support << " : " << toOutput_pattern->a_S << " : " << toOutput_pattern->p_value << " : " << toOutput_pattern->log_p_value << endl;
      #endif
      significant_itemsets++;


  		// update statistics
  		inOutput[toOutput_pattern->itemset[0]+cloemptyLenght]++;

      //cout << "increased stat of index "<< toOutput_pattern->itemset[0]+cloemptyLenght << endl;
      //cout << " now it is " << inOutput[toOutput_pattern->itemset[0]+cloemptyLenght] << endl;


      free(toOutput_pattern->itemset);
      free(toOutput_pattern);

      K_significant_patterns = K_significant_patterns - 1;
      //cout << "K_significant_patterns decreased to " << K_significant_patterns << endl;


		}

		#ifdef trace_state
		program_state = -55;
		#endif

    /* Top-K strategy */
    while(K_significant_patterns > 0 && q_res->getKElements() >= K_significant_patterns){
      int new_significance_level = q_res->getMaxIndex();
      new_significance_level++;
      if(new_significance_level > suppMin){
        suppMin = new_significance_level;
        suppMinCurr = new_significance_level;
      }

      cout << "Top-K strategy update! ";
      /*cout << "K_significant_patterns " << K_significant_patterns << endl;
      cout << "q_res->getKElements() returned " << q_res->getKElements() << endl;
      cout << "q_res->getMaxIndex() returned " << q_res->getMaxIndex() << endl;*/
      cout << " s_supp " << s_supp;
      cout << " new suppMin " << suppMin << endl;
      /*cout << "getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;
      cout << "getPsi(suppMin) " << getPsi(suppMin) << endl;*/

      delta = getPsi(suppMin);
      q_res->removeMax();
    }

		// remove from the queue all non-significant patterns
		//while(q_res->getMax() > -0.5 && q_res->getMax() > getPsi(suppMin) && q_res->getMin() <= getPsi(suppMin)){
		while(q_res->getInqueue() > 0 && q_res->getMax() > getPsi(suppMin)){
			//cout << "removing not more significant itemsets with p-value " << q_res->getMax() << endl;
			//q_res->printLowestHeapLevel();
			q_res->removeMax();
			//cout << " remove max done! " << endl;
			//q_res->printLowestHeapLevel();
      //cout << "new max p-value is now " << q_res->getMax() << endl;
		}

		//q_res->printLowestHeapLevel();



		/*cout << "here ok 8" << endl;*/



		/*if(max->support < 55){
			elapsed = get_cpu_time() - elapsed;
			cout << "int size: " << sizeof(tipoInt) << endl;
			cout << "int size: " << sizeof(short int) << endl;
			cout << "queue size: " << inqueue << endl;
			cout << "nodelist nodes in memory: " << nodelist_inqueue << endl;
			cout << "Elapsed total time: " << elapsed << endl;
			printTimes();
			cout << "Time for nodelist operations: " << time_keeper << endl;
			sleep(10);
			return 0;
		}*/

		/*cout << "here ok 8.5" << endl;*/

		/* p-value computation */
		//double* p_values = computePValues(max->support , max->info->a_S);
		//double p_value = p_values[0];
		#ifdef trace_state
		program_state = -6;
		#endif

		/*double min_p_value = getPsi(max->support);
		if(SOFT_DEBUG){
		cout << "support of itemset " << max->support << endl;
		cout << "a_S of itemset " << max->info->a_S[0] << endl;
		if(jp>0){
			cout << "a_S of first permutation " << max->info->a_S[1] << endl;
			cout << "psi of first permutation " << getPsi(max->info->a_S[1]) << endl;
		}
		//cout << "p-value of itemset " << p_value << endl;
		cout << "minimum attainable p-value of itemset " << min_p_value << endl;
		cout << "delta: " << delta << endl;
		cout << "minimum p-values: " << endl;}
		if(DEEP_DEBUG)
		printMinPvalues();*/
		/* SignificantMiner */
		// check if infrequent itemsets can be removed from queue

		/*cout << "here ok 8.7" << endl;*/


		#ifdef trace_state
		program_state = -7;
		#endif

		/* SignificantMiner */
			// check if infrequent itemsets can be removed from queue
			if(DEEP_DEBUG)
			cout << "checking to remove infrequent itemsets..." << endl;
			//Itemset* todel;
			//todel = qMin->getMin();
			//temp_time = get_cpu_time();
			//while(todel != NULL && todel->support < suppMin){
			while(q->getMin() > 0 && q->getMin() < suppMin && q->getMax() >= suppMin){
				//.cout << "free itemset (11) with adress "<< todel <<endl;
				//living_itemsets--;
				//free(todel);

				//todel=qMin->removeMin();
				//q->printLowestHeapLevel();
				q->removeMin();
				//free(todel->info->prefixCI);
				//free(todel->info->IDlist);
				//delete todel->info;
				//todel->info = NULL;
				//if(DEEP_DEBUG)
				//cout << "removed infrequent itemset with support " << todel->support << endl;

				//.cout << "free itemset (1) with adress "<< todel <<endl;
				//free(todel);
				//living_itemsets--;
				//todel = qMin->getMin();
			}
			//if(todel!=NULL){
				//.cout << "free itemset (12) with adress "<< todel <<endl;
				//living_itemsets--;
				//free(todel);
			//}
			//time_keeper += get_cpu_time() - temp_time;
			/* SignificantMiner */

		///update statistics

		//inOutput[toOutput[0]+cloemptyLenght+max->info->prefixCI[0]]++;

		#ifdef trace_state
		program_state = -8;
		#endif

		///now delete the startNodeList of max and the nodeIDlist
		free(tmpnidl);
		if (created)
			free(max->info->IDlist);
		max->info->IDlist = NULL;
		free(max->info->prefixCI);
		max->info->prefixCI = NULL;
		if(!new_allocation){
			max->info->nl_a_S.nl_a_S_t=NULL;
			max->info->nl_a_S_type = 0;
		}
		delete max->info;
		//.cout << "free itemset (5) with adress "<< max <<endl;
		free(max);
		living_itemsets--;

		//if (suppMinCurr==q->heap->support)	firstOfSuppMin=true;
		//if (suppMinCurr > q->heap->support)	firstOfSuppMin=true;

		/*cout << "here ok 9" << endl;*/


		//timekeepers[2]+= get_cpu_time() - temptime;


	}


  cout << "Number of patterns found during exploration: " << significant_itemsets << endl;




		#ifdef trace_state
		program_state = -541;
		#endif

		// output all itemset in the queue which are surely significant
		//while(q_res->getMin() > -0.5 && q_res->getMin() <= getPsi(s_supp)){
		while(q_res->getInqueue() > 0 && q_res->getLogMin() <= getLogPsi(suppMin)){

			/*cout << "getLogPsi(s_supp) " << getLogPsi(s_supp) << endl;
      cout << "getLogPsi(suppMin) " << getLogPsi(suppMin) << endl;
			cout << "get min returned " << q_res->getLogMin() << endl;*/

			//q_res->printLowestHeapLevel();

			ResultPattern* toOutput_pattern = q_res->removeMin();

			//cout << " remove min done! " << endl;
			//q_res->printLowestHeapLevel();

      if(SOFT_DEBUG){
        cout << "New significant pattern found: " << endl;
        cout << "getLogPsi(s_supp): " << getLogPsi(s_supp) << endl;
        for (int i=1; i<=toOutput_pattern->itemset[0];i++){
          cout << toOutput_pattern->itemset[i] << " ";
        }
        cout << ": " << toOutput_pattern->support << " : " << toOutput_pattern->a_S << " : " << toOutput_pattern->p_value << " : " << toOutput_pattern->log_p_value << endl;
      }

			//cout << ": " << toOutput_pattern->support << " : " << toOutput_pattern->a_S << " : " << toOutput_pattern->p_value << endl;

			///print the closed itemset in output file
      #ifdef write_results_to_file
			outfile << cloempty ;	///in every closed itemset there is Clo(emptyset)
			for (int i=1; i<=toOutput_pattern->itemset[0];i++){
				outfile << toOutput_pattern->itemset[i] << " ";
			}
			outfile << ": " << toOutput_pattern->support << " : " << toOutput_pattern->a_S << " : " << toOutput_pattern->p_value << " : " << toOutput_pattern->log_p_value << endl;
      #endif
      significant_itemsets++;
			free(toOutput_pattern->itemset);
			free(toOutput_pattern);
		}

		#ifdef trace_state
		program_state = -551;
		#endif

    /*cout << "queue of results size: " << q_res->getInqueue() << endl;
    cout << "q_res->getMax() " << q_res->getMax() << endl;
    cout << "q_res->getMin() " << q_res->getMin() << endl;
    cout << "q_res->getLogMin() " << q_res->getLogMin() << endl;
    cout << "getPsi(s_supp) " << getPsi(s_supp) << endl;

    if(q_res->getInqueue() > 0){
      for(int aj=0; aj<q_res->data.size(); aj++){
        if(q_res->data[aj].size() > 0){
          cout << " q_res->data["<<aj<<"].size() " << q_res->data[aj].size() << endl;
        }
      }
    }*/

	#ifdef trace_state
	program_state = -9;
	#endif

	///now we must extract all the itemset in the priority queue with the same support of the k-th extracted without generating their ppc-e; we must remember what are these items because

	notExpList toExpLater=NULL;	///list of itemset not expanded: it never contains an itemset with coreIndex=0
	notExpanded** nextfield=&toExpLater;	///pointer were the next notExpanded must be inserted

	if(SOFT_DEBUG)
	cout << "Extracting all itemsets with minimum support! " << endl;

	while(false){//q->heap->support == suppMinCurr){

		max=q->removeMax();	///remove the first itemset from the priority queue

		if(max == NULL)
			cout << "Extracted max IS NULL! " << endl;

		if (max->info->typeList==0){
			tmpnidl=changeNIDL(max->info->IDlist);
			created=true;
		}
		else{
			if (max->info->typeList==2){
				//cout << "reading from file 2" << endl;
				cout<<"RD " << max->info->fileName <<endl;
				FILE * pFile;

				for(int aj=0; aj<20; aj++) file_path[aj]=' ';
				sprintf(file_path, "%d.l", max->info->fileName);

				pFile = fopen ( file_path , "rb" );
				tipoInt length;

				fread(&length,sizeof(tipoInt),1,pFile);
				tmpnidl=(tipoInt*) malloc(length*sizeof(tipoInt));
				fread(tmpnidl,sizeof(tipoInt),length,pFile);

				if(max->info->nl_a_S_type==8){
  					fread (max->info->nl_a_S.nl_a_S_8, sizeof(uint8_t), max->info->nl_a_S_size*(jp+1), pFile );
					nlinqueue_8+=max->info->nl_a_S_size;
  				}
  				else
  					if(max->info->nl_a_S_type==16){
  						fread (max->info->nl_a_S.nl_a_S_16, sizeof(uint16_t), max->info->nl_a_S_size*(jp+1), pFile );
  						nlinqueue_16+=max->info->nl_a_S_size;
  					}
  					else{
  						fread (max->info->nl_a_S.nl_a_S_32, sizeof(uint32_t), max->info->nl_a_S_size*(jp+1), pFile );
  						nlinqueue_32+=max->info->nl_a_S_size;
  					}

				fclose (pFile);
				created=false;
			}
			else{
				tmpnidl=max->info->IDlist;
				created=false;
			}
		}

		#ifdef trace_state
		program_state = -10;
		#endif

		toOutput=findCloCI(tmpnidl, max);	///find the closure of max
		produced++;

		clearHT(max->info->minItem);		///now clear the HT for closure operation

		///now we replace the items in toOutput with their ID
		for (int i=1; i<=toOutput[0];i++){
			toOutput[i]=p->IL[toOutput[i]].item;
		}

		///update statistics
		inOutput[toOutput[0]+cloemptyLenght+max->info->prefixCI[0]]++;

		///if the support is >1 and coreIndex>0...
		if ((max->support>1)&&(max->info->coreIndex>0)){
			///... insert max in the notExpList
			(*nextfield)=new notExpanded();
			(*nextfield)->its=max;
			nextfield=&((*nextfield)->next);
		}

		///write on the output file
    #ifdef write_results_to_file
		outfile << cloempty ;	///in every closed itemset there is Clo(emptyset)
		for (int i=1; i<=max->info->prefixCI[0];i++){
			outfile << p->IL[max->info->prefixCI[i]].item << " ";
		}
		for (int i=1; i<=toOutput[0];i++){
			outfile << toOutput[i] << " ";
		}
		//outfile << ": " << max->support << endl;
		/* SignificantMiner */
		outfile << ": " << max->support << " : " << max->info->a_S[0] << endl;
    #endif

		///now delete the startNodeList of max

		if (created)
			free(tmpnidl);
	}

	total_time = get_cpu_time() - total_time;
  time_withoutpatricia = get_cpu_time() - time_withoutpatricia;


	//cout << "\nMinimum computed p-values: " <<endl;
	//printMinPvalues();

	#ifdef LAMP

	cout << "itemsets_LAMP_count_naive.size() = " << itemsets_LAMP_count_naive.size() << endl;
	cout << "itemsets_LAMP_count.size() = " << itemsets_LAMP_count.size() << endl;

	for(int aj=itemsets_LAMP_count_naive.size()-2; aj>0; aj--){
		itemsets_LAMP_count_naive[aj] = itemsets_LAMP_count_naive[aj] + itemsets_LAMP_count_naive[aj+1];
	}

	cout << "psi["<<(itemsets_LAMP_count.size()-1)<<"] = "<<getPsi(itemsets_LAMP_count.size()-1)<<" itemsets_LAMP_count["<<(itemsets_LAMP_count.size()-1)<<"] = " << itemsets_LAMP_count[itemsets_LAMP_count.size()-1] << " - itemsets_LAMP_count_naive["<<(itemsets_LAMP_count.size()-1)<<"] = " << itemsets_LAMP_count_naive[itemsets_LAMP_count.size()-1] << endl;
	for(int aj=itemsets_LAMP_count.size()-2; aj>=0; aj--){
		itemsets_LAMP_count[aj] = itemsets_LAMP_count[aj] + itemsets_LAMP_count[aj+1];
		cout << "psi["<<(aj)<<"] = "<<getPsi(aj)<<" itemsets_LAMP_count["<<aj<<"] = " << itemsets_LAMP_count[aj] << " - itemsets_LAMP_count_naive["<<aj<<"] = " << itemsets_LAMP_count_naive[aj] << endl;
	}

	int minSupp_LAMP = 2;
	while(itemsets_LAMP_count_naive[minSupp_LAMP+1] > (0.05/getPsi(minSupp_LAMP+1))){
		minSupp_LAMP++;
	}

	cout << "\nLAMP minimum support: " << minSupp_LAMP <<endl;

	minSupp_LAMP = 2;
	while(itemsets_LAMP_count[minSupp_LAMP+1] > (0.05/getPsi(minSupp_LAMP+1)) ){
		minSupp_LAMP++;
	}


	cout << "\nUpdated LAMP minimum support: " << minSupp_LAMP <<endl;

	#endif

	cout << "\nCurrent minimum support: " << suppMinCurr <<endl;
	/*cout << "Global minimum support: " << suppMin <<endl;
	cout << "\nNumber of itemset required: " ;
	if ( emptyIsNotClosed ) cout << (k+1) <<endl;
	else cout << k << endl;*/
	cout << "Number of itemset explored: " ;
	if ( emptyIsNotClosed ) cout << (produced+1) <<endl;
	else cout << produced <<endl;
	cout << "tested " << tested  << endl;
	cout << "underbound " << inferior_psi_bound_happened_neglect  << endl;
	cout << "\n------------- PRINTING STATISTICS ------------\n" << endl;
	cout << "S.I. Length\tNumber of S.I" << endl;
	for (int i=1; i<=max_trans_length; i++){
		if (inOutput[i]>0){
			cout <<i<<"\t\t\t" << inOutput[i] << endl;
		}
	}

	cout << "Number of significant itemsets: " << significant_itemsets << endl;

	cout << "Number of itemset that have been pruned: " << pruned_patterns << endl;

	freeToDel();

  #ifdef print_list_stats
	// compute list statistics
	double tot = 0.0;
	double avg = 0.0;
	for(int aj=0; aj < list_stat.size(); aj++){
		tot += list_stat[aj];
	}
	for(int aj=0; aj < list_stat.size(); aj++){
		avg += list_stat[aj] * aj / tot;
	}
	tot = 0.0;
	double avg_op = 0.0;
	for(int aj=0; aj < list_op_stat.size(); aj++){
		tot += list_op_stat[aj];
	}
	for(int aj=0; aj < list_op_stat.size(); aj++){
		avg_op += list_op_stat[aj] * aj / tot;
	}
  #endif

	/* SignificantMiner */
	// end here
	cout << "Total running time time: " << total_time << endl;
	cout << "Time to build patricia trie: " << time_pat << endl;
	printTimes();
  cout << "Time for exploration: " << time_withoutpatricia << endl;
	cout << "Time to process as values: " << as_processing_time << endl;
	//cout << "Time for a_S operations: " << time_keeper << endl;
	/*cout << "Time for ppce 1 operations: " << timekeepers[0] << endl;
	cout << "Time for ppce 0 operations: " << timekeepers[1] << endl;
	cout << "Time for cleanup operations: " << timekeepers[2] << endl;*/
	cout << "queue size: " << q->getInqueue() << endl;
	cout << "max queue size: " << maxinqueue << endl;
	cout << "alive itemsets in memory: " << living_itemsets << endl;
	cout << "alive nodePointers in memory: " << living_nodepointers << endl;
	cout << "alive infos in memory: " << living_infos << endl;
	cout << "deleted infos: " << info_deleted << endl;
	cout << "sizeof(Itemset) " << sizeof(Itemset) << endl;
  #ifdef print_list_stats
	cout << "avg nlists length " << avg << endl;
	cout << "avg nlists operations " << avg_op << endl;
  #endif
	//cout << "q->maxSize " << q->maxSize << endl;
	//cout << "qMin->maxSize " << qMin->maxSize << endl;
	//cout << "max number of nl_a_S in memory " << maxnlinqueue << endl;
	cout << "max number of info in memory " << maxlivinginfos << endl;
	cout << "memory for patricia trie: " << pat_tree_ram << " MByte" << endl;
	cout << "memory usage for permutation matrix: " << perm_matrix_space << " MByte" << endl;
	cout << "memory usage for count_a_S_element: " << count_a_S_element_space << " MByte" << endl;
	cout << "total peak memory usage: " << measurePeakMemory() /*memory_peak*/ << " MByte" << endl;
	cout << "nlasmalloc " << nlasmalloc <<endl;
	//sleep(10);

  // debug results queue
  q_res->printElements();

	std::ofstream outfile_;

  	outfile_.open("data_sm.csv", std::ios_base::app);
  	outfile_ << jp << ";" << memory_peak << ";" << total_time << ";" << suppMin << ";" << delta << ";" << (produced+1) << ";" << alpha << endl;

	outfile_.close();

  	#ifdef trace_state
	program_state = -91;
	#endif

	cout << "Program closed"<< endl;

  kill(getpid(), SIGKILL);

	return 0;


	//cout << "Finished extracting all itemsets with minimum support! " << endl;

	int otherIts;	///indicate if you want other itemset or not
	notExpanded* tmptoexp;	///pointer to a node that were not expanded
	notExpanded* prev=NULL;	///to delete the list
	char choice;

	while (produced<k_max && (/*(q->last!=-1) ||*/(toExpLater!=NULL)) && (!finish)) {

		cout << "\n************************************" << endl;
		cout << "\nDo you want more itemsets? [y/n]" << endl;
		while( choice != 'y' && choice != 'n' ){
			cin >> choice;
		}
//		cin >> otherIts;
//		if (otherIts==1){
		if (choice == 'y'){
			choice = 'a';
			cout << "\nHow many?" << endl;
			cin >> otherIts;

			//if too much itemsets desired, set a limit..
			if ( k+otherIts > k_max){
				cout << "\nTOO MANY ITEMSETS REQUIRED!\nPRODUCING ALL THE ITEMSETS TO REACH ";
				if ( emptyIsNotClosed ) cout << (k_max+1);
				else cout << k_max ;
				cout << " ( = k_global ) ITEMSETS IN OUTPUT\n" << endl;
				otherIts = k_max - k;
			}

			int j=0;	///number of itemset expanded NOW

			///we now take the list of itemset that were produced but not expanded and expand the first otherIts of them, so that at the end all these items were expanded

			while ((j<otherIts) && (toExpLater!=NULL)){
				tmptoexp=toExpLater;

				if (tmptoexp->its->info->typeList==0){
					tmpnidl=changeNIDL(tmptoexp->its->info->IDlist);
					created=true;
				}
				else{
					tmpnidl=tmptoexp->its->info->IDlist;
					created=false;
				}
				fillHTCI(tmpnidl, tmptoexp->its, q);		///fill the header table to generate the ppc-e
				toOutput=findPpceCI(tmptoexp->its, q, q_res);//, qMin);	///closure of tmptoexp
				j++;

				clearHT(tmptoexp->its->info->minItem);		///now clear the HT for closure operation

						///now delete the startNodeList of max and the nodeIDlist
				free(tmpnidl);
				if (created)
					free(tmptoexp->its->info->IDlist);
				free(tmptoexp->its->info->prefixCI);
				prev=tmptoexp;

				delete tmptoexp->its->info;
				delete prev;
				tmptoexp=tmptoexp->next;
				toExpLater=tmptoexp;

			}
			k+=otherIts;
			if (k<=produced){
				cout << "\nTHE REQUIRED (TOTAL) NUMBER OF ITEMSETS HAS ALREADY BEEN PRODUCED!\n" << endl;
			}
			else {
				///we expand ALL the itemset not previously expanded

				tmptoexp=toExpLater;
				while(tmptoexp!=NULL){

					if (tmptoexp->its->info->typeList==0){
						tmpnidl=changeNIDL(tmptoexp->its->info->IDlist);
						created=true;
					}
					else{
						tmpnidl=tmptoexp->its->info->IDlist;
						created=false;
					}
					fillHTCI(tmpnidl, tmptoexp->its, q);		///fill the header table to generate the ppc-e
					toOutput=findPpceCI(tmptoexp->its, q, q_res);//, qMin);	///closure of tmptoexp

					clearHT(tmptoexp->its->info->minItem);		///now clear the HT for closure operation

							///now delete the startNodeList of max and the nodeIDlist
					free(tmpnidl);
					if (created)
						free(tmptoexp->its->info->IDlist);
					free(tmptoexp->its->info->prefixCI);
					prev=tmptoexp;
					delete tmptoexp->its->info;
					delete prev;
					tmptoexp=tmptoexp->next;
					toExpLater=tmptoexp;
//					j++;
				}
				///now update the toExpLater list
				toExpLater=NULL;
				///first of all: use the trivial heuristic for suppMinCurr

				///trivial heuristic on suppMinCurr
				/*tipoInt suppMinCurrFut=suppMin;
				for (tipoInt l=suppMin; l<maximumSupport; l++){
					if (utilNodeArray[l]>=k){
						suppMinCurrFut=l;
					}
				}
				suppMinCurr=suppMinCurrFut;*/

				while((produced<k) /*&& (q->last!=-1)*/ && !(finish&&firstOfSuppMin)){

					max=q->removeMax();	///remove the first itemset from the priority queue


					if(max == NULL)
						cout << "Extracted max IS NULL! " << endl;



					if (max->info->typeList==0){
						tmpnidl=changeNIDL(max->info->IDlist);
						created=true;
					}
					else{
						tmpnidl=max->info->IDlist;
						created=false;
					}
					if (max->info->coreIndex==0){	///there is no need to generate the ppce
						toOutput=findCloCI(tmpnidl, max);	///find the closure of max
						produced++;
					}
					else{
						fillHTCI(tmpnidl, max, q);		///fill the header table to generate the ppc-e
						produced++;
						toOutput=findPpceCI(max, q, q_res);//, qMin);	///closure of max
					}

					clearHT(max->info->minItem);		///now clear the HT for next ppc-e generation
					/* SignificantMiner */
					/*for(int aj=0; aj < todelete->size(); aj++){
						//cout << "(2) freed memory at " << todelete->at(aj)->a_S << endl;
						free(todelete->at(aj)->a_S);
						living_nodepointers--;
					}
					todelete->clear();*/
					resetToDel();
					resetListIndex();
					/* SignificantMiner */
					nextfree=manualmem;	///free the manualmem

					for (int i=1; i<=toOutput[0];i++){
						toOutput[i]=p->IL[toOutput[i]].item;
					}

					///print the closed itemset in output file
          #ifdef write_results_to_file
					outfile << cloempty ;	///in every closed itemset there is Clo(emptyset)
					for (int i=1; i<=max->info->prefixCI[0];i++){
						outfile << p->IL[max->info->prefixCI[i]].item << " ";
					}
					for (int i=1; i<=toOutput[0];i++){
						outfile << toOutput[i] << " ";
					}
					//outfile << ": " << max->support << endl;
					/* SignificantMiner */
					outfile << ": " << max->support << " : " << max->info->a_S[0] << endl;
          #endif

					///update statistics

					inOutput[toOutput[0]+cloemptyLenght+max->info->prefixCI[0]]++;

					///now delete the startNodeList of max and the nodeIDlist

					free(tmpnidl);
					if (created)
						free(max->info->IDlist);
					free(max->info->prefixCI);
					delete max->info;
					//.cout << "free itemset (6) with adress "<< max <<endl;
					free(max);
					living_itemsets--;

					//if (suppMinCurr==q->heap->support)	firstOfSuppMin=true;
				}

				///now we must extract all the itemset in the priority queue with the same support of the k-th extracted without generating their ppc-e; we must remember what are these items because

				nextfield=&toExpLater;	///pointer were the next notExpanded must be inserted

				//while(q->heap->support == suppMinCurr){
				while(q->getMax() == suppMinCurr){

					max=q->removeMax();	///remove the first itemset from the priority queue


					if(max == NULL)
						cout << "Extracted max IS NULL! " << endl;

					if (max->info->typeList==0){
						tmpnidl=changeNIDL(max->info->IDlist);
						created=true;
					}
					else{
						tmpnidl=max->info->IDlist;
						created=false;
					}

					toOutput=findCloCI(tmpnidl, max);	///find the closure of max
					produced++;
					clearHT(max->info->minItem);		///now clear the HT for closure operation

					///now we replace the items in toOutput with their ID
					for (int i=1; i<=toOutput[0];i++){
						toOutput[i]=p->IL[toOutput[i]].item;
					}

					///update statistics
					inOutput[toOutput[0]+cloemptyLenght+max->info->prefixCI[0]]++;

					///if the support is >1 and coreIndex>0...
					if ((max->support>1)&&(max->info->coreIndex>0)){
						///... insert max in the notExpList
						(*nextfield)=new notExpanded();
						(*nextfield)->its=max;
						(*nextfield)->next=NULL;
						nextfield=&((*nextfield)->next);
					}

					///write on the output file
          #ifdef write_results_to_file
					outfile << cloempty ;	///in every closed itemset there is Clo(emptyset)
					for (int i=1; i<=max->info->prefixCI[0];i++){
						outfile << p->IL[max->info->prefixCI[i]].item << " ";
					}
					for (int i=1; i<=toOutput[0];i++){
						outfile << toOutput[i] << " ";
					}
					//outfile << ": " << max->support << endl;
					/* SignificantMiner */
					outfile << ": " << max->support << " : " << max->info->a_S[0] << endl;
          #endif

					///now delete the startNodeList of max

					if (created)
						free(tmpnidl);
				}
			}
		}
		else{
			return 0;
		}
	}

//	double ns=1000000;

	cout << "\nCurrent minimum support = " << suppMinCurr <<endl;
	cout << "Global minimum support = " << suppMin <<endl;
	cout << "\nNumber of itemset required: " ;
		if ( emptyIsNotClosed ) cout << (k+1) <<endl;
		else cout << k << endl;
		cout << "Number of itemset produced: " ;
		if ( emptyIsNotClosed ) cout << (produced+1) <<endl;
		else cout << produced <<endl;
	cout << "\n------------- PRINTING STATISTICS ------------\n" << endl;
	cout << "F.C.I. Length\tNumber of F.C.I\n";
	for (int i=1; i<=max_trans_length; i++){
		if (inOutput[i]>0){
			cout <<i<<"\t\t\t" << inOutput[i] << "\n";
		}
	}
	cout << endl;
}
