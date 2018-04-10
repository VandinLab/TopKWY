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

//#include <vector>
#include <vector>
#include <deque>

#define sizeIts sizeof(tipoInt)+sizeof(ItsInfo*)

int nodelist_inqueue;
int living_itemsets;
int living_infos;
int living_nodepointers;
int info_deleted;
//int nlinqueue;
int nlinqueue_8;
int nlinqueue_16;
int nlinqueue_32;
tipoInt nlasmalloc;

tipoInt s_supp;
tipoInt max_queue_size;
tipoInt max_queue_last;

int queue_program_state;

union nl_list_pointer{
	uint32_t *nl_a_S_32;
	uint16_t *nl_a_S_16;
	uint8_t *nl_a_S_8;
	tipoInt *nl_a_S_t;
};


/**class for the infos associated with a closed itemset memorized in a queue
 */

class ItsInfo{
	public:
		NodeIDList IDlist;	///list of ID of PatriciaNode in which the minimum item of the Itemset is contained
		tipoInt minItem;	///minimum item of the Itemset: it is contained in all the nodes in the IDlist
		tipoInt coreIndex;	///coreIndex of the itemset: it is memorized because we need only to generate ppc-e with items>than coreIndex
		short typeList;		///typeList is =0 iff the IDlist is referred to minItem nodes; typeList is =1 iff the IDlist is referred to coreIndex nodes; typeList is =2 if the list is memorized on file;
		tipoInt tagliamem;		///for statistic: the total of the memory used for this item  // deprecated
		tipoInt sizelista;		///for statistic: the total of the memory used for the list of ID  // deprecated
		Transaction prefixCI; ///here we memorize the coreIndex-prefix of the itemset, so there is no need to memorize them in the intersection of the IDlist
		tipoInt fileName;		///if typeList is =2, this is the name of the file that contains the list  // deprecated
		/* TopKWY */
		tipoInt* a_S; /// a_S of the itemset // deprecated
		tipoInt max_a_S;
		tipoInt min_a_S;
		tipoInt a_S_;
		double p_value;  // deprecated
		//tipoInt* nl_a_S; /// permutated a_S associated to typelist=1 nodes
		nl_list_pointer nl_a_S;  // deprecated
		uint8_t nl_a_S_type;  // deprecated

		//std::vector< std::vector<int> > nl_nodes;
		std::vector<int> minnodes_indexes;
		std::vector<int> minnodes;

		tipoInt nl_a_S_size; /// number of nodes for typelist=1 nodes  // deprecated
		~ItsInfo(void);
		/* TopKWY */
};

ItsInfo::~ItsInfo(void) {
	if(DEEP_DEBUG)
	cout << "info is being deleted!" << endl;
	//cout << "(1) freed memory at " << a_S << endl;
	//free(a_S);

	if(nl_a_S_type > 0){
		if(nl_a_S_type == 8){
			free(nl_a_S.nl_a_S_8);
			nl_a_S.nl_a_S_8 = NULL;
			nlinqueue_8-=nl_a_S_size;
		}
		else
			if(nl_a_S_type == 16){
				free(nl_a_S.nl_a_S_16);
				nl_a_S.nl_a_S_16 = NULL;
				nlinqueue_16-=nl_a_S_size;
			}
			else
				if(nl_a_S_type == 32){
					free(nl_a_S.nl_a_S_32);
					nl_a_S.nl_a_S_32 = NULL;
					nlinqueue_32-=nl_a_S_size;
				}
		//nlinqueue-=nl_a_S_size;
		nlasmalloc--;
		nl_a_S_type = 0;
	}

	a_S = NULL;
	nodelist_inqueue-=nl_a_S_size;
	living_infos--;
	info_deleted++;
	if(prefixCI)
	free(prefixCI);
	if(IDlist)
	free(IDlist);
	if(DEEP_DEBUG)
	cout << "info deleted!" << endl;
}

///class for the element of the priority queue: we memorize only the support of the itemset and a pointer to ItsInfo
class Itemset{
	public:
		tipoInt support;	///support of the itemset
		ItsInfo* info;		///pointer to the info associated with the itemset
		void print();		///print the information associate with the itemset
		~Itemset(void);
};

Itemset::~Itemset(void) {
	if(DEEP_DEBUG)
	cout << "itemset is being deleted!" << endl;
	if(info != NULL)
		delete info;
	info = NULL;
	living_itemsets--;
}

///class for the priority queue implemented as a max heap
class Queue{
	public:
		Itemset* heap;	///pointer to the array that implements the heap
		int last;	/// index of last element in the array; it is one minus than the number of element in the queue
		int maxSize; 	/// maximum number of element in the heap:it represent the number of element allocated in memory
		void insert(Itemset* its);	///insert the itemset its in the heap
		Itemset* removeMax();		///remove the max element from the queue
		Queue();			///default constructor
		void doubleQueue();			///double the space for the queue
		void trimQueue();			///trim the space for the queue

		tipoInt notFreqInLists(ItemList* IL);	///return the number of items that are not frequent in the lists
		tipoInt spaceSaved();
		long long spaceUsedByLists();
};

///class to remember what itemset were extracted from the queue BUT not generated ppc-e

class notExpanded{
	public:
		Itemset* its;
		notExpanded* next;
};

///list of notExpanded itemset

typedef notExpanded* notExpList;


///class for the min-based priority queue used to limit the number of itemsets in queue: it is the same of Queue but the Itemset are ordered by increasing value of support
class QueueMin{
		public:
		Itemset* heap;	///pointer to the array that implements the heap
		int last;	/// index of last element in the array; it is one minus than the number of element in the queue
		int maxSize; 	/// maximum number of element in the heap:it represent the number of element allocated in memory
		void insert(Itemset* its);	///insert the itemset its in the heap
		Itemset* removeMin();		///remove the max element from the queue
		Itemset* getMin();		///remove the max element from the queue
		QueueMin();			///default constructor
		void doubleQueueMin();			///double the space for the queue
		void trimQueueMin();			///trim the space for the queue
};


///class for the min-max-based priority queue
class QueueMinMax{
		public:
		Itemset* heap;	///pointer to the array that implements the heap
		int last;	/// index of last element in the array; it is one minus than the number of element in the queue
		int maxSize; 	/// maximum number of element in the heap:it represent the number of element allocated in memory
		void insert(Itemset* its);	///insert the itemset its in the heap
		void removeMin();		///remove the min element from the queue
		Itemset* removeMax();		///remove the max element from the queue
		tipoInt getMin();		///get the support of the min element from the queue
		tipoInt getMax();		///get the support of the max element from the queue
		QueueMinMax();			///default constructor
		void doubleQueueMaxMin();			///double the space for the queue
		void printLowestHeapLevel();
	private:
		tipoInt min_support;
		tipoInt min_index;
};

// class for an output pattern
class ResultPattern{
	public:
	Transaction itemset;
	tipoInt support;
	tipoInt a_S;
	double p_value;
	double log_p_value;

};

///class for the min-max-based priority queue for results
class QueueMinMaxResult{
	public:
	ResultPattern* heap;	///pointer to the array that implements the heap
	int last;	/// index of last element in the array; it is one minus than the number of element in the queue
	int maxSize; 	/// maximum number of element in the heap:it represent the number of element allocated in memory
	void insert(ResultPattern* res);	///insert the itemset its in the heap
	void removeMax();		///remove the min element from the queue
	ResultPattern* removeMin();		///remove the max element from the queue
	double getMin();		///get the support of the min element from the queue
	double getMax();		///get the support of the max element from the queue
	QueueMinMaxResult();			///default constructor
	void doubleQueueMinMaxResult();			///double the space for the queue
	void printLowestHeapLevel();
	int getQueueSize();
	private:
	double max_p_value;
	tipoInt max_index;
};


///class for the min-max-based priority queue
class QueueMinMax_test{
		public:
		std::vector< std::vector<Itemset*> > data;
		void insert(Itemset* its);	///insert the itemset its in the heap
		void removeMin();		///remove the min element from the queue
		Itemset* removeMax();		///remove the max element from the queue
		Itemset* removeLow(); // remove an item from the frontier
		tipoInt getMin();		///get the support of the min element from the queue
		tipoInt getMax();		///get the support of the max element from the queue
		QueueMinMax_test(tipoInt max_support);			///default constructor
		tipoInt getInqueue();   // returns number of elements in the queue
		double low_factor;
		double getQueueMemory();
	private:
		tipoInt min_index;
		tipoInt max_index;
		tipoInt elements_stored;
};


///class for the min-max-based priority queue for results
class QueueMinMax_Restest{
		public:
		std::vector< std::vector<ResultPattern*> > data;
		std::vector< int > observed;
		void insert(ResultPattern* its);	///insert the ResultPattern its in the heap
		void removeMax();		///remove the min element from the queue
		ResultPattern* removeMin();		///remove the max element from the queue
		double getMin();		///get the support of the min element from the queue
		double getLogMin();		///get the support of the min element from the queue
		double getMax();		///get the support of the max element from the queue
		int getMaxIndex();		///get the support of the max element from the queue
		QueueMinMax_Restest(tipoInt max_support);			///default constructor
		tipoInt getInqueue();   // returns number of elements in the queue
		tipoInt getKElements();
		void insert_observed(double log_observed_pval);
		void printElements();
	private:
		tipoInt min_index;
		double min_value;
		double min_log_value;
		tipoInt max_index;
		double max_value;
		tipoInt elements_stored;
		tipoInt elements_observed;
		tipoInt max_index_obs;
};
