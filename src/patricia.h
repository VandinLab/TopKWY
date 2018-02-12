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

/** \ file patricia.h
* \ definition of classes used by the PatriciaTrie
*/

#ifndef _patricia_h_
#define _patricia_h_

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#define DEEP_DEBUG 0
#define SOFT_DEBUG 0

using namespace std;

///the definition is the following: there are 4 tipoInt in the PatriciaNode and two tipoInt* in the PatriciaNode (if sizeof(HashTable*)=sizeof(tipoInt)) and a PatriciaNode*; so if you want to allocate space for a node with n item, use (sizeNode+(n-1)*sizeof(tipoInt))

///NB: the definition is ok IF AND ONLY IF sizeof(tipoInt)>=sizeof(HashTable*) AND sizeof(tipoInt)>=sizeof(NodePointer*)

///this definition used for malloc


#define maxSizeNode 7

#define maxListSize 5242880

///this definition is used for the manual allocation of NodePointer: if we allocate a NodePointer with an intersection of length l, we advance nextfree of sizePnt+l;

//#define sizePnt 8 // the evil is here
#define sizePnt (sizeof(NodePointer)/sizeof(tipoInt)) // the evil is here
#define sizeNode sizeof(PatriciaNode)

/// dataset file

ifstream dataset;
string fileinput;
ifstream dataset1;
string c_fileinput;
ifstream c_labels;

/// defines the type for the integers used

typedef int tipoInt;

/// a transaction is an array of item, and every item is a tipoInt

typedef tipoInt* Transaction;

///GLOBAL VARIABLES and CONSTANT

///ONLY FOR DEBUG

tipoInt inqueue;	///number of itemsets inserted in queue
tipoInt numPnt; 	///number of pointer memorized for the itemsets in the queue
tipoInt numResize=0;

int pat_tree_ram;
int perm_matrix_space;
int count_a_S_element_space;


///total number of itemset produced

//tipoInt produced;
long produced;
long tested;
tipoInt producedOld;	///debug
tipoInt per=0;		///debug
tipoInt nofold=0;	///debug

tipoInt max_ram;
tipoInt max_ram_l1, max_ram_l2;
tipoInt max_nl_size;

int max_trans_length;

///array for the number of nodes that contains an item

tipoInt* numNPI;

///memory for manual allocation

tipoInt* manualmem;

///pointer to the next free position in manualmem

tipoInt* nextfree;

///last valid position of manualmem: used to check if manualmem is full

tipoInt* lastvalid;

///number of itemset inserted in queue (with the extracted, too)

//tipoInt inserted;
long inserted;

///number of itemset produced as significant

long significant_itemsets;

///cumulative number of items in the nodes of PatriciaTrie (without in\tems in intersections)

tipoInt totItem;

///hash load mean value

tipoInt hash_limit=7;

///maximum number of f.c.i. that can be produced in different iterations

tipoInt k_max;

///number of f.c.i. that is required in this iterations

tipoInt k;

///minimum support for the maximum number of f.c.i. that can be produced in output in one or more iterations

tipoInt suppMin;

///minimum support for the top-K f.c.i. that will be produced in output in this iteration

tipoInt suppMinCurr;

///this indicates if we have already cumulated utilNodeArray for the bound on suppMinCurr

bool firstBound=false;

///total number of items in the dataset: it is used for the size of the HeaderTable

tipoInt num_item;

long long cached_a_S;

///max ID of a item: it is specified in the in input file

tipoInt maxID;

///maximum lenght of one transaction: it is specified in the input specification file

tipoInt max_tr_length;

///number of transaction: it is specified in the input specification file

tipoInt num_tr;

///number of transactions

tipoInt N;

///number of transactions with class label = 1 in the class label file

tipoInt n1;

///number of permutations

tipoInt jp;

///target FWER

double alpha;

///number of significant patterns to find
int K_significant_patterns;

///permutations matrix

bool *permutations;
double *outcomes;
bool *labels;

///number of transactions with length > 0 in the dataset file

tipoInt effect_num_tr;

///temp transaction used in most methods

Transaction temp;

///set of items that are common to all the transaction: Clo(\emptyset); they are represented as ID, not as index of the IL.

string cloempty;

///maximum support of a item: it is used for the util nodes heuristic

tipoInt maximumSupport;

/**array for the util nodes heuristic and for the trivial heuristic:

- in the building: if j>minimum support--> utilNodeArray[j]=number of util nodes with count >= j

-in the computation: if j>minimum support--> utilNodeArray[j]=number of closed itemset inserted in queue with support >= j
 */

tipoInt* utilNodeArray;

/**array for the number of items in queue that have the suppMin as support: it is used to choose when to delete ppce from the queue
*/

tipoInt* itsOfSupp;

/**number of nodes of the PatriciaTrie
*/

tipoInt numNodes;

///CLASSES

/**the ItemList is used to remember the id of the items and their support
*/

class ItemList{
	public:
		tipoInt item;	///< id of the item
		tipoInt count;	///< count of the item
		/* SignificantMiner */
		tipoInt a_S;	///< a_S of the item
		tipoInt *a_S_permutations; // permuted as of the item
		/* SignificantMiner */
};

/// temporary ItemList used for quicksort

ItemList tempIL;

/// the HashTable is described by its size and by the first bucket of the table

class PatriciaNode;

class NodePointer;

class HashTable{
	public:
		///this is the number of bucket of the hash table

		tipoInt size;

		/**this is the first bucket (bucket[0]) of the hash table:
		*  bucket[i] will indicate the (i+1)-th bucket, so we allocate (malloc)
		* the right space to have ''size'' bucket
		*/

		///first bucket

		PatriciaNode* bucket;
};

///this is definition of the node structure for the PatriciaTrie

class PatriciaNode{
	public:
		tipoInt count;		///< count of the node
		tipoInt* lastItem;	///< last item memorized in the node

		/* SignificantMiner */
		union {
			uint8_t* a_S;	///< a_S of the node (compressed into 1 byte)
			tipoInt* a_S_; 	///< a_S of the node if count > 127
		};
		/* SignificantMiner */

		/**an union is used because the two following field are used in a MUTUALLY EXCLUSIVE manner (one in the building, the other one in the computation)
		*/

		/**in the building: pointer to the hash table used to access the children of the node
		 *in the computation:pointer to the last item in the intersection
		 */

		union {

			HashTable* hash;	///< pointer to the hash table used ONLY during the building of the PatriciaTrie
			tipoInt* lastInter; 	///< pointer to the last item in the intersection of the node: it is used in the computation of the top-k f.c.i.
		};

		/**in the building: number of chldren of the node
		 *in the computation: ID of the node (it is assigned after the building)
		 */

		union {

			tipoInt numChildren;	///<in the building: number of chldren of the node

			tipoInt ID;		///<in the computation: ID of the node
		};

		/**in the building: pointer to the node following the current one in the bucket of the hash table used by the father
		 *in the computation: pointer to the father
		 */

		union{

			PatriciaNode* next;	///<in the building: pointer to the node following the current one in the bucket of the hash table used by the father

			PatriciaNode* father;	///<in the computation: pointer to the father
		};


		/**in the building: it is >0 iff the intersection is empty (else =0)
		 *in the computation: timestamp used to indicate if the node is already visited
		 */

		union{

			tipoInt isEmpty;	///<in the building: it is >0 iff the intersection is empty (else =0)
			NodePointer* pntToNode;	///<in the computation: the pntToNode is the node pointer that is used to visit the tree to generate ppc-e
		};

		tipoInt item; 	///<the first item memorized in the node; the (i+1)-th item is memorized in item[i];

		void print(int a); ///used only for debug: a=0 if during construction, a=1 if during computation

		void printHash(); ///used only for debug
};

/**definition of the PatriciaTrie
*/

class PatriciaTrie{
	public:

		PatriciaNode* root;	///<pointer to the root;

		/**in the building phase: count of the j-th item (j-th is the id of the item)
		 * in the computation: id of the j-th item sorted by decreasing count
		 */
		tipoInt* count_element;
		/* SignificantMiner */
		tipoInt* count_a_S_element;
		/* SignificantMiner */

		/**list for the frequent item (and their id)
		 */

		ItemList* IL;

		PatriciaTrie();	///< constructor of the trie: num_item is the number of items in the dataset

		void makeIL();	///< create the ItemList from the num_item array

		void ILRaise(); ///< raise minimum support (global and current) from IL

		void IL_sort(tipoInt left, tipoInt right); ///< quicksort of the IL using "count" as key

		Transaction visitNode(PatriciaNode* p, tipoInt* ID, PatriciaNode** pntInHash);	///recursive method used by assignIDAndIntersect: it is used to assign ID in post-order

		void assignIDAndIntersect();	///assign the ID and compute the intersection for the nodes of the patriciaTrie; to be used only after all the transactions were inserted

		void freeNodeSpace(PatriciaNode* pNode);	///free the memory previous used by pNode

		void secondScan();
		void secondScan_new();

		/**USED ONLY FOR DEBUG PURPOSE
		*/

		void printIL();	///< prints the IL

		void printCountElement(); ///< print the count_element

		void printAfter();	///<prints node after ID were assigned

		void printCountDistribution();

		void printBefore();	///< prints node before ID were assigned

		void traverse();
};

/** hash function: used to individuate a node in the hash table; the key is the first item memorized in the node;
 */

tipoInt f_hash(tipoInt item, tipoInt elm_hash){
	return (item % elm_hash);
}

/** create a node from a transaction; it is used when an insertion causes a node creation
 */


PatriciaNode* createNode(Transaction t, bool toSum){
	#ifdef debug_patricia_creation
	cout << "creating new node... "  << endl;
	#endif
	///t[0]=number of items in the transaction
	PatriciaNode* pnode=(PatriciaNode*)malloc(sizeNode+t[0]*sizeof(tipoInt));
	if(!pnode){
		cout << "problems in createnode! " << endl;
		abort();
	}
	memset(pnode,0,sizeNode+t[0]*sizeof(tipoInt));
	//pnode->count=1;
	pnode->count=1;
	/* SignificantMiner */
	pnode->a_S=(uint8_t *)malloc((jp + 1) * sizeof(uint8_t));
	if(!pnode->a_S){
		cout << "problems in createnode with as! " << endl;
		abort();
	}
	/*for (int j=0;j<jp+1;j++){
		pnode->a_S[j]=0;
	}*/
	/* SignificantMiner */
	pnode->lastItem=&((&pnode->item)[t[0]-1]);
	pnode->hash=NULL;
	pnode->numChildren=0;
	pnode->next=NULL;
	if (!toSum){
		for (int i=0;i<t[0];i++){
			(&pnode->item)[i]=t[i+1];
		}
	}
	else{
		for (int i=0;i<t[0];i++){
			(&pnode->item)[i]=t[i+1];
			numNPI[t[i+1]]++;
		}
		totItem+=t[0];
	}
	pnode->isEmpty=1;	///this will always have empty intersection

	#ifdef debug_patricia_creation
	cout << "created new node! "  << endl;
	#endif
	return pnode;
}

/** sorting for a transaction
 */

void tran_sort(Transaction t, tipoInt left, tipoInt right) {
	tipoInt pivot;
	tipoInt p, l, r, temp;
	if (left >= right) return;
	p=(left+right)/2; /* pivot index */
	pivot = t[p]; /* pivot value */
	/* swap pivot with last element */
	t[p] = t[right];
	t[right] = pivot;
	/* Partition elements between left and right-1 */
	l=left; r=right-1;
	while (l<=r) {
		while ((l<=r) && (t[l] <= pivot)) l++;
		while ((r>=l) && (t[r] >= pivot)) r--;
		if (l<r) {
			temp = t[l];
			t[l] = t[r];
			t[r] = temp;
		}
	}
	temp = t[l];
	t[l] = t[right];
	t[right] = temp;

	/* Recurse */
	tran_sort(t, left, l-1);
	tran_sort(t, l+1,right);
}

///insert a node in a hash table specifyng the number of the bucket

void insert_in_hash(HashTable *fg, PatriciaNode* f, tipoInt pos){
	PatriciaNode* c;
	PatriciaNode* tmp;
	PatriciaNode** prev;	///pointer to the field in which the node must be inserted
	tipoInt elem;
	tipoInt new_elem;
	c=(&fg->bucket)[pos]; /* Bucket for insertion */
	new_elem=f->item;
	prev=&((&fg->bucket)[pos]);	///at the beginning, the node must be inserted directly in  the bucket
	/// Find correct position in bucket
	while (c!=NULL){
		elem=c->item;
		if(elem > new_elem) break;
		prev=&(c->next);
		c=c->next;
	}
	if(c==NULL){ ///insert at the end of bucket ///
		(*prev)=f;
		f->next=NULL;
	}
	else{
		tmp=c;
		(*prev)=f;
		f->next=tmp;
	}
}

///create a new hash and copy the element in the new hash

void new_hash(PatriciaNode *p, tipoInt num_buckets){
	HashTable* fg;
	PatriciaNode* tmp;
	PatriciaNode* nextNode;
	int i=0;
	/// Create new hash table
	fg=(HashTable*) malloc(sizeof(HashTable)+((num_buckets-1)*sizeof(PatriciaNode*)));
	///set all the fields of new hash to 0
	fg=(HashTable*) memset(fg,0,sizeof(HashTable)+((num_buckets-1)*sizeof(PatriciaNode*)));
	///set the new size of hash table
	fg->size=num_buckets;
	///insert the element in the new hash table
	for(i=0;i<num_buckets/2;i++){
		tmp=(&p->hash->bucket)[i];
		while (tmp!=NULL){
			nextNode=tmp->next;
			insert_in_hash(fg,tmp,f_hash((&tmp->item)[0],num_buckets));
			tmp=nextNode;
		}
	}
	/// Free old hash table and assign new one to node
	free(p->hash);
	p->hash=fg;
}



void convertNodeAs(PatriciaNode* pNode , bool copy_values){

	#ifdef debug_patricia_creation
	cout << "converting node with count " << pNode->count << endl;
	#endif

	// create the new array with larger representation
	tipoInt* a_S_larger = (tipoInt *)malloc((jp + 1) * sizeof(tipoInt));

	if(copy_values){
		uint8_t *temp_ptr = pNode->a_S;
		for(int j=0; j < jp+1; j++){
				a_S_larger[j] = temp_ptr[j];
		}
	}

	free(pNode->a_S);

	pNode->a_S_ = a_S_larger;

	#ifdef debug_patricia_creation
	cout << "done! " << endl;
	#endif

}

/** create two node from one node: it is used when a transaction is the prefix of another one
*/

PatriciaNode* split_node(PatriciaNode* node, tipoInt ind1){
	//if(DEEP_DEBUG)
	#ifdef debug_patricia_creation
	cout << "split node call , node = " << node << endl;
	#endif
	tipoInt i;
	PatriciaNode* c;
	HashTable* fg;

	/// Save pointer to old hash table for node
	fg=node->hash;

	/// Create transaction with last segment of items

	int numItem=(tipoInt) 1+(node->lastItem-&(node->item))-ind1; ///number of item that must be inserted n the new node
	for(i=0;i<numItem;i++){
		temp[i+1]=(&node->item)[i+ind1];
	}
	temp[0]=numItem;


						/*if(node->count==0){
							cout << "count = 0" << endl;
						}*/


	/// Resize node
	PatriciaNode* new_node;
	//node=(PatriciaNode*)realloc(node,sizeNode+(ind1-1)*sizeof(tipoInt));
	new_node=(PatriciaNode*)realloc(node,sizeNode+(ind1-1)*sizeof(tipoInt));
	if(new_node){
		node = new_node;


							/*if(node->count==0){
								cout << "count = 0" << endl;
							}*/


	}
	else{
		cout << " failed realloc with size " << sizeNode+(ind1-1) << endl;
		abort();
	}

	///update lastItem pnt
	node->lastItem=&((&node->item)[ind1-1]);
	if(temp[0]>0){
		/// Create new hash table of size 1 for node
		node->hash=(HashTable*)calloc(1,sizeof(HashTable));
		node->hash->size=1;
		node->hash->bucket=NULL;
		/// Insert transaction temp as (unique) child of node
		PatriciaNode* tmpPnt=createNode(temp,false);	///pointer to the new node
		node->hash->bucket=tmpPnt;			///the node has only one child
		tmpPnt->count=node->count;			///set the count of the new node

		/* SignificantMiner */
		if(node->count >= 128){
			for(int j=0; j < jp+1; j++){
				tmpPnt->a_S_[j]=node->a_S_[j];		///set the a_S of the new node
			}
		}
		else{
			for(int j=0; j < jp+1; j++){
				tmpPnt->a_S[j]=node->a_S[j];		///set the a_S of the new node
			}
		}
		/* SignificantMiner */

		tmpPnt->isEmpty=node->isEmpty;			///if the old node was empty, so is the new one
		if (tmpPnt->isEmpty==0){
			utilNodeArray[tmpPnt->count]--;
		}

		tmpPnt->numChildren=node->numChildren; 		///the node will have the chilwood of the splitted node
		node->numChildren=1;				///the node splitted has only one child
		tmpPnt->hash=fg;
	}
	//if(DEEP_DEBUG)
	#ifdef debug_patricia_creation
	cout << "split node done , node = " << node << endl;
	#endif
	return node;
}

/// split a node into 2 part; it is used when a transaction shares a prefix with node's items

PatriciaNode* split_2node(PatriciaNode* node, tipoInt ind1, Transaction tran, bool* permuted_transaction){
	#ifdef debug_patricia_creation
	cout << "split 2 node call "  << endl;
	#endif
	/* SignificantMiner */
	/*long long offset = ((long long)t_index) * (jp_long+1);*/
	/* SignificantMiner */
	//if(DEEP_DEBUG)
	//cout << "split 2 node call "  << endl;
	int i;
	PatriciaNode* c;
	PatriciaNode* tp;
	HashTable* fg;
	/// Save pointer to old hash table for node
	fg=node->hash;
	/// Create transaction with last segment of items
	tipoInt numItem=(tipoInt) 1+(node->lastItem-&(node->item))-ind1; ///number of item that must be inserted n the new node
	for(i=0;i<numItem;i++){
		temp[i+1]=0;
		temp[i+1]=(&node->item)[i+ind1];
	}
	temp[0]=i;
	/// Resize node

	//cout << "split_2node: call to realloc with size " << sizeNode+(ind1-1) << endl;
	PatriciaNode* new_node;
	//node=(PatriciaNode*)realloc(node,sizeNode+(ind1-1)*sizeof(tipoInt));
	new_node=(PatriciaNode*)realloc(node,sizeNode+(ind1-1)*sizeof(tipoInt));
	if(new_node){
		node = new_node;
	}
	else{
		cout << " split_2node: failed realloc with size " << sizeNode+(ind1-1) << endl;
		abort();
	}

	///update lastItem pnt
	node->lastItem=&((&node->item)[ind1-1]);
	tipoInt proof=(tipoInt) 1+(node->lastItem-&(node->item));
	/// Create new hash table of size 1 for node: this implies that hash_limit is ALWAYS > 1
	node->hash=(HashTable*)calloc(1,sizeof(HashTable));
	node->hash->size=1;

	///Insert transaction temp as child of node */
	PatriciaNode* tmpPnt=createNode(temp,false);	///pointer to the new node
	node->hash->bucket=tmpPnt;			///the node has only one child at the moment
	tmpPnt->count=node->count;
	/* SignificantMiner */
	if(node->count >= 128){
		convertNodeAs(tmpPnt , false);
		for(int j=0; j < jp+1; j++){
			tmpPnt->a_S_[j]=node->a_S_[j];		///set the a_S of the new node
		}
	}
	else{
		for(int j=0; j < jp+1; j++){
			tmpPnt->a_S[j]=node->a_S[j];		///set the a_S of the new node
		}
	}
	/* SignificantMiner */
	///if the node splitted was empty, so it is its children
	if (node->isEmpty>0){
		tmpPnt->isEmpty=1;
		node->isEmpty=0;
	}
	else{
		tmpPnt->isEmpty=0;
	}
	tmpPnt->numChildren=node->numChildren;///the node will have the chilwood of the splitted node
	node->numChildren=2;	///the father has 2 children
	tmpPnt->hash=fg;
	/// Insert transaction tran as child of node */
	if ((tmpPnt->item)>tran[1]) {
		tp=tmpPnt;
		tmpPnt=createNode(tran,true);
		/* SignificantMiner */
		for(int j=0; j < jp+1; j++){
			//index = offset + j;
			/*if(permutations[index])
				tmpPnt->a_S[j] = 1;*/
			//tmpPnt->a_S[j] = permutations[index];
			tmpPnt->a_S[j] = permuted_transaction[j];
		}
		/* SignificantMiner */
		tmpPnt->next=tp;
		node->hash->bucket=tmpPnt;
	}
	else{
		tp=createNode(tran,true);
		/* SignificantMiner */
		for(int j=0; j < jp+1; j++){
			//index = offset + j;
			/*if(permutations[index])
				tp->a_S[j] = 1;*/
			//tp->a_S[j] = permutations[index];
			tp->a_S[j] = permuted_transaction[j];
		}
		/* SignificantMiner */
		tmpPnt->next=tp;
	}
	#ifdef debug_patricia_creation
	cout << "split 2 node done "  << endl;
	#endif
	return node;
}

int readSpecFile(int argc, char *argv[]){
	string line;
	if (argc<4){
		cout<< "\nUsage: ./topkminer file.spec k_current k_global [jp]\n" << endl;
		return 0;
	}
	ifstream specfile (argv[1]);
	k=10000000;
	k_max=10000000;
	max_ram=atoi(argv[2]);
	K_significant_patterns=atoi(argv[3]);
	jp=0;
	alpha=0.0;
	if(argc>=5)
		jp=atoi(argv[4]);
	if(argc>=6){
		alpha=atof(argv[5]);
		//cout << "alpha = " << alpha << endl;
	}
	if ( k_max < k ){
		cout << "k_global IS GREATER THAN k_current!\nSETTING k_current EQUAL TO k_global\n" << endl;
		k = k_max;
	}
	if (!specfile.is_open()){
		cout << "Unable to open specification file\n";
		return 0;
	}
	int i=0;
	while (! specfile.eof() && i<7){
		i++;
		getline (specfile,line);
		if(i==1){
			fileinput=line;
			dataset.open(fileinput.c_str(),ifstream::in); ///<dataset file;
			if(!dataset.is_open()){
				cout << "Unable to open dataset file\n";
				return 0;
			}
		}
		if(i==2){
			///maximum ID of an item in the dataset
			maxID=atoi(line.c_str());
		}
		if(i==3){
			///maximum length of a transaction
			max_tr_length=atoi(line.c_str());
		}
		if(i==4){
			///number of transactions in the dataset
			num_tr=atoi(line.c_str());
		}
		/* SignificantMiner */
		if(i==5){ // class label file
			c_fileinput=line;
			dataset.close();
			dataset.open(c_fileinput.c_str(),ifstream::in); ///<class labels file;
			if(!dataset.is_open()){
				cout << "Unable to open class label file\n";
				return 0;
			}
			dataset.close();
			dataset.open(fileinput.c_str(),ifstream::in); ///<reopen dataset file;
		}
		if(i==6){ // target FWER
			if(alpha < 0.000001)
			alpha=atof(line.c_str());
		}
		if(i==7){ // number of permutations
			if(jp==0)
			jp=atoi(line.c_str());
		}
		/* SignificantMiner */
	}
	if(i!=7){
		cout << "\nThere is some problem in the specification file!\n" << endl;
		return 0;
	}
	specfile.close();

	cout << "dataset = " << fileinput << endl;
	cout << "alpha = " << alpha << endl;
	cout << "jp = " << jp << endl;
	return 1;
}

void firstScan_new(PatriciaTrie* p){

	if(DEEP_DEBUG)
	cout << "call first scan " << endl;

	string line, line1;
	char* tmp;
	char* c_tmp;
	string end=""; /// to verify that a string is not null;
	const char *space = " "; ///delimiter
	char* charline; /// to memorize tokens

	/* SignificantMiner */
	/*long long jp_long = (long long)jp;
	long long count_a_s_size = ((long long)maxID+1)*(jp_long+1);
	count_a_S_element_space = sizeof(tipoInt) * ((double)count_a_s_size / 1000000.0);
	p->count_a_S_element=(tipoInt*)calloc(count_a_s_size,sizeof(tipoInt));
	if(!p->count_a_S_element){
		cout << "ERROR: can't allocate p->count_a_S_element with count_a_s_size = " << count_a_s_size << endl;
	}
	long long index = 0;*/
	p->count_a_S_element=(tipoInt*)calloc(maxID+1,sizeof(tipoInt));
	/* SignificantMiner */

	///initialize the count_element array
	p->count_element=(tipoInt*)calloc(maxID+1,sizeof(tipoInt));
	int num=0;
	int empty_lines=0;
	effect_num_tr=0;
	max_trans_length = 0;
	int trans_length = 0;
	long long i = 0;
	while(!dataset.eof()){
		getline (dataset,line);
		trans_length = 0;
		num++;
		//cout << line << endl;
		charline = ( char * )line.c_str();
		///control if the string is empty
		if (end.compare(charline)!=0){
			effect_num_tr++;
			trans_length++;
			///the first token
			tmp=strtok(charline,space);
			while(tmp!=NULL && !isspace(*tmp)){

				if(trans_length > max_trans_length)
					max_trans_length = trans_length;

				/*cout << "item = " << tmp << endl;
				cout << "atoi(tmp) = " << atoi(tmp) << endl;
				cout << "a_S = " << permutations[i * (jp + 1)] << endl;*/
				///increment the count of the item with ID tmp

				/*p->count_element[atoi(tmp)]++;
				int offset = atoi(tmp) * (jp+1);
				for(int j = 0; j < jp+1; j++){
					index = i * (jp + 1) + j;
					if(permutations[index]){
						p->count_a_S_element[offset + j]++;
					}
				}*/

				int element = atoi(tmp);
				p->count_element[element]++;
				/* SignificantMiner */
				/*long long offset1 = ((long long)element) * (jp_long+1);
				long long offset2 = i * (jp_long + 1);
				for(long long j = 0; j < jp_long+1; j++){
					index = offset2 + j;
					p->count_a_S_element[offset1 + j]+=permutations[index];
				}*/
				p->count_a_S_element[element]+=labels[i];
				/* SignificantMiner */

				tmp=strtok(NULL,space);
				trans_length++;
			}

		}
		else{
			empty_lines++;
		}
		i++;
	}

	if(N != effect_num_tr){
		cout << "Number of transactions in label and dataset files are different! " << endl;
		cout << "Number of labels = " << N << endl;
		cout << "Number of transactions = " << effect_num_tr << endl;
		cout << "Number of empty lines found in dataset = " << empty_lines << endl;
		abort();
	}

	/*for(int aj=0; aj < maxID+1; aj++)
		cout << aj << " - " << p->count_element[aj] << endl;*/

	if(DEEP_DEBUG)
	cout << "first scan ok " << endl;

	cout << "max_trans_length = " << max_trans_length << endl;
}

void firstScan(PatriciaTrie* p , bool* permutations){

	if(DEEP_DEBUG)
	cout << "call first scan " << endl;

	string line, line1;
	char* tmp;
	char* c_tmp;
	string end=""; /// to verify that a string is not null;
	const char *space = " "; ///delimiter
	char* charline; /// to memorize tokens

	/* SignificantMiner */
	p->count_a_S_element=(tipoInt*)calloc((maxID+1)*(jp+1),sizeof(tipoInt));
	tipoInt index = 0;
	/* SignificantMiner */

	///initialize the count_element array
	p->count_element=(tipoInt*)calloc(maxID+1,sizeof(tipoInt));
	int num=0;
	effect_num_tr=0;
	int i = 0;
	while(!dataset.eof()){
		getline (dataset,line);
		num++;
		//cout << line << endl;
		charline = ( char * )line.c_str();
		///control if the string is empty
		if (end.compare(charline)!=0){
			///the first token is the transaction's length
			tmp=strtok(charline,space);
			///to consider if there are some void transaction
			if (atoi(tmp)!=0) effect_num_tr++;
			///from the second token items begin
			tmp=strtok(NULL,space);

			while(tmp!=NULL && !isspace(*tmp)){
				cout << "item = " << tmp << endl;
				cout << "atoi(tmp) = " << atoi(tmp) << endl;
				cout << "a_S = " << permutations[i * (jp + 1)] << endl;
				///increment the count of the item with ID tmp
				p->count_element[atoi(tmp)]++;
				cout << "p->count_element[atoi(tmp)] = " << p->count_element[atoi(tmp)] << endl;
				/* SignificantMiner */
				int offset = atoi(tmp) * (jp+1);
				for(int j = 0; j < jp+1; j++){
					index = i * (jp + 1) + j;
					if(permutations[index]){
						p->count_a_S_element[offset + j]++;
					}
				}
				/* SignificantMiner */
				tmp=strtok(NULL,space);
			}

		}
		i++;
	}

	for(int aj=0; aj < maxID+1; aj++)
		cout << aj << " - " << p->count_element[aj] << endl;

	if(DEEP_DEBUG)
	cout << "first scan ok " << endl;
}

///function that converts an integer into a string: it is used to memorize the Clo(\emptyset) as a string

string itos(tipoInt i)	// convert int to string
{
	stringstream s;
	s << i;
	return s.str();
}

///prints the array for util node heuristic

void printUtilNodeArray(){
	cout << "\n************ utilNodeArray:BEGIN ****************\n" << endl;
	int showed=0;
	cout << maximumSupport << " | " << utilNodeArray[maximumSupport] << endl;
	for (int i=maximumSupport-1; i>= suppMin && showed<10; i--){
		if (utilNodeArray[i]!=utilNodeArray[i+1]){
			cout << i << " | " << utilNodeArray[i] << endl;
			showed++;
		}
	}
	cout << "\n************ utilNodeArray:END ****************\n" << endl;
}

///a NodeIDElem is a element of the list of pointer used to remember what nodes of PatriciaTrie contains the least element of a pattern P: used only in findCloEmpty

class NodeIDElem{
	public:
		tipoInt nodeID;			///id of the node pointed
		NodeIDElem* next;			///next NodeIDElem in the list
};

///a nodeIDList is an array with the IDs of the nodes coinvolted; in [0] there is the length

typedef tipoInt* NodeIDList;

///a StartNode is a pointer to a Node: it is used at the end to derive the final lists of ppc-e

class StartNode{
	public:
		tipoInt ID;
		short label;		///label at the end is 0 if the node with ID is not an extremal in the subtree, is 1 if the node is an extremal of a subtree that has 2 DISTINCT  extremal and is 2 if it is the only node in the subtree
};

/// a StartNodeList is an array of StartNode

typedef StartNode* StartNodeList;

/// a VisitedNodeElem is a kind of pointer that is used to memorize the information associated with a visited node for an item

class VisitedNodeElem{
	public:
		tipoInt offPntNode;	///offset of the pnToNode of the node in manualmem
		tipoInt IDnode;		///ID of the node visited
		tipoInt numInInter;	///number of items in the intersectoin that must be considered for the memorization of the pointer: this is the length of the prefix of the intersection that must be memorized
};

///a VisitedNodeList is an array of VisitedNodeElem

typedef VisitedNodeElem* VisitedNodeList;

///a NodePointer is a pointer to a PatriciaNode used in the computation of the ppc-e of an itemset

class NodePointer{
	public:
		tipoInt node; 	///ID of the node associated to the pointer
		tipoInt count;		///counter associated with the node
		/* SignificantMiner */
		tipoInt count_a_S;		///counter of a_S associated with the node
		//~NodePointer();
		tipoInt nl_list_index;
		// if you add something here, remind to increase size sizePnt, which is not dynamic
		// otherwise you will meet blasfemy and pain
		/* SignificantMiner */
		tipoInt minIDNodePnt;	///index of the node with minimum ID in the subtree rooted in the node: the index is referred to the startList associated
		tipoInt maxIDNodePnt;	///index of the node with maximum ID in the subtree rooted in the node: the index is referred to the startList associated
		tipoInt nextNodePointer; 	/// index in manualmem of next poiter of the current list: it is =-1 if next doesn't exist
		tipoInt intersection;	///intersection associated with the pointer
};

/*NodePointer::~NodePointer(){
	cout << "nodepointer deleted!"<<endl;
    free(a_S);
}*/

typedef NodePointer* NodePointerList;

/** header table: it is used in computation to memorize the count, the intersection and the node list of a ppc-e
*/

class HeaderTable{
	public:
		tipoInt count;			///count of the extension
		/* SignificantMiner */
		tipoInt a_S;			///a_s of the extension
		/* SignificantMiner*/
		Transaction intersection;	///intersection on the prefix
		VisitedNodeList visitedList;		///list of the nodes of the PatriciaTrie that contains the item
		tipoInt pointerList;	///list of pointer to the nodes touched but not visited: we identify a node with is offset in manualmem
		tipoInt* lastNodePnt;	///pointer to the field in which the next NodePointer must be added
		void print();			///print the information stored in the HT;
};

/**header table that will used in all the computation
*/

HeaderTable* HT;

/**array for the corrispondence between ID and address of a PatriciaNode: addressOf[ID]= address of PatriciaNode with ID
 */

PatriciaNode** addressOf;

/** array of list for the NodeIDList of the ppc-e of Clo(emptyset): they are treated in a different way because their construction is different from other ppc-e extension
*/

NodeIDElem** CloEmptyNodeIDList;

/**array of pointers to the right field in which the next NodeIDElem of a list in CloEmptyNodeIDList must be inserted
*/

NodeIDElem*** lastCloEmpty;

void resizeManualMem(){
	cout << "manual mem resize"<< endl;
	tipoInt* prev=manualmem;
	numResize++;
	tipoInt nextfreeOff=nextfree-manualmem;		///offset of nextfree in manualmem
	tipoInt sizeBefore=lastvalid-manualmem+1;	///size of manual memory allocated before

	cout << "current size "<< sizeBefore << endl;

	sizeBefore=2*sizeBefore;			///new size of manualmem

	cout << "new size "<< sizeBefore << endl;

	manualmem=(tipoInt*)realloc(manualmem, sizeBefore*sizeof(tipoInt));
	if(!manualmem){
		cout << "realloc failed!" << endl;
		abort();
	}
	nextfree=manualmem+nextfreeOff;	///in the case that realloc will move the block in memory
	lastvalid=&(manualmem[sizeBefore-1]);	///new last valid address for the manualmem
	cout << "manual mem done"<< endl;
}

#endif
