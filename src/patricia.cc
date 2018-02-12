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

#include "patricia.h"
#include "queueheap.cc"
#include "wy_.cc"

 //#define debug_patricia_creation 1

/**raise minimum support using the trivial heuristic AND inserting a closed itemset in queue of support supp; used in computation
NB: utilNodeArray is updated with the insertion of "supp" value in the execution of the method
 */

void trivialHeur(tipoInt supp){
	return;
	///trivial heuristic on suppMinCurr and suppMin
	for (tipoInt i=supp; i>suppMin; i--){
		if (utilNodeArray[i]<(k_max-1)){
			utilNodeArray[i]++;
			if (utilNodeArray[i]==k&&i>suppMinCurr){
				suppMinCurr=i;
			}
		}
		else{
			if (utilNodeArray[i]==(k_max-1)){
				suppMin=i;
				if (k==k_max){
					suppMinCurr=suppMin;
				}
			}
			break;
		}
	}
}

/**raise current minimum support using the trivial heuristic AND inserting a closed itemset in queue of support supp; used in computation
NB: utilNodeArray is update with the insertion of "supp" value in the execution of the method
 */

void trivialHeurCurr(tipoInt supp){
	return;
	///trivial heuristic on suppMinCurr
	for (tipoInt i=suppMin; i<=supp; i++){
		utilNodeArray[i]++;
		if (utilNodeArray[i]>=k&&i>suppMinCurr){
			suppMinCurr=i;
		}
	}
}


/**prints the informations stored in PatriciaNode: they are different if we are during the construction (a=0) or during the computation (a=1)
 */

void PatriciaNode::print(int a){
	cout << "\n---------PatriciaNode-----------\n" << endl;
	cout << "address= " << this << endl;
	cout << "count= " << count << endl;
	/* SignificantMiner */
	cout << "a_S= " << a_S << endl;
	/* SignificantMiner */
	cout << "address of pointer to lastItem= " << &lastItem << endl;
	cout << "lastItem address= " << lastItem <<endl;
	if (lastItem!=NULL)
	cout << "lastItem index= " << *(lastItem) << endl;
	if (a==0){
		cout << "hash address= " << hash << endl;
		cout << "numChildren= " << numChildren << endl;
		cout << "next address= " << next << endl;
		cout << "isEmpty= " <<  isEmpty << endl;
	}
	else{
		cout << "lastInter address= " << lastInter << endl;
		cout << "lastInter index= " << *(lastInter) << endl;
		cout << "ID= " << ID <<endl;
		cout << "father address= " << father << endl;
		cout << "nodePointer address= " <<  pntToNode << endl;
	}
	///calculate the number of item in the node
	tipoInt numItem=(tipoInt) 1+(lastItem-&(item));
	cout << "number of items= " << numItem << endl;
	cout << "items= ";
	for (int i=0; i<numItem; i++){
		cout << (&item)[i] << " ";
	}
	cout<<endl;
	if (a==1){///print the intersection
		numItem=(tipoInt) (lastInter-lastItem);
		cout << "number of items in intersection= " << numItem << endl;
		if (numItem>0)
			cout << "items= ";
		for (int i=0; i<numItem; i++){
			cout << (lastItem)[i+1] << " ";
		}
	}
	cout << "\n ------------------------------\n" << endl;
}

/**prints the HashTable that a patricia node has during the building
*/

void PatriciaNode::printHash(){
	cout<< "\n*************** Printing the Hash Table: BEGIN ************\n" << endl;
	HashTable* hashT=hash;
	tipoInt numBucket=hashT->size;
	PatriciaNode* tmp;
	for (int i=0;i<numBucket; i++){
		cout << "\n--------------------\n bucket= " << i << endl;
		tmp=(&hashT->bucket)[i];
		while(tmp!=NULL){
			tmp->print(0);
			tmp=tmp->next;
		}
	}
	cout<< "\n*************** Printing the Hash Table: END ************\n" << endl;
}


/** constructor of the PatriciaTrie; num_item is the number of items in the dataset
*/

PatriciaTrie::PatriciaTrie(){
	///create the root
	root=(PatriciaNode*)calloc(1,sizeNode);
	///assign the standard values for the root
	root->count = 0;
	/* SignificantMiner */
	root->a_S = (uint8_t*)malloc((jp+1)*sizeof(uint8_t));
	for(int aj=0; aj<jp+1; aj++)
		root->a_S[aj]=0;
	/* SignificantMiner */
	root->lastItem= NULL;
	///create a hash table initially NULL
	root->hash=NULL;
	root->numChildren=0;
	root->next=NULL;
	root->isEmpty=1;	///for us the root has always empty intersection
	///create a util node list with one element having count=1 and number=0
	suppMin=1;
	suppMinCurr=1;
	///initialize the number of items in the PatriciaTrie
	totItem=0;
}

/**make the ItemList: to use after evaluating num_item e count_element
*/

void PatriciaTrie::makeIL(){

	long long maxID_long = (long long)maxID;

	if(DEEP_DEBUG){
	cout << "count_a_S_element: " << endl;
	for(long long i=0; i < maxID_long+1; i++){
			cout << count_a_S_element[i] << " | ";
		cout << endl;
	}
	}

	IL=(ItemList*)calloc((maxID+1),sizeof(ItemList));
	for (long long j=0;j<=maxID_long;j++){
		IL[j].count=count_element[j];
		/*IL[j].a_S=(tipoInt*)calloc(jp+1,sizeof(tipoInt));
		long long offset = j * (jp_long + 1);
		for(long long aj=0; aj<jp_long+1; aj++)
			IL[j].a_S[aj]=count_a_S_element[offset+aj];*/
		IL[j].a_S=count_a_S_element[j];
		IL[j].a_S_permutations=NULL;
		IL[j].item=j;

		if(DEEP_DEBUG){
		cout << "makeIL: item: " << j << " has now support = " << IL[j].count << endl;
		cout << "makeIL: item: " << j << " has now a_S = " << IL[j].a_S << endl;}
	}
}

/**raise minimum support (global and current) using IL and refresh count_element so that it indicates the index of frequent items in the IL (it is -1 if an item is not frequent)
*/

void PatriciaTrie::ILRaise(){

	long long maxID_long = (long long)maxID;

	/*int numDiffSupp=1; ///<counter of different value of support for the singletons
	bool kOld=true;
	bool k_maxOld=true;
	if (k_max==1) {
		suppMin=IL[0].count;
		k_maxOld=false;
	}
	if (k==1) {
		suppMinCurr=IL[0].count;
		kOld=false;
	}
	for(int i=1; i<=maxID && (kOld || k_maxOld); i++){
		if (IL[i].count<IL[i-1].count) numDiffSupp++;
		if (numDiffSupp>=k_max && k_maxOld && IL[i].count> suppMin) {
			suppMin=IL[i].count;
			k_maxOld=false;
		}
		if (numDiffSupp>=k && kOld && IL[i].count> suppMinCurr) {
			suppMinCurr=IL[i].count;
			kOld=false;
		}
	}*/

	suppMin = 1;
	while(getPsi(suppMin) > alpha)
		suppMin++;
	suppMinCurr = suppMin;

	cout << "suppMin initialized to " << suppMin << endl;

	for(long long i=0; i<=maxID_long; i++){
		if (IL[i].count>=suppMin){
			count_element[IL[i].item]=i;
			/*long long offset = ((long long)IL[i].item) * (jp_long+1);
			for(long long aj=0; aj<jp_long+1; aj++)
				count_a_S_element[offset+aj]=i;*/
		}
		else{
			count_element[IL[i].item]=-1;
			//cout << "setting infrequent the item " << IL[i].item << endl;
			/*long long offset = ((long long)IL[i].item) * (jp_long+1);
			for(long long aj=0; aj<jp_long+1; aj++)
				count_a_S_element[offset+aj]=-1;*/
		}
	}
}

/**quicksort for the ItemList
*/

void PatriciaTrie::IL_sort(tipoInt left, tipoInt right) {
	if (left>=right) return;
	tipoInt pivot;
	tipoInt p, l, r;
	tipoInt x;

	p=right-left;
	x=p*(double)((double)rand()/((double)RAND_MAX));
	p=left+x;

	pivot = IL[p].count; /* pivot value */
	tipoInt numPiv=1; //number of elements equal to the pivot
	/* swap pivot with last element */
	tempIL.item = IL[p].item;
	tempIL.count = IL[p].count;
	/* SignificantMiner */
	tempIL.a_S = IL[p].a_S;
	/* SignificantMiner */
	IL[p].item = IL[right].item;
	IL[p].count = IL[right].count;
	/* SignificantMiner */
	IL[p].a_S = IL[right].a_S;
	/* SignificantMiner */
	IL[right].item=tempIL.item;
	IL[right].count=tempIL.count;
	/* SignificantMiner */
	IL[right].a_S=tempIL.a_S;
	/* SignificantMiner */
	tipoInt re=right-1; //re=rightmost element < pivot
	//set the right re
	while (re>=left && IL[re].count==pivot){
		numPiv++;
		re--;
	}
	if (re<left) return;	//all the elements are equal

	/* Partition elements between left and right-1 */
	l=left; r=re;
	while (l<=r) {
		while ((l<=r) && (IL[r].count <= pivot)){
			//if the element is equal to the pivot, swap immediately with IL[re]
			if(IL[r].count==pivot){
				tempIL.item = IL[r].item;
				tempIL.count = IL[r].count;
				/* SignificantMiner */
				tempIL.a_S = IL[r].a_S;
				IL[r].item = IL[re].item;
				IL[r].count = IL[re].count;
				/* SignificantMiner */
				IL[r].a_S = IL[re].a_S;
				IL[re].item=tempIL.item;
				IL[re].count=tempIL.count;
				/* SignificantMiner */
				IL[re].a_S=tempIL.a_S;
				re--;
				numPiv++;
			}
			r--;
		}
		while ((l<=r) && (IL[l].count > pivot)) l++;
		if (l<r) {
			//if IL[l] is equal to pivot...
			if (IL[l].count==pivot){
				// l<- r,re<-l,r<-re
				tempIL.count=IL[l].count;
				/* SignificantMiner */
				tempIL.a_S=IL[l].a_S;
				tempIL.item=IL[l].item;
				IL[l].count=IL[r].count;
				/* SignificantMiner */
				IL[l].a_S=IL[r].a_S;
				IL[l].item=IL[r].item;
				IL[r].count=IL[re].count;
				/* SignificantMiner */
				IL[r].a_S=IL[re].a_S;
				IL[r].item=IL[re].item;
				IL[re].count=tempIL.count;
				/* SignificantMiner */
				IL[re].a_S=tempIL.a_S;
				IL[re].item=tempIL.item;
				r--;
				re--;
				numPiv++;
			}
			else{
				tempIL.item = IL[l].item;
				tempIL.count = IL[l].count;
				/* SignificantMiner */
				tempIL.a_S = IL[l].a_S;
				IL[l].item = IL[r].item;
				IL[l].count = IL[r].count;
				/* SignificantMiner */
				IL[l].a_S = IL[r].a_S;
				IL[r].item=tempIL.item;
				IL[r].count=tempIL.count;
				/* SignificantMiner */
				IL[r].a_S=tempIL.a_S;
			}
		}
	}

	//copy the element equal to pivot in the middle, stopping when swap indices meet
	for (int j=0; j<numPiv && j<=((right-l)/2); j++){
		tempIL.count=IL[l+j].count;
		/* SignificantMiner */
		tempIL.a_S=IL[l+j].a_S;
		tempIL.item=IL[l+j].item;
		IL[l+j].count=IL[right-j].count;
		/* SignificantMiner */
		IL[l+j].a_S=IL[right-j].a_S;
		IL[l+j].item=IL[right-j].item;
		IL[right-j].count=tempIL.count;
		/* SignificantMiner */
		IL[right-j].a_S=tempIL.a_S;
		IL[right-j].item=tempIL.item;
	}
		/* Recurse */
	IL_sort(left,l-1);
	IL_sort(l+numPiv,right);
}

/**print the ItemList
*/

void PatriciaTrie::printIL(){
	cout<< "\n*************** Printing the ItemList: BEGIN ************\n" << endl;
	for (int i=0; i<maxID; i++){
		cout << "i: " << i << " ";
		cout << "Item ID= " << IL[i].item << " | count = " << IL[i].count <<endl;
	}
	cout<< "\n*************** Printing the ItemList: END ************\n" << endl;
}

/**insert transaction t in the PatriciaTrie at the node pNode;
*/

void insert(Transaction t, PatriciaNode* pNode, bool* permuted_transaction){ // SignificantMiner
//void insert(Transaction t, PatriciaNode* pNode){

	/*if(DEEP_DEBUG){
	cout << "call to insert transaction ";
	for(int ik = 1; ik <= t[0]; ik++)
		cout << t[ik] << " ";
	cout << " with label " << permutations[t_index * (jp+1)] << endl;
	}*/

	#ifdef debug_patricia_creation
	cout << "inserting t... "  << endl;
	#endif

	Transaction tnew; 	///used in recursive call of insert
	tipoInt num_bucket;	///number of bucket in the hash table of pNode; used for f_hash
	PatriciaNode* child; 	///child where to look for the insertion
	PatriciaNode** prev;	///pointer to the field in which the new node will be inserted
	tipoInt elem=0; 	///first item in a nodes that is in a bucket of the hash table
	/* SignificantMiner */
	//long long offset = ((long long)t_index) * (jp_long+1);
	/* SignificantMiner */

	if (t[0]==0) return;	/// base case: insertion of a null transaction

	///if the node has no children, it hasn't a hash table
	if (pNode->numChildren==0){
		pNode->hash=(HashTable*)calloc(1,sizeof(HashTable));
		pNode->hash->size=1;
		pNode->hash->bucket=NULL;
		num_bucket=1;
	}
	else{
		num_bucket=pNode->hash->size;
	}

	///find the bucket where to search for the correct child
	///THIS DEFINITION IS DIFFERENT FROM PATRICIA MINE!
	child=(&pNode->hash->bucket)[f_hash(t[1],num_bucket)];
	prev=&((&pNode->hash->bucket)[f_hash(t[1],num_bucket)]);
	///now look for the correct node in the patricia
	while(child!=NULL){
		elem=(&child->item)[0];	///first item in the node
		if (elem>=t[1]){	///HP:elements in a bucket of the hash table are ordered by decreasing indexes
			break;
		}
		else{	///consider the next node in the bucket
		prev=&(child->next);
		child=child->next;
		}
	}
	if(child==NULL){	///there are two possibilities: the new node will be the first in the bucket or the last (if it's not so, child!=NULL)

		///create a new node and insert him as child
		pNode->numChildren+=1;
		child=createNode(t,true);
		//cout << "   point 1"  << endl;
		/* SignificantMiner */
		uint8_t *temp_ptr = child->a_S;
		for(int j=0; j < jp+1; j++){
			//index = offset + j;
			/*if(permutations[index])
				child->a_S[j] = 1;*/
				temp_ptr[j] = permuted_transaction[j];
		}
		//cout << "   point 2"  << endl;
		/* SignificantMiner */
		(*prev)=child;
		///update the number of nodes of the PatriciaTrie
		numNodes++;
		//if(DEEP_DEBUG)
		//cout << "created node with count " << child->count << " has a_S " << child->a_S[0] << endl;
	}
	else {
		///the node is not the first in the bucket
		if (elem==t[1]){	///the node share at least the prefix with the transaction
			int i=1;	///i=number of item shared by transaction and node's items
			int numItem=(tipoInt) 1+(child->lastItem-&(child->item));
			while (i<numItem){

				/// Case 1: t is prefix of items of child
				if(i+2>t[0]&&(i!=(t[0]-1))) {

					/*if(child->count==0){
						cout << "count = 0" << endl;
					}*/

					#ifdef debug_patricia_creation
					cout << "before split node , child = " << child << endl;
					#endif
					child = split_node(child,i);///after split_node, child points to the father of the new node
					// ERROR HERE, CHILD IS NOT CORRECTLY UPDATED
					#ifdef debug_patricia_creation
					cout << "after split node , child = " << child << endl;
					#endif

					// update adress of child if changed
					(*prev)=child;

					child->count++;
					//cout << "   point 3"  << endl;
					if(child->count >= 128){
						if(child->count == 128) convertNodeAs(child , true);
						tipoInt *temp_ptr = child->a_S_;
						for(int j=0; j < jp+1; j++){
							//index = offset + j;
							//child->a_S_[j]+=permutations[index];
							temp_ptr[j]+=permuted_transaction[j];
						}
					}
					else{
						uint8_t *temp_ptr = child->a_S;
						for(int j=0; j < jp+1; j++){
							//index = offset + j;
							//child->a_S[j]+=permutations[index];
							temp_ptr[j]+=permuted_transaction[j];
						}
					}
					//cout << "   point 4"  << endl;
					/* SignificantMiner */

					///for the util node heuristic
					/*if (child->count>suppMin){
						trivialHeur(child->count);
					}*/

					///update the number of nodes of the PatriciaTrie
					numNodes++;
					i=t[0]+1;
					break;
				}

				/// Case 2: t and items of child share a prefix
				if( (&child->item)[i]!= t[i+1]){
					t[i]=t[0]-i;
					tnew=&t[i];
					child = split_2node(child,i,tnew,permuted_transaction); ///after split_2node, child points to the father of the new node
					child->count++;
					// update adress of child if changed
					(*prev)=child;
					//cout << "   point 5"  << endl;
					if(child->count >= 128){
						if(child->count == 128) convertNodeAs(child , true);
						tipoInt *temp_ptr = child->a_S_;
						for(int j=0; j < jp+1; j++){
							//index = offset + j;
							//child->a_S_[j]+=permutations[index];
							temp_ptr[j]+=permuted_transaction[j];
						}
					}
					else{
						uint8_t *temp_ptr = child->a_S;
						for(int j=0; j < jp+1; j++){
								//index = offset + j;
								//child->a_S[j]+=permutations[index];
								temp_ptr[j]+=permuted_transaction[j];
						}
					}
					//cout << "   point 6"  << endl;
					/* SignificantMiner */
					i=t[0]+1;
					///update the number of nodes of the PatriciaTrie
					numNodes=numNodes+2;
					break;
				}
				i++;
			}

			/// Case 3: items of child are prefix of t
			if(i<=t[0]){
				child->count++;
				//cout << "   point 7"  << endl;
				if(child->count >= 128){
					if(child->count == 128) convertNodeAs(child , true);
					tipoInt *temp_ptr = child->a_S_;
					for(int j=0; j < jp+1; j++){
						//index = offset + j;
						//child->a_S_[j]+=permutations[index];
						temp_ptr[j]+=permuted_transaction[j];
					}
				}
				else{
					uint8_t *temp_ptr = child->a_S;
					for(int j=0; j < jp+1; j++){
						//index = offset + j;
							//child->a_S[j]+=permutations[index];
							temp_ptr[j]+=permuted_transaction[j];
					}
				}
				//cout << "   point 8"  << endl;

				///for the util node heuristic
				if ((child->isEmpty>0) && child->count>suppMin){
					utilNodeArray[child->count-1]--;
					utilNodeArray[child->count]++;
					if (utilNodeArray[child->count]==k_max){
						suppMin=child->count;
					}
					if (utilNodeArray[child->count]==k){
						suppMinCurr=child->count;
					}
				}

				if(i<t[0]){
					t[i]=t[0]-i;
					tnew=&t[i];
					//insert(tnew,child); /** Recurse */
					/* SignificantMiner */
					//if(DEEP_DEBUG)
					#ifdef debug_patricia_creation
					cout << "recursion! " << endl;
					#endif
					insert(tnew,child,permuted_transaction); /** Recurse */
				}
			}
		}
		else{
			/// if elem>t[1] create new child
			pNode->numChildren++;
			PatriciaNode* tmpPoint=child;
			child=createNode(t,true);
			//cout << "   point 9"  << endl;
			uint8_t *temp_ptr = child->a_S;
			for(int j=0; j < jp+1; j++){
				//index = offset + j;
				//child->a_S[j] = permutations[index];
				temp_ptr[j] = permuted_transaction[j];
			}
			//cout << "   point 10"  << endl;
			/* SignificantMiner */
			(*prev)=child;
			child->next=tmpPoint;
			///update the number of nodes of the PatriciaTrie
			numNodes++;
			//if(DEEP_DEBUG)
			//cout << "created node with count " << child->count << " has a_S " << child->a_S << endl;
		}
	}
	/// If hash table is full do rehash
	if ((pNode->numChildren/num_bucket) > hash_limit){
		new_hash(pNode,(2*num_bucket));
	}

	#ifdef debug_patricia_creation
	cout << "inserted without problems "  << endl;
	#endif
}

///print the count_element array

void PatriciaTrie::printCountElement(){
	cout<< "\n*************** Printing the count_element array: BEGIN ************\n" << endl;
	cout << "\n num_item= " << maxID+1 << "\n" <<endl;
	cout << "\n count_element \n --------------------" << endl;
	for(int j=0; j<=maxID; j++){
		cout << "| "<< count_element[j] << "| ";
	}
	cout << "\n"<< endl;
	cout<< "\n*************** Printing the count_element array: END ************\n" << endl;
}

///recursively visit a node, memorizing its ID and finding the intersection that it must store

Transaction  PatriciaTrie::visitNode(PatriciaNode* p, tipoInt* IDpnt, PatriciaNode** pntInHash){


	if(DEEP_DEBUG)
	cout << "p before intersection: count and a_s = " << p->count << " " << p->a_S << endl;


	PatriciaNode* tmp;		///used to construct the new PatriciaNode with intersection
	PatriciaNode* prev;		///used to construct the new PatriciaNode with intersection
	Transaction inter;		///used to memorize the current intersection of the node
	Transaction interTmp;		///used to memorize the intersection that a child passes to his father
	tipoInt minID=(*IDpnt); 	///the minimum ID value of a node in the subtree rooted in p; we can assign it only at the end
	tipoInt numberChil=p->numChildren;
	if (numberChil==0){	///the node is a leaf
		///a leaf has empty intersections
		inter=(tipoInt*)calloc(1,sizeof(tipoInt));
		inter[0]=0;
	}
	else{
		int numVisited=0;	///number of children that were yet visited
		/// if the node has >=1 children, it has a HashTable
		tipoInt num_bucket=p->hash->size;
		if (p->isEmpty==0){///we must visit the node and calculate its intersection
			for (int i=0; i<num_bucket;i++){
				PatriciaNode** pnt=&(p->hash->bucket)+i;
				tmp=(&p->hash->bucket)[i];
				prev=NULL;
				while (tmp!=NULL){
					interTmp=visitNode(tmp,IDpnt,pnt);
					numVisited++;
					///we reallign the ponters
					if (prev==NULL){
						tmp=(&p->hash->bucket)[i];
					}
					else{
						tmp=prev->next;
					}
					pnt=&(tmp->next);
					prev=tmp;
					tmp=tmp->next;
					if(numVisited==1){
						inter=interTmp;
					}
					else{
						///we intersect the transaction and memorize the intersection in the "inter" transaction
						int i=1;	///pointer to elements of inter
						int j=1;	///pointer to element of interTmp
						int res=1;	///pointer to next element of inter that will contain an item in the intersection
						while(i<=inter[0]&&j<=interTmp[0]){
							if (inter[i]==interTmp[j]){
							///inter[i] is in the intersection
								inter[res]=inter[i];
								res++;
								i++;
								j++;
							}
							else{
								if(inter[i]<interTmp[j]){
								///inter[i] is not in the intersection
									i++;
								}
								else{
								///inter[i] can still be in the intersection
									j++;
								}
							}
						}
						inter[0]=res-1;
						free(interTmp);	///the child construct a new intersection (tipoInt array) and passes it to the father: so the father must kill it when it is not still necessary
					}
				}
			}
		}
		if (p->isEmpty>0){///we must ONLY visit the node: there is not need to find the intersection
			for (int i=0; i<num_bucket;i++){
				tmp=(&p->hash->bucket)[i];
				PatriciaNode** pnt=&(p->hash->bucket)+i;
				prev=NULL;
				while (tmp!=NULL){
					interTmp=visitNode(tmp,IDpnt,pnt);
					numVisited++;
					///we reallign the ponters
					if (prev==NULL){
						tmp=(&p->hash->bucket)[i];
					}
					else{
						tmp=prev->next;
					}
					pnt=&(tmp->next);
					prev=tmp;
					tmp=tmp->next;
					free(interTmp);
				}
			}
			///the node has empty intersection
			inter=(tipoInt*)calloc(1,sizeof(tipoInt));
			inter[0]=0;
		}
	}
	tipoInt numItem=(tipoInt) 1+(p->lastItem-&(p->item));
	tipoInt num_buckets;	///used to assign the father pointer
	PatriciaNode* pNext;	///used to assign the father pointer
	PatriciaNode* p2Next;	///used to assign the father pointer
	tipoInt numInter=inter[0];
	///now in inter there is the intersection that must be memorized in the node p
	if(inter[0]==0){

		if(p->hash!=NULL){
			///assign the rigth father pointer to children
			num_buckets=p->hash->size;
			for (int i=0; i< num_buckets; i++){
				pNext=(&p->hash->bucket)[i];
				while(pNext!=NULL){
					p2Next=pNext->next;
					pNext->father=p;
					pNext=p2Next;
				}
			}

			/// we must free the memory for its HashTable before overwrite it using the lastInter
			free(p->hash);
		}
		/// the node has empty intersection: there is no need to resize it
		p->lastInter=p->lastItem;
		///now we can free memory from old node
		///now we construct the intersection to pass it to the father
		free(inter);
		inter=(tipoInt*)calloc((numItem+1),sizeof(tipoInt));
		inter[0]=numItem;
		tipoInt prova;
		for (int i=0; i< numItem; i++){
			prova=(&p->item)[i];
			inter[i+1]=prova;
		}
	}
	else{
		///we must resize the node to insert in it the intersection
		numItem=(tipoInt) 1+(p->lastItem-&(p->item));
		///the new node will have numItem+inter[0] in it: numItem as normal items and inter[0] as intersection
		tmp=(PatriciaNode*)malloc(sizeNode+((numItem+inter[0]-1)*sizeof(tipoInt)));
		///now we can copy the old part of the node
		memcpy(tmp,p,sizeNode+((numItem-1)*sizeof(tipoInt)));
		///we insert the intersection in the new node
		for (int i=0;i< inter[0];i++){
			(&tmp->item)[numItem+i]=inter[i+1];
		}
		numInter=inter[0];

		///now we can free memory from old node
		free(p);
		p=tmp;
		*pntInHash=p;
		if(p->hash!=NULL){
			///assign the rigth father pointer to children
			num_buckets=p->hash->size;
			for (int i=0; i< num_buckets; i++){
				pNext=(&p->hash->bucket)[i];
				while(pNext!=NULL){
					p2Next=pNext->next;
					pNext->father=p;
					pNext=p2Next;
				}
			}

			/// we must free the memory for its HashTable before overwrite it using the lastInter
			free(p->hash);
		}
		///we update the ponters in the node
		tmp->lastItem=&((&tmp->item)[numItem-1]);
		tmp->lastInter=&((&tmp->item)[numItem+numInter-1]);
		///now we construct the intersection to pass it to the father
		free(inter);
		tipoInt sizesA=(numItem+numInter+1)*sizeof(tipoInt);
		inter=(tipoInt*)calloc((numItem+numInter+1),sizeof(tipoInt));
		inter[0]=numItem+numInter;
		for (int i=0; i< inter[0]; i++){
			inter[i+1]=(&p->item)[i];
		}
	}
	p->ID=*IDpnt;
	(*IDpnt)++;
	p->pntToNode=NULL;

	///insert the address of the node in the addressOf array
	addressOf[p->ID]=p;

	///insert the current node to the NodeIDList of every item it contains (not in intersection) and update numNodesPerItem
	NodeIDElem* nodepnt;

	///now update the intersection of every item in the node

	for (int s=numItem-1; s>=0; s--){
		if (IL[(&p->item)[s]].count>=suppMin){

			/// intersect the current suffix of the node with the intersection in HT
			if (HT[(&p->item)[s]].intersection!=NULL){

				///update intersection: we intersect the transactions and memorize the intersection in the "inter" transaction
				tipoInt i=1;	///pointer to elements of HT[l].intersection
				tipoInt j=numItem+numInter-1;	///pointer to element of the node
				tipoInt res=1;	///pointer to next element of intersection that will contain an item in the intersection
				while(i<=HT[(&p->item)[s]].intersection[0]&&j>=s+1){
					if (HT[(&p->item)[s]].intersection[i]==(&p->item)[j]){
						///intersection[i] is in the intersection
						HT[(&p->item)[s]].intersection[res]=HT[(&p->item)[s]].intersection[i];
						res++;
						i++;
						j--;
					}
					else{
						if(HT[(&p->item)[s]].intersection[i]>(&p->item)[j]){
							///HT.intersection[i] is not in the intersection
							i++;
						}
						else{
							///HT.intersection[i] can still be in he intersection
							j--;
						}
					}
				}
				HT[(&p->item)[s]].intersection[0]=res-1;
			}
			else{
				///set the intersection
				if (nextfree+numItem+numInter-s>lastvalid){
					cout << "call manual mem 21" << endl;
					resizeManualMem();
				}
				HT[(&p->item)[s]].intersection=nextfree;
				nextfree=nextfree+numItem+numInter-s;
				for (int o=numItem+numInter-1;o>=s+1; o--){
					HT[(&p->item)[s]].intersection[numItem+numInter-o]=(&p->item)[o];
				}
				HT[(&p->item)[s]].intersection[0]=numItem+numInter-s-1;
			}

			///now update the visitedList of the item
			tipoInt numElem=HT[(&p->item)[s]].visitedList[0].IDnode;
			HT[(&p->item)[s]].visitedList[numElem+1].offPntNode=-1;
			HT[(&p->item)[s]].visitedList[numElem+1].IDnode=p->ID;
			HT[(&p->item)[s]].visitedList[numElem+1].numInInter=-1;
			HT[(&p->item)[s]].visitedList[0].IDnode++;
		}
	}

	return inter;
}


/**visit the tree and assign the ID to nodes and computes the intersections for nodes
*/

void PatriciaTrie::assignIDAndIntersect(){
	///no work for the root
	tipoInt ID=0;
	if (root->hash==NULL){
		return;
	}
	tipoInt numBucket=root->hash->size;
	PatriciaNode* nextNode;
	for (int i=0; i< numBucket; i++){
		nextNode=(&root->hash->bucket)[i];
		PatriciaNode** pnt=&((&root->hash->bucket)[i]);
		PatriciaNode* prev=NULL;
		while(nextNode!=NULL){
			visitNode(nextNode,&ID,pnt);
			///reallign the pointers
			if (prev==NULL){
				nextNode=(&root->hash->bucket)[i];
			}
			else{
				nextNode=prev->next;
			}
			pnt=&(nextNode->next);
			prev=nextNode;
			nextNode=nextNode->next;
		}
	}
	///assign the right father to the children of the root, i.e. NULL
	PatriciaNode* next2Node;
	for (int i=0; i< numBucket; i++){
		nextNode=(&root->hash->bucket)[i];
		while(nextNode!=NULL){
			next2Node=nextNode->next;
			nextNode->father=NULL;
			nextNode=next2Node;
		}
	}
}

/**print the nodes ordered by their ID after the assignment of the PatriciaTrie
*/

void PatriciaTrie::printAfter(){
	for (int l=0; l< numNodes; l++){
		addressOf[l]->print(1);
	}
}

/**print the count distribution of nodes
*/

void PatriciaTrie::printCountDistribution(){
	std::vector<int> count_distribution(1);
	//cout << "sizeof(bool) = " << sizeof(bool) << endl;
	//cout << "sizeof(uint8_t) = " << sizeof(uint8_t) << endl;
	perm_matrix_space = sizeof(bool) * ((double)jp / 1000000.0) * ((long)effect_num_tr);
	cout << "Permutation matrix space = " << perm_matrix_space << endl;
	cout << "Number of nodes: = " << numNodes << endl;
	cout << "Patricia (naive) space = " << sizeof(int) * ((double)jp / 1000000.0) * ((long)numNodes) << endl;
	for (int l=0; l< numNodes; l++){
		if(addressOf[l]->count >=  count_distribution.size()){
			count_distribution.resize(addressOf[l]->count+1);
		}
		count_distribution[addressOf[l]->count]++;
	}
	int compressed_count = 0;
	for (int l=0; l< count_distribution.size() && l < 128; l++){
			compressed_count+=count_distribution[l];
			//cout << "count_distribution["<<l<<"] = " << count_distribution[l] << endl;
	}
	pat_tree_ram = (sizeof(uint8_t) * ((double)jp / 1000000.0) * ((long)compressed_count)) + sizeof(int) * ((double)jp / 1000000.0) * ((long)numNodes-compressed_count);
	cout << "Patricia (compressed) space = " << pat_tree_ram << endl;
}

static int rand_int_(int x){
	int rnd;
	int limit = RAND_MAX - RAND_MAX % x;

	do{
		rnd = rand();
	}while(rnd >= limit);
	return rnd % x;
}


void generatePermutations(int *sampling_counts , bool *permutations_){

	double left_count = (double)sampling_counts[0];
	for(int aj = 1; aj < jp+1; aj++){
		permutations_[aj] = ((double)rand_int_(1000000))*0.000001 <= ((double)sampling_counts[aj])/(left_count);
		sampling_counts[aj] -= permutations_[aj];
	}
	sampling_counts[0]--;
	//cout << endl;
	//cout << endl;

}

/**second scan: "transform" transactions and insert them in the PatriciaTrie
 */

void PatriciaTrie::secondScan_new(){

	double random_time = 0.0;

	cached_a_S = 1000000000.0 / ((double)(sizeof(tipoInt) * jp));

	cout << "cached_a_S = " << cached_a_S << endl;

	for(int aj = 0; aj<=cached_a_S && aj<num_item; aj++){
		IL[aj].a_S_permutations=(tipoInt*)calloc(jp+1,sizeof(tipoInt));
	}

	permutations = (bool*)malloc((jp+1)*sizeof(bool)); // permutation for one transaction
	int *sampling_counts = (int*)malloc((jp+1)*sizeof(int)); // counts of 1 elements to allocate to the i-th permutation
	outcomes = (double*)malloc((jp+1)*sizeof(double)); // outcomes buffer for permutations computation

	for(int aj=1; aj < (jp+1); aj++)
		sampling_counts[aj]=n1;
	sampling_counts[0]=N;

	string line;
	char* tmp;
	string end=""; /// to verify that a string is not null;
	const char *space = " "; ///delimiter
	char* charline; /// to memorize tokens

	///allocate the memory for temp transaction: its size is bounded by the maximum transaction length (max_tr_length) and its size

	temp=(tipoInt*)calloc((max_trans_length+1),sizeof(tipoInt));

	dataset1.open(fileinput.c_str(),ifstream::in); ///<dataset file;

	///allocate the memory for the transaction and its multiplicity
	Transaction t=(tipoInt*)calloc((max_trans_length+1),sizeof(tipoInt)) ;
	///counter for the position  in the transaction
	tipoInt cnt,t_i;
	tipoInt index_tmp;
	numNodes=0;
	///read transactions from dataset
	t_i=0;
	int inserted=0;
	while(!dataset1.eof()){
		cnt=0;
		t[cnt++]=0;
		getline (dataset1,line);
		charline = ( char * )line.c_str();
		///control if the string is empty
		if (end.compare(charline)!=0){
			///first token
			//t[0]++;
			tmp=strtok(charline,space);
			//t[cnt++]=atoi(tmp);
			///second token
			//tmp=strtok(NULL,space);
			while(tmp!=NULL && !isspace(*tmp)){
				t[0]++;
				///find the index of the item
				index_tmp=count_element[atoi(tmp)];
				if(index_tmp!=-1){
					///if the item is frequent, substitute its index in the transaction
					t[cnt++]=index_tmp;
				}
				else{
					///update transaction length
					t[0]--;
					//cout << "infrequent item:" << atoi(tmp) << " " << index_tmp << endl;
				}
				///next item
				tmp=strtok(NULL,space);
			}

			///print the transaction
			/*for(int i=0; i <= t[0]; i++){
				cout << t[i] << " ";
			}
			cout << endl;*/
			///sort it
			tran_sort(t,1,t[0]);

			// generate the permutations for this transaction

			double temp_time_ = get_cpu_time();
			generatePermutations(sampling_counts , permutations);
			random_time += get_cpu_time() - temp_time_;
			permutations[0] = labels[t_i];

			// scan the transaction and update the count_as of singletons using bounded space
			int *temp_prt;
			for(int aj = 1; aj<=t[0]; aj++){
				if(t[aj] <= cached_a_S){
					temp_prt = IL[t[aj]].a_S_permutations;
					for(int h = 0; h<(jp+1); h++){
						temp_prt[h] += permutations[h];
					}
				}
			}

			///insert the transaction and its permutations in the Patricia
			//insert(t,root);
			insert(t,root,permutations);
			inserted++;
		}
		t_i++;
	}

	/*bool aborting = false;
	for(int aj=0; aj < (jp+1); aj++){
		if(sampling_counts[aj] > 0){
			cout << "ERROR on sampling_counts at index aj = " << aj << " sampling_counts[aj] =  " << sampling_counts[aj] << endl;
			aborting = true;
		}
	}
	if(aborting){
		abort();
	}*/

	free(sampling_counts);
	free(permutations);
	free(outcomes);
	free(t);


	if(DEEP_DEBUG)
	cout << "second scan ok"<< endl;

	cout << " random time = " << random_time << endl;
}

/*void PatriciaTrie::secondScan(){
	string line;
	char* tmp;
	string end=""; /// to verify that a string is not null;
	const char *space = " "; ///delimiter
	char* charline; /// to memorize tokens

	///allocate the memory for temp transaction: its size is bounded by the maximum transaction length (max_tr_length) plus one (for the size)

	temp=(tipoInt*)calloc((max_tr_length+1),sizeof(tipoInt));

	dataset1.open(fileinput.c_str(),ifstream::in); ///<dataset file;

	///allocate the memory for the transaction and its multiplicity
	Transaction t=(tipoInt*)calloc((max_tr_length+1),sizeof(tipoInt)) ;
	///counter for the position  in the transaction
	tipoInt cnt,t_i;
	tipoInt index_tmp;
	numNodes=0;
	///read transactions from dataset
	t_i=0;
	int inserted=0;
	while(!dataset1.eof()){
		cnt=0;
		getline (dataset1,line);
		charline = ( char * )line.c_str();
		///control if the string is empty
		if (end.compare(charline)!=0){
			///the first token is the transaction's length
			tmp=strtok(charline,space);
			///in t[0] there is the length of the transaction
			t[cnt++]=atoi(tmp);
			///from the second token items begin
			tmp=strtok(NULL,space);
			while(tmp!=NULL){
				///find the index of the item
				index_tmp=count_element[atoi(tmp)];
				if(index_tmp!=-1){
					///if the item is frequent, substitute its index in the transaction
					t[cnt++]=index_tmp;
				}
				else{
					///update transaction length
					t[0]--;
					cout << "infrequent item:" << atoi(tmp) << " " << index_tmp << endl;
				}
				///next item
				tmp=strtok(NULL,space);
			}

			///print the transaction
			for(int i=0; i <= t[0]; i++){
				cout << t[i] << " ";
			}
			cout << endl;
			///sort it
			tran_sort(t,1,t[0]);
			///insert the transaction in the Patricia
			//insert(t,root);
			insert(t,root,t_i);
			inserted++;
		}
		t_i++;
	}

	if(DEEP_DEBUG)
	cout << "second scan ok"<< endl;
}*/

void visitToPrint(PatriciaNode* pNode,int* prof){
	if(pNode->numChildren==0){
		cout << "\ndepth= " << *prof <<  "\n "<< endl;
		pNode->print(0);
		(*prof)--;
		return;
	}
	tipoInt num_bucket=pNode->hash->size;
	PatriciaNode* pNext;
	for (int i=0; i<num_bucket;i++){
		pNext=(&pNode->hash->bucket)[i];
		while(pNext!=NULL){
			(*prof)++;
			visitToPrint(pNext,prof);
			pNext=pNext->next;
		}
	}

	cout << "\ndepth= " << *prof <<  "\n "<< endl;
	pNode->print(0);
	(*prof)--;
}

void PatriciaTrie::printBefore(){

	cout<< "\n*************** Printing the PatriciaTrie before the assignment of ID and intersection: BEGIN ************\n" << endl;
	int prof=1;
	if (root->hash==NULL){
		cout << "\n The Patricia is empty!!!\n" << endl;
	}
	else{
		tipoInt num_bucket=root->hash->size;
		PatriciaNode* pNext;
		for (int i=0; i<num_bucket;i++){
			pNext=(&root->hash->bucket)[i];
			while(pNext!=NULL){
				prof=1;
				visitToPrint(pNext,&prof);
				pNext=pNext->next;
			}
		}
	}
	cout<< "\n*************** Printing the PatriciaTrie before the assignment of ID and intersection: END  ************\n" << endl;
}


//this method free the memory used by a NodeIDList: the method that invokes it must be sure that nodepnt is not NULL

void deleteNodeIDList(NodeIDElem* nodepnt){
	NodeIDElem* prev=nodepnt;
	nodepnt=nodepnt->next;
	while (nodepnt!=NULL){
		delete prev;
		prev=nodepnt;
		nodepnt=nodepnt->next;
		numPnt--;
	}
	delete prev;
}

/** prints the informations stored in the HT
*/

void HeaderTable::print(){
	cout << "\n***************** HeaderTable: BEGIN ****************\n" << endl;
	VisitedNodeElem* nextVis;
	for (int i=0; i<=maxID; i++){
		cout << "\n---------\n" << endl;
		cout << "\nitem= " << i << "\n"<< endl;
		cout << "\ncount= " << HT[i].count << "\n" << endl;
		if (HT[i].intersection!=NULL){
			tipoInt lengthTrans=HT[i].intersection[0];
			cout << "intersezione=";
			for (int j=0; j< lengthTrans; j++){
				cout << "ok" << endl;
				cout << " " << HT[i].intersection[j+1];
			}
			cout << endl;
			cout << "ok" << endl;
		}
	}
	cout << "\n***************** HeaderTable: END ****************\n" << endl;
}

void printCloEmptyList(){
	NodeIDElem* nextN;
	cout << "\n------------- LISTA DEI NODI PER CLO(EMPTYSET) ----------------------\n" << endl;
	for (int i=0; i<=num_item; i++){
		if (CloEmptyNodeIDList[i]!=NULL){
			cout << "\n------------- ITEM = " << i << " ---------------\n" << endl;
			cout << "lista dei nodi:" ;
			nextN=CloEmptyNodeIDList[i];
			while (nextN!=NULL){
				cout << " " << nextN->nodeID;
				nextN=nextN->next;
			}
			cout << endl;
		}
	}
}

/**return a StartNodeList from a NodeIDList
*/

StartNodeList giveStartNodeList(NodeIDList nodelist ){
	StartNodeList stnl;
	tipoInt toalloc=nodelist[0]*sizeof(StartNode);
	if ((stnl=(StartNode*)malloc((nodelist[0]+1)*sizeof(StartNode)))==NULL){
		cout << "error in allocation!" << endl;
		exit(0);
	}
	stnl[0].ID=nodelist[0];
	for(int i=1;i<=nodelist[0];i++){
		stnl[i].ID=nodelist[i];
		stnl[i].label=0;
	}
	return stnl;
}
