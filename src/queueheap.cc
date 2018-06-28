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
#include "wy_.cc"
#include <queue>


///default constructor: it allocate the memory for the heap and initialize the descriptor variables

Queue::Queue(){
	///allocate the space for a empty heap
	heap=(Itemset*)calloc(1024,sizeof(Itemset));
	last=-1;
	maxSize=1024;
}

///doubles the memory space for the queue; it is used before insert an itemset in the queue

void Queue::doubleQueue(){
	///allocate the new space
	Itemset* tmpHeap=(Itemset*) calloc(2*maxSize+1,sizeof(Itemset));
	///copy the past heap
	memcpy(tmpHeap,heap,(last+1)*sizeof(Itemset));
	///free the memory used by last heap
	free(heap);
	///assign the new heap
	heap=tmpHeap;
	maxSize=2*maxSize+1;
	cout << "max queue doubled to " << maxSize << endl;
}

///insert a element in the queue;

void Queue::insert(Itemset* its){
	///if the heap is full, double it
	if ((last+1)==maxSize){
		doubleQueue();
	}
	///insert the new element in last position
	heap[last+1].support=its->support;
	heap[last+1].info=its->info;
	int current=last+1;	///current position of the new element
	///up-heap bubbling
	while(current!=0){	///while the new element is != root
		if(heap[(current-1)/2].support>=its->support){	///comparison with the father
			break;
		}
		///swap
		heap[current].support=heap[(current-1)/2].support;
		heap[current].info=heap[(current-1)/2].info;
		heap[(current-1)/2].support=its->support;
		heap[(current-1)/2].info=its->info;
		current=(current-1)/2;
	}
	last++;
	max_queue_last=last;
	max_queue_size=maxSize;
}

Itemset* Queue::removeMax(){
	Itemset* max=(Itemset*)calloc(1,sizeof(Itemset));	///itemset to return
	//.cout << "created itemset (6) with adress "<<max<<endl;
	living_itemsets++;
	if (last==-1) return NULL;	///the queue is empty
	max->support=heap[0].support;
	max->info=heap[0].info;
	Itemset* tmp=(Itemset*)calloc(1,sizeof(Itemset));
	//.cout << "created itemset (7) with adress "<<tmp<<endl;
	living_itemsets++;
	if (last==0){	///the heap has only one element
		heap[0].support=0;
		heap[0].info=NULL;
		///update last
		last--;
	}
	else{	///down-heap bubbling
		heap[0].support=heap[last].support;
		heap[0].info=heap[last].info;
		heap[last].support=0;
		heap[last].info=NULL;
		///update last
		last--;
		int current=0;	///current position of the element to swap
		while(current<((last+1)/2)){	///current is not leaf
			if ((((2*current)+2)>last)||((heap[(2*current)+1].support)>=(heap[(2*current)+2].support))){	///choose the left child
				if ((heap[(2*current)+1].support)>heap[current].support){	///swap
					tmp->support=heap[current].support;
					tmp->info=heap[current].info;
					heap[current].support=heap[(2*current)+1].support;
					heap[current].info=heap[(2*current)+1].info;
					heap[(2*current)+1].support=tmp->support;
					heap[(2*current)+1].info=tmp->info;
					///update current
					current=(2*current)+1;
				}
				else{
					break;
				}
			}
			else{	///choose the right child
				if ((heap[(2*current)+2].support)>heap[current].support){
					tmp->support=heap[current].support;
					tmp->info=heap[current].info;
					heap[current].support=heap[(2*current)+2].support;
					heap[current].info=heap[(2*current)+2].info;
					heap[(2*current)+2].support=tmp->support;
					heap[(2*current)+2].info=tmp->info;
						///update current
					current=(2*current)+2;
				}
				else{
					break;
				}
			}
		}
	}
	//.cout << "free itemset (7) with adress "<< tmp <<endl;
	free(tmp);
	living_itemsets--;
	max_queue_last=last;
	///return the max
	return max;
}


/** print the infos associated with the itemset
*/

void Itemset::print(){
	tipoInt no;
	cout << "\n------------ Itemset infos: BEGIN ----------------\n" << endl;
	cout << "support: " << support << endl;
	/* TopKWY */
	cout << "target class support (a_S): " << info->a_S[0] << endl;
	/* TopKWY */
	cout << "minItem: " << info->minItem << endl;
	cout << "coreIndex: " << info->coreIndex << endl;
	cout << "typelist: " << info->typeList << endl;
	cout << "NodeIDList: " << endl;
	if (info->typeList==0){
		cout << "number of nodes: " << info->IDlist[0] << endl;
		for (int i=1; i<=info->IDlist[0];i++){
			cout << " " << info->IDlist[i];
		}
		cout << endl;
	}
	else{
		tipoInt off=1;
		cout << "number of nodes: " << info->IDlist[0] << endl;
		for (int i=1; i<=info->IDlist[0];i++){
			cout << "\nnode " << i <<" : " << info->IDlist[off] << ", count: " << info->IDlist[off+1] << endl;
			cout << "intersection: " ;
			if (info->IDlist[off+2]!=-1){
				for (int j=1;j<=info->IDlist[off+2];j++){
					cout <<info->IDlist[off+2+j] << " " ;
				}
				cout << endl;
				off=off+info->IDlist[off+2]+3;
			}
			else{
				cout << "in the node!" << endl;
				off=off+3;
			}
		}
		cout << endl;
	}
	cout << "prefixCI: " << endl;
	for (int i=1; i<=info->prefixCI[0];i++){
		cout << " " << info->prefixCI[i];
	}
	cout << endl;
	cout << "\n------------ Itemset infos: END ----------------\n" << endl;
}

/** print the notExtList
*/

void printNotExpandedList(notExpList notexlist ){
	notExpanded* tmp=notexlist;
	cout << "\n--------- print notExtList: BEGIN ---------------\n" << endl;
	while(tmp!=NULL){
		tmp->its->print();
		tmp=tmp->next;
	}
	cout << "\n--------- print notExtList: END ---------------\n" << endl;
}

///default constructor of QueueMin: it allocate the memory for the heap and initialize the descriptor variables

QueueMin::QueueMin(){
	///allocate the space for a empty heap
	heap=(Itemset*)calloc(1024,sizeof(Itemset));
	last=-1;
	maxSize=1024;
}

///doubles the memory space for the queue; it is used before insert an itemset in the queue

void QueueMin::doubleQueueMin(){

	// the duplication needs to be performed only if the queue is actually full
	// we check the size of the max queue to be sure
	if(max_queue_size >= maxSize){
		///allocate the new space
		Itemset* tmpHeap=(Itemset*) calloc(2*maxSize+1,sizeof(Itemset));
		///copy the past heap
		memcpy(tmpHeap,heap,(last+1)*sizeof(Itemset));
		///free the memory used by last heap
		free(heap);
		///assign the new heap
		heap=tmpHeap;
		maxSize=2*maxSize+1;

		cout << "min queue doubled to " << maxSize << endl;
	}
	else{
		// we can free space by keeping in queue only itemsets having support lower than the current support upper bound
		int j = 0;
		for(int i=0; i <= last; i++){
			if(heap[i].support<s_supp){
				j++;
			}
		}
		maxSize = (j+1)*2+1;
		Itemset* tmpHeap=(Itemset*) calloc(maxSize,sizeof(Itemset));

		j = 0;
		for(int i=0; i <= last; i++){
			if(heap[i].support<s_supp){
				tmpHeap[j]=heap[i];
				j++;
			}
		}

		last = j;
		free(heap);
		heap=tmpHeap;
		cout << "no min queue duplication: only copied " << (last+1) << " items" << endl;
		cout << "min queue extended to " << maxSize << endl;
	}
}

///insert a element in the queue;

void QueueMin::insert(Itemset* its){

	queue_program_state = -1;

	///if the heap is full, double it
	if ((last+1)==maxSize){
		queue_program_state = -2;
		doubleQueueMin();
		queue_program_state = -3;
	}
	///insert the new element in last position
	heap[last+1].support=its->support;
	heap[last+1].info=its->info;
	int current=last+1;	///current position of the new element

	queue_program_state = -4;


	///up-heap bubbling
	while(current!=0){	///while the new element is != root
		if(heap[(current-1)/2].support<=its->support){	///comparison with the father
			break;
		}
		///swap
		heap[current].support=heap[(current-1)/2].support;
		heap[current].info=heap[(current-1)/2].info;
		heap[(current-1)/2].support=its->support;
		heap[(current-1)/2].info=its->info;
		current=(current-1)/2;
	}
	last++;

	queue_program_state = 0;

}

Itemset* QueueMin::removeMin(){
	if (last==-1) return NULL;	///the queue is empty
	Itemset* min=(Itemset*)calloc(1,sizeof(Itemset));	///itemset to return
	//.cout << "created itemset (8) with adress "<< min <<endl;
	living_itemsets++;
	min->support=heap[0].support;
	min->info=heap[0].info;
	Itemset* tmp=(Itemset*)calloc(1,sizeof(Itemset));
	//.cout << "created itemset (9) with adress "<< tmp <<endl;
	living_itemsets++;
	if (last==0){	///the heap has only one element
		heap[0].support=0;
		heap[0].info=NULL;
		///update last
		last--;
	}
	else{	///down-heap bubbling
		heap[0].support=heap[last].support;
		heap[0].info=heap[last].info;
		heap[last].support=0;
		heap[last].info=NULL;
		///update last
		last--;
		int current=0;	///current position of the element to swap
		while(current<((last+1)/2)){	///current is not leaf
			if ((((2*current)+2)>last)||((heap[(2*current)+1].support)<=(heap[(2*current)+2].support))){	///choose the left child
				if ((heap[(2*current)+1].support)<heap[current].support){	///swap
					tmp->support=heap[current].support;
					tmp->info=heap[current].info;
					heap[current].support=heap[(2*current)+1].support;
					heap[current].info=heap[(2*current)+1].info;
					heap[(2*current)+1].support=tmp->support;
					heap[(2*current)+1].info=tmp->info;
					///update current
					current=(2*current)+1;
				}
				else{
					break;
				}
			}
			else{	///choose the right child
				if ((heap[(2*current)+2].support)<heap[current].support){
					tmp->support=heap[current].support;
					tmp->info=heap[current].info;
					heap[current].support=heap[(2*current)+2].support;
					heap[current].info=heap[(2*current)+2].info;
					heap[(2*current)+2].support=tmp->support;
					heap[(2*current)+2].info=tmp->info;
						///update current
					current=(2*current)+2;
				}
				else{
					break;
				}
			}
		}
	}
	//.cout << "free itemset (8) with adress "<< tmp <<endl;
	free(tmp);
	living_itemsets--;
	///return the min
	return min;
}

Itemset* QueueMin::getMin(){
	if (last==-1) return NULL;	///the queue is empty
	Itemset* min=(Itemset*)calloc(1,sizeof(Itemset));	///itemset to return
	//.cout << "created itemset (10) with adress "<< min <<endl;
	living_itemsets++;
	min->support=heap[0].support;
	min->info=heap[0].info;
	///return the min
	return min;
}

///returns the number of items that are not frequent and are in the intersection associated with an itemset

tipoInt Queue::notFreqInLists(ItemList* IL){
	tipoInt notFreq=0;
	for (int i=0; i<= last; i++){
		if (heap[i].support>=suppMin&&heap[i].info->typeList==1){
			tipoInt curr=1;
			for (int j=0; j<heap[i].info->IDlist[0];j++){
				if (heap[i].info->IDlist[curr+2]==-1){
					curr+=3;
				}
				else{
					int leng=heap[i].info->IDlist[curr+2];
					for (int k=0; k<leng; k++){
						if(IL[heap[i].info->IDlist[curr+3+k]].count<suppMin){
							notFreq++;
						}
						if (k>0&&heap[i].info->IDlist[curr+3+k]>=heap[i].info->IDlist[curr+3+k-1]){
							cout << "Error!" << endl;
						}
					}
					curr=curr+leng+3;
				}
			}
		}
	}
	return notFreq;
}

///returns the space that can be saved using a different coding for the intersection list
tipoInt Queue::spaceSaved(){
	tipoInt saved=0;
	for (int i=0; i<= last; i++){
		if (heap[i].support>=suppMin&&heap[i].info->typeList==1){
			tipoInt curr=1;
			for (int j=0; j<heap[i].info->IDlist[0];j++){
				if (heap[i].info->IDlist[curr+2]==-1){
					curr+=3;
					saved++;
				}
				else{
					int leng=heap[i].info->IDlist[curr+2];
					curr=curr+leng+3;
				}
			}
		}
	}
	return saved;
}

///returns the total space used by lists
long long Queue::spaceUsedByLists(){
	long long used=0;
	long long max=0;
	for (int i=0; i<= last; i++){
		if (heap[i].support>=suppMin&&heap[i].info->typeList==1){
			tipoInt curr=1;
			for (int j=0; j<heap[i].info->IDlist[0];j++){
				if (heap[i].info->IDlist[curr+2]==-1){
					curr+=3;
				}
				else{
					int leng=heap[i].info->IDlist[curr+2];
					curr=curr+leng+3;
				}
			}
			used+=curr;
			if (curr>max) max=curr;
		}
		else{
			if (heap[i].support>=suppMin&&heap[i].info->typeList==0){
				used=used+heap[i].info->IDlist[0]+1;
			}
		}
	}
	double MB=1048576;
	cout << "MAX SIZE OF A LIST (MB): " << ((double) max)/MB*sizeof(tipoInt) << endl;
	return used;
}




// Min-Max queue implementation for TopKWY


///default constructor of QueueMinMax: it allocate the memory for the heap and initialize the descriptor variables

QueueMinMax::QueueMinMax(){
	///allocate the space for a empty heap
	heap=(Itemset*)calloc(1024,sizeof(Itemset));
	last=-1;
	maxSize=1024;
	min_support=-1;
	min_index=-1;
}

///doubles the memory space for the queue; it is used before insert an itemset in the queue

void QueueMinMax::doubleQueueMaxMin(){

	#ifdef trace_state
	queue_program_state = 1;
	#endif

	///allocate the new space
	Itemset* tmpHeap=(Itemset*) calloc(2*maxSize+1,sizeof(Itemset));
	///copy the past heap
	memcpy(tmpHeap,heap,(last+1)*sizeof(Itemset));
	///free the memory used by last heap
	free(heap);
	///assign the new heap
	heap=tmpHeap;
	maxSize=2*maxSize+1;
	cout << "minmax queue doubled to " << maxSize << endl;

	#ifdef trace_state
	queue_program_state = 0;
	#endif

}

///insert a element in the queue;

void QueueMinMax::insert(Itemset* its){

	#ifdef trace_state
	queue_program_state = 2;
	#endif

	///if the heap is full, double it
	if ((last+1)==maxSize){
		doubleQueueMaxMin();
	}
	///insert the new element in last position
	heap[last+1].support=its->support;
	heap[last+1].info=its->info;
	int current=last+1;	///current position of the new element
	///up-heap bubbling
	while(current!=0){	///while the new element is != root
		if(heap[(current-1)/2].support>=its->support){	///comparison with the father
			break;
		}
		///swap
		heap[current].support=heap[(current-1)/2].support;
		heap[current].info=heap[(current-1)/2].info;
		heap[(current-1)/2].support=its->support;
		heap[(current-1)/2].info=its->info;
		current=(current-1)/2;
	}
	// update minimum element if needed
	if(min_support==-1 || its->support < min_support){
		min_support=its->support;
		min_index=current;
		/*cout << "-------------" << endl;
		cout << "insert:" << endl;
		cout << "new min element has support " << min_support << endl;
		cout << "stored in position " << min_index << endl;
		cout << "heap contains " << heap[min_index].support << endl;
		cout << "-------------" << endl;*/
	}
	last++;

	#ifdef trace_state
	queue_program_state = 0;
	#endif

}

Itemset* QueueMinMax::removeMax(){

	#ifdef trace_state
	queue_program_state = 3;
	#endif


	Itemset* max=(Itemset*)calloc(1,sizeof(Itemset));	///itemset to return
	//.cout << "created itemset (6) with adress "<<max<<endl;
	living_itemsets++;
	if (last==-1) return NULL;	///the queue is empty
	max->support=heap[0].support;
	max->info=heap[0].info;
	Itemset* tmp=(Itemset*)calloc(1,sizeof(Itemset));
	//.cout << "created itemset (7) with adress "<<tmp<<endl;
	living_itemsets++;
	if (last==0){	///the heap has only one element
		heap[0].support=0;
		heap[0].info=NULL;
		///update last
		last--;
	}
	else{	///down-heap bubbling
		heap[0].support=heap[last].support;
		heap[0].info=heap[last].info;
		heap[last].support=0;
		heap[last].info=NULL;

		///update last
		last--;
		int current=0;	///current position of the element to swap
		while(current<((last+1)/2)){	///current is not leaf
			if ((((2*current)+2)>last)||((heap[(2*current)+1].support)>=(heap[(2*current)+2].support))){	///choose the left child
				if ((heap[(2*current)+1].support)>heap[current].support){	///swap
					tmp->support=heap[current].support;
					tmp->info=heap[current].info;
					heap[current].support=heap[(2*current)+1].support;
					heap[current].info=heap[(2*current)+1].info;
					heap[(2*current)+1].support=tmp->support;
					heap[(2*current)+1].info=tmp->info;
					///update current
					current=(2*current)+1;
				}
				else{
					break;
				}
			}
			else{	///choose the right child
				if ((heap[(2*current)+2].support)>heap[current].support){
					tmp->support=heap[current].support;
					tmp->info=heap[current].info;
					heap[current].support=heap[(2*current)+2].support;
					heap[current].info=heap[(2*current)+2].info;
					heap[(2*current)+2].support=tmp->support;
					heap[(2*current)+2].info=tmp->info;
						///update current
					current=(2*current)+2;
				}
				else{
					break;
				}
			}
		}
		if(heap[current].support == min_support){
			min_index = current;
			//cout << "updated min_index to " << min_index << endl;
		}
	}
	//.cout << "free itemset (7) with adress "<< tmp <<endl;
	free(tmp);
	living_itemsets--;
	max_queue_last=last;

	#ifdef trace_state
	queue_program_state = 0;
	#endif

	///return the max
	return max;
}

void QueueMinMax::printLowestHeapLevel(){

	cout << "-------------" << endl;
	cout << "lowest heap level:" << endl;

	// test: print the lowest heap level
	int j = (last-1)/2;
	j++;
	cout << "last = " << last << endl;
	cout << "j = " << j << endl;
	// scan last heap level
	while(j <= last){
		cout << "heap[" << j << "].support = " << heap[j].support << endl;
		j++;
	}
	cout << "-------------" << endl;

}

void QueueMinMax::removeMin(){

	#ifdef trace_state
	queue_program_state = 4;
	#endif


	//cout << "-------------" << endl;
	//cout << "remove min:" << endl;
	if (last==-1) return;	///the queue is empty

	ItsInfo* todel;
	Itemset* tmp=(Itemset*)calloc(1,sizeof(Itemset));
	//.cout << "created itemset (7) with adress "<<tmp<<endl;
	living_itemsets++;

	if (last==0){	///the heap has only one element
		todel = heap[0].info;
		heap[0].support=0;
		heap[0].info=NULL;
		///update last
		last--;
		//update min index
		min_index=-1;
		min_support=-1;
	}
	else{
		cout << "removing min value " << min_support << endl;
		cout << "queue size " << (last+1) << endl;
		/*cout << "position of min value " << min_index << endl;
		cout << "heap contains " << heap[min_index].support << endl;
		cout << "-------------" << endl;
		cout << "operations:" << endl;*/
		if(heap[min_index].support != min_support){
			cout << "ERROR: MIN ITEM IS IN WRONG LOCATION" << endl;
			int j = (last-1)/2;
			j++;
			min_support=heap[0].support;
			min_index=0;
			while(j <= last){
				if(heap[j].support < min_support){
					min_support = heap[j].support;
					min_index = j;
				}
				j++;
			}
		}
		todel=heap[min_index].info;
		// swap last with removed node
		heap[min_index].support=heap[last].support;
		heap[min_index].info=heap[last].info;
		heap[last].support=0;
		heap[last].info=NULL;
		///update last
		last--;
		//cout << "heap[" << min_index << "].support = "<< heap[min_index].support << endl;
		int current=min_index;	///current position of the element to swap
		/// we need to performs an up-heap since we are removing a leaf node
		///up-heap bubbling
		while(current!=0){	///while the new element is != root
			if(heap[(current-1)/2].support>=heap[current].support){	///comparison with the father
				break;
			}
			///swap
			/*cout << "swap operation! " <<  endl;
			cout << "indexes: " << current << " -> " << (current-1)/2 <<  endl;
			cout << "values: " << heap[current].support << " -> " << heap[(current-1)/2].support <<  endl;*/
			tmp->support=heap[current].support;
			tmp->info=heap[current].info;
			heap[current].support=heap[(current-1)/2].support;
			heap[current].info=heap[(current-1)/2].info;
			heap[(current-1)/2].support=tmp->support;
			heap[(current-1)/2].info=tmp->info;
			current=(current-1)/2;
		}

		//cout << "-------------" << endl;

		// now, update the minimum value and its index
		// we scan only the last heap level
		int j = (last-1)/2;
		j++;
		/*cout << "last = " << last << endl;
		cout << "j = " << j << endl;*/
		//reset min
		min_support=heap[0].support;
		min_index=0;
		//cout << "min reset to max = " << min_support << endl;
		// scan last heap level
		while(j <= last){
			//cout << "heap[" << j << "].support = " << heap[j].support << endl;
			if(heap[j].support < min_support){
				min_support = heap[j].support;
				min_index = j;
				//cout << "min set to = " << min_support << endl;
			}
			j++;
		}
		//cout << "new min value set to " << min_support << endl;
	}
	//.cout << "free itemset (7) with adress "<< tmp <<endl;
	free(tmp);
	living_itemsets--;
	// free memory or space for the minimum itemset
	delete todel;
	//cout << "-------------" << endl;

	#ifdef trace_state
	queue_program_state = 0;
	#endif


}

tipoInt QueueMinMax::getMin(){
	///return the cached minimum support
	if(min_index > 0)
		return min_support;
	return -1;
}

tipoInt QueueMinMax::getMax(){
	///return the support of max element, which is stored as first element on the heap
	if(last > 0)
		return heap[0].support;
	return -1;
}













// QueueMinMax for Results for TopKWY


// Min-Max queue implementation for TopKWY


///default constructor of QueueMinMax: it allocate the memory for the heap and initialize the descriptor variables

QueueMinMaxResult::QueueMinMaxResult(){
	///allocate the space for a empty heap
	heap=(ResultPattern*)calloc(1024,sizeof(ResultPattern));
	last=-1;
	maxSize=1024;
	max_p_value=-1.0;
	max_index=-1;
}

///doubles the memory space for the queue; it is used before insert an itemset in the queue

void QueueMinMaxResult::doubleQueueMinMaxResult(){

	#ifdef trace_state
	queue_program_state = 5;
	#endif


	///allocate the new space
	ResultPattern* tmpHeap=(ResultPattern*) calloc(2*maxSize+1,sizeof(ResultPattern));
	///copy the past heap
	memcpy(tmpHeap,heap,(last+1)*sizeof(ResultPattern));
	///free the memory used by last heap
	free(heap);
	///assign the new heap
	heap=tmpHeap;
	maxSize=2*maxSize+1;
	cout << "minmax result queue doubled to " << maxSize << endl;

	#ifdef trace_state
	queue_program_state = 0;
	#endif


}


int QueueMinMaxResult::getQueueSize(){
	return last;
}

///insert a element in the queue;

void QueueMinMaxResult::insert(ResultPattern* res){

	#ifdef trace_state
	queue_program_state = 6;
	#endif


	///if the heap is full, double it
	if ((last+1)==maxSize){
		doubleQueueMinMaxResult();
	}

	if(res->support==0){
		cout << "ERROR! INSERTED WITH 0 SUPPORT" << endl;
	}

	///insert the new element in last position
	heap[last+1].itemset=res->itemset;
	heap[last+1].support=res->support;
	heap[last+1].a_S=res->a_S;
	heap[last+1].p_value=res->p_value;
	int current=last+1;	///current position of the new element
	///up-heap bubbling
	while(current!=0){	///while the new element is != root
		if(heap[(current-1)/2].p_value<=res->p_value){	///comparison with the father
			break;
		}
		///swap
		heap[current].itemset=heap[(current-1)/2].itemset;
		heap[current].support=heap[(current-1)/2].support;
		heap[current].a_S=heap[(current-1)/2].a_S;
		heap[current].p_value=heap[(current-1)/2].p_value;

		heap[(current-1)/2].itemset=res->itemset;
		heap[(current-1)/2].support=res->support;
		heap[(current-1)/2].a_S=res->a_S;
		heap[(current-1)/2].p_value=res->p_value;

		current=(current-1)/2;
	}
	// update minimum element if needed
	if(max_p_value < -0.5 || res->p_value > max_p_value){
		max_p_value=res->p_value;
		max_index=current;
		/*cout << "-------------" << endl;
		cout << "result insert:" << endl;
		cout << "new max element has p_value " << max_p_value << endl;
		cout << "stored in position " << max_index << endl;
		cout << "heap contains " << heap[max_index].p_value << endl;
		cout << "-------------" << endl;*/
	}
	else{
		/*cout << "-------------" << endl;
		cout << "result insert:" << endl;
		cout << "new element has p_value " << res->p_value << endl;
		cout << "stored in position " << current << endl;
		cout << "heap contains " << heap[current].p_value << endl;
		cout << "-------------" << endl;*/
	}
	last++;

	#ifdef trace_state
	queue_program_state = 0;
	#endif


}

ResultPattern* QueueMinMaxResult::removeMin(){

	#ifdef trace_state
	queue_program_state = 7;
	#endif


	if (last==-1) return NULL;	///the queue is empty
	ResultPattern* min=(ResultPattern*)calloc(1,sizeof(ResultPattern));	///itemset to return

	min->itemset=heap[0].itemset;
	min->support=heap[0].support;
	min->a_S=heap[0].a_S;
	min->p_value=heap[0].p_value;

	ResultPattern* tmp=(ResultPattern*)calloc(1,sizeof(ResultPattern));

	if (last==0){	///the heap has only one element
		heap[0].support=0;
		heap[0].a_S=0;
		heap[0].p_value=1.0;
		heap[0].itemset=NULL;
		///update last
		last--;
	}
	else{	///down-heap bubbling
		heap[0].support=heap[last].support;
		heap[0].a_S=heap[last].a_S;
		heap[0].p_value=heap[last].p_value;
		heap[0].itemset=heap[last].itemset;

		heap[last].support=0;
		heap[last].a_S=0;
		heap[last].p_value=1.0;
		heap[last].itemset=NULL;


		///update last
		last--;
		int current=0;	///current position of the element to swap
		while(current<((last+1)/2)){	///current is not leaf
			if ((((2*current)+2)>last)||((heap[(2*current)+1].p_value)<=(heap[(2*current)+2].p_value))){	///choose the left child
				if ((heap[(2*current)+1].p_value)<heap[current].p_value){	///swap
					tmp->support=heap[current].support;
					tmp->a_S=heap[current].a_S;
					tmp->p_value=heap[current].p_value;
					tmp->itemset=heap[current].itemset;

					heap[current].support=heap[(2*current)+1].support;
					heap[current].a_S=heap[(2*current)+1].a_S;
					heap[current].p_value=heap[(2*current)+1].p_value;
					heap[current].itemset=heap[(2*current)+1].itemset;

					heap[(2*current)+1].support=tmp->support;
					heap[(2*current)+1].a_S=tmp->a_S;
					heap[(2*current)+1].p_value=tmp->p_value;
					heap[(2*current)+1].itemset=tmp->itemset;

					///update current
					current=(2*current)+1;
				}
				else{
					break;
				}
			}
			else{	///choose the right child
				if ((heap[(2*current)+2].p_value)<heap[current].p_value){
					tmp->support=heap[current].support;
					tmp->a_S=heap[current].a_S;
					tmp->p_value=heap[current].p_value;
					tmp->itemset=heap[current].itemset;

					heap[current].support=heap[(2*current)+2].support;
					heap[current].a_S=heap[(2*current)+2].a_S;
					heap[current].p_value=heap[(2*current)+2].p_value;
					heap[current].itemset=heap[(2*current)+2].itemset;

					heap[(2*current)+2].support=tmp->support;
					heap[(2*current)+2].a_S=tmp->a_S;
					heap[(2*current)+2].p_value=tmp->p_value;
					heap[(2*current)+2].itemset=tmp->itemset;
						///update current
					current=(2*current)+2;
				}
				else{
					break;
				}
			}
		}
		/*if(heap[current].p_value == max_p_value){
			max_index = current;
			//cout << "updated max_index to " << min_index << endl;
		}*/
	}
	//cout << "free tmp" << endl;
	free(tmp);

	#ifdef trace_state
	queue_program_state = 0;
	#endif


	///return the min
	return min;
}

void QueueMinMaxResult::printLowestHeapLevel(){

	if(last < 0) return;

	cout << "-------------" << endl;
	cout << "min item: " << heap[0].p_value << endl;
	cout << "lowest heap level:" << endl;

	// test: print the lowest heap level
	int j = 0;//(last-1)/2;
	//j++;
	cout << "last = " << last << endl;
	cout << "j = " << j << endl;
	// scan last heap level
	while(j <= last){
		cout << "heap[" << j << "].support = " << heap[j].support << endl;
		cout << "heap[" << j << "].p_value = " << heap[j].p_value << endl;
		j++;
	}
	cout << "-------------" << endl;

}

void QueueMinMaxResult::removeMax(){

	#ifdef trace_state
	queue_program_state = 8;
	#endif


	//cout << "remove max call" << endl;
	//cout << "-------------" << endl;
	//cout << "remove min:" << endl;
	if (last==-1) return;	///the queue is empty

	Transaction todel;
	ResultPattern* tmp=(ResultPattern*)calloc(1,sizeof(ResultPattern));

	if (last==0){	///the heap has only one element
		todel = heap[0].itemset;
		heap[0].support=0;
		heap[0].a_S=0;
		heap[0].p_value=1.0;
		heap[0].itemset=NULL;
		///update last
		last--;
		//update min index
		max_index=-1;
		max_p_value=-1.0;
	}
	else{
		/*cout << "removing max result value " << max_p_value << endl;
		cout << "queue size " << (last+1) << endl;*/
		/*cout << "position of min value " << min_index << endl;
		cout << "heap contains " << heap[min_index].support << endl;
		cout << "-------------" << endl;
		cout << "operations:" << endl;*/
		//if(heap[max_index].p_value != max_p_value){
			//cout << "ERROR: MAX RESULT IS IN WRONG LOCATION" << endl;
			int j = (last-1)/2;
			j++;
			max_p_value=heap[last].p_value;
			max_index=last;
			while(j <= last){
				if(heap[j].p_value > max_p_value){
					max_p_value = heap[j].p_value;
					max_index = j;
				}
				j++;
			}

			//cout << "removing p-value " << max_p_value << " at position " << max_index << endl;
		//}
		todel=heap[max_index].itemset;
		// swap last with removed node
		heap[max_index].support=heap[last].support;
		heap[max_index].a_S=heap[last].a_S;
		heap[max_index].p_value=heap[last].p_value;
		heap[max_index].itemset=heap[last].itemset;

		heap[last].support=0;
		heap[last].a_S=0;
		heap[last].p_value=1.0;
		heap[last].itemset=NULL;
		///update last
		last--;
		//cout << "heap[" << max_index << "].support = "<< heap[max_index].support << endl;
		int current=max_index;	///current position of the element to swap
		/// we need to performs an up-heap since we are removing a leaf node
		///up-heap bubbling
		while(current!=0){	///while the new element is != root
			if(heap[(current-1)/2].p_value<=heap[current].p_value){	///comparison with the father
				break;
			}
			///swap
			/*cout << "swap operation! " <<  endl;
			cout << "indexes: " << current << " -> " << (current-1)/2 <<  endl;
			cout << "values: " << heap[current].support << " -> " << heap[(current-1)/2].support <<  endl;*/
			tmp->support=heap[current].support;
			tmp->a_S=heap[current].a_S;
			tmp->p_value=heap[current].p_value;
			tmp->itemset=heap[current].itemset;

			heap[current].support=heap[(current-1)/2].support;
			heap[current].a_S=heap[(current-1)/2].a_S;
			heap[current].p_value=heap[(current-1)/2].p_value;
			heap[current].itemset=heap[(current-1)/2].itemset;

			heap[(current-1)/2].support=tmp->support;
			heap[(current-1)/2].a_S=tmp->a_S;
			heap[(current-1)/2].p_value=tmp->p_value;
			heap[(current-1)/2].itemset=tmp->itemset;

			current=(current-1)/2;
		}

		//cout << "-------------" << endl;

		// now, update the minimum value and its index
		// we scan only the last heap level
		/*int */j = (last-1)/2;
		j++;
		/*cout << "last = " << last << endl;
		cout << "j = " << j << endl;*/
		//reset min
		max_p_value=heap[0].p_value;
		max_index=0;
		//cout << "min reset to max = " << min_support << endl;
		// scan last heap level
		while(j <= last){
			//cout << "heap[" << j << "].support = " << heap[j].support << endl;
			if(heap[j].p_value > max_p_value){
				max_p_value = heap[j].p_value;
				max_index = j;
				//cout << "min set to = " << min_support << endl;
			}
			j++;
		}
		//cout << "new min value set to " << min_support << endl;
	}
	// free memory or space for the minimum itemset
	//cout << "free todel" << endl;
	free(todel);
	//cout << "free tmp 2" << endl;
	free(tmp);
	//cout << "-------------" << endl;

	#ifdef trace_state
	queue_program_state = 0;
	#endif


}

double QueueMinMaxResult::getMax(){
	///return the cached maximum support
	if(max_index > 0)
		return max_p_value;
	return -1.0;
}

double QueueMinMaxResult::getMin(){
	///return the support of max element, which is stored as first element on the heap
	if(last > 0)
		return heap[0].p_value;
	return -1.0;
}




// test: novel data structure for itemset ordering

void QueueMinMax_test::insert(Itemset* its){

	#ifdef trace_state
	queue_program_state = 9;
	#endif

	Itemset* its_=(Itemset*)malloc(sizeof(Itemset));
	its_->support = its->support;
	its_->info = its->info;

	/*if(data[its_->support].capacity() < 1024){
		data[its_->support].reserve(1024);
	}*/

	//data[its_->support].push_front(its_);
	data[its_->support].push_back(its_);

	if(elements_stored == 0){
		min_index = its_->support;
		max_index = its_->support;
	}
	else{
		// update min index
		if(its_->support < min_index)
			min_index = its_->support;
		// update max index
		if(its_->support > max_index)
			max_index = its_->support;
	}
	elements_stored++;

	#ifdef trace_state
	queue_program_state = 0;
	#endif

}

void QueueMinMax_test::removeMin(){

	#ifdef trace_state
	queue_program_state = 10;
	#endif

	if(elements_stored > 0){
		elements_stored = elements_stored - data[min_index].size();
		for(int i=0; i < data[min_index].size(); i++){
			delete data[min_index].at(i)->info;
			free(data[min_index].at(i));
		}
		//cout << "clear min level " << min_index << " size " << data[min_index].size() << " cap " << data[min_index].capacity() << endl;
		std::vector<Itemset*> foo (0);
		foo.swap(data[min_index]);
		//while(data[min_index].size() > 0){
			//delete data[min_index].front()->info;
			//free(data[min_index].front());
			//data[max_index].pop_front();
		//}
		if(elements_stored > 0){
			min_index++;
			while(data[min_index].size()==0)
				min_index++;
			//cout << "new min index after removeMin = " << min_index << endl;
		}
		else{
			min_index = -1;
			max_index = -1;
		}
	}

	#ifdef trace_state
	queue_program_state = 0;
	#endif


}

Itemset* QueueMinMax_test::removeMax(){

	#ifdef trace_state
	queue_program_state = 11;
	#endif


	if(elements_stored > 0){
		Itemset* max;
		elements_stored = elements_stored - 1;
		//max = data[max_index].front();
		//data[max_index].pop_front();
		max = data[max_index].back();
		data[max_index].pop_back();
		if(elements_stored > 0){
			if(data[max_index].size() == 0){

				//cout << "clear level " << max_index << " cap " << data[max_index].capacity() << endl;
				if(data[max_index].capacity() > 0){
					std::vector<Itemset*> foo (0);
					foo.swap(data[max_index]);
					//cout << "after level " << max_index << " cap " << data[max_index].capacity() << endl;
				}

				max_index--;
				while(data[max_index].size()==0 && max_index > -1){
					//cout << "clear another level " << max_index << " cap " << data[max_index].capacity() << endl;
					if(data[max_index].capacity() > 0){
						std::vector<Itemset*> foo_ (0);
						foo_.swap(data[max_index]);
						//cout << "after level " << max_index << " cap " << data[max_index].capacity() << endl;
					}
					max_index--;
				}
			}
		}
		else{
			min_index = -1;
			max_index = -1;
		}

		#ifdef trace_state
		queue_program_state = 0;
		#endif


		return max;
	}

	#ifdef trace_state
	queue_program_state = 0;
	#endif


	return NULL;
}


Itemset* QueueMinMax_test::removeLow(){

	#ifdef trace_state
	queue_program_state = 12;
	#endif


	if(elements_stored > 0){
		Itemset* max;
		elements_stored = elements_stored - 1;
		int low_index = (int)(((double)suppMin) * low_factor);
		if(low_index < min_index)
			low_index = min_index;
		while(data[low_index].size() == 0 && low_index > min_index)
			low_index--;

		max = data[low_index].back();
		data[low_index].pop_back();
		if(elements_stored > 0){
			while(data[min_index].size()==0)
				min_index++;
			//cout << "new min index after removeLow = " << min_index << endl;
		}
		else{
			min_index = -1;
			max_index = -1;
		}

		#ifdef trace_state
		queue_program_state = 0;
		#endif


		return max;
	}

	#ifdef trace_state
	queue_program_state = 0;
	#endif


	return NULL;
}

tipoInt QueueMinMax_test::getMin(){
	return min_index;
}

tipoInt QueueMinMax_test::getMax(){
	return max_index;
}

tipoInt QueueMinMax_test::getInqueue(){
	return elements_stored;
}

QueueMinMax_test::QueueMinMax_test(tipoInt max_support){

	// allocate data
	data.resize(max_support);
	/*for(int i=0; i < data.size(); i++){
		data[i].reserve(1024);
	}*/
	min_index = -1;
	max_index = -1;
	elements_stored = 0;
	low_factor = 1.5;

}


double QueueMinMax_test::getQueueMemory(){

	double queue_memory = 0.0;

	int total_capacity = 0;

	for(int i=0; i < data.size(); i++){
		total_capacity += data[i].capacity();
	}

	queue_memory += (double)data.capacity() * (sizeof(data[0]));
	queue_memory += (double)total_capacity * (sizeof(Itemset));
	queue_memory += (double)elements_stored * (sizeof(ItsInfo));

	ItsInfo* temp_info;

	for(int i=min_index; i <= max_index; i++){
		for(int j = 0; j < data[i].size(); j++){
			temp_info = data[i][j]->info;

			//queue_memory += temp_info->nl_nodes.capacity() * (sizeof(tipoInt) + sizeof(temp_info->nl_nodes[0]));
			queue_memory += (temp_info->minnodes.capacity() + temp_info->minnodes_indexes.capacity()) * (sizeof(tipoInt));
			queue_memory += temp_info->sizelista * (sizeof(tipoInt));

			/*for(int aj = 0; aj < temp_info->nl_nodes.size(); aj++){
				queue_memory += (sizeof(tipoInt)) * temp_info->nl_nodes[aj].capacity();
			}*/
		}
	}

	return queue_memory / 1000000.0;

}






// test: novel data structure for results ordering

void QueueMinMax_Restest::insert(ResultPattern* res){

	#ifdef trace_state
	queue_program_state = 9;
	#endif

	//int index = (int)(-res->log_p_value);
	//cout << "index for insertion = " << index << endl;
	/*int index_inefficient = 0;
	while(res->log_p_value < getLogPsi(index_inefficient) && index_inefficient < data.size()){
		index_inefficient++;
	}
	index_inefficient--;*/

	int max_index_bsearch = data.size() - 1;
	int min_index_bsearch = 0;
	int index = (int)(((double)( max_index_bsearch + min_index_bsearch )) / 2.0 + 0.5);
	bool found = false;
	while(!found){
	/*cout << "min_index_bsearch = " << min_index_bsearch << " max_index_bsearch = " << max_index_bsearch << " index = " << index << endl;*/
		if(res->log_p_value < getLogPsi(index)){
			min_index_bsearch = index;
			if(max_index_bsearch - min_index_bsearch < 6){
				// find it linearly
				while(res->log_p_value < getLogPsi(index) && index < data.size()){
					index++;
				}
				index--;
				found = true;
			}
			else{
				index = (int)(((double)( max_index_bsearch + min_index_bsearch )) / 2.0 + 0.5);
			}
		}
		else{
			max_index_bsearch = index;
			if(max_index_bsearch - min_index_bsearch < 6){
				// find it linearly
				index = min_index_bsearch;
				while(res->log_p_value < getLogPsi(index) && index < data.size()){
					index++;
				}
				index--;
				found = true;
			}
			else{
				index = (int)(((double)( max_index_bsearch + min_index_bsearch )) / 2.0 + 0.5);
			}
		}
	}
	/*cout << "min_index_bsearch = " << min_index_bsearch << " max_index_bsearch = " << max_index_bsearch << " index = " << index << endl;
	cout << "res->log_p_value = " << res->log_p_value << " getLogPsi(index) = " << getLogPsi(index) << " index = " << index << endl;

	if(index_inefficient != index){
		cout << "ERROR in index finding: index_inefficient " << index_inefficient << " index " << index << endl;
		abort();
	}*/
	/*if(data[its_->support].capacity() < 1024){
		data[its_->support].reserve(1024);
	}*/

	//data[its_->support].push_front(its_);
	data[index].push_back(res);

	if(elements_stored == 0){
		min_index = index;
		max_index = index;
		min_value = res->p_value;
		max_value = res->p_value;
		min_log_value = res->log_p_value;
		/*cout << "min_log_value initialized to = " << min_log_value << endl;
		cout << "min_value initialized to = " << min_value << endl;
		cout << "min_index initialized to = " << min_index << endl;*/
	}
	else{
		// update min index and min value
		/*if(res->p_value < min_value){
			min_value = res->p_value;
			min_index = index;
			cout << "min_value updated to = " << min_value << endl;
			cout << "min_index updated to = " << min_index << endl;
		}*/
		// update max index and max value
		if(res->p_value > max_value || index < max_index){
			max_value = res->p_value;
			max_index = index;
		}
		// update min index and min log value
		if(res->log_p_value < min_log_value || index > min_index){
			min_log_value = res->log_p_value;
			min_value = res->p_value;
			min_index = index;
			/*cout << "min_value updated to = " << min_value << endl;
			cout << "min_log_value updated to = " << min_log_value << endl;
			cout << "min_index updated to = " << min_index << endl;*/
		}
	}
	elements_stored++;

	//cout << "new element with p-val " << res->p_value << " inserted at position " << index << endl;
	//cout << "min_value " << min_value << " max_value " << max_value << endl;

	#ifdef trace_state
	queue_program_state = 0;
	#endif

}

void QueueMinMax_Restest::insert_observed(double log_observed_pval){

	// find the position for this observed pattern
	/*int index_inefficient = 0;
	while(log_observed_pval < getLogPsi(index_inefficient) && index_inefficient < observed.size()){
		index_inefficient++;
	}
	index_inefficient--;*/

	observed_p_values.push(log_observed_pval);

	// binary search of index of observed pattern
	int max_index_bsearch = observed.size() - 1;
	int min_index_bsearch = 0;
	int index = (int)(((double)( max_index_bsearch + min_index_bsearch )) / 2.0 + 0.5);
	bool found = false;
	while(!found){
		if(log_observed_pval < getLogPsi(index)){
			min_index_bsearch = index;
			if(max_index_bsearch - min_index_bsearch < 6){
				// find it linearly
				while(log_observed_pval < getLogPsi(index) && index < observed.size()){
					index++;
				}
				index--;
				found = true;
			}
			else{
				index = (int)(((double)( max_index_bsearch + min_index_bsearch )) / 2.0 + 0.5);
			}
		}
		else{
			max_index_bsearch = index;
			if(max_index_bsearch - min_index_bsearch < 6){
				// find it linearly
				index = min_index_bsearch;
				while(log_observed_pval < getLogPsi(index) && index < observed.size()){
					index++;
				}
				index--;
				found = true;
			}
			else{
				index = (int)(((double)( max_index_bsearch + min_index_bsearch )) / 2.0 + 0.5);
			}
		}
	}

	/*if(index_inefficient != index){
		cout << "ERROR in index finding: index_inefficient " << index_inefficient << " index " << index << endl;
		abort();
	}*/
	//cout << "observed.size() " << observed.size() << " index " << index << endl;
	observed[index]++;
	if(index < max_index_obs || elements_observed == 0){
		max_index_obs = index;
		//cout << "new observed inserted, max index set to " << max_index_obs << endl;
	}
	elements_observed++;
}

// removal of non-significant patterns
void QueueMinMax_Restest::removeMax(){

	#ifdef trace_state
	queue_program_state = 10;
	#endif

	int current_max_index = max_index_obs;

	// observed elements
	if(elements_observed > 0){
		elements_observed = elements_observed - observed[max_index_obs];
		observed[max_index_obs] = 0;
		if(elements_observed > 0){
			while(max_index_obs < observed.size() && observed[max_index_obs] == 0){
				max_index_obs++;
			}
		}
	}

	if(elements_observed == 0){
		max_index_obs = -1;
	}

	//cout << "after remove max, max index set to " << max_index_obs << endl;

	// stored elements
	if(elements_stored > 0 && max_index == current_max_index){
		elements_stored = elements_stored - data[max_index].size();
		elements_observed = elements_observed - observed[max_index];
		for(int i=0; i < data[max_index].size(); i++){
			//cout << " log-p-value of deleted pattern " << data[max_index].at(i)->log_p_value << endl;
			free(data[max_index].at(i)->itemset);
			free(data[max_index].at(i));
		}
		//cout << "clear res max level " << max_index << " size " << data[max_index].size() << " cap " << data[max_index].capacity() << endl;
		std::vector<ResultPattern*> foo (0);
		foo.swap(data[max_index]);
		//while(data[min_index].size() > 0){
			//delete data[min_index].front()->info;
			//free(data[min_index].front());
			//data[max_index].pop_front();
		//}
		if(elements_stored > 0){
			max_index++;
			while(data[max_index].size()==0)
				max_index++;
			max_value = getPsi(max_index+1);
			//cout << "new max index after removeMax = " << max_index << endl;
			//cout << "new max value after removeMax = " << max_value << endl;
		}
		else{
			min_index = -1;
			max_index = -1;
			min_value = -1.0;
			max_value = -1.0;
			min_log_value = 0;
		}
	}

	#ifdef trace_state
	queue_program_state = 0;
	#endif


}

// output a significant pattern
ResultPattern* QueueMinMax_Restest::removeMin(){

	#ifdef trace_state
	queue_program_state = 11;
	#endif

	// observed elements
	if(elements_observed > 0){
		elements_observed = elements_observed - 1;
		observed[min_index]--;
	}

	if(elements_observed == 0){
		max_index_obs = -1;
	}

	//cout << "after remove min, max index set to " << max_index_obs << endl;

  // stored elements
	if(elements_stored > 0){
		ResultPattern* max;
		elements_stored = elements_stored - 1;
		//max = data[max_index].front();
		//data[max_index].pop_front();
		max = data[min_index].back();
		data[min_index].pop_back();
		if(elements_stored > 0){
			if(data[min_index].size() == 0){

				//cout << "clear level " << max_index << " cap " << data[max_index].capacity() << endl;
				if(data[min_index].capacity() > 0){
					std::vector<ResultPattern*> foo (0);
					foo.swap(data[min_index]);
					//cout << "after level " << max_index << " cap " << data[max_index].capacity() << endl;
				}

				min_index--;
				while(data[min_index].size()==0 && min_index > 0){
					//cout << "clear another level " << min_index << " cap " << data[min_index].capacity() << endl;
					if(data[min_index].capacity() > 0){
						std::vector<ResultPattern*> foo_ (0);
						foo_.swap(data[min_index]);
						//cout << "after level " << min_index << " cap " << data[min_index].capacity() << endl;
					}
					min_index--;
				}
				min_value = getPsi(min_index);
				min_log_value = getLogPsi(min_index);
			}
		}
		else{
			min_index = -1;
			max_index = -1;
			min_value = -1;
			max_value = -1;
			min_log_value = 0;
		}

		/*cout << "after removeMin min_index = " << min_index << endl;
		cout << "min_value = " << min_value << endl;
		cout << "min_log_value = " << min_log_value << endl;*/

		#ifdef trace_state
		queue_program_state = 0;
		#endif


		return max;
	}

	#ifdef trace_state
	queue_program_state = 0;
	#endif


	return NULL;
}


double QueueMinMax_Restest::getMin(){
	return min_value;
}


double QueueMinMax_Restest::getLogMin(){
	return min_log_value;
}

double QueueMinMax_Restest::getMax(){
	return max_value;
}

int QueueMinMax_Restest::getMaxIndex(){
	return max_index_obs;
	//return max_index;
}

tipoInt QueueMinMax_Restest::getInqueue(){
	return elements_stored;
}

tipoInt QueueMinMax_Restest::getKElements(){
	/*cout << "getKElements call!" << endl;
	cout << "   max_index = " << max_index << endl;
	cout << "   min_index = " << min_index << endl;
	cout << "   elements_stored = " << elements_stored << endl;
	cout << "   data[max_index].size() = " << data[max_index].size() << endl;
	for(int aj=0; aj<data.size(); aj++)
		if(data[aj].size() > 0)
			cout << "   position " << aj << " contains " << data[aj].size() << endl;*/
	//if(elements_stored == 0) return 0;
	//return (elements_stored - data[max_index].size());
	/*cout << "getKElements call! (continue from above)" << endl;
	cout << "   elements_observed " << elements_observed << endl;
	cout << "   max_index_obs " << max_index_obs << endl;
	cout << "   observed[max_index_obs] " << observed[max_index_obs] << endl;
	for(int aj=0; aj<observed.size(); aj++)
		if(observed[aj] > 0)
			cout << "   position " << aj << " contains " << observed[aj] << endl;*/
	if(elements_observed == 0) return 0;
	return (elements_observed - observed[max_index_obs]);
}

void QueueMinMax_Restest::printElements(){
	cout << "printElements call!" << endl;
	cout << "   max_index = " << max_index << endl;
	cout << "   min_index = " << min_index << endl;
	cout << "   elements_stored = " << elements_stored << endl;
	if(max_index >= 0)
	cout << "   data[max_index].size() = " << data[max_index].size() << endl;
	cout << "data vector:" << endl;
	for(int aj=0; aj<data.size(); aj++)
		if(data[aj].size() > 0)
			cout << "   position " << aj << " contains " << data[aj].size() << endl;
	//if(elements_stored == 0) return 0;
	//return (elements_stored - data[max_index].size());
	cout << "printElements call! (continue from above)" << endl;
	cout << "   elements_observed " << elements_observed << endl;
	cout << "   max_index_obs " << max_index_obs << endl;
	if(max_index_obs >= 0)
	cout << "   observed[max_index_obs] " << observed[max_index_obs] << endl;
	cout << "observed vector:" << endl;
	for(int aj=0; aj<observed.size(); aj++)
		if(observed[aj] > 0)
			cout << "   position " << aj << " contains " << observed[aj] << endl;
	cout << "printElements call end!" << endl;
}

QueueMinMax_Restest::QueueMinMax_Restest(tipoInt max_support){

	// allocate data
	data.resize(max_support+1);
	observed.resize(max_support+1);
	/*for(int i=0; i < data.size(); i++){
		data[i].reserve(1024);
	}*/
	min_index = -1;
	max_index = -1;
	elements_stored = 0;
	elements_observed = 0;
	max_index_obs = -1;

}
