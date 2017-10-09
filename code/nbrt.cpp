/******************************************************************************
* MIT License
*
* Copyright (c) 2017 Siemens AG
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
* 
******************************************************************************/

#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <thread>
#include <random>

#define LABEL 			uint64_t
#define SIZE 			uint8_t
#define KEY				uint64_t
#define KEY_MIN			0	
#define KEY_MAX			0xFFFFFFFFFFFFFFFF	
#define KEY_SIZE		64
#define ALPHABET_CARDINALITY	2
#define ALPHABET_SIZE		4
#define NULL_BIT 		0x01
#define FREEZE			0x02
#define FLAG			0x04

int numThreads;
int searchPer;
int insertPer;
int deletePer;
int duration;
int range;
int seed;


// These variables control the loop for the time specified
volatile bool start = false, steadyState = false, stop = false;

/*===============================================================================================*/

struct Node {
	LABEL label;
	SIZE size;
	bool isLeaf;

	Node(LABEL _label, SIZE _size, bool _isLeaf) {
		this->label = _label;
		this->size = _size;
		this->isLeaf = _isLeaf;
	}
};

struct iNode : public Node{
	uintptr_t child[ALPHABET_SIZE];

	iNode(LABEL _label, SIZE _size, bool _isLeaf) : 
		Node(_label, _size, _isLeaf) {
		for (int i = 0; i < ALPHABET_SIZE; i++)
			child[i] = 0x01;
	}
};

iNode *root;
Node *lMin;
Node *lMax;


// SeekRecord stores the result of traversal The first five members
// gPar, par, curr, currIndex and parIndex are the result of traversal.
// The other members are for determining the throughput of the implementation
struct SeekRecord {
	uintptr_t gPar;
	uintptr_t par;
	uintptr_t curr;

	SIZE currIndex;
	SIZE parIndex;

	uint16_t tid;
	uint64_t casCounter;
	uint64_t successfulInserts;
	uint64_t successfulRemoves;
	uint64_t insertCount;
	uint64_t removeCount;
	uint64_t searchCount;
};

/*===============================================================================================*/
// Initialize the root node of the tree
void init() {
	//std::cout<<"Initializing tree"<<std::endl;
	root = new iNode(0, 0, 0);
	lMin = new Node(KEY_MIN, KEY_SIZE, true);
	lMax = new Node(KEY_MAX, KEY_SIZE, true);
	root->child[0] = (uintptr_t)lMin;
	root->child[ALPHABET_SIZE - 1] = (uintptr_t)lMax;
}

/*===============================================================================================*/

// Returns "true" if the null bit of a pointer is set
bool isNull(uintptr_t node) {
	return (node & NULL_BIT);
}

// Returns "true" if the "freeze" bit of a poiner is set
bool isFrozen(uintptr_t node) {
	return ((node & FREEZE) == FREEZE);
}

// Returns "true" if the "flag" bit of a poiner is set
bool isFlagged(uintptr_t node) {
	return ((node & FLAG) == FLAG);
}

// Returns the address of the node with null, freeze and flag bits set to '0'
uintptr_t getAddr(uintptr_t node) {
	return (node & ~0x07);
}

// Returns the index of the key in the array of child pointers
SIZE getIndex(KEY key, SIZE size) {
	SIZE index = (key >> size);
	return (index & (ALPHABET_SIZE - 1));
}

// Match the prefix of target key and label of the node. Return "true" if prefix matches
bool matchBits(KEY key, LABEL label, SIZE size) {
	return ((key & ((((LABEL)1) << size) - 1)) == label);
}

// Finds the common prefix between target key and node's label. Updates the result in
// commonLabel and commonSize structures
void findCommonBits(KEY key, LABEL label, SIZE size, LABEL *commonLabel, SIZE *commonSize) {
	LABEL temp = (key ^ label);
	*commonSize = 0;
	while (*commonSize < size) {
		if (temp & (ALPHABET_SIZE - 1))
			break;
		temp = (temp >> ALPHABET_CARDINALITY);
		(*commonSize) += ALPHABET_CARDINALITY;
	}
	*commonLabel = key & ((((LABEL)1) << *commonSize) - 1);
	return;
}
// Does compare and swap atomic operation
bool CAS(uintptr_t par, SIZE currIndex, uintptr_t curr, bool currFlag, bool currFreeze, bool currNull, uintptr_t newNode, bool newFlag, bool newFreeze, bool newNull ) {
	uintptr_t oldV = (curr & ~0x07) | (currFlag ? FLAG : 0x00) | (currFreeze ? FREEZE : 0x00) | (currNull ? NULL_BIT : 0);
	uintptr_t newV = (newNode & ~0x07) | (newFlag ? FLAG : 0x00) | (newFreeze ? FREEZE : 0x00) | (newNull ? NULL_BIT : 0);
	return __sync_bool_compare_and_swap(&(((iNode *)par)->child[currIndex]), oldV, newV);
}

// Counts for the number of non-null child pointers of a node
int checkChildCount(uintptr_t node) {
	iNode *curr = (iNode *)node;
	SIZE childCount = 0;
	for (int i = 0; i < ALPHABET_SIZE; i++)
		if (isNull(curr->child[i]) == false)
			childCount++;
	return childCount;
}

/*===============================================================================================*/


// Searches for the presence of a key in the tree
void search(KEY key, SeekRecord* record) {
	uintptr_t gPar = (uintptr_t)root;
	uintptr_t par = (uintptr_t)root;
	
	int parIndex = -1;
	int currIndex = getIndex(key, root->size);
	uintptr_t curr = root->child[currIndex];
	while ( (isNull(curr) == false) && 
		(((iNode *)getAddr(curr))->isLeaf == false) &&
		(matchBits(key, ((iNode *)getAddr(curr))->label, ((iNode *)getAddr(curr))->size)) ) {

		gPar = par;
		parIndex = currIndex;
		par = getAddr(curr);
		currIndex = getIndex(key, ((iNode *)getAddr(curr))->size);
		curr = ((iNode *)getAddr(curr))->child[currIndex];
	}

	record->gPar = gPar;
	record->par = par;
	record->curr = curr;
	record->parIndex = parIndex;
	record->currIndex = currIndex;
	return;
}

// Similar to Search method but looks for the presence of a node in the tree
void searchNode(uintptr_t node, SeekRecord *record) {
	uintptr_t gPar = (uintptr_t)root;
	uintptr_t par = (uintptr_t)root;
	LABEL label = ((Node *)node)->label;
	
	int parIndex = -1;
	int currIndex = getIndex(label, root->size);
	uintptr_t curr = root->child[currIndex];
	while ( (isNull(curr) == false) && 
		(getAddr(curr) != node) &&
		(((iNode *)getAddr(curr))->isLeaf == false) &&
		(matchBits(label, ((iNode *)getAddr(curr))->label, ((iNode *)getAddr(curr))->size)) ) {

		gPar = par;
		parIndex = currIndex;
		par = getAddr(curr);
		currIndex = getIndex(label, ((iNode *)getAddr(curr))->size);
		curr = ((iNode *)getAddr(curr))->child[currIndex];
	}

	record->gPar = gPar;
	record->par = par;
	record->curr = curr;
	record->parIndex = parIndex;
	record->currIndex = currIndex;
	return;
}

/*===============================================================================================*/

bool contains(KEY key, SeekRecord *record) {
	int currIndex = getIndex(key, root->size);
	uintptr_t curr = root->child[currIndex];
	
	while ( (isNull(curr) == false) && 
		(((iNode *)getAddr(curr))->isLeaf == false) &&
		(matchBits(key, ((iNode *)getAddr(curr))->label, ((iNode *)getAddr(curr))->size)) ) {
		currIndex = getIndex(key, ((iNode *)getAddr(curr))->size);
		curr = ((iNode *)getAddr(curr))->child[currIndex];
	}

	return (isNull(curr) ? false : ((((Node *)curr)->label == key) && (((Node *)curr)->size == KEY_SIZE)));
}

/*===============================================================================================*/


// Removes an internal node from the tree
bool removeINode(uintptr_t par, uintptr_t curr, SIZE currIndex) {
	int childCount = checkChildCount(curr);
	uintptr_t ptr;
	bool null = false;
	if (childCount == 0)
		ptr = 0, null = 1;
	else if (childCount == 1) {
		for (SIZE i = 0; i < ALPHABET_SIZE;i++) {
			if (isNull(((iNode *)curr)->child[i]) == false) {
				ptr = getAddr(((iNode *)curr)->child[i]);
				break;
			}
		}
	}
	else {
		iNode *newNode = new iNode(((iNode *)curr)->label, ((iNode *)curr)->size, false);
		ptr = (uintptr_t)newNode;
		for (SIZE i = 0; i < ALPHABET_SIZE; i++) {
			if (isNull(((iNode *)curr)->child[i]) == false) {
				((iNode *)ptr)->child[i] = getAddr(((iNode *)curr)->child[i]);
			}
		}
	}
	return CAS(par, currIndex, curr, 1, 0, 0, ptr, 0, 0, null); 
}

/*===============================================================================================*/

// Updates the "freeze" bit of all the child pointers to '1'
void helpFreeze(uintptr_t node) {
	for (SIZE i = 0; i < ALPHABET_SIZE; i++) {
		while (true) {
			uintptr_t child = ((iNode *)node)->child[i];
			bool flagged = isFlagged(child);
			bool freeze = isFrozen(child);
			if (flagged) {
				helpFreeze(getAddr(child));
				removeINode(node, getAddr(child), i);
			}
			else if (freeze) {
				break;
			}
			else if (CAS(node, i, child, 0, 0, isNull(child), child, 0, 1, isNull(child)))
				break;
		}
	}
}

void checkConsistencyAndRemove(uintptr_t gPar, uintptr_t par, uintptr_t curr, SIZE currIndex, SeekRecord* record) {
	
	int childCount = checkChildCount(curr);
	if (childCount > 1)
		return;

	SIZE parIndex = -1;
	if (par == NULL) {
		searchNode(curr, record);
		if (getAddr(record->curr) == curr) {
			gPar = record->gPar;
			par = record->par;
			currIndex = record->currIndex;
			parIndex = record->parIndex;
		}
		else 
			return;
	}

	bool res = CAS(par, currIndex, curr, 0, 0, 0, curr, 1, 0, 0);
	if (res) {
		helpFreeze(curr);
		removeINode(par, curr, currIndex);
		return checkConsistencyAndRemove(NULL, gPar, par, parIndex, record);
	}
	else {
		uintptr_t newCurr = ((iNode *)par)->child[currIndex];
		if (isNull(newCurr)) 
			return checkConsistencyAndRemove(NULL, gPar, par, parIndex, record);
		else if (isFlagged(newCurr))
			return;
		else if (isFrozen(newCurr)) {
			helpFreeze(par);
			searchNode(par, record);
			if (getAddr(record->curr) == par) {
				gPar = record->par;
				par = getAddr(record->curr);
				parIndex = record->currIndex;
				removeINode(gPar, par, parIndex);
				checkConsistencyAndRemove(NULL, record->gPar, record->par, record->parIndex, record);
			}
			return checkConsistencyAndRemove(NULL, NULL, curr, currIndex, record);
		}
		else if (((iNode *)getAddr(par))->child[currIndex] != curr)
			return checkConsistencyAndRemove(NULL, NULL, curr, currIndex, record);
	}	
}


// Removes a node containing the target key from the tree
bool removeT(KEY key, SeekRecord* record) {
	search(key, record);
	
	uintptr_t gPar = record->gPar;
	uintptr_t par = record->par;
	uintptr_t curr = record->curr;
	SIZE currIndex = record->currIndex;
	SIZE parIndex = record->parIndex;

	bool null = isNull(curr);
	bool flagged = isFlagged(curr);
	bool frozen = isFrozen(curr);
	if (null || flagged)
		return false;
	else {
		iNode *node = (iNode *)getAddr(curr);
		if ((node->label != key) || (node->size != KEY_SIZE))
			return false;

		if (frozen) {
			helpFreeze(par);
			removeINode(gPar,par, parIndex);
			return removeT(key, record);
		}
		else if (CAS(par, currIndex, curr, 0, 0, 0, NULL, 0, 0, 1)) {
			record->successfulRemoves++;
			checkConsistencyAndRemove(NULL, gPar, par, parIndex, record);
			return true;
		}
		else
			return removeT(key, record);
	}
}
/*===============================================================================================*/

// Inserts a node containing the target key into the tree
bool insert(KEY key, SeekRecord* record) {
	search(key, record);

	uintptr_t gPar = record->gPar;
	uintptr_t par = record->par;
	uintptr_t curr = record->curr;
	SIZE currIndex = record->currIndex;
	SIZE parIndex = record->parIndex;

	iNode *newNode;
	bool null = isNull(curr);
	bool flagged = isFlagged(curr);
	bool frozen = isFrozen(curr);
	LABEL commonLabel;
	SIZE commonSize;


	if (isNull(curr)) {
		commonLabel = key;
		commonSize = KEY_SIZE;
		newNode = new iNode(key, KEY_SIZE, true);
		null = 1;
	}
	else {
		curr = getAddr(curr);
		if ( (((iNode *)curr)->label == key) && (((iNode *)curr)->size == KEY_SIZE))
			return false;
		findCommonBits(key, ((iNode *)curr)->label, ((iNode *)curr)->size, &commonLabel, &commonSize);
	
		newNode = new iNode(commonLabel, commonSize, false);
		
		iNode *newLeafNode = new iNode(key, KEY_SIZE, true);
		SIZE newLeafIndex = getIndex(key, commonSize);
		newNode->child[newLeafIndex] = (uintptr_t)newLeafNode;
	
		SIZE newCurrIndex = getIndex(((iNode *)curr)->label, commonSize);
		newNode->child[newCurrIndex] = curr;
		null = 0;
	}
	
	if (CAS(par, currIndex, curr, 0, 0, null, (uintptr_t)newNode, 0, 0, 0)) {
		record->successfulInserts++;
		return true;
	}
	else {
		uintptr_t newCurr  = ((iNode *)par)->child[currIndex];
		bool flagged = isFlagged(newCurr);
		bool frozen = isFrozen(newCurr);
		if (flagged) {
			helpFreeze(getAddr(newCurr));
			removeINode(par, getAddr(newCurr), currIndex);
			return insert(key, record);
		}
		else if (frozen) {
			helpFreeze(par);
			removeINode(gPar, par, parIndex);
			return insert(key, record);
		}
		else 
			return insert(key, record);
	}
}

/*===============================================================================================*/

// Counts the number of nodes in the trie
void trieSize(uintptr_t node, uint64_t *size) {
	if (isNull(node))
		return;
	iNode *newNode = (iNode *)node;
	if (newNode->size == KEY_SIZE) {
		if (newNode == lMin || newNode == lMax)
			return;
		else
			*size = (*size)+1;
		return;
	}
	else {
		int childCount = 0;
		for (int i = 0; i < ALPHABET_SIZE; i++) {
			trieSize(newNode->child[i], size);
			if (isNull((uintptr_t)(newNode->child[i])))
				childCount++;
		}
		if (childCount > (ALPHABET_SIZE - 2)) {
			std::cout<<"ERROR: Inconsistent nodes found"<<std::endl;
			exit(1);
		}
	}
}
/*===============================================================================================*/

// Checks if the trie is valid
void isValidTrie(iNode *root, uint64_t succIns, uint64_t succRem, int initialSize) {
	std::cout<<"Validating trie..."<<std::endl;
	uint64_t size = 0;
	trieSize((uintptr_t)root, &size);	
	if (size != (initialSize + succIns - succRem)) {
		std::cout<<"Invalid trie sizes: "<<std::endl; 
		std::cout<<"Final Size: "<<size<<std::endl;
		std::cout<<"Successful inserts: "<<succIns<<std::endl;
		std::cout<<"Successful removes: "<<succRem<<std::endl;
		std::cout<<"Inital size: "<<initialSize<<std::endl;
		exit(1);
	}
}

void printTrie_Normal(uintptr_t node) {
	if (isNull(node))
		return;
	//node = getAddr(node);
	iNode *newNode = (iNode *)node;
	if (newNode->size == KEY_SIZE) {
		if (newNode == lMin || newNode == lMax)
			return;
		else
			std::cout<<newNode->label<<std::endl;
		return;
	}
	else {
		for (int i = 0; i < ALPHABET_SIZE; i++)
			printTrie_Normal(newNode->child[i]);
	}
}

void printTrie(uintptr_t node) {
	if (isNull(node))
		return;
	//node = getAddr(node);
	iNode *newNode = (iNode *)node;
	if (newNode->size == KEY_SIZE) {
		if (newNode == lMin || newNode == lMax)
			return;
		else
			std::cout<<"VARUN: "<<newNode->label<<std::endl;
		return;
	}
	else {
		for (int i = 0; i < ALPHABET_SIZE; i++)
			printTrie(newNode->child[i]);
	}
}
/*===============================================================================================*/

void test_insertRemoveParallel() {
	const int NUM_THREADS = 1000;
	srand(time(NULL));
	int arr[NUM_THREADS];
	for (int i = 0; i < NUM_THREADS; i++)
		arr[i]= rand() + 1;

	SeekRecord *ins[NUM_THREADS];
	for (int i = 0; i < NUM_THREADS; i++)
		ins[i] = new SeekRecord();
	SeekRecord *rem[NUM_THREADS];
	for (int i = 0; i < NUM_THREADS; i++)
		rem[i] = new SeekRecord();

	std::cout<<"Inserting elements"<<std::endl;	
	std::vector<std::thread> insT(NUM_THREADS);
	for (int i = 0; i < NUM_THREADS; i++) {
		insT[i] = std::thread(&insert, arr[i], ins[i]);
	}
	for (int i = 0; i < NUM_THREADS; i++)
		insT[i].join();
	printTrie_Normal((uintptr_t)root);
//	
	uint64_t totalSucc = 0;
	for (int i = 0; i < NUM_THREADS; i++)
		totalSucc += ins[i]->successfulInserts;

	if (totalSucc != NUM_THREADS) {
		std::cout<<"ERROR"<<std::endl;
		exit(1);
	}


	std::cout<<"Removing elements..."<<std::endl;
	std::vector<std::thread> remT(NUM_THREADS);
	for (int i = 0; i < NUM_THREADS; i++) {
		remT[i] = std::thread(&removeT, arr[i], rem[i]);
	}
	for (int i = 0; i < NUM_THREADS; i++)
		remT[i].join();
	//testTrie();
	printTrie((uintptr_t)root);
	std::cout<<"Printing after removing elements..."<<std::endl;
}

int *arr;

void *operateOnTree(void *arg) {
	SeekRecord *record = (SeekRecord*)arg;
	
	//std::mt19937 myrand(seed + record->tid);
	std::mt19937 myrand(seed);
	std::uniform_int_distribution<int> operDist(0, 100);
	auto operationGenerator = std::bind(operDist, myrand);	
	int mySearchPer = searchPer;
	int myInsertPer = insertPer + searchPer;	

	int chooseOperation;
	while (!start) {
	}

	record->successfulInserts = 0;	
	record->successfulRemoves = 0;	
	while (!steadyState) {
		chooseOperation = operationGenerator();
		uint32_t key = myrand() % range + 1;

		if (chooseOperation < mySearchPer) {
			contains(key, record);
		}
		else if (chooseOperation < (myInsertPer))  {
			insert(key, record);
		}
		else {
			removeT(key, record);
		}
	}

	while (!stop) {
		chooseOperation = operationGenerator();
		uint32_t key = myrand() % range + 1;

		if (chooseOperation < mySearchPer) {
			contains(key, record);
			record->searchCount++;
		}
		else if (chooseOperation < (myInsertPer))  {
			insert(key, record);
			record->insertCount++;
		}
		else {
			removeT(key, record);
			record->removeCount++;
		}
	}

}

void test_simultaneousInsertRemove() {
	struct timespec runTime, transientTime;
	transientTime.tv_sec = 0;
	transientTime.tv_nsec = 2000000;
	runTime.tv_sec = duration;
	runTime.tv_nsec = 0;
	SeekRecord* record[numThreads];
	for (int i = 0; i < numThreads; i++) {
		record[i] = new SeekRecord();
		record[i]->tid = i;
		record[i]->successfulInserts = 0;
		record[i]->successfulRemoves = 0;
		record[i]->insertCount = 0;
		record[i]->removeCount = 0;
		record[i]->searchCount = 0;
		record[i]->casCounter = 0;
	}

	std::mt19937 myRand(seed);
	for (int i = 0; i < range; i++)
		arr[i] = (myRand()) + 1; 
	int i = 0;
	while (i < range/2) {
		if (insert(myRand()%range + 1,  record[0]))
			i++;
	}

	uint64_t size = 0;
	trieSize((uintptr_t)root, &size);
	std::cout<<"Initial size of trie is : "<<size<<std::endl;

	pthread_t threadArray[numThreads];
	for(int i=0;i<numThreads;i++)
        {
                pthread_create(&threadArray[i], NULL, operateOnTree, (void*) record[i] );
        }

        start=true; //start operations
        nanosleep(&transientTime,NULL); //warmup
        steadyState=true;
        nanosleep(&runTime,NULL);
        stop=true;      //stop operations

        for(int i=0;i<numThreads;i++)
        {
                pthread_join(threadArray[i], NULL);
        }


	uint64_t totalInsertCount = 0;
	uint64_t totalRemoveCount = 0;
	uint64_t totalSearchCount = 0;
	uint64_t succIns = 0;
	uint64_t succRem = 0;

	for (int i = 0; i < numThreads; i++) {
		totalInsertCount += record[i]->insertCount;
		totalRemoveCount += record[i]->removeCount;
		totalSearchCount += record[i]->searchCount;
		succIns += record[i]->successfulInserts;
		succRem += record[i]->successfulRemoves;
	}

	uint64_t totalOperations =  totalInsertCount + totalRemoveCount + totalSearchCount;
	double mops = ((totalOperations)/(duration * 1000000.0));
	std::cout<<"Throughput is : "<<mops<<" MOPS"<<std::endl;
	std::cout<<"SearchPer: "<<searchPer<<std::endl;
	std::cout<<"InsertPer: "<<insertPer<<std::endl;
	std::cout<<"RemovePer: "<<deletePer<<std::endl;
	std::cout<<"Duration: "<<duration<<std::endl;
	std::cout<<"TotalInsertCount: "<<totalInsertCount<<std::endl;
	std::cout<<"TotalRemoveCount: "<<totalRemoveCount<<std::endl;
	std::cout<<"TotalSearchcount: "<<totalSearchCount<<std::endl;
	std::cout<<"Total operations: "<<totalOperations<<std::endl;
	isValidTrie(root, succIns, succRem, size);

		
}


int main(int argc, char **argv) {
	std::cout<<"Size: "<<sizeof(iNode)<<std::endl;
	std::cout<<"Recoed Size: "<<sizeof(SeekRecord)<<std::endl;
	numThreads = atoi(argv[1]);
	searchPer = atoi(argv[2]);
	insertPer = atoi(argv[3]);
	deletePer = atoi(argv[4]);
	duration = atoi(argv[5]);
	range = atoi(argv[6]);
	seed = atoi(argv[7]);
	arr = (int *)malloc(range * sizeof(int));
	init();
	// For testing the getIndex function
	//test_getIndex();

	// For testing the findCommonAlphabets function
	//test_findCommonBits();

	// For testing sequential version of insert
	//test_insertSequential();
	
	// For testing parallel removes only
	//test_removeParallel();

	// Both insert and remove in parallel but seperately
//	test_insertRemoveParallel();

	test_simultaneousInsertRemove();
	
	return 0;
}

