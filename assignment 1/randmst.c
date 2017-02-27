/*
 *
 *	randmst.c
 *
 *	Implements Kruskal's algorithm for finding Minimum Spanning Tree over
 *	a randomly-generated complete undirected graph. 
 *
 *	Edge weights are calculated as Euclidean distance between vertex locations,
 *	which are randomly selected in a unit object of variable dimension. 	
 *
 *	Usage: compile with "make randmst"
 *	Execute as "./randmst flag numpoints numtrials dimension" (all integers)
 *
 *	Output is "average numpoints numtrials dimension" where average is the mean
 *	weight of the MST found over numtrials in a graph with numpoints vertices
 *	where locations are drawn from a unit object of dimension = dimension.
 *
 */

#import <stdio.h>
#import <stdlib.h>
#import <time.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

// a set has a parent and a rank and an 'included' tag for testing 
typedef struct set{
	struct set* parent;
	int rank;
	bool included;
} set;

typedef struct edge{
	int v1;
	int v2;
	float weight;
} edge;

// generates a randomized adjacency matrix over n points
edge* generate(int n, int dimension);
// finds euclidean distance between two points input as lists of coordinates
float euclideanDist(float* p1, float* p2, int dimension);
// makes a new set
set* makeSet(void);
// finds the parent of a set
set* find(set* i);
// links two sets
set* link(set* s1, set* s2);
// finds union of two sets
bool U(set* s1, set* s2);
// compare 2 edges (for sorting)
int compare(const void* a, const void* b);

// implements Kruskal's Algorithm
int main(int argc, char* argv[]){

	// parse args
	if(argc != 5){
		fprintf(stderr, "error - incorrect arguments\n");
		printf("usage: ./randmst flag numpoints numtrials dimension\n");
		return(-1.);
	}
	int flag, numpoints, numtrials, dimension;
	if (sscanf (argv[1], "%i", &flag) != 1) {
	    fprintf(stderr, "error - not an integer\n");
	    printf("usage: ./randmst flag numpoints numtrials dimension\n");
	    return(-1.);
	}
	if (sscanf (argv[2], "%i", &numpoints) != 1) {
	    fprintf(stderr, "error - not an integer\n");
	    printf("usage: ./randmst flag numpoints numtrials dimension\n");
	    return(-1.);
	}
	if (sscanf (argv[3], "%i", &numtrials) != 1) {
	    fprintf(stderr, "error - not an integer\n");
	    printf("usage: ./randmst flag numpoints numtrials dimension\n");
	    return(-1.);
	}
	if (sscanf (argv[4], "%i", &dimension) != 1) {
	    fprintf(stderr, "error - not an integer\n");
	    printf("usage: ./randmst flag numpoints numtrials dimension\n");
	    return(-1.);
	}
	clock_t begin = clock();
	int k;
	float total = 0;
	// for numtrials iterations...
	for(k = 0; k < numtrials; k++){
		printf("Running trial # %d\r",k + 1);
		fflush(stdout);

		edge* edgeList = generate(numpoints, dimension);
		int nedges = numpoints*(numpoints-1)/2;
		qsort(edgeList, nedges, sizeof(edge), compare);
		
		int nIncluded = 0, i, j;
		float minWeight;
		set* a;
		set* b;
		int index = 0; 
		// initialize every vertex as a singleton set
		set** sets = (set**) malloc(numpoints * sizeof(set*));
		for(i = 0; i < numpoints; i++){
			set* s = makeSet();
			if (s == NULL)
				return(-1);
			sets[i] = s;
		}
		// we need to find numpoints - 1 best edges
		edge ei;
		while (nIncluded < numpoints - 1){
			// printf("index = %d\n", index);
			ei = edgeList[index];
			a = sets[ei.v1];
			b = sets[ei.v2];
			// if a and c are already linked, find will find the same root
			a = find(a);
			b = find(b);
			// and union will return false
			if(U(a,b)){
				// otherwise we take this edge in the MST
				nIncluded++;
				total += ei.weight;
				sets[ei.v1]->included = sets[ei.v2]->included = true;
			}
			index ++;
		}

		for(i = 0; i < numpoints; i++){
			// this was for testing
			assert(sets[i]->included == true);
			free(sets[i]);
		}
		free(sets);
		free(edgeList);
	}
	total /= numtrials;
	clock_t end = clock();
	printf("Time spent: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("%f %d %d %d\n", total, numpoints, numtrials, dimension);
	return 0;
}

// randomly generates a new adjacency matrix
edge* generate(int n, int dimension){
	// seed the RNG with current time -- always new
	srand(time(NULL));
	int i, j;
	// allocate an n * dimensions matrix to assign a location to every vertex
	float** locations = (float**)malloc(n * sizeof(float*));
	if (locations == NULL)
		return NULL;
	float* l;
    for (i = 0; i < n; i++){
        l = (float*)malloc(dimension * sizeof(float));
    	if (l == NULL)
			return NULL;
		locations[i] = l;
	}
	// randomly assign coordinates
	for(i = 0; i < n; i++){
		for(j = 0; j < dimension; j++){
			locations[i][j] = (float) rand() / (float) RAND_MAX;
		}
	}
	// turn the coordinates into an adjacency list with distance values 
	int nedges = n*(n-1)/2;
	edge* edgeList = (edge*) malloc(nedges * sizeof(edge));
	edge* e;
	for(i = 0; i < n; i++){
		for(j = 0; j < i; j++){
			e = &edgeList[i*(i-1)/2+j]; // triangular numbers plus offset
			e->v1 = i;
			e->v2 = j;
			e->weight = euclideanDist(locations[i], locations[j], dimension);
		}
	}
	// we don't need the locations anymore!
	for(i = 0; i < n; i++)
		free(locations[i]);
	free(locations);
	return (edgeList);
}

// calculates distance in 'dimension'-dimensional space via Pythagorean theorem
float euclideanDist(float* p1, float* p2, int dimension){
	float sum = 0.;
	int i;
	for(i = 0; i < dimension; i++)
		sum += pow((p1[i] - p2[i]), 2);
	return(sqrt(sum));
}

// allocates and initializes a singleton set
set* makeSet(){
	set* s = (set *)malloc(sizeof(set));
	if(s == NULL)
		return NULL;
	s->parent = s;
	s->rank = 0;
	s->included = false;
	return s;
}

// find implements path compression
set* find(set* i){
	if (i->parent != i)
		i->parent = find(i->parent);
	return i->parent;
}

// link implements union by rank
set* link(set* s1, set* s2){
	if(s1->rank == s2->rank){
		s1->parent = s2;
		s1->rank++;
		return s2;
	}
	else if (s1->rank < s2->rank){
		s1->parent = s2;
		return s2;
	}
	else{
		s2->parent = s1;
		return s1;
	}
}

// calls find on the top of set tree
bool U(set* s1, set* s2){
	if(s1 != s2){
		link(find(s1), find(s2));
		return true;
	}
	return false;
}

int compare(const void* a, const void* b){
	edge* e1 = (edge*) a;
	edge* e2 = (edge*) b;
	if(e1->weight > e2->weight)
		return 1;
	if(e1->weight < e2->weight)
		return -1;
	return 0;
}