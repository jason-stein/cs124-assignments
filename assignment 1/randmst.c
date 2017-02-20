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
	// bool included;
} set;

// generates a randomized adjacency matrix over n points
float** generate(int n, int dimension);
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
		// make a new graph
		float** adjMat = generate(numpoints, dimension);
		if (adjMat == NULL)
			return(-1);
		int nIncluded = 0, i = 0, j = 0;
		float minWeight;
		set* a; 
		int b; 
		set* c; 
		int d;
		// initialize every vertex as a singleton set
		set** sets = (set**) malloc(numpoints * sizeof(set*));
		for(i = 0; i < numpoints; i++){
			set* s = makeSet();
			if (s == NULL)
				return(-1);
			sets[i] = s;
		}
		// we need to find numpoints - 1 best edges
		while (nIncluded < numpoints - 1){
			// find the min-weight edge
			for(i = 0, minWeight = 2.0; i < numpoints; i++){
				// adjMat is symmetric so we only need to iterate half (up to i)
				for(j = 0; j < i; j++){
					if(adjMat[i][j] < minWeight){
						minWeight = adjMat[i][j];
						a = sets[i];
						b = i;
						c = sets[j];
						d = j;
					}
				}
			}
			// if a and c are already linked, find will find the same root
			a = find(a);
			c = find(c);
			// and union will return false
			if(U(a,c)){
				// otherwise we take this edge in the MST
				nIncluded++;
				total += minWeight;
				// sets[b]->included = sets[d]->included = true;
			}
			// printf("Taking edge (%d,%d) for cost of %f\n",b,d,adjMat[b][d]);
			
			// throw out this edge by maxing its cost
			adjMat[b][d] = adjMat[d][b] = 2.0;
		}

		for(i = 0; i < numpoints; i++){
			// this was for testing
			// assert(sets[i]->included == true);
			free(sets[i]);
			free(adjMat[i]);
		}
		free(sets);
		free(adjMat);
	}
	total /= numtrials;
	clock_t end = clock();
	printf("Time spent: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("%f %d %d %d\n", total, numpoints, numtrials, dimension);
	return 0;
}

// randomly generates a new adjacency matrix
float** generate(int n, int dimension){
	// seed the RNG with current time -- always new
	srand(time(NULL));
	int i, j;
	// allocate an n * dimensions matrix to assign a location to every vertex
	float** locations = (float**)malloc(n * sizeof(float*));
	if (locations == NULL)
		return NULL;
    for (i = 0; i < n; i++){
        float* l = (float*)malloc(dimension * sizeof(float));
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
	// turn the coordinates into an n * n adjacency matrix with distance values 
	float** adjMat = (float **)malloc(n * sizeof(float *));
	if (adjMat == NULL)
		return NULL;
    for (i = 0; i < n; i++){
    	float* a = (float *)malloc(n * sizeof(float));
    	if (a == NULL)
    		return NULL;
    	adjMat[i] = a;
    }
	for(i = 0; i < n; i++){
		for(j = 0; j < i; j++)
			adjMat[i][j] = adjMat[j][i] = 
				euclideanDist(locations[i], locations[j], dimension);
	}
	// we don't need the locations anymore!
	for(i = 0; i < n; i++)
		free(locations[i]);
	free(locations);
	return (adjMat);
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
	// s->included = false;
	return s;
}

// find implements path compression
set* find(set* i){
	// printf("%p, %p\n", i, i->parent);
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
// nice
