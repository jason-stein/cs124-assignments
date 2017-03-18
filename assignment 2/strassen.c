#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int** squareMatMult(int dim, int** M1, int** M2);
int printMat(int dim, int** M);
int printDiag(int dim, int** M);

int main(int argc, char* argv[]){
	if(argc != 4){
		printf("usage: ./strassen flag dimension infile\n");
		return(-1.);
	}
	int dim, flag;
	char* fileName;
	if (sscanf (argv[1], "%i", &flag) != 1) {
	    printf("usage: ./strassen flag dimension infile\n");
	    return(-1);
	}
	if (sscanf (argv[2], "%i", &dim) != 1) {
	    printf("usage: ./strassen flag dimension infile\n");
	    return(-1);
	}
	FILE* infile = fopen(argv[3], "r");
	if (infile == NULL){
		printf("Couldn't open file.");
		return(-1);
	}
	// M1 * M2 = M3
	int i, j;
	int** M1 = (int**) malloc(dim * sizeof(int*));
	int** M2 = (int**) malloc(dim * sizeof(int*));
	for(i = 0; i < dim; i++){
		M1[i] = (int*) malloc(dim * sizeof(int));
		M2[i] = (int*) malloc(dim * sizeof(int));
	}
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			fscanf(infile, "%d", &M1[i][j]);
		}
	}
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			fscanf(infile, "%d", &M2[i][j]);
		}
	}
	clock_t start = clock();
	int** M3 = squareMatMult(dim, M1, M2);
	clock_t stop = clock();
	printMat(dim, M3);
	printDiag(dim, M3);
	printf("Time spent: %f seconds\n", (double)(stop - start) / CLOCKS_PER_SEC);
	return 0;
}

int** squareMatMult(int dim, int** M1, int** M2){
	int i, j, k, sum;
	int** M3 = malloc(dim * sizeof(int*));
	for(i = 0; i < dim; i++)
		M3[i] = malloc(dim * sizeof(int));
	for (i = 0; i < dim; i++) {
    	for (j = 0; j < dim; j++) {
    		sum = 0;
			for (k = 0; k < dim; k++) {
				sum = sum + M1[i][k] * M2[k][j];
			}
			M3[i][j] = sum;
		}
    }
    return M3;
}

int printMat(int dim, int** M){
	int i, j;
	for (i = 0; i < dim; i++){
		for (j = 0; j < dim; j++){
			printf("%d ", M[i][j]);
		}
		printf("\n");
	}
	return 0;
}

int printDiag(int dim, int** M){
	int i;
	for (i = 0; i < dim; i++){
		printf("%d\n", M[i][i]);
	}
	return 0;
}