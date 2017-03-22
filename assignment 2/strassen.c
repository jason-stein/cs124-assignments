/*
 *
 *	strassen.c: implements standard (naive) matrix multiplication, as well as 
 *	Strassen's algorithm. Usage is: ./strassen flag dimension infile, 
 *	where flag is an integer, dimension is the size of the matrices (two 
 * 	square matrices of size dimension * dimension), and infile is the source 
 *	of the matrices: a text file with 2 * dimension ^ 2 integers, one per line.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define n0 2048

int** squareMatMult(int dim, int** M1, int** M2);
int** strassen(int dim, int** M1, int** M2);
int** squareMatAdd(int flag, int i0, int j0, int dim, int** M1, int** M2);
int** subMat(int** M, int i0, int j0, int i1, int j1);
int freeMat(int** M, int dim);
int printMat(int dim, int** M);
int printDiag(int dim, int** M);

int main(int argc, char* argv[]){
	// parse args
	if(argc != 4){
		printf("usage: ./strassen flag dimension infile\n");
		return(-1.);
	}
	int dim, flag;
	if (sscanf (argv[1], "%i", &flag) != 1) {
	    printf("usage: ./strassen flag dimension infile\n");
	    return(-1);
	}
	if (sscanf (argv[2], "%i", &dim) != 1) {
	    printf("usage: ./strassen flag dimension infile\n");
	    return(-1);
	}
	
	// open the infile
	FILE* infile = fopen(argv[3], "r");
	if (infile == NULL){
		printf("Couldn't open file.");
		return(-1);
	}
	
	// allocate space for input matrices
	// M1 * M2 = M3
	int i, j;
	int** M1 = (int**) malloc(dim * sizeof(int*));
	int** M2 = (int**) malloc(dim * sizeof(int*));
	for(i = 0; i < dim; i++){
		M1[i] = (int*) malloc(dim * sizeof(int));
		M2[i] = (int*) malloc(dim * sizeof(int));
	}
	
	// read infile in two parts - each dim^2 integers
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

	// execute multiplications and time them
	clock_t start, stop;

	start = clock();
	int** M3 = squareMatMult(dim, M1, M2);
	stop = clock();
	printMat(dim, M3);
	// printDiag(dim, M3);
	printf("Time spent: %f seconds\n", (double)(stop - start) / CLOCKS_PER_SEC);
	
	start = clock();
	int** M4 = strassen(dim, M1, M2);
	stop = clock();
	printMat(dim, M4);
	printf("Time spent: %f seconds\n", (double)(stop - start) / CLOCKS_PER_SEC);
	
	// memory freeing
	freeMat(M1, dim);
	freeMat(M2, dim);
	freeMat(M3, dim);
	freeMat(M4, dim);
	return 0;
}

// multiplies two square matrices--naive method
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

int** strassen(int dim, int** M1, int** M2){
	if(dim <= n0)
		return squareMatMult(dim, M1, M2);
	return M1;
}

// adds (or subtracts) two square matrices
// if flag == 0, adds M1 to M2, otherwise subtracts M1 from M2
// M2 can be larger than M1, addition / subtraction starts from (i0, j0) in M2
int** squareMatAdd(int flag, int i0, int j0, int dim, int** M1, int** M2){
	int i, j;
	for(i = i0; i < dim + i0; i++){
		for(j = j0; j < dim + j0; j++){
			// indexing [i][j] is cache-optimal
			if(flag == 0)
				M2[i + i0][j + j0] += M1[i][j];
			else
				M2[i + i0][j + j0] -= M1[i][j];
		}
	}
	return M2;
}

int** subMat(int** M, int i0, int j0, int i1, int j1){
	int i, j, m = i1 - i0, n = j1 - j0;
	int** sub = (int**) malloc(m * sizeof(int*));
	for(i = 0; i < m; i++){
		sub[i] = (int*) malloc(n * sizeof(int));
		for(j = 0; j < n; j++)
			sub[i][j] = M[i + i0][j + j0];
	}
	return sub;
}

int freeMat(int** M, int dim){
	int i;
	for(i = 0; i < dim; i++)
		free(M[i]);
	free(M);
	return 0;
}

// prints an entire matrix
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

// prints the diagonal elements of a matrix
int printDiag(int dim, int** M){
	int i;
	for (i = 0; i < dim; i++){
		printf("%d\n", M[i][i]);
	}
	return 0;
}