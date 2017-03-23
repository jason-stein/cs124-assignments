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
#include <assert.h>

// #define n0 128

int** squareMatMult(int dim, int** M1, int** M2);
int** strassen(int dim, int n0, int** A, int** B);
int** matAdd(int flag, int dim1, int dim2, int i1, int j1, int i2, int j2, int** M1, int** M2);
int zeroMat(int dim, int** M);
int freeMat(int dim, int** M);
int printMat(int dim, int** M);
int printDiag(int dim, int** M);
void assertEqual(int dim, int** M1, int** M2);

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
		printf("Couldn't open file.\n");
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

	// start = clock();
	// int** M3 = squareMatMult(dim, M1, M2);
	// stop = clock();
	// printf("Time spent: %f seconds\n", (double)(stop - start) / CLOCKS_PER_SEC);
	
	// experimenting with values of n0
	// for(i = 0; i < 10; i ++){
	// 	int n0 = 20 * (i + 4);
	// 	start = clock();
	// 	int** M4 = strassen(dim, n0, M1, M2);
	// 	stop = clock();
	// 	// assertEqual(dim, M3, M4);
	// 	printf("n0 = %d, Time spent: %f seconds\n", n0, (double)(stop - start) / CLOCKS_PER_SEC);
	// 	freeMat(dim, M4);
	// }

	
	// for submission:
	int** M3 = strassen(dim, 110, M1, M2);
	printDiag(dim, M3);
	
	// memory freeing
	freeMat(dim, M1);
	freeMat(dim, M2);
	freeMat(dim, M3);
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

int** strassen(int dim, int n0, int** A, int** B){
	if(dim <= n0 || dim == 1)
		return squareMatMult(dim, A, B);
	int i, div1, div2;
	// div1 is larger subdivision, div2 is smaller (equal if even)
	div1 = dim / 2;
	div2 = dim / 2;
	if (dim % 2 == 1)
		div1 += 1;
	// allocate two temp div1 * div1 matrices and one final dim * dim matrix
	int** tmp1 = (int**) malloc(div1 * sizeof(int*));
	int** tmp2 = (int**) malloc(div1 * sizeof(int*));
	int** C = (int**) malloc(dim * sizeof(int*));
	for(i = 0; i < div1; i++){
		tmp1[i] = calloc(div1, sizeof(int));
		tmp2[i] = calloc(div1, sizeof(int));
	}
	for(i = 0; i < dim; i++)
		C[i] = calloc(dim, sizeof(int));

	// now for the fun part
	// using matAdd, place the right values in tmp1 and tmp2 and recursively
	// multiply them. A11 is size dim1 * dim1, A12 is dim1 * dim2, etc...
	// A11 starts at (0, 0) in A, A12 starts at (0, dim1), etc...
	// zero padding is automatic by adding smaller matrix to larger 
	matAdd(0, div1, div1, 0, 0, 0, 0, A, tmp1);			// tmp1 = A11
	matAdd(0, div2, div2, div1, div1, 0, 0, A, tmp1);	// tmp1 = (A11 + A22)
	matAdd(0, div1, div1, 0, 0, 0, 0, B, tmp2);			// tmp2 = B11
	matAdd(0, div2, div2, div1, div1, 0, 0, B, tmp2);	// tmp2 = (B11 + B22)
	int** M1 = strassen(div1, n0, tmp1, tmp2); 			// M1 = (A11 + A22)(B11 + B22)
	zeroMat(div1, tmp1);
	zeroMat(div1, tmp2);
	matAdd(0, div2, div1, div1, 0, 0, 0, A, tmp1);		// tmp1 = A21
	matAdd(0, div2, div2, div1, div1, 0, 0, A, tmp1);	// tmp1 = (A21 + A22)
	matAdd(0, div1, div1, 0, 0, 0, 0, B, tmp2);			// tmp2 = B11
	int** M2 = strassen(div1, n0, tmp1, tmp2);			// M2 = (A21 + A22)B11
	zeroMat(div1, tmp1);
	zeroMat(div1, tmp2);
	matAdd(0, div1, div1, 0, 0, 0, 0, A, tmp1);			// tmp1 = A11
	matAdd(0, div1, div2, 0, div1, 0, 0, B, tmp2);		// tmp2 = B12
	matAdd(1, div2, div2, div1, div1, 0, 0, B, tmp2);	// tmp2 = (B12 - B22)
	int** M3 = strassen(div1, n0, tmp1, tmp2);			// M3 = A11(B12 - B22)
	zeroMat(div1, tmp1);
	zeroMat(div1, tmp2);
	matAdd(0, div2, div2, div1, div1, 0, 0, A, tmp1);	// tmp1 = A22
	matAdd(0, div2, div1, div1, 0, 0, 0, B, tmp2);		// tmp2 = B21
	matAdd(1, div1, div1, 0, 0, 0, 0, B, tmp2);			// tmp2 = (B21 - B11)
	int** M4 = strassen(div1, n0, tmp1, tmp2);			// M4 = A22(B21 - B11)
	zeroMat(div1, tmp1);
	zeroMat(div1, tmp2);
	matAdd(0, div1, div1, 0, 0, 0, 0, A, tmp1);			// tmp1 = A11
	matAdd(0, div1, div2, 0, div1, 0, 0, A, tmp1);		// tmp1 = (A11 + A12)
	matAdd(0, div2, div2, div1, div1, 0, 0, B, tmp2);	// tmp2 = B22
	int** M5 = strassen(div1, n0, tmp1, tmp2);			// M5 = (A11 + A12)B22
	zeroMat(div1, tmp1);
	zeroMat(div1, tmp2);
	matAdd(0, div2, div1, div1, 0, 0, 0, A, tmp1);		// tmp1 = A21
	matAdd(1, div1, div1, 0, 0, 0, 0, A, tmp1);			// tmp1 = (A21 - A11)
	matAdd(0, div1, div1, 0, 0, 0, 0, B, tmp2);			// tmp2 = B11
	matAdd(0, div1, div2, 0, div1, 0, 0, B, tmp2);		// tmp2 = (B11 + B12)
	int** M6 = strassen(div1, n0, tmp1, tmp2);			// M6 = (A21 - A11)(B11 + B12)
	zeroMat(div1, tmp1);
	zeroMat(div1, tmp2);
	matAdd(0, div1, div2, 0, div1, 0, 0, A, tmp1);		// tmp1 = A12
	matAdd(1, div2, div2, div1, div1, 0, 0, A, tmp1);	// tmp1 = (A12 - A22)
	matAdd(0, div2, div1, div1, 0, 0, 0, B, tmp2);		// tmp2 = B21
	matAdd(0, div2, div2, div1, div1, 0, 0, B, tmp2);	// tmp2 = (B21 + B22)
	int** M7 = strassen(div1, n0, tmp1, tmp2);			// M7 = (A12 - A22)(B21 + B22)
	freeMat(div1, tmp1);
	freeMat(div1, tmp2);

	// now put M1-7 into C where necessary. 
	// trim extra zeroes by taking smaller chunks of matrices (e.g. C22 is dim2 * dim2)
	matAdd(0, div1, div1, 0, 0, 0, 0, M1, C);			// C11 = M1
	matAdd(0, div1, div1, 0, 0, 0, 0, M4, C);			// C11 = M1 + M4
	matAdd(1, div1, div1, 0, 0, 0, 0, M5, C);			// C11 = M1 + M4 - M5
	matAdd(0, div1, div1, 0, 0, 0, 0, M7, C);			// C11 = M1 + M4 - M5 + M7
	matAdd(0, div1, div2, 0, 0, 0, div1, M3, C);		// C12 = M3
	matAdd(0, div1, div2, 0, 0, 0, div1, M5, C);		// C12 = M3 + M5
	matAdd(0, div2, div1, 0, 0, div1, 0, M2, C);		// C21 = M2
	matAdd(0, div2, div1, 0, 0, div1, 0, M4, C);		// C21 = M2 + M4
	matAdd(0, div2, div2, 0, 0, div1, div1, M1, C);		// C22 = M1
	matAdd(1, div2, div2, 0, 0, div1, div1, M2, C);		// C22 = M1 - M2
	matAdd(0, div2, div2, 0, 0, div1, div1, M3, C);		// C22 = M1 - M2 + M3
	matAdd(0, div2, div2, 0, 0, div1, div1, M6, C);		// C22 = M1 - M2 + M3 + M6

	freeMat(div1, M1);
	freeMat(div1, M2);
	freeMat(div1, M3);
	freeMat(div1, M4);
	freeMat(div1, M5);
	freeMat(div1, M6);
	freeMat(div1, M7);

	return C;
}

/* 
 * adds (or subtracts) two matrices
 * if flag == 0, adds M1 to M2, otherwise subtracts M1 from M2
 * adds / subtracts subsections of size dim1 * dim2
 * starting at (i1, j1) in M1 and (i2, j2) in M2
 * this allows me to pad with zeros by adding a smaller matrix to a larger calloc'ed matrix,
 * and remove zeroes by setting dim1 and dim2 to only take the appropriate dimensions
 * it also prevents me from needing to allocate submatrices A11..C22 and instead manipulate in place
 * this implementation is great and i am proud of it
 */
int** matAdd(int flag, int dim1, int dim2, int i1, int j1, int i2, int j2, int** M1, int** M2){
	int i, j;
	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			if(flag == 0)
				M2[i + i2][j + j2] += M1[i + i1][j + j1];
			else
				M2[i + i2][j + j2] -= M1[i + i1][j + j1];
		}
	}
	return M2;
}

// sets all entries of matrix to 0
int zeroMat(int dim, int** M){
	int i, j;
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			M[i][j] = 0;
		}
	}
	return 0;
}

// frees all pointers in a matrix
int freeMat(int dim, int** M){
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

// asserts that all entries are equal
void assertEqual(int dim, int** M1, int** M2){
	int i, j;
	for(i = 0; i < dim; i++){
		for (int j = 0; j < dim; j++)
		{
			assert(M1[i][j] == M2[i][j]);
		}
	}
}
// nice 