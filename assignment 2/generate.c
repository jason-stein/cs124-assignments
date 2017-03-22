/*
 *
 *	generate.c: generates randomized matrices for testing
 *	Usage: ./generate dimension filename, where dimension is an integer that
 *	defines the dimension of the two matrices (generate creates two dimension *
 *	dimension matrices), and filename is the name of the file to write to.
 *
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]){
	int n, i, m;
	float p;
	char* filename;
	if (sscanf (argv[1], "%i", &n) != 1) {
	    return(-1);
	}
	FILE* outfile = fopen(argv[2], "w");
	n = 2 * n * n;
	for(i = 0; i < n; i++){
		p = (float) rand() / RAND_MAX;
		if(p < 0.33)
			m = 0;
		else if(p < 0.67)
			m = 1;
		else
			m = 2;
		fprintf(outfile, "%i\n", m);
	}
}