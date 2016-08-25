#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void sifter(float *s, unsigned int *r, unsigned int NUM, float cutoff);
void mask(unsigned int ***maskArray, unsigned int NUM);
void partition(unsigned int *s, unsigned int ***array3D, unsigned int NUM);
void componentMult(unsigned int ***s, unsigned int ***r, unsigned int ***z,
		   unsigned int NUM);
unsigned int ***tidyZ(unsigned int ***s, int *count, unsigned int NUM);
unsigned int ***tidyY(unsigned int ***s, int *count, unsigned int NUM);
unsigned int ***tidyX(unsigned int ***s, int *count, unsigned int NUM);
void flatten(unsigned int ***s, unsigned int *r, unsigned int NUM);
int compare(const void *a, const void *b);
unsigned int *lower(unsigned int *s, unsigned int *r, unsigned int NUM);

void partitionfloat(float *s, float ***r, unsigned int NUM);
void oscillonstatistic(unsigned int *flatsiftedmask, float *arrayFlat,
		       unsigned int NUM, float aval, unsigned int oscnum,
		       char **argv);

int main(int argc, char *argv[])
{
	/*argv[1]=input file, argv[2]=file to write to, argv[3]=# of lattice points, argv[4]=(4or8)#of bytes per number, argv[5]=cutoff, argv[6]=a, argv[7]= statistics file */

	unsigned int i, j, l, m;

	unsigned int NUM;
	NUM = (unsigned int)atoi(argv[3]);
	int size;
	size = (int)atoi(argv[4]);
	float cutoff;
	cutoff = (float)atof(argv[5]);
	float aval;
	aval = (float)atof(argv[6]);

	float *arrayFlat;	/* 1D array to store field values */
	arrayFlat = (float *)malloc(NUM * NUM * NUM * sizeof(float));
	unsigned int *arrayflatsifted;
	arrayflatsifted =
	    (unsigned int *)malloc(NUM * NUM * NUM * sizeof(unsigned int));

	FILE *fp;
	fp = fopen(argv[1], "rb");

	if (fp == NULL) {
		printf("\n file did not open, try again \n");
		return 0;
	}
	if (size != 4 && size != 8) {
		printf("invalid input for argument 4, terminating");
		return 0;
	}

	double *real64;
	real64 = (double *)malloc(sizeof(double));

	if (size == 8) {
		for (i = 0; i < (NUM * NUM * NUM); ++i) {
			fread(real64, size, 1, fp);
			arrayFlat[i] = (float)*real64;
		}
	}
	free(real64);

	if (size == 4) {
		for (i = 0; i < NUM * NUM * NUM; ++i) {
			fread(&arrayFlat[i], size, 1, fp);
		}

	}

	sifter(arrayFlat, arrayflatsifted, NUM, cutoff);	/*fitler out small field values, yields arrayflatsifted with 0s and 1s */

	unsigned int ***maskArray;
	maskArray = (unsigned int ***)malloc(NUM * sizeof(unsigned int **));
	for (l = 0; l < NUM; l++) {
		maskArray[l] =
		    (unsigned int **)malloc(NUM * sizeof(unsigned int *));
		for (m = 0; m < NUM; m++) {
			maskArray[l][m] =
			    (unsigned int *)malloc(NUM *
						   sizeof(unsigned int *));
		}
	}

	mask(maskArray, NUM);	/*create the 3D mask array with integer labels at each point */

	unsigned int ***array3D;
	array3D = (unsigned int ***)malloc(NUM * sizeof(unsigned int **));
	for (l = 0; l < NUM; l++) {
		array3D[l] =
		    (unsigned int **)malloc(NUM * sizeof(unsigned int *));
		for (m = 0; m < NUM; m++) {
			array3D[l][m] =
			    (unsigned int *)malloc(NUM *
						   sizeof(unsigned int *));
		}
	}

	partition(arrayflatsifted, array3D, NUM);	/*Turn arrayflatsifted into array3D, each point ijk assigned the field value 0 or 1 */
	free(arrayflatsifted);

	unsigned int ***siftedmask;
	siftedmask = (unsigned int ***)malloc(NUM * sizeof(unsigned int **));
	for (l = 0; l < NUM; l++) {
		siftedmask[l] =
		    (unsigned int **)malloc(NUM * sizeof(unsigned int *));
		for (m = 0; m < NUM; m++) {
			siftedmask[l][m] =
			    (unsigned int *)malloc(NUM *
						   sizeof(unsigned int *));
		}
	}

	componentMult(maskArray, array3D, siftedmask, NUM);	/*multiply integer labels onto sifted field values, integers only assigned to points above the CUTOFF */

	for (i = 0; i < NUM; ++i) {
		for (j = 0; j < NUM; ++j) {
			free(maskArray[i][j]);
			free(array3D[i][j]);
		}
	}
	free(maskArray);
	free(array3D);

	int *count;
	count = (int *)malloc(sizeof(int));
	*count = 1;

	while ((*count) > 0) {	/*run through each dimension of array shrinking adjacent integer labels in a fixed point iteration */
		tidyZ(siftedmask, count, NUM);
		tidyY(siftedmask, count, NUM);
		tidyX(siftedmask, count, NUM);
		//              printf("%d\t",*count); 
	}

	printf("tidied up\n");

	unsigned int *interimflat, *flatsiftedmask;	/*flatten the tidied array, then order using quicksort */
	interimflat =
	    (unsigned int *)malloc(NUM * NUM * NUM * sizeof(unsigned int));
	flatsiftedmask =
	    (unsigned int *)malloc(NUM * NUM * NUM * sizeof(unsigned int));

	flatten(siftedmask, interimflat, NUM);
	flatten(siftedmask, flatsiftedmask, NUM);
	qsort(interimflat, NUM * NUM * NUM, sizeof(unsigned int), compare);

	for (i = 0; i < NUM; ++i) {
		for (j = 0; j < NUM; ++j) {
			free(siftedmask[i][j]);
		}
	}
	free(siftedmask);

	flatsiftedmask = lower(interimflat, flatsiftedmask, NUM);	/*replace each large integer label with a small integer ~ 1-100 */
	free(interimflat);

	FILE *fp2;
	fp2 = fopen(argv[2], "wb");
	fwrite(flatsiftedmask, 4, NUM * NUM * NUM, fp2);

	unsigned int max;
	max = flatsiftedmask[0];
	for (i = 1; i < NUM * NUM * NUM; ++i) {	/*find largest value in the reduced mask */
		if (flatsiftedmask[i] > max) {
			max = flatsiftedmask[i];
		}
	}
	printf("wrote the file, %d oscillons found\n", max);

//      FILE* fp3;
//      fp3=fopen(argv[7], "w");

	oscillonstatistic(flatsiftedmask, arrayFlat, NUM, aval, max, argv);

	printf("wrote the stats, done\n");
	return 0;
}

void sifter(float *s, unsigned int *r, unsigned int NUM, float cutoff)
{				/*sifts through field values, replacing values below CUTOFF with 0 and values above CUTOFF with 1 */

	unsigned int i;

	for (i = 0; i < NUM * NUM * NUM; i++) {
		if ((s[i] - cutoff) < 0.0) {
			r[i] = 0;
		} else {
			if ((s[i] - cutoff) > 0) {
				r[i] = 1;
			}
		}
	}
}

void mask(unsigned int ***maskArray, unsigned int NUM)
{				/*create a lattice with a unique integer label for each point */
	unsigned int i, j, k;

	for (i = 0; i < NUM; i++) {
		for (j = 0; j < NUM; j++) {
			for (k = 0; k < NUM; k++) {
				maskArray[i][j][k] =
				    1 + i + (NUM + 1) * j + (NUM + 1) * (NUM +
									 1) * k;
			}
		}
	}
}

void partition(unsigned int *s, unsigned int ***array3D, unsigned int NUM)
{				/*partition takes in a flat array of values, and then splits up the values into a 3D array of field values */

	unsigned int i, j, k, d;

	d = 0;
	for (i = 0; i < NUM; i++) {
		for (j = 0; j < NUM; j++) {
			for (k = 0; k < NUM; k++) {
				array3D[i][j][k] = s[d];
				++d;
			}
		}
	}
}

void componentMult(unsigned int ***s, unsigned int ***r, unsigned int ***z,
		   unsigned int NUM)
{				/*multiply two arrays compoonent by component */

	unsigned int i, j, k;

	for (i = 0; i < NUM; ++i) {
		for (j = 0; j < NUM; ++j) {
			for (k = 0; k < NUM; ++k) {
				z[i][j][k] = (s[i][j][k]) * (r[i][j][k]);
			}
		}
	}
}

unsigned int ***tidyZ(unsigned int ***s, int *count, unsigned int NUM)
{				/*bring the integer clumps into single valued clumps */

	int i, j, k;
	*count = 0;

	for (i = 0; i < (NUM); ++i) {
		for (j = 0; j < (NUM); ++j) {
			if (s[i][j][NUM - 1] > 0 && s[i][j][0] > 0) {
				if (s[i][j][NUM - 1] > s[i][j][0]) {
					s[i][j][NUM - 1] = s[i][j][0];
					(*count)++;
				} else if (s[i][j][NUM - 1] < s[i][j][0]) {
					s[i][j][0] = s[i][j][NUM - 1];
					(*count)++;
				}
			}
		}
	}

	for (i = 0; i < (NUM); ++i) {
		for (j = 0; j < (NUM); ++j) {
			for (k = 0; k < (NUM - 1); ++k) {
				if (s[i][j][k] > 0 && s[i][j][k + 1] > 0) {
					if (s[i][j][k] > s[i][j][k + 1]) {
						s[i][j][k] = s[i][j][k + 1];
						(*count)++;
					} else if (s[i][j][k] < s[i][j][k + 1]) {
						s[i][j][k + 1] = s[i][j][k];
						(*count)++;
					}

				}
			}
		}
	}

	return s;
}

unsigned int ***tidyY(unsigned int ***s, int *count, unsigned int NUM)
{

	int i, j, k;

	for (i = 0; i < (NUM); ++i) {
		for (k = 0; k < (NUM); k++) {
			if (s[i][NUM - 1][k] > 0 && s[i][0][k] > 0) {
				if (s[i][NUM - 1][k] > s[i][0][k]) {
					s[i][NUM - 1][k] = s[i][0][k];
					(*count)++;
				} else if (s[i][NUM - 1][k] < s[i][0][k]) {
					s[i][0][k] = s[i][NUM - 1][k];
					(*count)++;
				}
			}
		}
	}

	for (i = 0; i < (NUM); ++i) {
		for (k = 0; k < (NUM); ++k) {
			for (j = 0; j < (NUM - 1); ++j) {
				if (s[i][j][k] > 0 && s[i][j + 1][k] > 0) {
					if (s[i][j][k] > s[i][j + 1][k]) {
						s[i][j][k] = s[i][j + 1][k];
						(*count)++;
					} else if (s[i][j][k] < s[i][j + 1][k]) {
						s[i][j + 1][k] = s[i][j][k];
						(*count)++;
					}

				}
			}
		}
	}
	return s;
}

unsigned int ***tidyX(unsigned int ***s, int *count, unsigned int NUM)
{

	int i, j, k;

	for (j = 0; j < (NUM); ++j) {
		for (k = 0; k < (NUM); ++k) {
			if (s[NUM - 1][j][k] > 0 && s[0][j][k] > 0) {
				if (s[NUM - 1][j][k] > s[0][j][k]) {
					s[NUM - 1][j][k] = s[0][j][k];
					(*count)++;
				} else if (s[NUM - 1][j][k] < s[0][j][k]) {
					s[0][j][k] = s[NUM - 1][j][k];
					(*count)++;
				}
			}
		}
	}

	for (j = 0; j < (NUM); ++j) {
		for (k = 0; k < (NUM); ++k) {
			for (i = 0; i < (NUM - 1); ++i) {
				if (s[i][j][k] > 0 && s[i + 1][j][k] > 0) {
					if (s[i][j][k] > s[i + 1][j][k]) {
						s[i][j][k] = s[i + 1][j][k];
						(*count)++;
					} else if (s[i][j][k] < s[i + 1][j][k]) {
						s[i + 1][j][k] = s[i][j][k];
						(*count)++;
					}

				}
			}
		}
	}
//      printf("%d\t",*count);
	return s;

}

void flatten(unsigned int ***s, unsigned int *r, unsigned int NUM)
{

	int i, j, k, d;
	d = 0;

	for (i = 0; i < NUM; ++i) {
		for (j = 0; j < NUM; ++j) {
			for (k = 0; k < NUM; ++k) {
				r[d] = s[i][j][k];
				++d;
			}
		}
	}
}

int compare(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}

unsigned int *lower(unsigned int *s, unsigned int *r, unsigned int NUM)
{				/*takes the sorted array of integers, sets all but one copy of each to 0, reads nonzero integers into label array, then lowers the labels */

	unsigned int *reduced, i, j, d;
	reduced =
	    (unsigned int *)malloc(NUM * NUM * NUM * sizeof(unsigned int));
	d = 0;
	for (i = (NUM * NUM * NUM) - 1; i > 0; --i) {
		if (s[i] != 0) {
			if (s[i] == s[i - 1]) {
				s[i] = 0;
			}
		}
	}
	s[0] = 0;

	for (i = 0; i < NUM * NUM * NUM; ++i) {
		reduced[i] = 0;
	}

	for (i = 0; i < NUM * NUM * NUM; ++i) {
		if (s[i] != 0) {
			reduced[d] = s[i];
			++d;
		}
	}

	for (i = 0; i < NUM * NUM * NUM; ++i) {
		if (r[i] != 0) {
			for (j = 0; reduced[j] != 0; ++j) {
				if (r[i] == reduced[j]) {
					r[i] = j + 1;
				}
			}
		}
	}

	return r;
}

void partitionfloat(float *s, float ***r, unsigned int NUM)
{				/*partition takes in a flat array of values, and then splits up the values into a 3D array of field values */

	unsigned int i, j, k, d;

	d = 0;
	for (i = 0; i < NUM; i++) {
		for (j = 0; j < NUM; j++) {
			for (k = 0; k < NUM; k++) {
				r[i][j][k] = s[d];
				++d;
			}
		}
	}
}

void oscillonstatistic(unsigned int *flatsiftedmask, float *arrayFlat,
		       unsigned int NUM, float aval, unsigned int max,
		       char **argv)
{
	unsigned int d, i, j, k, l, m;

	FILE *fp3;
	fp3 = fopen(argv[7], "w");
	fprintf(fp3, "{");

	unsigned int ***reducedmask3D;
	reducedmask3D = (unsigned int ***)malloc(NUM * sizeof(unsigned int **));
	for (l = 0; l < NUM; l++) {
		reducedmask3D[l] =
		    (unsigned int **)malloc(NUM * sizeof(unsigned int *));
		for (m = 0; m < NUM; m++) {
			reducedmask3D[l][m] =
			    (unsigned int *)malloc(NUM *
						   sizeof(unsigned int *));
		}
	}

	partition(flatsiftedmask, reducedmask3D, NUM);

	float ***energygrid;
	energygrid = (float ***)malloc(NUM * sizeof(float **));
	for (l = 0; l < NUM; l++) {
		energygrid[l] = (float **)malloc(NUM * sizeof(float *));
		for (m = 0; m < NUM; m++) {
			energygrid[l][m] =
			    (float *)malloc(NUM * sizeof(float *));
		}
	}

	partitionfloat(arrayFlat, energygrid, NUM);

	for (d = 1; d <= max; ++d) {

		float peak = 0.0;

		for (i = 0; i < NUM; ++i) {	/*find max field value of oscillon d */
			for (j = 0; j < NUM; ++j) {
				for (k = 0; k < NUM; ++k) {
					if (reducedmask3D[i][j][k] == d) {
						if (energygrid[i][j][k] - peak >
						    0.0) {
							peak =
							    energygrid[i][j][k];
						}
					}
				}
			}
		}

		/*now we have peak */

		float energy = 0.0;
		double widthcubed = 0.0;

		for (i = 0; i < NUM; ++i) {
			for (j = 0; j < NUM; ++j) {
				for (k = 0; k < NUM; ++k) {
					if (reducedmask3D[i][j][k] == d) {
						if (energygrid[i][j][k] -
						    (peak / (3.0 * 2.718281)) >
						    0) {
							energy =
							    energy +
							    energygrid[i][j][k];
							widthcubed =
							    widthcubed + 1.0;
						}
					}
				}
			}
		}

		double width = cbrt(widthcubed) * aval;

		if (d == max) {
			fprintf(fp3, "{%f,%f,%f,%lf}}\n", aval, peak, energy,
				width);
		}

		else {
			fprintf(fp3, "{%f,%f,%f,%lf},\n", aval, peak, energy,
				width);
		}

//              printf("{%f, %f, %f, %lf}\n", aval, peak, energy, width);

	}

	fclose(fp3);
	for (i = 0; i < NUM; ++i) {
		for (j = 0; j < NUM; ++j) {
			free(energygrid[i][j]);
			free(reducedmask3D[i][j]);
		}
	}
	free(energygrid);
	free(reducedmask3D);
}
