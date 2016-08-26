#include <stdio.h>
#include <stdlib.h>

enum { N = 128 };	      /* # of lattice points per dimension */
typedef double val_t;	      /* source data type (float or double) */

struct point {
	int x, y, z;
};

static inline int to_linear_index(struct point p) {
	return p.z + N * (p.y + N * p.x);
}

static inline struct point from_linear_index(int idx) {
	struct point p;
	p.x = idx / (N * N);
	p.y = (idx % (N * N)) / N;
	p.z = idx % N;
	return p;
}

struct field {
	struct {
		struct point pt;
		val_t val;
	} *array;
	int npoints;
};

/* Isolate oscillon profile from mask file
 * phi --- field value file from simulation
 * idx --- index of oscillon to extract (from stats file)
 * mask --- mask array from extraction tool
 */
struct field clip_oscillon(val_t *phi, int idx, int *mask)
{
	struct field f = {0};
	for (int i = 0; i < N * N * N; i++) {
		if (mask[i] == idx) {
			f.npoints++;
		}
	}
	f.array = malloc(f.npoints * sizeof(*f.array));
	for (int i = 0, c = 0; i < N * N * N; i++) {
		if (mask[i] == idx) {
			f.array[c].pt = from_linear_index(i);
			f.array[c].val = phi[i];
			c++;
		}
	}
	return f;
}

void print_field(struct field f)
{
	for (int i = 0; i < f.npoints; i++) {
		printf("%d %d %d %f\n", f.array[i].pt.x,
		       f.array[i].pt.y, f.array[i].pt.y,
		       f.array[i].val);
	}
}

void *load_file(const char *fname, int elem_size)
{
	void *array = malloc(elem_size * N * N * N);
	if (array == NULL) {
		fprintf(stderr, "malloc() failed.");
		exit(1);
	}
	FILE *f = fopen(fname, "rb");
	if (fread(array, elem_size, N * N * N, f) < N * N * N) {
		fprintf(stderr, "fread failed: %s\n", fname);
		exit(1);
	}
	return array;
}

val_t *load_phi(const char *fname)
{
	return (val_t *) load_file(fname, sizeof(val_t));
}

int *load_mask(const char *fname)
{
	return (int *) load_file(fname, sizeof(int));
}

int main(int argc, char **argv)
{
	if (argc != 4) {
		fprintf(stderr, "usage: extract_profile PHI INDEX MASK\n");
		return 1;
	}

	const char *phi_file = argv[1];
	int index = atoi(argv[2]);
	const char *mask_file = argv[3];
	val_t *phi = load_phi(phi_file);
	int *mask = load_mask(mask_file);
	print_field(clip_oscillon(phi, index, mask));

	return 0;
}
