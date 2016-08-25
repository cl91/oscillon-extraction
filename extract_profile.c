#include <stdio.h>

enum { N = 128 };	      /* # of lattice points per dimension */
typedef val_t double;	      /* source data type (float or double) */

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

/* Isola oscillon from mask file
 * phi --- field value file from simulation
 * idx --- index of oscillon to extract (from stats file)
 * mask --- mask array from 
 */
struct field clip_oscillon(val_t *phi, int idx, 
