/* Converts input x into DELTAT*linenr x/PHYFIELD */

#include <stdio.h>

#define DELTAT ((double) 0.125)
#define PHYFIELD ((double) 0.003989422804014327)

int main()
{
	int nline = 0;
	while (!feof(stdin)) {
		double field = 0;
		scanf("%lf", &field);
		printf("%f %f\n", nline * DELTAT, field / PHYFIELD);
		nline++;
	}
	return 0;
}
