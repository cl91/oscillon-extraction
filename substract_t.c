#include <stdio.h>

int main()
{
	while (!feof(stdin)) {
		double x, y;
		scanf("%lf %lf", &x, &y);
		printf("%f %f\n", x - 390, y);
	}
	return 0;
}
