#include <cmath>

void Filter1(double *u, double *y, int N, double T, double *U_1, double *Y_1 );
void Clip(double *u, double *y, int N);
void Filter2(double *u, double *y, int N, double vol, double T, double *U_1, double *Y_1 );
