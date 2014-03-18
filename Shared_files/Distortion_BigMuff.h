#include <cmath>

void Filter1(double *u, double *y, int N, double T, double *U_1, double *Y_1 );
void Filter2(double *u, double *y, int N, double T, double *U_1, double *Y_1, double *U_2, double *Y_2, double *U_3, double *Y_3 );
void Filter3(double *u, double *y, int N, double T, double *U_1, double *Y_1, double *U_2, double *Y_2, double x, double x_1 );
void Clip(double *u, double *y, int N, double T, double *U_1, double *Y_1 );
void Filter4(double *u, double *y, int N, double T, double *U_1, double *Y_1, double *U_2, double *Y_2, double *U_3, double *Y_3, double tone, double vol );
