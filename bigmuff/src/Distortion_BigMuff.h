#include <cmath>

void Filter1(float *u, float *y, int N, double T, float *U_1, float *Y_1 );
void Filter2(float *u, float *y, int N, double T, float *U_1, float *Y_1, float *U_2, float *Y_2, float*U_3, float *Y_3 );
void Filter3(float *u, float *y, int N, double T, float *U_1, float *Y_1, float *U_2, float *Y_2, double x, double x_1 );
void Clip(float *u, float *y, int N, double T, float *U_1, float *Y_1 );
void Filter4(float *u, float *y, int N, double T, float *U_1, float *Y_1, float *U_2, float *Y_2, float *U_3, float *Y_3, double tone, double vol );
