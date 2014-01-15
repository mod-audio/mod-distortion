#include <cmath>

void LPF1(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 );
void HPF1(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 );
void LPF2(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 );
void HPF2(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 );
void LPF3(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 );
void HPF3(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 );
void BPF1(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3, double *Y_4, double *U_4, double *Y_5, double *U_5, double *Y_6, double *U_6 );
void BPF2(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3, double *Y_4, double *U_4, double *Y_5, double *U_5, double *Y_6, double *U_6 );
void BPF3(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3, double *Y_4, double *U_4, double *Y_5, double *U_5, double *Y_6, double *U_6 );

