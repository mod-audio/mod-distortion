#include <cmath>

void Filter1(double *u, double *y, int N, double T, double *U_1, double *Y_1 );
void Filter2(double *u, double *y, int N, double T, double *U_1, double *Y_1, double *U_2, double *Y_2, double *U_3, double *Y_3, double *U_4, double *Y_4 );
void FilterGain(double *u, double *y, int N, double Dist, double T, double *U_1, double *Y_1, double *U_2, double *Y_2 );
void DS1_Clip_Tone(double *u, double *y, double *v1, double *v2, double *v3,  int N, double T, double *U_1, double *Y_1, double *V1_1, double *V2_1, double *V3_1, double t, double vol);
void Filter1_48000(double *u, double *y, int N, double *U_1, double *Y_1 );
void Filter2_48000(double *u, double *y, int N, double *U_1, double *Y_1, double *U_2, double *Y_2, double *U_3, double *Y_3, double *U_4, double *Y_4 );
void FilterGain_48000(double *u, double *y, int N, double Dist, double *U_1, double *Y_1, double *U_2, double *Y_2 );
void DS1_Clip_Tone_48000(double *u, double *y, double *v1, double *v2, double *v3,  int N, double *U_1, double *Y_1, double *V1_1, double *V2_1, double *V3_1, double t, double vol);
