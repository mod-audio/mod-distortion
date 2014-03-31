#include <cmath>

void Filter1(float *u, float *y, int N, float SampleRate, float *U_1, float *Y_1 );
void Filter2(float *u, float *y, int N, float SampleRate, float *U_1, float *Y_1, float *U_2, float *Y_2, float *U_3, float *Y_3, float *U_4, float *Y_4 );
void FilterGain(float *u, float *y, int N, float Dist, float SampleRate, float *U_1, float *Y_1, float *U_2, float *Y_2 );
void DS1_Clip_Tone(float *u, float *y, float *v1, float *v2, float *v3,  int N, float T, float *U_1, float *Y_1, float *V1_1, float *V2_1, float *V3_1, float t, float vol);
void Filter1_48000(float *u, float *y, int N, float *U_1, float *Y_1 );
void Filter2_48000(float *u, float *y, int N, float *U_1, float *Y_1, float *U_2, float *Y_2, float *U_3, float *Y_3, float *U_4, float *Y_4 );
void FilterGain_48000(float *u, float *y, int N, float Dist, float *U_1, float *Y_1, float *U_2, float *Y_2 );
void DS1_Clip_Tone_48000(float *u, float *y, float *v1, float *v2, float *v3,  int N, float *U_1, float *Y_1, float *V1_1, float *V2_1, float *V3_1, float t, float vol);
