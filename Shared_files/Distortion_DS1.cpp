#include <cmath>
#include "Distortion_DS1.h"
#include "HyperbolicTables.h"

void Filter1(float *u, float *y, int N, float SampleRate, float *U_1, float *Y_1 )
{
	const float R1 = 1000;
	const float R2 = 470e3;

	const float C1 = 47e-9;
	
	float c = 2*SampleRate;
	
	float y_1 = Y_1[0];
	float u_1 = U_1[0];
	
	/*
	             R2C1s       b1s + b0
	 G(s) = -------------- = --------
	        (R1+R2)C1s + 1   a1s + a0

	         B0 + B1z⁻¹
	 G(z) =  ----------
             A0 + A1z⁻¹
     
     y[k] = (-A1*y[k-1] + B0*u[k] + B1*u[k-1] )/A0
	*/
	
	const float b0 = 0;
	const float b1 = R2*C1;
	const float a0 = 1;
	const float a1 = (R1+R2)*C1;
	
	const float B0 = b0 + b1*c;
	const float B1 = b0 - b1*c;
	const float A0 = a0 + a1*c;
	const float A1 = a0 - a1*c;	
	
	y[0] = (-A1*y_1 + B0*u[0] + B1*u_1)/A0;
		
	for (int i=1; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] + B0*u[i] + B1*u[i-1] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	
}

void Filter2(float *u, float *y, int N, float SampleRate, float *U_1, float *Y_1, float *U_2, float *Y_2, float *U_3, float *Y_3, float *U_4, float *Y_4 )
{
	const double R45 = 90.909090909090907e3;
	const double R6 = 100e3;
	const double R7 = 470e3;
	const double R8 = 10e3;
	const double R9 = 22;
	const double R10 = 100e3;
	
	const double G45 = 1/R45;
	const double G6 = 1/R6;
	const double G7 = 1/R7;
	const double G8 = 1/R8;
	const double G9 = 1/R9;
	const double G10 = 1/R10;

	const double C2 = 470e-9;
	const double C3 = 47e-9;
	const double C4 = 250e-12;
	const double C5 = 68e-9;
	
	double c = 2*SampleRate;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	double y_3 = Y_3[0];
	double u_3 = U_3[0];
	double y_4 = Y_4[0];
	double u_4 = U_4[0];
	
	
	
	double const b4 = C2*C3*C4*C5;
	double const b3 = C2*C3*C5*(G7 - G9);
	double const b2 = 0;
	double const b1 = 0;
	double const b0 = 0;
	
	double const a4 = C2*C3*C4*C5;
	double const a3 = (C2*C3*C4 + C2*C3*C5 + C2*C4*C5 + C3*C4*C5)*G10 + C2*C3*C5*G7 + C2*C4*C5*G6 + C2*C3*C5*G8 + C3*C4*C5*G6 + C2*C4*C5*G8 + C2*C4*C5*G9 + C3*C4*C5*G8 + C3*C4*C5*G9 + C3*C4*C5*G45;
	double const a2 = (C2*C3*G7 + C2*C4*G6 + C2*C3*G8 + C2*C5*G6 + C3*C4*G6 + C2*C4*G8 + C2*C5*G7 + C3*C5*G6 + C2*C4*G9 + C3*C4*G8 + C3*C5*G7 + C3*C4*G9 + C3*C4*G45 + C3*C5*G45 + C4*C5*G45)*G10 + C2*C5*G6*G7 + C2*C5*G6*G8 + C3*C5*G6*G7 + C2*C5*G7*G8 + C3*C5*G6*G8 + C2*C5*G7*G9 + C3*C5*G7*G8 + C3*C5*G7*G9 + C3*C5*G7*G45 + C4*C5*G6*G45 + C3*C5*G8*G45 + C4*C5*G8*G45 + C4*C5*G9*G45;
	double const a1 = (C2*G6*G7 + C2*G6*G8 + C3*G6*G7 + C2*G7*G8 + C3*G6*G8 + C2*G7*G9 + C3*G7*G8 + C3*G7*G9 + C3*G7*G45 + C4*G6*G45 + C3*G8*G45 + C5*G6*G45 + C4*G8*G45 + C5*G7*G45 + C4*G9*G45)*G10 + C5*G6*G7*G45 + C5*G6*G8*G45 + C5*G7*G8*G45 + C5*G7*G9*G45;
	double const a0 = G10*G45*(G6*G7 + G6*G8 + G7*G8 + G7*G9);

	
	const double B0 = b4*pow(c,4) + b3*pow(c,3) + b2*c*c + b1*c + b0;
	const double B1 = -4*b4*pow(c,4) - 2*b3*pow(c,3) + 2*b1*c + 4*b0;
	const double B2 = 6*b4*pow(c,4) - 2*b2*c*c + 6*b0;
	const double B3 = -4*b4*pow(c,4) + 2*b3*pow(c,3) - 2*b1*c + 4*b0;
	const double B4 = b4*pow(c,4) - b3*pow(c,3) + b2*c*c - b1*c + b0;
	const double A0 = a4*pow(c,4) + a3*pow(c,3) + a2*c*c + a1*c + a0;
	const double A1 = -4*a4*pow(c,4) - 2*a3*pow(c,3) + 2*a1*c + 4*a0;
	const double A2 = 6*a4*pow(c,4) - 2*a2*c*c + 6*a0;
	const double A3 = -4*a4*pow(c,4) + 2*a3*pow(c,3) - 2*a1*c + 4*a0;
	const double A4 = a4*pow(c,4) - a3*pow(c,3) + a2*c*c - a1*c + a0;
	
	y[0] = (-A1*y_1 -A2*y_2 - A3*y_3 - A4*y_4 + B0*u[0] + B1*u_1 + B2*u_2 + B3*u_3 + B4*u_4 )/A0;
	y[1] = (-A1*y[0] -A2*y_1 - A3*y_2 - A4*y_3 + B0*u[1] + B1*u[0] + B2*u_1 + B3*u_2 + B4*u_3 )/A0;
	y[2] = (-A1*y[1] -A2*y[0] - A3*y_1 - A4*y_2 + B0*u[2] + B1*u[1] + B2*u[0] + B3*u_1 + B4*u_2 )/A0;
	y[3] = (-A1*y[2] -A2*y[1] - A3*y[0] - A4*y_1 + B0*u[3] + B1*u[2] + B2*u[1] + B3*u[0] + B4*u_1 )/A0;
		
	for (int i=4; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] -A2*y[i-2] - A3*y[i-3] - A4*y[i-4] + B0*u[i] + B1*u[i-1] + B2*u[i-2] + B3*u[i-3] + B4*u[i-4] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	U_4[0] = u[N-4];
	Y_4[0] = y[N-4];
	
}

void FilterGain(float *u, float *y, int N, float Dist, float SampleRate, float *U_1, float *Y_1, float *U_2, float *Y_2 )
{
	
	const float Rd = 100e3;
	const float R13 = 4.7e3;
	const float Cc = 10e-12;
	const float Cz = 470e-9;
	
	//const float Cc = 250e-12;
	//const float Cz = 1000e-9;
	
	float Rt = Dist*Rd;
	float Rb = (1-Dist)*Rd + R13;
	
	float Gt = 1/Rt;
	float Gb = 1/Rb;
	
	float b0;
	float b1;
	float b2;
	float a0;
	float a1;
	float a2;
	
	float B0;
	float B1;
	float B2;
	float A0;
	float A1;
	float A2;
	
	float c = 2*SampleRate;
	
	float y_1 = Y_1[0];
	float u_1 = U_1[0];
	float y_2 = Y_2[0];
	float u_2 = U_2[0];
	
	/*
	  
	          Cz Cc s² + (Cc Gb + Cz Gb + Cz Gt) s + Gb Gt    b2s² + b1s + b0
	 G(s) =  --------------------------------------------- =  ----------------
	          Cz Cc s² + (Cc Gb + Cz Gt) s + Gb Gt            a2s² + a1s + a0
	           
	Rt = DRd
	Rb = (1-D)Rd + R13 
	                   
	         B0 + B1z⁻¹ + B2z⁻²
	 G(z) =  -------------------
             A0 + A1z⁻¹ + A2z⁻²
     
     y[k] = (-A1*y[k-1] -A2*y[k-2] + B0*u[k] + B1*u[k-1] + B2*u[k-2] )/A0
     
	*/
	
	a0 = Gb*Gt;
	a1 = Cc*Gb + Cz*Gt;
	a2 = Cz*Cc;
	b0 = a0;
	b1 = a1 + Cz*Gb;
	b2 = a2;
	
	//Bilinear transform
	
	B0 = b0 + b1*c + b2*c*c;
	B1 = 2*b0 - 2*b2*c*c;
	B2 = b0 - b1*c + b2*c*c;
	A0 = a0 + a1*c + a2*c*c;
	A1 = 2*a0 - 2*a2*c*c;
	A2 = a0 - a1*c + a2*c*c;
	
	y[0] = (-A1*y_1 -A2*y_2 + B0*u[0] + B1*u_1 + B2*u_2 )/A0;
	y[1] = (-A1*y[0] -A2*y_1 + B0*u[1] + B1*u[0] + B2*u_1 )/A0;
		
	for (int i=2; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] -A2*y[i-2] + B0*u[i] + B1*u[i-1] + B2*u[i-2] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	
}

void DS1_Clip_Tone(float *u, float *y, float *v1, float *v2, float *v3,  int N, float T, float *U_1, float *Y_1, float *V1_1, float *V2_1, float *V3_1, float t, float vol)
{
	
	const float R1 = 2.2e3;
	const float R2 = 6.8e3;
	const float R3 = 2.2e3;
	const float R4 = 6.8e3;
	const float Rt = 20e3;
	const float Rv = 100e3;
	
	const float C1 = 470e-9;
	const float C2 = 10e-9;
	const float C3 = 100e-9;
	const float C4 = 22e-9;


	const float Is = 2.52e-9;
	const float Vt = 45.3e-3;
	
	float c1 = (1 + R3*( 1/R4 + 1/Rt ) )/( t*vol ) + (1-t)*(Rt/Rv)*(1 + R3*(1/R4 + 1/( (1-t)*Rt) ) )/vol;
	float c2 = ( (t-1)/t )*(1 + R3*(1/R4 + 1/( (1-t)*Rt) ) );
	float c3 = ( 1/R4 + 1/Rt )/( t*vol ) + (1-t)*(Rt/Rv)*(1/R4 + 1/( (1-t)*Rt) )/vol;
	float c4 = ( (t-1)/t )*(1/R4 + 1/( (1-t)*Rt) );
	
	float y_1 = Y_1[0];
	float u_1 = U_1[0];
	float v1_1 = V1_1[0];
	float v2_1 = V2_1[0];
	float v3_1 = V3_1[0];
	
	const float E1_1 = 1;
	//float const E1_2 = -1;
	//float const E1_3 = 0;
	//float const E1_4 = 0;
	
	const float E2_1 = 0;
	float E2_2 = 1 + ((T*Is)/(C2*Vt))*COSH( v2_1/Vt ); //Muda
	const float E2_3 = 0;
	const float E2_4 = 0;
	
	//float const E3_1 = 0;
	const float E3_2 = 0;
	const float E3_3 = 1;
	const float E3_4 = 0;
	
	//float const E4_1 = 0;
	//float const E4_2 = 1;
	float E4_3 = -c2;
	float E4_4 = -c1;
	
	const float F1_1 = 1/( 2*R1*C1 );
	//float const F1_2 = 0;
	//float const F1_3 = 0;
	//float const F1_4 = 0;
	
	const float F2_1 = 1/( 2*R1*C2 );
	const float F2_2 = 1/( 2*R2*C2 );
	float F2_3 = (c4 - 1/R2)/( 2*C2 );
	float F2_4 = c3/(2*C2);
	
	//float const F3_1 = 0;
	float const F3_2 = -1/( 2*R2*C3 );
	float F3_3 = ( 1/R2 + 1/(t*Rt) )/( 2*C3 );
	float F3_4 = -1/( 2*vol*t*Rt*C3 );
	
	//float const F4_1 = 0;
	//float const F4_2 = 0;
	float F4_3 = -c4/(2*C4);
	float F4_4 = -c3/(2*C4);
	
	
	const float A1_1 = E1_1 + T*F1_1;
	//float A1_2 = -1;
	//float A1_3 = 0;
	//float A1_4 = 0;
	
	const float A2_1 = E2_1 + T*F2_1;
	float A2_2 = E2_2 + T*F2_2; //Muda
	float A2_3 = E2_3 + T*F2_3;
	float A2_4 = E2_4 + T*F2_4;
	
	//float A3_1 = 0;
	const float A3_2 = E3_2 + T*F3_2;
	float A3_3 = E3_3 + T*F3_3;
	float A3_4 = E3_4 + T*F3_4;
	
	//float A4_1 = 0;
	//float A4_2 = 1;
	float A4_3 = E4_3 + T*F4_3;
	float A4_4 = E4_4 + T*F4_4;

	
	const float A_1_1 = E1_1 - T*F1_1;
	const float A_1_2 = -1;
	const float A_1_3 = 0;
	const float A_1_4 = 0;
	
	const float A_2_1 = E2_1 - T*F2_1;
	float A_2_2 = E2_2 - T*F2_2; //Muda
	float A_2_3 = E2_3 - T*F2_3;
	float A_2_4 = E2_4 - T*F2_4;
	
	const float A_3_1 = 0;
	const float A_3_2 = E3_2 - T*F3_2;
	float A_3_3 = E3_3 - T*F3_3;
	float A_3_4 = E3_4 - T*F3_4;
	
	const float A_4_1 = 0;
	const float A_4_2 = 1;
	float A_4_3 = E4_3 - T*F4_3;
	float A_4_4 = E4_4 - T*F4_4;
	
	
	float B1 = (T/( 2*R1*C1 ))*( u[0] + u_1) + A_1_1*v1_1 + A_1_2*v2_1 + A_1_3*v3_1 + A_1_4*y_1 ; //Muda
	float B2 = (T/( 2*R1*C2 ))*( u[0] + u_1) + A_2_1*v1_1 + A_2_2*v2_1 + A_2_3*v3_1 + A_2_4*y_1 -(2*T*Is/C2)*SINH( v2_1/Vt) ; //Muda
	float B3 = A_3_1*v1_1 + A_3_2*v2_1 + A_3_3*v3_1 + A_3_4*y_1; //Muda
	float B4 = A_4_1*v1_1 + A_4_2*v2_1 + A_4_3*v3_1 + A_4_4*y_1; //Muda
	
	/*
		
		A*X[k] = B*X[k-1] + Y
		
		X[k] = [V1[k] V2[k] V3[k] Vout[k] ]'
		X[k-1] = [V1[k-1] V2[k-1] V3[k-1] Vout[k-1] ]'
		
	*/
	
	v1[0] =  (A2_3*A3_4*B1 - A2_4*A3_3*B1 + A2_3*A3_4*B4 - A2_4*A3_3*B4 - A2_3*A4_4*B3 + A2_4*A4_3*B3 + A3_3*A4_4*B2 - A3_4*A4_3*B2 + A2_2*A3_3*A4_4*B1 - A2_2*A3_4*A4_3*B1 - A2_3*A3_2*A4_4*B1 + A2_4*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	v2[0] =  (A1_1*A2_3*A3_4*B4 - A1_1*A2_4*A3_3*B4 - A1_1*A2_3*A4_4*B3 + A1_1*A2_4*A4_3*B3 + A1_1*A3_3*A4_4*B2 - A1_1*A3_4*A4_3*B2 - A2_1*A3_3*A4_4*B1 + A2_1*A3_4*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	v3[0] = -(A1_1*A2_4*B3 - A1_1*A3_4*B2 + A2_1*A3_4*B1 + A2_1*A3_4*B4 - A2_1*A4_4*B3 + A1_1*A2_2*A3_4*B4 - A1_1*A2_4*A3_2*B4 - A1_1*A2_2*A4_4*B3 + A1_1*A3_2*A4_4*B2 - A2_1*A3_2*A4_4*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	y[0] =  (A1_1*A2_3*B3 - A1_1*A3_3*B2 + A2_1*A3_3*B1 + A2_1*A3_3*B4 - A2_1*A4_3*B3 + A1_1*A2_2*A3_3*B4 - A1_1*A2_3*A3_2*B4 - A1_1*A2_2*A4_3*B3 + A1_1*A3_2*A4_3*B2 - A2_1*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	
	for (int i=1; i<=N-1; i++)
	{
		E2_2 = 1 + ((T*Is)/(C2*Vt))*COSH( v2[i-1]/Vt );
		A2_2 = E2_2 + T*F2_2;
		A_2_2 = E2_2 - T*F2_2;
		
		B1 = (T/( 2*R1*C1 ))*( u[i] + u[i-1]) + A_1_1*v1[i-1] + A_1_2*v2[i-1] + A_1_3*v3[i-1] + A_1_4*y[i-1] ; //Muda
		B2 = (T/( 2*R1*C2 ))*( u[i] + u[i-1]) + A_2_1*v1[i-1] + A_2_2*v2[i-1] + A_2_3*v3[i-1] + A_2_4*y[i-1] -(2*T*Is/C2)*SINH( v2[i-1]/Vt) ; //Muda
		B3 = A_3_1*v1[i-1] + A_3_2*v2[i-1] + A_3_3*v3[i-1] + A_3_4*y[i-1]; //Muda
		B4 = A_4_1*v1[i-1] + A_4_2*v2[i-1] + A_4_3*v3[i-1] + A_4_4*y[i-1];
		
		v1[i] =  (A2_3*A3_4*B1 - A2_4*A3_3*B1 + A2_3*A3_4*B4 - A2_4*A3_3*B4 - A2_3*A4_4*B3 + A2_4*A4_3*B3 + A3_3*A4_4*B2 - A3_4*A4_3*B2 + A2_2*A3_3*A4_4*B1 - A2_2*A3_4*A4_3*B1 - A2_3*A3_2*A4_4*B1 + A2_4*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
		v2[i] =  (A1_1*A2_3*A3_4*B4 - A1_1*A2_4*A3_3*B4 - A1_1*A2_3*A4_4*B3 + A1_1*A2_4*A4_3*B3 + A1_1*A3_3*A4_4*B2 - A1_1*A3_4*A4_3*B2 - A2_1*A3_3*A4_4*B1 + A2_1*A3_4*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
		v3[i] = -(A1_1*A2_4*B3 - A1_1*A3_4*B2 + A2_1*A3_4*B1 + A2_1*A3_4*B4 - A2_1*A4_4*B3 + A1_1*A2_2*A3_4*B4 - A1_1*A2_4*A3_2*B4 - A1_1*A2_2*A4_4*B3 + A1_1*A3_2*A4_4*B2 - A2_1*A3_2*A4_4*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
		y[i] =  (A1_1*A2_3*B3 - A1_1*A3_3*B2 + A2_1*A3_3*B1 + A2_1*A3_3*B4 - A2_1*A4_3*B3 + A1_1*A2_2*A3_3*B4 - A1_1*A2_3*A3_2*B4 - A1_1*A2_2*A4_3*B3 + A1_1*A3_2*A4_3*B2 - A2_1*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	V1_1[0] = v1[N-1];
	V2_1[0] = v2[N-1];
	V3_1[0] = v3[N-1];
	
}

void Filter1_48000(float *u, float *y, int N, float *U_1, float *Y_1 )
{
	const float R1 = 1000;
	const float R2 = 470e3;

	const float C1 = 47e-9;
	
	const float c = 2*2*48000;
	
	float y_1 = Y_1[0];
	float u_1 = U_1[0];
	
	/*
	             R2C1s       b1s + b0
	 G(s) = -------------- = --------
	        (R1+R2)C1s + 1   a1s + a0

	         B0 + B1z⁻¹
	 G(z) =  ----------
             A0 + A1z⁻¹
     
     y[k] = (-A1*y[k-1] + B0*u[k] + B1*u[k-1] )/A0
	*/
	
	const float b0 = 0;
	const float b1 = R2*C1;
	const float a0 = 1;
	const float a1 = (R1+R2)*C1;
	
	const float B0 = b0 + b1*c;
	const float B1 = b0 - b1*c;
	const float A0 = a0 + a1*c;
	const float A1 = a0 - a1*c;	
	
	y[0] = (-A1*y_1 + B0*u[0] + B1*u_1)/A0;
		
	for (int i=1; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] + B0*u[i] + B1*u[i-1] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	
}

void Filter2_48000(float *u, float *y, int N, float *U_1, float *Y_1, float *U_2, float *Y_2, float *U_3, float *Y_3, float *U_4, float *Y_4 )
{
	const double R45 = 90.909090909090907e3;
	const double R6 = 100e3;
	const double R7 = 470e3;
	const double R8 = 10e3;
	const double R9 = 22;
	const double R10 = 100e3;
	
	const double G45 = 1/R45;
	const double G6 = 1/R6;
	const double G7 = 1/R7;
	const double G8 = 1/R8;
	const double G9 = 1/R9;
	const double G10 = 1/R10;

	const double C2 = 470e-9;
	const double C3 = 47e-9;
	const double C4 = 250e-12;
	const double C5 = 68e-9;
	
	const double c = 2*2*48000;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	double y_3 = Y_3[0];
	double u_3 = U_3[0];
	double y_4 = Y_4[0];
	double u_4 = U_4[0];
	
	
	
	double const b4 = C2*C3*C4*C5;
	double const b3 = C2*C3*C5*(G7 - G9);
	double const b2 = 0;
	double const b1 = 0;
	double const b0 = 0;
	
	double const a4 = C2*C3*C4*C5;
	double const a3 = (C2*C3*C4 + C2*C3*C5 + C2*C4*C5 + C3*C4*C5)*G10 + C2*C3*C5*G7 + C2*C4*C5*G6 + C2*C3*C5*G8 + C3*C4*C5*G6 + C2*C4*C5*G8 + C2*C4*C5*G9 + C3*C4*C5*G8 + C3*C4*C5*G9 + C3*C4*C5*G45;
	double const a2 = (C2*C3*G7 + C2*C4*G6 + C2*C3*G8 + C2*C5*G6 + C3*C4*G6 + C2*C4*G8 + C2*C5*G7 + C3*C5*G6 + C2*C4*G9 + C3*C4*G8 + C3*C5*G7 + C3*C4*G9 + C3*C4*G45 + C3*C5*G45 + C4*C5*G45)*G10 + C2*C5*G6*G7 + C2*C5*G6*G8 + C3*C5*G6*G7 + C2*C5*G7*G8 + C3*C5*G6*G8 + C2*C5*G7*G9 + C3*C5*G7*G8 + C3*C5*G7*G9 + C3*C5*G7*G45 + C4*C5*G6*G45 + C3*C5*G8*G45 + C4*C5*G8*G45 + C4*C5*G9*G45;
	double const a1 = (C2*G6*G7 + C2*G6*G8 + C3*G6*G7 + C2*G7*G8 + C3*G6*G8 + C2*G7*G9 + C3*G7*G8 + C3*G7*G9 + C3*G7*G45 + C4*G6*G45 + C3*G8*G45 + C5*G6*G45 + C4*G8*G45 + C5*G7*G45 + C4*G9*G45)*G10 + C5*G6*G7*G45 + C5*G6*G8*G45 + C5*G7*G8*G45 + C5*G7*G9*G45;
	double const a0 = G10*G45*(G6*G7 + G6*G8 + G7*G8 + G7*G9);

	
	const double B0 = b4*pow(c,4) + b3*pow(c,3) + b2*c*c + b1*c + b0;
	const double B1 = -4*b4*pow(c,4) - 2*b3*pow(c,3) + 2*b1*c + 4*b0;
	const double B2 = 6*b4*pow(c,4) - 2*b2*c*c + 6*b0;
	const double B3 = -4*b4*pow(c,4) + 2*b3*pow(c,3) - 2*b1*c + 4*b0;
	const double B4 = b4*pow(c,4) - b3*pow(c,3) + b2*c*c - b1*c + b0;
	const double A0 = a4*pow(c,4) + a3*pow(c,3) + a2*c*c + a1*c + a0;
	const double A1 = -4*a4*pow(c,4) - 2*a3*pow(c,3) + 2*a1*c + 4*a0;
	const double A2 = 6*a4*pow(c,4) - 2*a2*c*c + 6*a0;
	const double A3 = -4*a4*pow(c,4) + 2*a3*pow(c,3) - 2*a1*c + 4*a0;
	const double A4 = a4*pow(c,4) - a3*pow(c,3) + a2*c*c - a1*c + a0;
	
	y[0] = (-A1*y_1 -A2*y_2 - A3*y_3 - A4*y_4 + B0*u[0] + B1*u_1 + B2*u_2 + B3*u_3 + B4*u_4 )/A0;
	y[1] = (-A1*y[0] -A2*y_1 - A3*y_2 - A4*y_3 + B0*u[1] + B1*u[0] + B2*u_1 + B3*u_2 + B4*u_3 )/A0;
	y[2] = (-A1*y[1] -A2*y[0] - A3*y_1 - A4*y_2 + B0*u[2] + B1*u[1] + B2*u[0] + B3*u_1 + B4*u_2 )/A0;
	y[3] = (-A1*y[2] -A2*y[1] - A3*y[0] - A4*y_1 + B0*u[3] + B1*u[2] + B2*u[1] + B3*u[0] + B4*u_1 )/A0;
		
	for (int i=4; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] -A2*y[i-2] - A3*y[i-3] - A4*y[i-4] + B0*u[i] + B1*u[i-1] + B2*u[i-2] + B3*u[i-3] + B4*u[i-4] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	U_4[0] = u[N-4];
	Y_4[0] = y[N-4];
	
}

void FilterGain_48000(float *u, float *y, int N, float Dist, float *U_1, float *Y_1, float *U_2, float *Y_2 )
{
	const float Rd = 100e3;
	const float R13 = 4.7e3;
	const float Cc = 10e-12;
	const float Cz = 470e-9;
	
	//const float Cc = 250e-12;
	//const float Cz = 1000e-9;
	
	float Rt = Dist*Rd;
	float Rb = (1-Dist)*Rd + R13;
	
	float Gt = 1/Rt;
	float Gb = 1/Rb;
	
	float b0;
	float b1;
	float b2;
	float a0;
	float a1;
	float a2;
	
	float B0;
	float B1;
	float B2;
	float A0;
	float A1;
	float A2;
	
	const float c = 2*2*48000;
	
	float y_1 = Y_1[0];
	float u_1 = U_1[0];
	float y_2 = Y_2[0];
	float u_2 = U_2[0];
	
	/*
	  
	          Cz Cc s² + (Cc Gb + Cz Gb + Cz Gt) s + Gb Gt    b2s² + b1s + b0
	 G(s) =  --------------------------------------------- =  ----------------
	          Cz Cc s² + (Cc Gb + Cz Gt) s + Gb Gt            a2s² + a1s + a0
	           
	Rt = DRd
	Rb = (1-D)Rd + R13 
	                   
	         B0 + B1z⁻¹ + B2z⁻²
	 G(z) =  -------------------
             A0 + A1z⁻¹ + A2z⁻²
     
     y[k] = (-A1*y[k-1] -A2*y[k-2] + B0*u[k] + B1*u[k-1] + B2*u[k-2] )/A0
     
	*/
	
	a0 = Gb*Gt;
	a1 = Cc*Gb + Cz*Gt;
	a2 = Cz*Cc;
	b0 = a0;
	b1 = a1 + Cz*Gb;
	b2 = a2;
	
	//Bilinear transform
	
	B0 = b0 + b1*c + b2*c*c;
	B1 = 2*b0 - 2*b2*c*c;
	B2 = b0 - b1*c + b2*c*c;
	A0 = a0 + a1*c + a2*c*c;
	A1 = 2*a0 - 2*a2*c*c;
	A2 = a0 - a1*c + a2*c*c;
	
	y[0] = (-A1*y_1 -A2*y_2 + B0*u[0] + B1*u_1 + B2*u_2 )/A0;
	y[1] = (-A1*y[0] -A2*y_1 + B0*u[1] + B1*u[0] + B2*u_1 )/A0;
		
	for (int i=2; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] -A2*y[i-2] + B0*u[i] + B1*u[i-1] + B2*u[i-2] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	
}

void DS1_Clip_Tone_48000(float *u, float *y, float *v1, float *v2, float *v3,  int N, float *U_1, float *Y_1, float *V1_1, float *V2_1, float *V3_1, float t, float vol)
{
	const float SampleRate = 8*48000;
	const float T = 1/SampleRate;
	
	const float R1 = 2.2e3;
	const float R2 = 6.8e3;
	const float R3 = 2.2e3;
	const float R4 = 6.8e3;
	const float Rt = 20e3;
	const float Rv = 100e3;
	
	const float C1 = 470e-9;
	const float C2 = 10e-9;
	const float C3 = 100e-9;
	const float C4 = 22e-9;


	const float Is = 2.52e-9;
	const float Vt = 45.3e-3;
	
	float c1 = (1 + R3*( 1/R4 + 1/Rt ) )/( t*vol ) + (1-t)*(Rt/Rv)*(1 + R3*(1/R4 + 1/( (1-t)*Rt) ) )/vol;
	float c2 = ( (t-1)/t )*(1 + R3*(1/R4 + 1/( (1-t)*Rt) ) );
	float c3 = ( 1/R4 + 1/Rt )/( t*vol ) + (1-t)*(Rt/Rv)*(1/R4 + 1/( (1-t)*Rt) )/vol;
	float c4 = ( (t-1)/t )*(1/R4 + 1/( (1-t)*Rt) );
	
	float y_1 = Y_1[0];
	float u_1 = U_1[0];
	float v1_1 = V1_1[0];
	float v2_1 = V2_1[0];
	float v3_1 = V3_1[0];
	
	const float E1_1 = 1;
	//float const E1_2 = -1;
	//float const E1_3 = 0;
	//float const E1_4 = 0;
	
	const float E2_1 = 0;
	float E2_2 = 1 + ((T*Is)/(C2*Vt))*COSH( v2_1/Vt ); //Muda
	const float E2_3 = 0;
	const float E2_4 = 0;
	
	//float const E3_1 = 0;
	const float E3_2 = 0;
	const float E3_3 = 1;
	const float E3_4 = 0;
	
	//float const E4_1 = 0;
	//float const E4_2 = 1;
	float E4_3 = -c2;
	float E4_4 = -c1;
	
	const float F1_1 = 1/( 2*R1*C1 );
	//float const F1_2 = 0;
	//float const F1_3 = 0;
	//float const F1_4 = 0;
	
	const float F2_1 = 1/( 2*R1*C2 );
	const float F2_2 = 1/( 2*R2*C2 );
	float F2_3 = (c4 - 1/R2)/( 2*C2 );
	float F2_4 = c3/(2*C2);
	
	//float const F3_1 = 0;
	float const F3_2 = -1/( 2*R2*C3 );
	float F3_3 = ( 1/R2 + 1/(t*Rt) )/( 2*C3 );
	float F3_4 = -1/( 2*vol*t*Rt*C3 );
	
	//float const F4_1 = 0;
	//float const F4_2 = 0;
	float F4_3 = -c4/(2*C4);
	float F4_4 = -c3/(2*C4);
	
	
	float A1_1 = E1_1 + T*F1_1;
	//float A1_2 = -1;
	//float A1_3 = 0;
	//float A1_4 = 0;
	
	float A2_1 = E2_1 + T*F2_1;
	float A2_2 = E2_2 + T*F2_2; //Muda
	float A2_3 = E2_3 + T*F2_3;
	float A2_4 = E2_4 + T*F2_4;
	
	//float A3_1 = 0;
	float A3_2 = E3_2 + T*F3_2;
	float A3_3 = E3_3 + T*F3_3;
	float A3_4 = E3_4 + T*F3_4;
	
	//float A4_1 = 0;
	//float A4_2 = 1;
	float A4_3 = E4_3 + T*F4_3;
	float A4_4 = E4_4 + T*F4_4;

	
	float A_1_1 = E1_1 - T*F1_1;
	const float A_1_2 = -1;
	const float A_1_3 = 0;
	const float A_1_4 = 0;
	
	float A_2_1 = E2_1 - T*F2_1;
	float A_2_2 = E2_2 - T*F2_2; //Muda
	float A_2_3 = E2_3 - T*F2_3;
	float A_2_4 = E2_4 - T*F2_4;
	
	const float A_3_1 = 0;
	float A_3_2 = E3_2 - T*F3_2;
	float A_3_3 = E3_3 - T*F3_3;
	float A_3_4 = E3_4 - T*F3_4;
	
	const float A_4_1 = 0;
	const float A_4_2 = 1;
	float A_4_3 = E4_3 - T*F4_3;
	float A_4_4 = E4_4 - T*F4_4;
	
	
	float B1 = (T/( 2*R1*C1 ))*( u[0] + u_1) + A_1_1*v1_1 + A_1_2*v2_1 + A_1_3*v3_1 + A_1_4*y_1 ; //Muda
	float B2 = (T/( 2*R1*C2 ))*( u[0] + u_1) + A_2_1*v1_1 + A_2_2*v2_1 + A_2_3*v3_1 + A_2_4*y_1 -(2*T*Is/C2)*SINH( v2_1/Vt) ; //Muda
	float B3 = A_3_1*v1_1 + A_3_2*v2_1 + A_3_3*v3_1 + A_3_4*y_1; //Muda
	float B4 = A_4_1*v1_1 + A_4_2*v2_1 + A_4_3*v3_1 + A_4_4*y_1; //Muda
	
	/*
		
		A*X[k] = B*X[k-1] + Y
		
		X[k] = [V1[k] V2[k] V3[k] Vout[k] ]'
		X[k-1] = [V1[k-1] V2[k-1] V3[k-1] Vout[k-1] ]'
		
	*/
	
	v1[0] =  (A2_3*A3_4*B1 - A2_4*A3_3*B1 + A2_3*A3_4*B4 - A2_4*A3_3*B4 - A2_3*A4_4*B3 + A2_4*A4_3*B3 + A3_3*A4_4*B2 - A3_4*A4_3*B2 + A2_2*A3_3*A4_4*B1 - A2_2*A3_4*A4_3*B1 - A2_3*A3_2*A4_4*B1 + A2_4*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	v2[0] =  (A1_1*A2_3*A3_4*B4 - A1_1*A2_4*A3_3*B4 - A1_1*A2_3*A4_4*B3 + A1_1*A2_4*A4_3*B3 + A1_1*A3_3*A4_4*B2 - A1_1*A3_4*A4_3*B2 - A2_1*A3_3*A4_4*B1 + A2_1*A3_4*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	v3[0] = -(A1_1*A2_4*B3 - A1_1*A3_4*B2 + A2_1*A3_4*B1 + A2_1*A3_4*B4 - A2_1*A4_4*B3 + A1_1*A2_2*A3_4*B4 - A1_1*A2_4*A3_2*B4 - A1_1*A2_2*A4_4*B3 + A1_1*A3_2*A4_4*B2 - A2_1*A3_2*A4_4*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	y[0] =  (A1_1*A2_3*B3 - A1_1*A3_3*B2 + A2_1*A3_3*B1 + A2_1*A3_3*B4 - A2_1*A4_3*B3 + A1_1*A2_2*A3_3*B4 - A1_1*A2_3*A3_2*B4 - A1_1*A2_2*A4_3*B3 + A1_1*A3_2*A4_3*B2 - A2_1*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	
	for (int i=1; i<=N-1; i++)
	{
		E2_2 = 1 + ((T*Is)/(C2*Vt))*COSH( v2[i-1]/Vt );
		A2_2 = E2_2 + T*F2_2;
		A_2_2 = E2_2 - T*F2_2;
		
		B1 = (T/( 2*R1*C1 ))*( u[i] + u[i-1]) + A_1_1*v1[i-1] + A_1_2*v2[i-1] + A_1_3*v3[i-1] + A_1_4*y[i-1] ; //Muda
		B2 = (T/( 2*R1*C2 ))*( u[i] + u[i-1]) + A_2_1*v1[i-1] + A_2_2*v2[i-1] + A_2_3*v3[i-1] + A_2_4*y[i-1] -(2*T*Is/C2)*SINH( v2[i-1]/Vt) ; //Muda
		B3 = A_3_1*v1[i-1] + A_3_2*v2[i-1] + A_3_3*v3[i-1] + A_3_4*y[i-1]; //Muda
		B4 = A_4_1*v1[i-1] + A_4_2*v2[i-1] + A_4_3*v3[i-1] + A_4_4*y[i-1];
		
		v1[i] =  (A2_3*A3_4*B1 - A2_4*A3_3*B1 + A2_3*A3_4*B4 - A2_4*A3_3*B4 - A2_3*A4_4*B3 + A2_4*A4_3*B3 + A3_3*A4_4*B2 - A3_4*A4_3*B2 + A2_2*A3_3*A4_4*B1 - A2_2*A3_4*A4_3*B1 - A2_3*A3_2*A4_4*B1 + A2_4*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
		v2[i] =  (A1_1*A2_3*A3_4*B4 - A1_1*A2_4*A3_3*B4 - A1_1*A2_3*A4_4*B3 + A1_1*A2_4*A4_3*B3 + A1_1*A3_3*A4_4*B2 - A1_1*A3_4*A4_3*B2 - A2_1*A3_3*A4_4*B1 + A2_1*A3_4*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
		v3[i] = -(A1_1*A2_4*B3 - A1_1*A3_4*B2 + A2_1*A3_4*B1 + A2_1*A3_4*B4 - A2_1*A4_4*B3 + A1_1*A2_2*A3_4*B4 - A1_1*A2_4*A3_2*B4 - A1_1*A2_2*A4_4*B3 + A1_1*A3_2*A4_4*B2 - A2_1*A3_2*A4_4*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
		y[i] =  (A1_1*A2_3*B3 - A1_1*A3_3*B2 + A2_1*A3_3*B1 + A2_1*A3_3*B4 - A2_1*A4_3*B3 + A1_1*A2_2*A3_3*B4 - A1_1*A2_3*A3_2*B4 - A1_1*A2_2*A4_3*B3 + A1_1*A3_2*A4_3*B2 - A2_1*A3_2*A4_3*B1)/(A1_1*A2_3*A3_4 - A1_1*A2_4*A3_3 + A2_1*A3_3*A4_4 - A2_1*A3_4*A4_3 + A1_1*A2_2*A3_3*A4_4 - A1_1*A2_2*A3_4*A4_3 - A1_1*A2_3*A3_2*A4_4 + A1_1*A2_4*A3_2*A4_3);
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	V1_1[0] = v1[N-1];
	V2_1[0] = v2[N-1];
	V3_1[0] = v3[N-1];
	
}
