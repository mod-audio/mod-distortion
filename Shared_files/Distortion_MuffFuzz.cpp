#include <cmath>
#include "Distortion_MuffFuzz.h"
#include "HyperbolicTables.h"

void Filter1(double *u, double *y, int N, double T, double *U_1, double *Y_1 )
{
	const double R1 = 100e3;
	const double R2 = 470e3;

	const double C1 = 0.01e-6;
	
	double c = 2/T;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	
	/*
	             R2C1s       b1s + b0
	 G(s) = -------------- = --------
	           R1C1s + 1     a1s + a0

	         B0 + B1z⁻¹
	 G(z) =  ----------
             A0 + A1z⁻¹
     
     y[k] = (-A1*y[k-1] + B0*u[k] + B1*u[k-1] )/A0
	*/
	
	const double b0 = 0;
	const double b1 = R2*C1;
	const double a0 = 1;
	const double a1 = R1*C1;
	
	double B0 = b0 + b1*c;
	double B1 = b0 - b1*c;
	double A0 = a0 + a1*c;
	double A1 = a0 - a1*c;	
	
	y[0] = (-A1*y_1 + B0*u[0] + B1*u_1)/A0;
		
	for (int i=1; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] + B0*u[i] + B1*u[i-1] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	
}

void Clip(double *u, double *y, int N)
{
	const double R3 = 10e3;
	
	const double Is = 2.52e-9;
	const double Vt = 45.3e-3;
	
	const double V_ = 2*Is*R3;
	
	for (int i=0; i<=N-1; i++)
	{
		y[i] = Vt*asinh(u[i]/V_);
	}
}

void Filter2(double *u, double *y, int N, double vol, double T, double *U_1, double *Y_1 )
{
	const double Rv = 100e3;

	const double C2 = 0.1e-6;
	
	double c = 2/T;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	
	/*
	             volRvC2s     b1s + b0
	 G(s) = --------------- = --------
	           volRvC2s + 1   a1s + a0

	         B0 + B1z⁻¹
	 G(z) =  ----------
             A0 + A1z⁻¹
     
     y[k] = (-A1*y[k-1] + B0*u[k] + B1*u[k-1] )/A0
	*/
	
	const double b0 = 0;
	double b1 = vol*Rv*C2;
	const double a0 = 1;
	double a1 = b1;
	
	double B0 = b0 + b1*c;
	double B1 = b0 - b1*c;
	double A0 = a0 + a1*c;
	double A1 = a0 - a1*c;	
	
	y[0] = (-A1*y_1 + B0*u[0] + B1*u_1)/A0;
		
	for (int i=1; i<=N-1; i++)
	{
		y[i] = (-A1*y[i-1] + B0*u[i] + B1*u[i-1] )/A0;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	
}
