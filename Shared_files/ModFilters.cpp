#include <cmath>
#include "ModFilters.h"

void LPF1(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 )
{
	double K;
	double A;
	double B;
	double f;
	double wc;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	
	/*
	            z + 1
	 G(z) = K ----------
	           Az + B 
	*/         
	
	f = f_before;
	wc = 2*M_PI*f;
	K = wc*T;
	A = wc*T+2;
	B = wc*T-2;
	
	y[0] = (-B*y_1 + K*(u[0] + u_1) )/A;
	
	for (int i=1; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		K = wc*T;
		A = wc*T+2;
		B = wc*T-2;
		y[i] = (-B*y[i-1] + K*(u[i] + u[i-1]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	
}

void HPF1(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 )
{
	double K;
	double A;
	double B;
	double f;
	double wc;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	
	/*
	            z - 1
	 G(z) = K ----------
	           Az + B 
	*/         
	
	f = f_before;
	wc = 2*M_PI*f;
	K = 2;
	A = wc*T+2;
	B = wc*T-2;
	
	y[0] = (-B*y_1 + K*(u[0] - u_1) )/A;
	
	for (int i=1; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		K = 2;
		A = wc*T+2;
		B = wc*T-2;
		y[i] = (-B*y[i-1] + K*(u[i] - u[i-1]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	
}

void LPF2(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 )
{
	double K;
	double A;
	double B;
	double C;
	double f;
	double wc;
	double wcT;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	
	/*
	            z² + 2z + 1
	 G(z) = K --------------
	           Az² + Bz + C 
	*/         
	
	f = f_before;
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = wcT*wcT;
	A = 4 + 2*sqrt(2)*wcT + wcT*wcT;
	B = 2*wcT*wcT - 8;
	C = 4 - 2*sqrt(2)*wcT + wcT*wcT;
	
	y[0] = (-B*y_1 - C*y_2 + K*(u[0] +  2*u_1 + u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = wcT*wcT;
	A = 4 + 2*sqrt(2)*wcT + wcT*wcT;
	B = 2*wcT*wcT - 8;
	C = 4 - 2*sqrt(2)*wcT + wcT*wcT;
	
	y[1] = (-B*y[0] -C*y_1 + K*(u[1] +  2*u[0] + u_1) )/A;
		
	for (int i=2; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		wcT = wc*T;
		K = wcT*wcT;
		A = 4 + 2*sqrt(2)*wcT + wcT*wcT;
		B = 2*wcT*wcT - 8;
		C = 4 - 2*sqrt(2)*wcT + wcT*wcT;
		y[i] = (-B*y[i-1] -C*y[i-2] + K*(u[i] +  2*u[i-1] + u[i-2]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	
}

void HPF2(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 )
{
	double K;
	double A;
	double B;
	double C;
	double f;
	double wc;
	double wcT;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	
	/*
	            z² - 2z + 1
	 G(z) = K --------------
	           Az² + Bz + C 
	*/         
	
	f = f_before;
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = 4;
	A = 4 + 2*sqrt(2)*wcT + wcT*wcT;
	B = 2*wcT*wcT - 8;
	C = 4 - 2*sqrt(2)*wcT + wcT*wcT;
	
	y[0] = (-B*y_1 - C*y_2 + K*(u[0] -  2*u_1 + u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = 4;
	A = 4 + 2*sqrt(2)*wcT + wcT*wcT;
	B = 2*wcT*wcT - 8;
	C = 4 - 2*sqrt(2)*wcT + wcT*wcT;
	
	y[1] = (-B*y[0] -C*y_1 + K*(u[1] -  2*u[0] + u_1) )/A;
		
	for (int i=2; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		wcT = wc*T;
		K = 4;
		A = 4 + 2*sqrt(2)*wcT + wcT*wcT;
		B = 2*wcT*wcT - 8;
		C = 4 - 2*sqrt(2)*wcT + wcT*wcT;
		y[i] = (-B*y[i-1] -C*y[i-2] + K*(u[i] -  2*u[i-1] + u[i-2]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	
}

void LPF3(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 )
{
	double K;
	double A;
	double B;
	double C;
	double D;
	double f;
	double wc;
	double wcT;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	double y_3 = Y_3[0];
	double u_3 = U_3[0];
	
	/*
	            z³ + 3z² + 3z + 1
	 G(z) = K --------------------
	           Az³ + Bz² + Cz + D 
	*/         
	
	f = f_before;
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = wcT*wcT*wcT;
	A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
	B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
	C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
	D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
	
	y[0] = (-B*y_1 - C*y_2 - D*y_3 + K*(u[0] +  3*u_1 + 3*u_2 + u_3) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = wcT*wcT*wcT;
	A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
	B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
	C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
	D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
	
	y[1] = (-B*y[0] - C*y_1 - D*y_2 + K*(u[1] +  3*u[0] + 3*u_1 + u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(2);
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = wcT*wcT*wcT;
	A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
	B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
	C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
	D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
	
	y[2] = (-B*y[1] - C*y[0] - D*y_1 + K*(u[2] +  3*u[1] + 3*u[0] + u_1) )/A;
	
	
		
	for (int i=3; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		wcT = wc*T;
		K = wcT*wcT*wcT;
		A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
		B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
		C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
		D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
		y[i] = (-B*y[i-1] - C*y[i-2] - D*y[i-3] + K*(u[i] +  3*u[i-1] + 3*u[i-2] + u[i-3]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	
}

void HPF3(double *u, double *y, int N, double f_before, double f_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3 )
{
	double K;
	double A;
	double B;
	double C;
	double D;
	double f;
	double wc;
	double wcT;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	double y_3 = Y_3[0];
	double u_3 = U_3[0];
	
	/*
	            z³ - 3z² + 3z - 1
	 G(z) = K --------------------
	           Az³ + Bz² + Cz + D 
	*/         
	
	f = f_before;
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = 8;
	A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
	B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
	C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
	D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
	
	y[0] = (-B*y_1 - C*y_2 - D*y_3 + K*(u[0] -  3*u_1 + 3*u_2 - u_3) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = 8;
	A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
	B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
	C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
	D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
	
	y[1] = (-B*y[0] - C*y_1 - D*y_2 + K*(u[1] -  3*u[0] + 3*u_1 - u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(2);
	wc = 2*M_PI*f;
	wcT = wc*T;
	K = 8;
	A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
	B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
	C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
	D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
	
	y[2] = (-B*y[1] - C*y[0] - D*y_1 + K*(u[2] -  3*u[1] + 3*u[0] - u_1) )/A;
	
	
		
	for (int i=3; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		wcT = wc*T;
		K = 8;
		A = wcT*wcT*wcT + 4*wcT*wcT  + 8*wcT + 8;
		B = 3*wcT*wcT*wcT + 4*wcT*wcT  - 8*wcT -24;
		C = 3*wcT*wcT*wcT - 4*wcT*wcT  - 8*wcT +24;
		D = wcT*wcT*wcT - 4*wcT*wcT  + 8*wcT - 8;
		y[i] = (-B*y[i-1] - C*y[i-2] - D*y[i-3] + K*(u[i] -  3*u[i-1] + 3*u[i-2] - u[i-3]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	
}

void BPF1(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3, double *Y_4, double *U_4, double *Y_5, double *U_5, double *Y_6, double *U_6 )
{
	double K;
	double A;
	double B;
	double C;
	double f;
	double BW;
	double wc;
	double dw;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	
	/*
	              z² - 1
	 G(z) = K --------------
	           Az² + Bz + C 
	*/         
	
	f = f_before;
	BW = BW_before;
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = 2*dw*T;
	A = 4 + 2*dw*T + (wc*T)*(wc*T);
	B = 2*(wc*T)*(wc*T) - 8;
	C = 4 - 2*dw*T + (wc*T)*(wc*T);
	
	y[0] = (-B*y_1 - C*y_2 + K*(u[0] - u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = 2*dw*T;
	A = 4 + 2*dw*T + (wc*T)*(wc*T);
	B = 2*(wc*T)*(wc*T) - 8;
	C = 4 - 2*dw*T + (wc*T)*(wc*T);
	
	y[1] = (-B*y[0] -C*y_1 + K*(u[1] - u_1) )/A;
		
	for (int i=2; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(i);
		BW = BW_before + ((BW_now - BW_before)/(N-1))*(i);
		wc = 2*M_PI*f;
		dw = 2*M_PI*BW;
		K = 2*dw*T;
		A = 4 + 2*dw*T + (wc*T)*(wc*T);
		B = 2*(wc*T)*(wc*T) - 8;
		C = 4 - 2*dw*T + (wc*T)*(wc*T);
		y[i] = (-B*y[i-1] -C*y[i-2] + K*(u[i] - u[i-2]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	U_4[0] = u[N-4];
	Y_4[0] = y[N-4];
	U_5[0] = u[N-5];
	Y_5[0] = y[N-5];
	U_6[0] = u[N-6];
	Y_6[0] = y[N-6];
	
}

void BPF2(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3, double *Y_4, double *U_4, double *Y_5, double *U_5, double *Y_6, double *U_6 )
{
	double K;
	double A;
	double B;
	double C;
	double D;
	double E;
	double f;
	double BW;
	double wc;
	double dw;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	double y_3 = Y_3[0];
	double u_3 = U_3[0];
	double y_4 = Y_4[0];
	double u_4 = U_4[0];
	
	/*
	                  z⁴ - 2z² + 1
	 G(z) = K --------------------------
	           Az⁴ + Bz³ + Cz² + Dz + E 
	*/         
	
	f = f_before;
	BW = BW_before;
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = 4*(dw*T)*(dw*T);
	A = 16 + 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T + 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	B = -64 - 16*sqrt(2)*dw*T + 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	C = 96 - 8*( dw*dw + 2*wc*wc)*T*T + 6*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	D = -64 + 16*sqrt(2)*dw*T - 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	E = 16 - 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T - 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	
	y[0] = (-B*y_1 - C*y_2 - D*y_3 - E*y_4 + K*(u[0] - 2*u_2 + u_4) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = 4*(dw*T)*(dw*T);
	A = 16 + 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T + 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	B = -64 - 16*sqrt(2)*dw*T + 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	C = 96 - 8*( dw*dw + 2*wc*wc)*T*T + 6*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	D = -64 + 16*sqrt(2)*dw*T - 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	E = 16 - 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T - 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	
	y[1] = (-B*y[0] - C*y_1 - D*y_2 - E*y_3 + K*(u[1] - 2*u_1 + u_3) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(2);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(2);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = 4*(dw*T)*(dw*T);
	A = 16 + 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T + 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	B = -64 - 16*sqrt(2)*dw*T + 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	C = 96 - 8*( dw*dw + 2*wc*wc)*T*T + 6*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	D = -64 + 16*sqrt(2)*dw*T - 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	E = 16 - 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T - 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	
	y[2] = (-B*y[1] - C*y[0] - D*y_1 - E*y_2 + K*(u[2] - 2*u[0] + u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(3);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(3);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = 4*(dw*T)*(dw*T);
	A = 16 + 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T + 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	B = -64 - 16*sqrt(2)*dw*T + 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	C = 96 - 8*( dw*dw + 2*wc*wc)*T*T + 6*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	D = -64 + 16*sqrt(2)*dw*T - 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
	E = 16 - 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T - 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
	
	y[3] = (-B*y[2] - C*y[1] - D*y[0] - E*y_1 + K*(u[3] - 2*u[1] + u_1) )/A;
		
	for (int i=4; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(3);
		BW = BW_before + ((BW_now - BW_before)/(N-1))*(3);
		wc = 2*M_PI*f;
		dw = 2*M_PI*BW;
		K = 4*(dw*T)*(dw*T);
		A = 16 + 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T + 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
		B = -64 - 16*sqrt(2)*dw*T + 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
		C = 96 - 8*( dw*dw + 2*wc*wc)*T*T + 6*(wc*T)*(wc*T)*(wc*T)*(wc*T);
		D = -64 + 16*sqrt(2)*dw*T - 4*sqrt(2)*dw*wc*wc*T*T*T + 4*(wc*T)*(wc*T)*(wc*T)*(wc*T);
		E = 16 - 8*sqrt(2)*dw*T + 4*( dw*dw + 2*wc*wc)*T*T - 2*sqrt(2)*dw*wc*wc*T*T*T + (wc*T)*(wc*T)*(wc*T)*(wc*T);
		
		y[i] = (-B*y[i-1] - C*y[i-2] - D*y[i-3] - E*y[i-4] + K*(u[i] - 2*u[i-2] + u[i-4]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	U_4[0] = u[N-4];
	Y_4[0] = y[N-4];
	U_5[0] = u[N-5];
	Y_5[0] = y[N-5];
	U_6[0] = u[N-6];
	Y_6[0] = y[N-6];
	
}

void BPF3(double *u, double *y, int N, double f_before, double f_now, double BW_before, double BW_now, double T, double *Y_1, double *U_1, double *Y_2, double *U_2, double *Y_3, double *U_3, double *Y_4, double *U_4, double *Y_5, double *U_5, double *Y_6, double *U_6 )
{
	double K;
	double A;
	double B;
	double C;
	double D;
	double E;
	double F;
	double G;
	double f;
	double BW;
	double wc;
	double dw;
	
	double y_1 = Y_1[0];
	double u_1 = U_1[0];
	double y_2 = Y_2[0];
	double u_2 = U_2[0];
	double y_3 = Y_3[0];
	double u_3 = U_3[0];
	double y_4 = Y_4[0];
	double u_4 = U_4[0];
	double y_5 = Y_5[0];
	double u_5 = U_5[0];
	double y_6 = Y_6[0];
	double u_6 = U_6[0];
	
	/*
	                    z⁶ - 3z⁴ + 3z² - 1
	 G(z) = K --------------------------------------
	           Az⁶ + Bz⁵ + Cz⁴ + Dz³ + Ez² + Fz + G 
	*/         
	
	f = f_before;
	BW = BW_before;
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = pow(2*dw*T,3);
	A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
	B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
	C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
	D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
	E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
	F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
	G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

	
	y[0] = (-B*y_1 - C*y_2 - D*y_3 - E*y_4 - F*y_5 - G*y_6 + K*(u[0] - 3*u_2 + 3*u_4 - u_6) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(1);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(1);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = pow(2*dw*T,3);
	A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
	B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
	C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
	D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
	E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
	F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
	G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

	
	y[1] = (-B*y[0] - C*y_1 - D*y_2 - E*y_3 - F*y_4 - G*y_5 + K*(u[1] - 3*u_1 + 3*u_3 - u_5) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(2);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(2);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = pow(2*dw*T,3);
	A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
	B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
	C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
	D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
	E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
	F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
	G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

	y[2] = (-B*y[1] - C*y[0] - D*y_1 - E*y_2 - F*y_3 - G*y_4 + K*(u[2] - 3*u[0] + 3*u_2 - u_4) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(3);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(3);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = pow(2*dw*T,3);
	A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
	B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
	C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
	D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
	E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
	F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
	G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

	
	y[3] = (-B*y[2] - C*y[1] - D*y[0] - E*y_1 - F*y_2 - G*y_3 + K*(u[3] - 3*u[1] + 3*u_1 - u_3) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(4);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(4);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = pow(2*dw*T,3);
	A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
	B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
	C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
	D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
	E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
	F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
	G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

	
	y[4] = (-B*y[3] - C*y[2] - D*y[1] - E*y[0] - F*y_1 - G*y_2 + K*(u[4] - 3*u[2] + 3*u[0] - u_2) )/A;
	
	f = f_before + ((f_now - f_before)/(N-1))*(5);
	BW = BW_before + ((BW_now - BW_before)/(N-1))*(5);
	wc = 2*M_PI*f;
	dw = 2*M_PI*BW;
	K = pow(2*dw*T,3);
	A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
	B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
	C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
	D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
	E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
	F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
	G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

	
	y[5] = (-B*y[4] - C*y[3] - D*y[2] - E*y[1] - F*y[0] - G*y_1 + K*(u[5] - 3*u[3] + 3*u[1] - u_1) )/A;
		
	for (int i=6; i<=N-1; i++)
	{
		f = f_before + ((f_now - f_before)/(N-1))*(3);
		BW = BW_before + ((BW_now - BW_before)/(N-1))*(3);
		wc = 2*M_PI*f;
		dw = 2*M_PI*BW;
		K = pow(2*dw*T,3);
		A = pow(wc*T,6) + (4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (64*dw)*T + 64;
		B = 6*pow(wc*T,6) + (16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (-256*dw)*T - 384;
		C = 15*pow(wc*T,6) + (20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (320*dw)*T + 960;
		D = 20*pow(wc*T,6) + (-16*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (128*dw*dw + 192*wc*wc)*T*T - 1280;
		E = 15*pow(wc*T,6) + (-20*dw*pow(wc,4))*pow(T,5) + (-4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (24*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (- 32*dw*dw - 48*wc*wc)*T*T + (-320*dw)*T + 960;
		F = 6*pow(wc*T,6) + (-16*dw*pow(wc,4))*pow(T,5) + (8*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (- 64*dw*dw - 96*wc*wc)*T*T + (256*dw)*T - 384;
		G = pow(wc*T,6) + (-4*dw*pow(wc,4))*pow(T,5) + (4*wc*wc*(2*dw*dw + 3*wc*wc))*pow(T,4) + (-8*dw*(dw*dw + 4*wc*wc))*pow(T,3) + (32*dw*dw + 48*wc*wc)*T*T + (-64*dw)*T + 64;

		
		y[i] = (-B*y[i-1] - C*y[i-2] - D*y[i-3] - E*y[i-4] - F*y[i-5] - G*y[i-6] + K*(u[i] - 3*u[i-2] + 3*u[i-4] - u[i-6]) )/A;
	}
	
	U_1[0] = u[N-1];
	Y_1[0] = y[N-1];
	U_2[0] = u[N-2];
	Y_2[0] = y[N-2];
	U_3[0] = u[N-3];
	Y_3[0] = y[N-3];
	U_4[0] = u[N-4];
	Y_4[0] = y[N-4];
	U_5[0] = u[N-5];
	Y_5[0] = y[N-5];
	U_6[0] = u[N-6];
	Y_6[0] = y[N-6];
	
}

