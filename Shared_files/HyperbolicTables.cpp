#include <cmath>
#include <iostream>
#include "Sinh.h"
#include "Cosh.h"
#include "HyperbolicTables.h"


using namespace std;

double SINH( double x)
{
	int N = 250000;
	double dx = 0.0001;
	double inicio = 0;
	double fim = (N-1)*dx;
	int flag = 1;
	
	if (x < 0)
	{
		flag = -1;
		x = -x;
	}
	
	double n;
	n = ((x-inicio)/(fim-inicio))*(N-1);
	int n1 = floor(n);
	int n2 = ceil(n);
	
	double SinH;
	double SinH1;
	double SinH2;
	
	if( x > fim)
	{
		SinH = (Sinh[N-1] + Cosh[N-1]*(x-fim))*flag;
	}
	else
	{
		SinH1 = Sinh[n1];
		SinH2 = Sinh[n2];
		SinH = (SinH1 + (SinH2-SinH1)*(n-n1))*flag;
	}
	
	return SinH;
}

double COSH( double x)
{
	int N = 250000;
	double dx = 0.0001;
	double inicio = 0;
	double fim = (N-1)*dx;
	
	if (x < 0)
	{
		x = -x;
	}
	
	double n;
	n = ((x-inicio)/(fim-inicio))*(N-1);
	int n1 = floor(n);
	int n2 = ceil(n);
	
	double CosH;
	double CosH1;
	double CosH2;
	
	if( x > fim)
	{
		CosH = Cosh[N-1] + Sinh[N-1]*(x-fim);
	}
	else
	{
		CosH1 = Cosh[n1];
		CosH2 = Cosh[n2];
		CosH = CosH1 + (CosH2-CosH1)*(n-n1);
	}
	
	return CosH;
}
