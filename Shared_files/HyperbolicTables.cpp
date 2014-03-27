#include <cmath>
#include <iostream>
#include "Sinh.h"
#include "ASinh.h"
#include "Cosh.h"
#include "HyperbolicTables.h"


using namespace std;

double SINH( double x)
{
	int flag = 1;
	
	if (x < 0)
	{
		flag = -1;
		x = -x;
	}
	
	double SinH;
	
	if( x > SINH_fim)
	{
		SinH = (Sinh[SINH_N-1])*flag;
	}
	else
	{
		double naux = x*SINH_Idx;
		int n = round(naux);
		SinH = Sinh[n]*flag;
	}
	
	return SinH;
}

double COSH( double x)
{
	if (x < 0)
	{
		x = -x;
	}
	
	double CosH;
	
	if( x > COSH_fim)
	{
		CosH = Cosh[COSH_N-1];
	}
	else
	{
		double naux = x*COSH_Idx;
		int n = round(naux);
		CosH = Cosh[n];
	}
	
	return CosH;
}

double ASINH( double x)
{
	int flag = 1;
	
	if (x < 0)
	{
		flag = -1;
		x = -x;
	}
	
	double ASinH;
	
	if( x > ASINH_fim)
	{
		ASinH = (ASinh[ASINH_N-1])*flag;
	}
	else
	{
		double naux = x*ASINH_Idx;
		int n = round(naux);
		ASinH = ASinh[n]*flag;
	}
	
	return ASinH;
}
