#!/usr/bin/env python

import sys

from math import *
from decimal import *
getcontext().prec = 50

if len(sys.argv) > 1:
    N = int(sys.argv[1])
else:
    N = Decimal(8928)

dx = Decimal(89.28)/N

inicio = 0
fim = (N-1)*dx
Idx = Decimal(1/dx)

f = open('Cosh.h', 'w')
f2 = open('Sinh.h', 'w')
f3 = open('ASinh.h', 'w')

f.write('using namespace std; \n \n')
f2.write('using namespace std; \n \n')
f3.write('using namespace std; \n \n')

f.write('#define COSH_N ' + str(N) + '\n')
f.write('#define COSH_Idx ' + str(Idx) + '\n')
f.write('#define COSH_fim ' + str(fim) + '\n\n')

f2.write('#define SINH_N ' + str(N) + '\n')
f2.write('#define SINH_Idx ' + str(Idx) + '\n')
f3.write('#define SINH_fim ' + str(fim) + '\n\n')

f3.write('#define ASINH_N ' + str(N) + '\n')
f3.write('#define ASINH_Idx ' + str(Idx) + '\n')
f3.write('#define ASINH_fim ' + str(fim) + '\n\n')


f.write('const float Cosh[] = {')
f2.write('const float Sinh[] = {')
f3.write('const float ASinh[] = {')

for i in range(N):
    value = Decimal ( cosh( Decimal ( i*dx ) ) )
    value2 = Decimal ( sinh( Decimal ( i*dx ) ) )
    value3 = Decimal ( asinh( Decimal ( i*dx ) ) )
    v = "%.40e" % value
    v2 = "%.40e" % value2
    v3 = "%.40e" % value3
    f.write(v + ',\n')	
    f2.write(v2 + ',\n')
    f3.write(v3 + ',\n')

f.write('\n};')
f.close()

f2.write('\n};')
f2.close()

f3.write('\n};')
f3.close()
