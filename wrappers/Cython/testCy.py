from time import clock
import numpy as np

N=350000
T=np.linspace(250,340,N)

from CoolPropCy import nothing,PropsV
nothing(4)
print 'ter'
PropsV('H','T',np.array([300.0]),'D',np.array([1.5]),'R410A')
print 'ter'

t1=clock()
PropsV('H','T',T,'D',np.array([1.5]),'R410A')
elap=clock()-t1
print "PropsV: Elapsed time for {0:d} calls is {1:g} s or {2:g} us/call".format(N,elap,elap/N*1e6)

from CoolPropCy import PropsVempty
t1=clock()
PropsVempty('H','T',T,'D',np.array([1.5]),'R410A')
elap=clock()-t1
print "PropsVempty: Elapsed time for {0:d} calls is {1:g} s or {2:g} us/call".format(N,elap,elap/N*1e6)

from CoolProp.CoolProp import Props,UseSinglePhaseLUT

UseSinglePhaseLUT(1)
Props('H','T',300,'P',100,'R410A')

t1=clock()
H=np.zeros_like(T)
for i in xrange(N):
    H[i]=Props('H','T',T[i],'D',1.5,'R410A')
elap=clock()-t1
print "Elapsed time for {0:d} calls is {1:g} s or {2:g} us/call".format(N,elap,elap/N*1e6)

t1=clock()
H=[Props('H','T',T_,'D',1.5,'R410A') for T_ in T]
elap=clock()-t1
print "Elapsed time for {0:d} calls is {1:g} s or {2:g} us/call".format(N,elap,elap/N*1e6)