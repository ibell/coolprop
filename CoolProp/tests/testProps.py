import CoolProp.State as CPS
import timeit

N = 1000000
print timeit.Timer("CP.Props('C','T',1200,'D',1.1,'Air')","import CoolProp.CoolProp as CP").timeit(N)*1000000/N
print timeit.Timer("CP.Props('O','T',1200,'D',1.1,'Air')","import CoolProp.CoolProp as CP").timeit(N)*1000000/N

print 'Enthalpy'
print timeit.Timer("CP.Props('H','T',1200,'D',1.1,'Air')","import CoolProp.CoolProp as CP").timeit(N)*1000000/N
print timeit.Timer("S.h","import CoolProp.State as CPS; S = CPS.State('Air',dict(T = 300, D = 1.1))").timeit(N)*1000000/N
print timeit.Timer("S.get_h()","import CoolProp.State as CPS; S = CPS.State('Air',dict(T = 300, D = 1.1))").timeit(N)*1000000/N

print 'Pressure'
print timeit.Timer("CP.Props('P','T',1200,'D',1.1,'Air')","import CoolProp.CoolProp as CP").timeit(N)*1000000/N
print timeit.Timer("S.p","import CoolProp.State as CPS; S = CPS.State('Air',dict(T = 300, D = 1.1))").timeit(N)*1000000/N
print timeit.Timer("S.get_p()","import CoolProp.State as CPS; S = CPS.State('Air',dict(T = 300, D = 1.1))").timeit(N)*1000000/N

S = CPS.State('Air',dict(T = 300, D = 10.1))

S.speed_test(N)