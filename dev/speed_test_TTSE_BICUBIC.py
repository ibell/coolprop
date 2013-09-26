import timeit
import CoolProp.CoolProp as CP

def time_check(N, h, p, TTSE = False, mode = 'TTSE'):

    if TTSE:
        if mode =='TTSE':
            setup = "import CoolProp.CoolProp as CP; CP.enable_TTSE_LUT('R245fa'); CP.set_TTSE_mode('R245fa','TTSE'); CP.Props('T','H',250,'P',1000,'R245fa')"
        elif mode =='BICUBIC':
            setup = "import CoolProp.CoolProp as CP; CP.enable_TTSE_LUT('R245fa'); CP.set_TTSE_mode('R245fa','BICUBIC'); CP.Props('T','H',250,'P',1000,'R245fa')"
        else:
            raise ValueError()
    else:
        setup = "import CoolProp.CoolProp as CP"
    
    time = timeit.Timer("CP.Props('D','H',"+str(h)+",'P',"+str(p)+",'R245fa')",setup).timeit(N)/N*1e6
    value = CP.Props('D','H',h,'P',p,'R245fa')
    
    return time, value
    
print 'EOS'
print time_check(10000,250,1000)
print time_check(10000,350,1000)
print time_check(10000,520,1000)
print 'TTSE'
print time_check(10000,250,1000,TTSE = True)
print time_check(10000,350,1000,TTSE = True)
print time_check(10000,520,1000,TTSE = True)
print 'BICUBIC'
print time_check(10000,250,1000,TTSE = True, mode='BICUBIC')
print time_check(10000,350,1000,TTSE = True, mode='BICUBIC')
print time_check(10000,520,1000,TTSE = True, mode='BICUBIC')
