import CoolProp.CoolProp as CP
import time
        
def runLine(in1=[]):
    out = []
    for i in range(len(in1)):
        out.append(CP.PropsSI(in1[i][0],in1[i][1],in1[i][2],in1[i][3],in1[i][4],in1[i][5]))
    print out 
    
def runNumVal(val=[],num=[]):
    for i in range(len(num[0])):
        out = []
        out.append(CP.PropsSI(val[0],val[1],num[1][i],val[2],num[2][i],fluid))
        time.sleep(0.1)
        out.append(CP.PropsSI(val[1],val[2],num[2][i],val[0],num[0][i],fluid))
        print out 

fluid = "n-Pentane"                                                                                                                                                                                
CP.enable_TTSE_LUT(fluid)                                                                                                                                                                                
CP.PropsSI('H','T',477.8,'D',31.77945717536664,fluid)
CP.disable_TTSE_LUT(fluid)


val = ['H','T','P']
H   = [320683.8718103034, -66337.56629141961, -384680.8070781441]
T   = [477.8            , 418.67414597586514,  300]
P   = [14.7e+05         , 14.7e+05          ,  1.47e+05]
num = [H,T,P]

print 
print "Run with EOS"
CP.disable_TTSE_LUT(fluid)
#CP.set_TTSE_mode(fluid, "TTSE")
CP.set_debug_level(0)
runNumVal(val,num)

print 
print "And now with TTSE"
CP.enable_TTSE_LUT(fluid)
CP.set_TTSE_mode(fluid, "TTSE")
CP.set_debug_level(4)
runNumVal(val,num)

print 
print "And now with BICUBIC"
CP.enable_TTSE_LUT(fluid)
CP.set_TTSE_mode(fluid, "BICUBIC")
CP.set_debug_level(4)
runNumVal(val,num)