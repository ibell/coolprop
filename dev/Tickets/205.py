import CoolProp.CoolProp as CP   

#def runArray(in1=[],val1=[],in2=[],val2=[],fluids=[]):
    #for i in range(len(in1)):
        #print CP.PropsSI('D',in1[i],val1[i],in2[i],val2[i],fluids[i])
        
def runLine(in1=[]):
    out = []
    for i in range(len(in1)):
        out.append(CP.PropsSI(in1[i][0],in1[i][1],in1[i][2],in1[i][3],in1[i][4],in1[i][5]))
    print out 


fluid = "n-Pentane"                                                                                                                                                                                
CP.enable_TTSE_LUT(fluid)                                                                                                                                                                                
CP.PropsSI('H','T',477.8,'D',31.77945717536664,fluid)
CP.disable_TTSE_LUT(fluid)


#val=[]
#val.append(['H','T',477.8,'P',14.7e+05,fluid])
#val.append(['H','Q',  0.0,'P',14.7e+05,fluid])
#val.append(['T','Q',  0.0,'P',14.7e+05,fluid])
#val.append(['H','P',1.4e+05,'T',3.0e+02,fluid])


#CP.disable_TTSE_LUT(fluid)
#CP.set_debug_level(0)
#runLine(val)

#CP.enable_TTSE_LUT(fluid)
#CP.set_TTSE_mode(fluid, "TTSE")
#CP.set_debug_level(0)
#runLine(val)

#CP.enable_TTSE_LUT(fluid)
#CP.set_TTSE_mode(fluid, "BICUBIC")
#CP.set_debug_level(0)
#runLine(val)



val = ['H','T','P']
H   = [320683.8718103034, -66337.56629141961, -384680.8070781441]
T   = [477.8            , 418.67414597586514,  300]
P   = [14.7e+05         , 14.7e+05          ,  1.47e+05]
num = [H,T,P]

for i in range(len(H)):
    out = []
    out.append(CP.PropsSI(val[0],val[1],num[1][i],val[2],num[2][i],fluid))
    out.append(CP.PropsSI(val[1],val[2],num[2][i],val[0],num[0][i],fluid))
    print out 

print 
print "And now with TTSE"
CP.enable_TTSE_LUT(fluid)
CP.set_TTSE_mode(fluid, "TTSE")
CP.set_debug_level(0)
for i in range(len(H)):
    out = []
    out.append(CP.PropsSI(val[0],val[1],num[1][i],val[2],num[2][i],fluid))
    out.append(CP.PropsSI(val[1],val[2],num[2][i],val[0],num[0][i],fluid))
    print out 

print 
print "And now with BICUBIC"
CP.enable_TTSE_LUT(fluid)
CP.set_TTSE_mode(fluid, "BICUBIC")
CP.set_debug_level(0)
for i in range(len(H)):
    out = []
    out.append(CP.PropsSI(val[0],val[1],num[1][i],val[2],num[2][i],fluid))
    out.append(CP.PropsSI(val[1],val[2],num[2][i],val[0],num[0][i],fluid))
    print out 