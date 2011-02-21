
from CoolProp import Props

errCode=0
errString=''

print Props('H','T',300,'P',500,"R134a")
print Props('S','T',300,'P',500,"R134a")
print Props('D','T',300,'P',500,"R134a")
