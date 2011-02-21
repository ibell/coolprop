
from CoolProp import Props
import matplotlib as mpl
import matplotlib.pyplot as pyplot

import numpy as np

errCode=0
errString=''

print Props('H','T',300,'P',500,"R134a")
print Props('S','T',300,'P',500,"R134a")
print Props('D','T',300,'P',500,"R134a")


