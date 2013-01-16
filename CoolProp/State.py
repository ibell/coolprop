
#We made everything build into one module for simplicity as it makes the code much nicer to compile.  But we want to hide __State in the compiled CoolProp module
from CoolProp import __State as State, set_1phase_LUT_params, debug, LUT, Props

