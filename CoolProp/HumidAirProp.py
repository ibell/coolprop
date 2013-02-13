from __future__ import absolute_import
from .CoolProp import HAProps as _HAProps, HAProps_Aux as _HAProps_Aux
#from matplotlib import docstring

#@docstring.copy_dedent(_HAProps) #Use docs from HAProps
def HAProps(*args):
    return _HAProps(*args)

#@docstring.copy_dedent(_HAProps_Aux) #Use docs from HAProps_Aux
def HAProps_Aux(*args):
    return _HAProps_Aux(*args)