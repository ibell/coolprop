from CoolProp import __HAProps as _HAProps, __HAProps_Aux as _HAProps_Aux
from matplotlib import docstring

@docstring.copy_dedent(_HAProps) #Use docs from HAProps
def HAProps(*args):
    return _HAProps(*args)

@docstring.copy_dedent(_HAProps_Aux) #Use docs from HAProps_Aux
def HAProps_Aux(*args):
    return _HAProps_Aux(*args)