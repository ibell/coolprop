#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension

import subprocess

subprocess.call(['swig','-python','CoolProp.i'])

CoolProp_module = Extension('_CoolProp',
                           sources=['CoolProp_wrap.c', 'CoolProp.c',
                           'R134a.c','R744.c','R290.c','R410A.c',
                           'Brine.c','R32.c','R717.c','R404A.c','R407C.c',
                           'R507A.c','Argon.c','Nitrogen.c'],
                           )

setup (name = 'CoolProps',
       version = '0.0',
       author      = "Ian Bell",
       description = """Properties of R134a, R744, R410A, R290, R717, R32, R404A, R407C, R507A, Argon, Nitrogen and brines""",
       ext_modules = [CoolProp_module],
       py_modules = ["CoolProp"],
       )



subprocess.call(['swig','-python','FloodProp.i'])

FloodProp_module = Extension('_FloodProp',
                           sources=['FloodProp_wrap.c', 'FloodProp.c',
                           'R134a.c','R744.c','R290.c','R410A.c',
                           'Brine.c','R32.c','R717.c','R404A.c','Nitrogen.c','Argon.c'],
                           )

setup (name = 'FloodProp',
       version = '0.0',
       author      = "Ian Bell",
       description = """Properties of R134a, R744, R410A, R290, R717, R32, R404A and brines""",
       ext_modules = [FloodProp_module],
       py_modules = ["FloodProp"],
       )
