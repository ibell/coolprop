# -*- coding: utf-8 -*-
"""
Created on Fri May 31 09:21:27 2013

@author: Belli
"""

import CoolProp.CoolProp as CP
import CoolProp

for fluid in CoolProp.__fluids__:
    try:
        lines = open('C:\\Program Files (x86)\\REFPROP\\fluids\\' + CP.get_REFPROPname(fluid) + '.fld','r').readlines()
        
        for line in lines:
            if line.find('CAS number') > -1:
                CAS_number = line.split('!')[0].strip()
                break
        
        if not CP.get_CAS_code(fluid) == CAS_number:
            print fluid, CP.get_CAS_code(fluid), CAS_number
        else:
            print 'CAS_number is ok for',fluid
        
    except IOError:
        pass#print "didn't find",fluid
    
    