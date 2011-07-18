from CoolProp.CoolProp import Props

print Props('D','T',300,'P',100,'REFPROP-MIX:R32[0.697615]&R125[0.302385]')
print Props('D','T',300,'P',100,'R410A')

Q=0.5
print Props('H','T',280,'Q',Q,'R290')
D=Props('D','T',280,'Q',Q,'R290')
print D
print Props('H','T',280,'D',D,'R290')

print Props('H','T',280,'Q',Q,'REFPROP-Propane')
D=Props('D','T',280,'Q',Q,'REFPROP-Propane')
print D
print Props('H','T',280,'D',D,'REFPROP-Propane')