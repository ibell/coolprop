from CoolProp.CoolProp import UseSinglePhaseLUT,Props

UseSinglePhaseLUT(1)
print Props('D','T',300,'P',100,'R410A')
UseSinglePhaseLUT(0)
print Props('D','T',300,'P',100,'R410A')