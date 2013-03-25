import unittest
from CoolProp.CoolProp import Props
import CoolProp

class Props1BadInputParameters(unittest.TestCase):
    """ All fluids, all parameters """
    def testEmptyFluid(self):
        self.assertRaises(ValueError,Props,'','Tcrit')
    def testIntegerFluid(self):
        self.assertRaises(TypeError,Props,1,'Tcrit')
    def testFloatFluid(self):
        self.assertRaises(TypeError,Props,1.0,'Tcrit')
    
    def testEmptyParam(self):
        self.assertRaises(ValueError,Props,'R134a','')
    def testBadParam(self):
        self.assertRaises(ValueError,Props,'R134a','R134a')
        
class Props1CheckInputs(unittest.TestCase):
    def testAllCoolPropPairs(self):
        for fluid in CoolProp.__fluids__:
            for param in ["Ttriple","Tcrit","pcrit","Tmin","molemass","rhocrit","accentric"]:
                val = Props(fluid,param)
                self.assertGreater(val, -1)
                self.assertLess(val, 1e9)
    
if __name__=='__main__':
    unittest.main()