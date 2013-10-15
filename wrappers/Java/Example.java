// runme.java

public class Example {
    static {
        System.loadLibrary("CoolProp");
    }
    
    public static void main(String argv[]){
        System.out.println(CoolProp.Props("P",'T',300,'Q',0,"R134a"));
        
        double T, h, p, D;
        int SI,KSI;
        System.out.println("CoolProp version: " + CoolProp.get_global_param_string("version") + "\n");
        System.out.println("CoolProp gitrevision: " + CoolProp.get_global_param_string("gitrevision") + "\n");
        System.out.println("CoolProp fluids: " + CoolProp.get_global_param_string("FluidsList") + "\n");

        System.out.println(" " + "\n");
        System.out.println("************ USING EOS *************" + "\n");
        System.out.println(" " + "\n");
        System.out.println("FLUID STATE INDEPENDENT INPUTS" + "\n");
        System.out.println("Critical Density Propane: " + CoolProp.Props1("Propane", "rhocrit") + "kg/m^3" + "\n");
        System.out.println("TWO PHASE INPUTS (Pressure)" + "\n");
        System.out.println("Density of saturated liquid Propane at 101.325 kPa: " + CoolProp.Props("D", 'P', 101.325, 'Q', 0, "Propane") + " kg/m^3" + "\n");
        System.out.println("Density of saturated vapor R290 at 101.325 kPa: " + CoolProp.Props("D", 'P', 101.325, 'Q', 1, "R290") + " kg/m^3" + "\n");
        System.out.println("TWO PHASE INPUTS (Temperature)" + "\n");
        System.out.println("Density of saturated liquid Propane at 300 K: " + CoolProp.Props("D", 'T', 300, 'Q', 0, "Propane") + " kg/m^3" + "\n");
        System.out.println("Density of saturated vapor R290 at 300 K: " + CoolProp.Props("D", 'T', 300, 'Q', 1, "R290") + " kg/m^3" + "\n");
        System.out.println("SINGLE PHASE CYCLE (propane)" + "\n");
        p = CoolProp.Props("P", 'T', 300, 'D', 1, "Propane");
        h = CoolProp.Props("H", 'T', 300, 'D', 1, "Propane");
        System.out.println("T,D -> P,H " + 300 + "," + 1 + " --> " + p + ',' + h + "\n");
        T = CoolProp.Props("T", 'P', p, 'H', h, "Propane");
        D = CoolProp.Props("D", 'P', p, 'H', h, "Propane");
        System.out.println("P,H -> T,D " + p + ',' + h + " --> " + T + ',' + D + "\n");

        System.out.println(" " + "\n");
        System.out.println("************ USING TTSE ***************" + "\n");
        System.out.println(" " + "\n");
        CoolProp.enable_TTSE_LUT("Propane");
        System.out.println("TWO PHASE INPUTS (Pressure)" + "\n");
        System.out.println("Density of saturated liquid Propane at 101.325 kPa: " + CoolProp.Props("D", 'P', 101.325, 'Q', 0, "Propane") + " kg/m^3" + "\n");
        System.out.println("Density of saturated vapor R290 at 101.325 kPa: " + CoolProp.Props("D", 'P', 101.325, 'Q', 1, "R290") + " kg/m^3" + "\n");
        System.out.println("TWO PHASE INPUTS (Temperature)" + "\n");
        System.out.println("Density of saturated liquid Propane at 300 K: " + CoolProp.Props("D", 'T', 300, 'Q', 0, "Propane") + " kg/m^3" + "\n");
        System.out.println("Density of saturated vapor R290 at 300 K: " + CoolProp.Props("D", 'T', 300, 'Q', 1, "R290") + " kg/m^3" + "\n");
        System.out.println("SINGLE PHASE CYCLE (propane)" + "\n");
        p = CoolProp.Props("P", 'T', 300, 'D', 1, "Propane");
        h = CoolProp.Props("H", 'T', 300, 'D', 1, "Propane");
        System.out.println("T,D -> P,H " + 300 + ","+ 1+ " --> " + p + ',' + h + "\n");
        T = CoolProp.Props("T", 'P', p, 'H', h, "Propane");
        D = CoolProp.Props("D", 'P', p, 'H', h, "Propane");
        System.out.println("P,H -> T,D " + p + ',' + h + " --> " + T + ',' + D + "\n");
        CoolProp.disable_TTSE_LUT("Propane");

        try
        {
            System.out.println(" " + "\n");
            System.out.println("************ USING REFPROP ***************" + "\n");
            System.out.println(" " + "\n");
            System.out.println("TWO PHASE INPUTS (Pressure)" + "\n");
            System.out.println("Density of saturated liquid Propane at 101.325 kPa: " + CoolProp.Props("D", 'P', 101.325, 'Q', 0, "Propane") + " kg/m^3" + "\n");
            System.out.println("Density of saturated vapor R290 at 101.325 kPa: " + CoolProp.Props("D", 'P', 101.325, 'Q', 1, "R290") + " kg/m^3" + "\n");
            System.out.println("TWO PHASE INPUTS (Temperature)" + "\n");
            System.out.println("Density of saturated liquid Propane at 300 K: " + CoolProp.Props("D", 'T', 300, 'Q', 0, "Propane") + " kg/m^3" + "\n");
            System.out.println("Density of saturated vapor R290 at 300 K: " + CoolProp.Props("D", 'T', 300, 'Q', 1, "R290") + " kg/m^3" + "\n");
            System.out.println("SINGLE PHASE CYCLE (propane)" + "\n");
            p = CoolProp.Props("P",'T',300,'D',1,"Propane"); 
            h = CoolProp.Props("H",'T',300,'D',1,"Propane");
            System.out.println("T,D -> P,H " + 300 + "," + 1 + " --> " + p + ',' + h + "\n");
            T = CoolProp.Props("T",'P',p,'H',h,"Propane"); 
            D = CoolProp.Props("D",'P',p,'H',h,"Propane");
            System.out.println("P,H -> T,D " + p + ',' + h + " --> " + T + ',' + D + "\n");
        }
        catch(Exception e)
        {
            System.out.println(" " + "\n");
            System.out.println("************ CANT USE REFPROP ************" + "\n");
            System.out.println(" " + "\n");
        }

        System.out.println(" " + "\n");
        System.out.println("************ CHANGE UNIT SYSTEM (default is kSI) *************" + "\n");
        System.out.println(" " + "\n");
        CoolProp.set_standard_unit_system(unit_systems.UNIT_SYSTEM_SI.swigValue());
        System.out.println("Vapor pressure of water at 373.15 K in SI units (Pa): " + CoolProp.Props("P", 'T', 373.15, 'Q', 0, "Water") + "\n");
        CoolProp.set_standard_unit_system(unit_systems.UNIT_SYSTEM_KSI.swigValue());
        System.out.println("Vapor pressure of water at 373.15 K in kSI units (kPa): " + CoolProp.Props("P", 'T', 373.15, 'Q', 0, "Water") + "\n");

        System.out.println(" " + "\n");
        System.out.println("************ BRINES AND SECONDARY WORKING FLUIDS *************" + "\n");
        System.out.println(" " + "\n");
        System.out.println("Density of 50% (mass) ethylene glycol/water at 300 K, 101.325 kPa: " + CoolProp.Props("D", 'T', 300, 'P', 101.325, "EG-50%") + "kg/m^3" + "\n");
        System.out.println("Viscosity of Therminol D12 at 350 K, 101.325 kPa: " + CoolProp.Props("V", 'T', 350, 'P', 101.325, "TD12") + "Pa-s" + "\n");

        System.out.println(" " + "\n");
        System.out.println("************ HUMID AIR PROPERTIES *************" + "\n");
        System.out.println(" " + "\n");
        System.out.println("Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa: " + CoolProp.HAProps("W", "T", 300, "P", 101.325, "R", 0.5) + " kg_w/kg_da" + "\n");
        System.out.println("Relative humidity from last calculation: " + CoolProp.HAProps("R", "T", 300, "P", 101.325, "W", CoolProp.HAProps("W", "T", 300, "P", 101.325, "R", 0.5)) + "(fractional)" + "\n");

        System.out.println("Enter to quit");
    }
}