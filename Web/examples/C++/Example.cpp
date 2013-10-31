#include "CoolProp.h"
#include "HumidAirProp.h"
#include "CPState.h"
#include <iostream>
#include <stdlib.h>

int main()
{
    double T, h, p, D;
    std::cout << "CoolProp version: " << get_global_param_string("version") << "\n";
    std::cout << "CoolProp gitrevision: " << get_global_param_string("gitrevision") << "\n";
    std::cout << "CoolProp fluids: " << get_global_param_string("FluidsList") << "\n";

    std::cout <<" " << "\n";
    std::cout <<"************ USING EOS *************" << "\n";
    std::cout <<" " << "\n";
    std::cout <<"FLUID STATE INDEPENDENT INPUTS" << "\n";
    std::cout <<"Critical Density Propane: " << Props1("Propane", "rhocrit") << " kg/m^3\n";
    std::cout <<"TWO PHASE INPUTS (Pressure)" << "\n";
    std::cout <<"Density of saturated liquid Propane at 101.325 kPa: " << Props("D", 'P', 101.325, 'Q', 0, "Propane") << " kg/m^3\n";
    std::cout <<"Density of saturated vapor R290 at 101.325 kPa: " << Props("D", 'P', 101.325, 'Q', 1, "R290") << " kg/m^3\n";
    std::cout <<"TWO PHASE INPUTS (Temperature)" << "\n";
    std::cout <<"Density of saturated liquid Propane at 300 K: " << Props("D", 'T', 300, 'Q', 0, "Propane") << " kg/m^3\n";
    std::cout <<"Density of saturated vapor R290 at 300 K: " << Props("D", 'T', 300, 'Q', 1, "R290") << " kg/m^3\n";
    std::cout <<"SINGLE PHASE CYCLE (propane)\n";
    p = Props("P", 'T', 300, 'D', 1, "Propane");
    h = Props("H", 'T', 300, 'D', 1, "Propane");
    std::cout << "T,D -> P,H " << 300 << "," << 1 << " --> " << p << ',' << h << "\n";
    T = Props("T", 'P', p, 'H', h, "Propane");
    D = Props("D", 'P', p, 'H', h, "Propane");
    std::cout << "P,H -> T,D " << p + ',' << h << " --> " << T << ',' << D << "\n";

    std::cout <<" \n";
    std::cout <<"************ USING TTSE ***************\n";
    std::cout <<" \n";
    enable_TTSE_LUT("Propane");
    std::cout <<"TWO PHASE INPUTS (Pressure)" << "\n";
    std::cout <<"Density of saturated liquid Propane at 101.325 kPa: " << Props("D", 'P', 101.325, 'Q', 0, "Propane") << " kg/m^3\n";
    std::cout <<"Density of saturated vapor R290 at 101.325 kPa: " << Props("D", 'P', 101.325, 'Q', 1, "R290") << " kg/m^3\n";
    std::cout <<"TWO PHASE INPUTS (Temperature)" << "\n";
    std::cout <<"Density of saturated liquid Propane at 300 K: " << Props("D", 'T', 300, 'Q', 0, "Propane") << " kg/m^3\n";
    std::cout <<"Density of saturated vapor R290 at 300 K: " << Props("D", 'T', 300, 'Q', 1, "R290") << " kg/m^3\n";
    std::cout <<"SINGLE PHASE CYCLE (propane)\n";
    p = Props("P", 'T', 300, 'D', 1, "Propane");
    h = Props("H", 'T', 300, 'D', 1, "Propane");
    std::cout << "T,D -> P,H " << 300 << "," << 1 << " --> " << p << ',' << h << "\n";
    T = Props("T", 'P', p, 'H', h, "Propane");
    D = Props("D", 'P', p, 'H', h, "Propane");
    std::cout << "P,H -> T,D " << p + ',' << h << " --> " << T << ',' << D << "\n";
    disable_TTSE_LUT("Propane");

    try
    {
        std::cout <<" \n";
        std::cout <<"************ USING REFPROP ***************\n";
        std::cout <<" \n";
        std::cout <<"TWO PHASE INPUTS (Pressure)" << "\n";
        std::cout <<"Density of saturated liquid Propane at 101.325 kPa: " << Props("D", 'P', 101.325, 'Q', 0, "Propane") << " kg/m^3\n";
        std::cout <<"Density of saturated vapor R290 at 101.325 kPa: " << Props("D", 'P', 101.325, 'Q', 1, "R290") << " kg/m^3\n";
        std::cout <<"TWO PHASE INPUTS (Temperature)" << "\n";
        std::cout <<"Density of saturated liquid Propane at 300 K: " << Props("D", 'T', 300, 'Q', 0, "Propane") << " kg/m^3\n";
        std::cout <<"Density of saturated vapor R290 at 300 K: " << Props("D", 'T', 300, 'Q', 1, "R290") << " kg/m^3\n";
        std::cout <<"SINGLE PHASE CYCLE (propane)\n";
        p = Props("P", 'T', 300, 'D', 1, "Propane");
        h = Props("H", 'T', 300, 'D', 1, "Propane");
        std::cout << "T,D -> P,H " << 300 << "," << 1 << " --> " << p << ',' << h << "\n";
        T = Props("T", 'P', p, 'H', h, "Propane");
        D = Props("D", 'P', p, 'H', h, "Propane");
        std::cout << "P,H -> T,D " << p + ',' << h << " --> " << T << ',' << D << "\n";
    }
    catch (std::exception)
    {
        std::cout <<" \n";
        std::cout <<"************ CANT USE REFPROP ************\n";
        std::cout <<" \n";
    }

    std::cout <<" \n";
    std::cout <<"************ CHANGE UNIT SYSTEM (default is kSI) *************\n";
    std::cout <<" \n";
    set_standard_unit_system(UNIT_SYSTEM_SI);
    std::cout << "Vapor pressure of water at 373.15 K in SI units (Pa): " << Props("P", 'T', 373.15, 'Q', 0, "Water") << "\n";
    set_standard_unit_system(UNIT_SYSTEM_KSI);
    std::cout << "Vapor pressure of water at 373.15 K in kSI units (kPa): " << Props("P", 'T', 373.15, 'Q', 0, "Water") << "\n";

    std::cout <<" \n";
    std::cout <<"************ BRINES AND SECONDARY WORKING FLUIDS *************\n";
    std::cout <<" \n";
    std::cout <<"Density of 50% (mass) ethylene glycol/water at 300 K, 101.325 kPa: " << Props("D", 'T', 300, 'P', 101.325, "EG-50%") << "kg/m^3\n";
    std::cout <<"Viscosity of Therminol D12 at 350 K, 101.325 kPa: " << Props("V", 'T', 350, 'P', 101.325, "TD12") << "Pa-s\n";

    std::cout <<" \n";
    std::cout <<"************ HUMID AIR PROPERTIES *************\n";
    std::cout <<" \n";
    std::cout <<"Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa: " << HAProps("W", "T", 300, "P", 101.325, "R", 0.5) << " kg_w/kg_da\n";
    std::cout <<"Relative humidity from last calculation: " << HAProps("R", "T", 300, "P", 101.325, "W", HAProps("W", "T", 300, "P", 101.325, "R", 0.5)) << "(fractional)\n";
    return 0;
}