within ;
package CoolProp2Modelica 
import SI = Modelica.SIunits;


  package Common "Package with common definitions"
    type InputChoice = enumeration(
      dT "(d,T) as inputs",
      ph "(p,h) as inputs",
      ps "(p,s) as inputs",
      pT "(p,T) as inputs");
    type InputChoiceMixture = enumeration(
      dTX "(d,T,X) as inputs",
      phX "(p,h,X) as inputs",
      psX "(p,s,X) as inputs",
      pTX "(p,T,X) as inputs");
  end Common;


  package Interfaces
  "Interface definitions for compatibility to other packages like ExternalMedia"
  extends Modelica.Icons.InterfacesPackage;
    partial package ExternalTwoPhaseMedium
    "Generic external two phase medium package"
      extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(
        mediumName = "ExternalMedium",
        singleState = false,
        onePhase = false,
        smoothModel = false,
        fluidConstants = {externalFluidConstants});
    import CoolProp2Modelica.Common.InputChoice;
      constant String libraryName = "UnusableExternalMedium"
      "Name of the external fluid property computation library";
      constant String substanceName = substanceNames[1]
      "Only one substance can be specified";
      redeclare record extends FluidConstants "external fluid constants"
        MolarMass molarMass "molecular mass";
        Temperature criticalTemperature "critical temperature";
        AbsolutePressure criticalPressure "critical pressure";
        MolarVolume criticalMolarVolume "critical molar Volume";
      end FluidConstants;
      constant FluidConstants externalFluidConstants = FluidConstants(
        iupacName=            "unknown",
        casRegistryNumber=    "unknown",
        chemicalFormula=      "unknown",
        structureFormula=     "unknown",
        molarMass=            getMolarMass(),
        criticalTemperature=  getCriticalTemperature(),
        criticalPressure=     getCriticalPressure(),
        criticalMolarVolume=  getCriticalMolarVolume(),
        acentricFactor=       0,
        triplePointTemperature=  280.0,
        triplePointPressure=  500.0,
        meltingPoint=         280,
        normalBoilingPoint=   380.0,
        dipoleMoment=         2.0);
      constant InputChoice inputChoice=InputChoice.ph
      "Default choice of input variables for property computations";
      redeclare replaceable record extends ThermodynamicState
        PrandtlNumber Pr "prandtl number";
        Temperature T "temperature";
        VelocityOfSound a "velocity of sound";
        Modelica.SIunits.CubicExpansionCoefficient beta
        "isobaric expansion coefficient";
        SpecificHeatCapacity cp "specific heat capacity cp";
        SpecificHeatCapacity cv "specific heat capacity cv";
        Density d "density";
        DerDensityByEnthalpy ddhp
        "derivative of density wrt enthalpy at constant pressure";
        DerDensityByPressure ddph
        "derivative of density wrt pressure at constant enthalpy";
        DynamicViscosity eta "dynamic viscosity";
        SpecificEnthalpy h "specific enthalpy";
        Modelica.SIunits.Compressibility kappa "compressibility";
        ThermalConductivity lambda "thermal conductivity";
        AbsolutePressure p "pressure";
        SpecificEntropy s "specific entropy";
      end ThermodynamicState;

      redeclare record extends SaturationProperties
        Temperature Tsat "saturation temperature";
        Real dTp "derivative of Ts wrt pressure";
        DerDensityByPressure ddldp "derivative of dls wrt pressure";
        DerDensityByPressure ddvdp "derivative of dvs wrt pressure";
        DerEnthalpyByPressure dhldp "derivative of hls wrt pressure";
        DerEnthalpyByPressure dhvdp "derivative of hvs wrt pressure";
        Density dl "density at bubble line (for pressure ps)";
        Density dv "density at dew line (for pressure ps)";
        SpecificEnthalpy hl
        "specific enthalpy at bubble line (for pressure ps)";
        SpecificEnthalpy hv "specific enthalpy at dew line (for pressure ps)";
        AbsolutePressure psat "saturation pressure";
        SurfaceTension sigma "surface tension";
        SpecificEntropy sl "specific entropy at bubble line (for pressure ps)";
        SpecificEntropy sv "specific entropy at dew line (for pressure ps)";
      end SaturationProperties;

      redeclare replaceable model extends BaseProperties(
        p(stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoice.ph or
                            basePropertiesInputChoice == InputChoice.pT or
                            basePropertiesInputChoice == InputChoice.ps) then
                                StateSelect.prefer else StateSelect.default),
        T(stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoice.pT or
                           basePropertiesInputChoice == InputChoice.dT) then
                             StateSelect.prefer else StateSelect.default),
        h(stateSelect = if preferredMediumStates and
                           basePropertiesInputChoice == InputChoice.ph then
                             StateSelect.prefer else StateSelect.default),
        d(stateSelect = if preferredMediumStates and
                           basePropertiesInputChoice == InputChoice.dT then
                             StateSelect.prefer else StateSelect.default))
      import CoolProp2Modelica.Common.InputChoice;
        parameter InputChoice basePropertiesInputChoice=inputChoice
        "Choice of input variables for property computations";
        FixedPhase phaseInput
        "Phase input for property computation functions, 2 for two-phase, 1 for one-phase, 0 if not known";
        Integer phaseOutput
        "Phase output for medium, 2 for two-phase, 1 for one-phase";
        SpecificEntropy s(
          stateSelect = if basePropertiesInputChoice == InputChoice.ps then
                           StateSelect.prefer else StateSelect.default)
        "Specific entropy";
        SaturationProperties sat "saturation property record";
      equation
        MM = externalFluidConstants.molarMass;
        R = Modelica.Constants.R/MM;
        if (onePhase or (basePropertiesInputChoice == InputChoice.pT)) then
          phaseInput = 1 "Force one-phase property computation";
        else
          phaseInput = 0 "Unknown phase";
        end if;
        if (basePropertiesInputChoice == InputChoice.ph) then
          // Compute the state record (including the unique ID)
          state = setState_ph(p, h, phaseInput);
          // Modification of the ExternalMedia code to reduce the number of calls:
          // SQ, January 2013:
          d = density(state);
          T = temperature(state);
          s = specificEntropy(state);
        elseif (basePropertiesInputChoice == InputChoice.dT) then
          state = setState_dT(d, T, phaseInput);
          h = specificEnthalpy(state);
          p = pressure(state);
          s = specificEntropy(state);
        elseif (basePropertiesInputChoice == InputChoice.pT) then
          state = setState_pT(p, T, phaseInput);
          d = density(state);
          h = specificEnthalpy(state);
          s = specificEntropy(state);
        elseif (basePropertiesInputChoice == InputChoice.ps) then
          state = setState_ps(p, s, phaseInput);
          d = density(state);
          h = specificEnthalpy(state);
          T = temperature(state);
        end if;
        // Compute the internal energy
        u = h - p/d;
        // Compute the saturation properties record
        sat = setSat_p_state(state);
        // Event generation for phase boundary crossing
         if smoothModel then
          // No event generation
          phaseOutput = state.phase;
         else
           // Event generation at phase boundary crossing
           if basePropertiesInputChoice == InputChoice.ph then
             phaseOutput = if ((h > bubbleEnthalpy(sat) and h < dewEnthalpy(sat)) and
                                p < fluidConstants[1].criticalPressure) then 2 else 1;
           elseif basePropertiesInputChoice == InputChoice.dT then
             phaseOutput = if  ((d < bubbleDensity(sat) and d > dewDensity(sat)) and
                                 T < fluidConstants[1].criticalTemperature) then 2 else 1;
           elseif basePropertiesInputChoice == InputChoice.ps then
             phaseOutput = if ((s > bubbleEntropy(sat) and s < dewEntropy(sat)) and
                                p < fluidConstants[1].criticalPressure) then 2 else 1;
           else
             // basePropertiesInputChoice == pT
             phaseOutput = 1;
           end if;
         end if;
      end BaseProperties;

      redeclare function molarMass "Return the molar mass of the medium"
          input ThermodynamicState state;
          output MolarMass MM "Mixture molar mass";
      algorithm
        MM := fluidConstants[1].molarMass;
      end molarMass;

      replaceable partial function getMolarMass
        output MolarMass MM "molar mass";
        external "C" MM = TwoPhaseMedium_getMolarMass_(mediumName, libraryName, substanceName)
          annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end getMolarMass;

      replaceable partial function getCriticalTemperature
        output Temperature Tc "Critical temperature";
        external "C" Tc = TwoPhaseMedium_getCriticalTemperature_(mediumName, libraryName, substanceName)
          annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end getCriticalTemperature;

      replaceable partial function getCriticalPressure
        output AbsolutePressure pc "Critical temperature";
        external "C" pc = TwoPhaseMedium_getCriticalPressure_(mediumName, libraryName, substanceName)
          annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end getCriticalPressure;

      replaceable partial function getCriticalMolarVolume
        output MolarVolume vc "Critical molar volume";
        external "C" vc = TwoPhaseMedium_getCriticalMolarVolume_(mediumName, libraryName, substanceName)
          annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end getCriticalMolarVolume;

      redeclare replaceable function setState_ph
      "Return thermodynamic state record from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input SpecificEnthalpy h "specific enthalpy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
      external "C" TwoPhaseMedium_setState_ph_(p, h, phase, state, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end setState_ph;

      redeclare replaceable function setState_pT
      "Return thermodynamic state record from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input Temperature T "temperature";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
      external "C" TwoPhaseMedium_setState_pT_(p, T, state, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end setState_pT;

      redeclare replaceable function setState_dT
      "Return thermodynamic state record from d and T"
        extends Modelica.Icons.Function;
        input Density d "density";
        input Temperature T "temperature";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
      external "C" TwoPhaseMedium_setState_dT_(d, T, phase, state, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end setState_dT;

      redeclare replaceable function setState_ps
      "Return thermodynamic state record from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input SpecificEntropy s "specific entropy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
      external "C" TwoPhaseMedium_setState_ps_(p, s, phase, state, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end setState_ps;

      replaceable partial function setSat_p_state_dummy
      "Return dummy values for the saturation properties from the state"
        extends Modelica.Icons.Function;
        input ThermodynamicState state;
        output SaturationProperties sat "saturation property record";
        // Standard definition
      algorithm
        sat.Tsat :=0;
        sat.dTp  :=0;
        sat.ddldp :=0;
        sat.ddvdp  :=0;
        sat.dhldp  :=0;
        sat.dhvdp  :=0;
        sat.dl  :=0;
        sat.dv  :=0;
        sat.hl  :=0;
        sat.hv  :=0;
        sat.psat  :=0;
        sat.sigma  :=0;
        sat.sl  :=0;
        sat.sv  :=0;
        annotation(Inline = true);
      end setSat_p_state_dummy;

      replaceable partial function setSat_p_state
      "Return saturation properties from the state"
        extends Modelica.Icons.Function;
        input ThermodynamicState state;
        output SaturationProperties sat "saturation property record";
        // Standard definition
      algorithm
        sat:=setSat_p(state.p);
        //Redeclare this function for more efficient implementations avoiding the repeated computation of saturation properties
      /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setSat_p_state_(state, sat)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end setSat_p_state;

      redeclare function extends setState_phX
      algorithm
        // The composition is an empty vector
        state :=setState_ph(p, h, phase);
      end setState_phX;

      redeclare function extends setState_pTX
      algorithm
        // The composition is an empty vector
        state :=setState_pT(p, T, phase);
      end setState_pTX;

      redeclare function extends setState_dTX
      algorithm
        // The composition is an empty vector
        state :=setState_dT(d, T, phase);
      end setState_dTX;

      redeclare function extends setState_psX
      algorithm
        // The composition is an empty vector
        state :=setState_ps(p, s, phase);
      end setState_psX;

      redeclare replaceable function density_ph "Return density from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_ph(p, h, phase));
        //removal of the reference to the density derivative.
        // SQ, January 2013:
        //     annotation(derivative(noDerivative = phase) = density_ph_der,
        //                Inline = true);
      annotation(Inline = true);
      end density_ph;

      redeclare replaceable function temperature_ph
      "Return temperature from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_ph(p, h, phase));
        annotation(Inline = true);
      end temperature_ph;

      replaceable function specificEntropy_ph
      "Return specific entropy from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEntropy s "specific entropy";
      algorithm
        s := specificEntropy(setState_ph(p, h, phase));
        annotation(Inline = true);
      end specificEntropy_ph;

      redeclare replaceable function density_pT "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_pT(p, T, phase));
        //  To be implemented:
        //     annotation(derivative(noDerivative = phase) = density_pT_der,
        //                Inline = true);
      end density_pT;

      replaceable partial function density_pT_der
      "Total derivative of density_pT"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase
        "2 for two-phase, 1 for one-phase, 0 if not known";
        input Real p_der;
        input Real T_der;
        output Real d_der;
        // To be implemented
        annotation(Inline = true);
      end density_pT_der;

      redeclare replaceable function specificEnthalpy_pT
      "Return specific enthalpy from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_pT(p, T, phase));
        annotation(Inline = true);
      end specificEnthalpy_pT;

      redeclare replaceable function pressure_dT "Return pressure from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output AbsolutePressure p "Pressure";
      algorithm
        p := pressure(setState_dT(d, T, phase));
        annotation(Inline = true);
      end pressure_dT;

      redeclare replaceable function specificEnthalpy_dT
      "Return specific enthalpy from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_dT(d, T, phase));
        annotation(Inline = true);
      end specificEnthalpy_dT;

      redeclare replaceable function density_ps "Return density from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_ps(p, s, phase));
        // To be implemented:
        //     annotation(derivative(noDerivative = phase) = density_ps_der,
        //                Inline = true);
      end density_ps;

      replaceable partial function density_ps_der
      "Total derivative of density_ps"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase
        "2 for two-phase, 1 for one-phase, 0 if not known";
        input Real p_der;
        input Real h_der;
        output Real d_der;
        // To be implemented
        annotation(Inline = true);
      end density_ps_der;

      redeclare replaceable function temperature_ps
      "Return temperature from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_ps(p, s, phase));
        annotation(Inline = true);
      end temperature_ps;

      redeclare replaceable function specificEnthalpy_ps
      "Return specific enthalpy from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase = 0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_ps(p,s, phase));
        annotation(Inline = true);
      end specificEnthalpy_ps;

      redeclare function extends prandtlNumber
        /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_prandtlNumber_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end prandtlNumber;

      redeclare replaceable function extends temperature
      "Return temperature from state"
        // Standard definition
      algorithm
        T := state.T;
        /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_temperature_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end temperature;

      redeclare replaceable function extends velocityOfSound
      "Return velocity of sound from state"
        // Standard definition
      algorithm
        a := state.a;
        /*  // If special definition in "C"
  external "C" a=  TwoPhaseMedium_velocityOfSound_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end velocityOfSound;

      redeclare replaceable function extends isobaricExpansionCoefficient
      "Return isobaric expansion coefficient from state"
        // Standard definition
      algorithm
        beta := state.beta;
        /*  // If special definition in "C"
  external "C" beta=  TwoPhaseMedium_isobaricExpansionCoefficient_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end isobaricExpansionCoefficient;

      redeclare replaceable function extends specificHeatCapacityCp
      "Return specific heat capacity cp from state"
        // Standard definition
      algorithm
        cp := state.cp;
        /*  // If special definition in "C"
  external "C" cp=  TwoPhaseMedium_specificHeatCapacityCp_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end specificHeatCapacityCp;

      redeclare replaceable function extends specificHeatCapacityCv
      "Return specific heat capacity cv from state"
        // Standard definition
      algorithm
        cv := state.cv;
        /*  // If special definition in "C"
  external "C" cv=  TwoPhaseMedium_specificHeatCapacityCv_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end specificHeatCapacityCv;

      redeclare replaceable function extends density
      "Return density from state"
        // Standard definition
      algorithm
        d := state.d;
        /*  // If special definition in "C"
  external "C" d=  TwoPhaseMedium_density_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end density;

      redeclare replaceable function extends density_derh_p
      "Return derivative of density wrt enthalpy at constant pressure from state"
        // Standard definition
      algorithm
        ddhp := state.ddhp;
        /*  // If special definition in "C"
  external "C" ddhp=  TwoPhaseMedium_density_derh_p_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end density_derh_p;

      redeclare replaceable function extends density_derp_h
      "Return derivative of density wrt pressure at constant enthalpy from state"
        // Standard definition
      algorithm
        ddph := state.ddph;
        /*  // If special definition in "C"
  external "C" ddph=  TwoPhaseMedium_density_derp_h_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end density_derp_h;

      redeclare replaceable function extends dynamicViscosity
      "Return dynamic viscosity from state"
        // Standard definition
      algorithm
        eta := state.eta;
        /*  // If special definition in "C"
  external "C" eta=  TwoPhaseMedium_dynamicViscosity_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dynamicViscosity;

      redeclare replaceable function extends specificEnthalpy
      "Return specific enthalpy from state"
        // Standard definition
      algorithm
        h := state.h;
        /*  // If special definition in "C"
  external "C" h=  TwoPhaseMedium_specificEnthalpy_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end specificEnthalpy;

      redeclare replaceable function extends isothermalCompressibility
      "Return isothermal compressibility from state"
        // Standard definition
      algorithm
        kappa := state.kappa;
        /*  // If special definition in "C"
  external "C" kappa=  TwoPhaseMedium_isothermalCompressibility_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end isothermalCompressibility;

      redeclare replaceable function extends thermalConductivity
      "Return thermal conductivity from state"
        // Standard definition
      algorithm
        lambda := state.lambda;
        /*  // If special definition in "C"
  external "C" lambda=  TwoPhaseMedium_thermalConductivity_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end thermalConductivity;

      redeclare replaceable function extends pressure
      "Return pressure from state"
        // Standard definition
      algorithm
        p := state.p;
        /*  // If special definition in "C"
  external "C" p=  TwoPhaseMedium_pressure_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end pressure;

      redeclare replaceable function extends specificEntropy
      "Return specific entropy from state"
        // Standard definition
      algorithm
        s := state.s;
        /*  // If special definition in "C"
    external "C" s=  TwoPhaseMedium_specificEntropy_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end specificEntropy;

      replaceable function density_ph_der "Total derivative of density_ph"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase
        "2 for two-phase, 1 for one-phase, 0 if not known";
        input Real p_der;
        input Real h_der;
        output Real d_der;
        // Standard definition
      algorithm
        d_der:=density_derp_h(setState_ph(p, h))*p_der + density_derh_p(setState_ph(
          p, h))*h_der;
        /*  // If special definition in "C"
  external "C" d_der=  TwoPhaseMedium_density_ph_der_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end density_ph_der;

      redeclare replaceable function extends isentropicEnthalpy
      external "C" h_is = TwoPhaseMedium_isentropicEnthalpy_(p_downstream, refState,
       mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end isentropicEnthalpy;

      redeclare replaceable function setSat_p
      "Return saturation properties from p"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        output SaturationProperties sat "saturation property record";
      external "C" TwoPhaseMedium_setSat_p_(p, sat, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end setSat_p;

      redeclare replaceable function setSat_T
      "Return saturation properties from p"
        extends Modelica.Icons.Function;
        input Temperature T "temperature";
        output SaturationProperties sat "saturation property record";
      external "C" TwoPhaseMedium_setSat_T_(T, sat, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end setSat_T;

      redeclare replaceable function extends setBubbleState
      "set the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase =  1 "phase: default is one phase";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
      algorithm
        state :=setState_ph(sat.psat, sat.hl, phase);
        /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setBubbleState_(sat, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end setBubbleState;

      redeclare replaceable function extends setDewState
      "set the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase =  1 "phase: default is one phase";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
      algorithm
        state :=setState_ph(sat.psat, sat.hv, phase);
        /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setDewState_(sat, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end setDewState;

      redeclare replaceable function extends saturationTemperature
        // Standard definition
      algorithm
        T :=saturationTemperature_sat(setSat_p(p));
        /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_saturationTemperature_(p, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end saturationTemperature;

      redeclare function extends saturationTemperature_sat
        annotation(Inline = true);
      end saturationTemperature_sat;

      redeclare replaceable function extends saturationTemperature_derp "Returns derivative of saturation temperature w.r.t.. pressureBeing this function inefficient, it is strongly recommended to use saturationTemperature_derp_sat
     and never use saturationTemperature_derp directly"
      external "C" dTp = TwoPhaseMedium_saturationTemperature_derp_(p, mediumName, libraryName, substanceName)
        annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end saturationTemperature_derp;

      redeclare replaceable function saturationTemperature_derp_sat
      "Returns derivative of saturation temperature w.r.t.. pressure"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output Real dTp "derivative of saturation temperature w.r.t. pressure";
        // Standard definition
      algorithm
        dTp := sat.dTp;
        /*  // If special definition in "C"
  external "C" dTp=  TwoPhaseMedium_saturationTemperature_derp_sat_(sat.psat, sat.Tsat, sat.uniqueID, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end saturationTemperature_derp_sat;

      redeclare replaceable function extends dBubbleDensity_dPressure
      "Returns bubble point density derivative"
        // Standard definition
      algorithm
        ddldp := sat.ddldp;
        /*  // If special definition in "C"
  external "C" ddldp=  TwoPhaseMedium_dBubbleDensity_dPressure_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dBubbleDensity_dPressure;

      redeclare replaceable function extends dDewDensity_dPressure
      "Returns dew point density derivative"
        // Standard definition
      algorithm
        ddvdp := sat.ddvdp;
        /*  // If special definition in "C"
  external "C" ddvdp=  TwoPhaseMedium_dDewDensity_dPressure_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dDewDensity_dPressure;

      redeclare replaceable function extends dBubbleEnthalpy_dPressure
      "Returns bubble point specific enthalpy derivative"
        // Standard definition
      algorithm
        dhldp := sat.dhldp;
        /*  // If special definition in "C"
  external "C" dhldp=  TwoPhaseMedium_dBubbleEnthalpy_dPressure_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dBubbleEnthalpy_dPressure;

      redeclare replaceable function extends dDewEnthalpy_dPressure
      "Returns dew point specific enthalpy derivative"
        // Standard definition
      algorithm
        dhvdp := sat.dhvdp;
        /*  // If special definition in "C"
  external "C" dhvdp=  TwoPhaseMedium_dDewEnthalpy_dPressure_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dDewEnthalpy_dPressure;

      redeclare replaceable function extends bubbleDensity
      "Returns bubble point density"
        // Standard definition
      algorithm
        dl := sat.dl;
        /*  // If special definition in "C"
  external "C" dl=  TwoPhaseMedium_bubbleDensity_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end bubbleDensity;

      redeclare replaceable function extends dewDensity
      "Returns dew point density"
        // Standard definition
      algorithm
        dv := sat.dv;
        /*  // If special definition in "C"
  external "C" dv=  TwoPhaseMedium_dewDensity_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dewDensity;

      redeclare replaceable function extends bubbleEnthalpy
      "Returns bubble point specific enthalpy"
        // Standard definition
      algorithm
        hl := sat.hl;
        /*  // If special definition in "C"
  external "C" hl=  TwoPhaseMedium_bubbleEnthalpy_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end bubbleEnthalpy;

      redeclare replaceable function extends dewEnthalpy
      "Returns dew point specific enthalpy"
        // Standard definition
      algorithm
        hv := sat.hv;
        /*  // If special definition in "C"
  external "C" hv=  TwoPhaseMedium_dewEnthalpy_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dewEnthalpy;

      redeclare replaceable function extends saturationPressure
        // Standard definition
      algorithm
        p :=saturationPressure_sat(setSat_T(T));
        /*  // If special definition in "C"
  external "C" p=  TwoPhaseMedium_saturationPressure_(T, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = false,
                   LateInline = true,
                   derivative = saturationPressure_der);
      end saturationPressure;

      redeclare function extends saturationPressure_sat
        annotation(Inline = true);
      end saturationPressure_sat;

      redeclare replaceable function extends surfaceTension
      "Returns surface tension sigma in the two phase region"
        //Standard definition
      algorithm
        sigma := sat.sigma;
        /*  //If special definition in "C"
  external "C" sigma=  TwoPhaseMedium_surfaceTension_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end surfaceTension;

      redeclare replaceable function extends bubbleEntropy
      "Returns bubble point specific entropy"
        //Standard definition
      algorithm
        sl := specificEntropy(setBubbleState(sat));
        /*  //If special definition in "C"
  external "C" sl=  TwoPhaseMedium_bubbleEntropy_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end bubbleEntropy;

      redeclare replaceable function extends dewEntropy
      "Returns dew point specific entropy"
        //Standard definition
      algorithm
        sv := specificEntropy(setDewState(sat));
        /*  //If special definition in "C"
  external "C" sv=  TwoPhaseMedium_dewEntropy_(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end dewEntropy;

      function saturationPressure_der
      "Return saturation pressure time derivative"
        extends Modelica.Icons.Function;
        input Temperature T "temperature";
        input Real T_der "Temperature derivative";
        output Real p_der "saturation pressure derivative";
        // Standard definition
      algorithm
        p_der :=T_der/saturationTemperature_derp_sat(setSat_T(T));
        annotation(Inline = true);
      end saturationPressure_der;
    end ExternalTwoPhaseMedium;

    partial package ExternalTwoPhaseMixture
    "Template class for two phase medium of a mixture of substances"
      extends Modelica.Media.Interfaces.PartialMixtureMedium(ThermoStates=inputChoice);
      constant Boolean smoothModel = false
      "true if the (derived) model should not generate state events";
      constant Boolean onePhase =    false
      "true if the (derived) model should never be called with two-phase inputs";

       //constant Real minValue = 1e-6;

       constant FluidConstants mixtureFluidConstants(
         iupacName =              "mixture",
         casRegistryNumber =      "mixture",
         chemicalFormula =        "mixture",
         structureFormula =       "mixture",
         molarMass =              0.001,
         criticalTemperature =    0,
         criticalPressure =       0,
         criticalMolarVolume =    0,
         acentricFactor =         0,
         triplePointTemperature = 0,
         triplePointPressure =    0,
         meltingPoint =           0,
         normalBoilingPoint =     0,
         dipoleMoment =           0);
       // These values are not constant for mixtures.
       // We set them to dummy values.

      import InputChoice = CoolProp2Modelica.Common.InputChoiceMixture;
      constant InputChoice inputChoice=InputChoice.phX
      "Default choice of input variables for property computations";

      type FixedPhase = Integer(min=0,max=2)
      "phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g. interactive use";
      record FluidLimits "validity limits for fluid model"
        extends Modelica.Icons.Record;
        Temperature TMIN "minimum temperature";
        Temperature TMAX "maximum temperature";
        Density DMIN "minimum density";
        Density DMAX "maximum density";
        AbsolutePressure PMIN "minimum pressure";
        AbsolutePressure PMAX "maximum pressure";
        SpecificEnthalpy HMIN "minimum enthalpy";
        SpecificEnthalpy HMAX "maximum enthalpy";
        SpecificEntropy SMIN "minimum entropy";
        SpecificEntropy SMAX "maximum entropy";
        annotation(Documentation(
            info="<html>
          <p>The minimum pressure mostly applies to the liquid state only.
          The minimum density is also arbitrary, but is reasonable for techical
          applications to limit iterations in non-linear systems. The limits in
          enthalpy and entropy are used as safeguards in inverse iterations.</p>
          </html>"));
      end FluidLimits;

    redeclare replaceable record extends ThermodynamicState
      "Thermodynamic state of two phase medium"
        MolarMass molarMass "Molar mass of bulk mixture";
        FixedPhase phase(min=0, max=2)
        "phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g. interactive use";
        PrandtlNumber Pr "prandtl number";
        Temperature T "temperature";
        VelocityOfSound a(min=1e-8) "velocity of sound";
        Modelica.SIunits.CubicExpansionCoefficient beta
        "isobaric expansion coefficient";
        SpecificHeatCapacity cp(min=1e-8) "specific heat capacity cp";
        SpecificHeatCapacity cv(min=1e-8) "specific heat capacity cv";
        Density d(min=1e-8) "density";
        DerDensityByEnthalpy ddhp
        "derivative of density wrt enthalpy at constant pressure and composition";
        DerDensityByPressure ddph
        "derivative of density wrt pressure at constant enthalpy and composition";
        DynamicViscosity eta "dynamic viscosity";
        SpecificEnthalpy h "specific enthalpy";
        Modelica.SIunits.Compressibility kappa "compressibility";
        ThermalConductivity lambda "thermal conductivity";
        AbsolutePressure p "pressure";
        SpecificEntropy s "specific entropy";
        MassFraction X[nX] "Mass fraction of components in kg/kg";
    end ThermodynamicState;

      replaceable record SaturationProperties
      "Saturation properties of two phase medium"
        extends Modelica.Icons.Record;
        Temperature Tsat(min=1e-8) "saturation temperature";
        Real dTp "derivative of Ts wrt pressure";
        DerDensityByPressure ddldp "derivative of dls wrt pressure";
        DerDensityByPressure ddvdp "derivative of dvs wrt pressure";
        DerEnthalpyByPressure dhldp "derivative of hls wrt pressure";
        DerEnthalpyByPressure dhvdp "derivative of hvs wrt pressure";
        Density dl "density at bubble line (for pressure ps)";
        Density dv "density at dew line (for pressure ps)";
        SpecificEnthalpy hl
        "specific enthalpy at bubble line (for pressure ps)";
        SpecificEnthalpy hv "specific enthalpy at dew line (for pressure ps)";
        AbsolutePressure psat(min=1e-8) "saturation pressure";
        SurfaceTension sigma "surface tension";
        SpecificEntropy sl "specific entropy at bubble line (for pressure ps)";
        SpecificEntropy sv "specific entropy at dew line (for pressure ps)";
        MassFraction X[nX] "Bulk mass fractions";
        MassFraction Xl[nX] "Mass fractions of liquid phase";
        MassFraction Xv[nX] "Mass fractions of gaseous phase";
      end SaturationProperties;

    redeclare replaceable model extends BaseProperties(
        p(min=1e-8,stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoice.phX or
                            basePropertiesInputChoice == InputChoice.pTX or
                            basePropertiesInputChoice == InputChoice.psX) then
                                StateSelect.prefer else StateSelect.default),
        T(min=1e-8,stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoice.pTX or
                           basePropertiesInputChoice == InputChoice.dTX) then
                             StateSelect.prefer else StateSelect.default),
        h(stateSelect = if preferredMediumStates and
                           basePropertiesInputChoice == InputChoice.phX then
                             StateSelect.prefer else StateSelect.default),
        d(min=1e-8,stateSelect = if preferredMediumStates and
                           basePropertiesInputChoice == InputChoice.dTX then
                             StateSelect.prefer else StateSelect.default))
        import InputChoice = CoolProp2Modelica.Common.InputChoiceMixture;
        parameter InputChoice basePropertiesInputChoice = inputChoice
        "Choice of input variables for property computations";
        FixedPhase phaseInput
        "Phase input for property computation functions, 2 for two-phase, 1 for one-phase, 0 if not known";
        Integer phaseOutput
        "Phase output for medium, 2 for two-phase, 1 for one-phase";
        SpecificEntropy s(
          stateSelect = if basePropertiesInputChoice == InputChoice.psX then
                           StateSelect.prefer else StateSelect.default)
        "Specific entropy";
        SaturationProperties sat "saturation property record";
    equation
        MM = state.molarMass;
        R = Modelica.Constants.R/max(1e-8,MM);
        if (onePhase or (basePropertiesInputChoice == InputChoice.pTX)) then
          phaseInput = 1 "Force one-phase property computation";
        else
          phaseInput = 0 "Unknown phase";
        end if;
        if (basePropertiesInputChoice == InputChoice.phX) then
          // Compute the state record (including the unique ID)
          state =
            setState_phX(
            p,
            h,
            X,
            phaseInput);
          d = density(state);
          s = specificEntropy(state);
          T = temperature(state);
        elseif (basePropertiesInputChoice == InputChoice.dTX) then
          state =
            setState_dTX(
            d,
            T,
            X,
            phaseInput);
          h = specificEnthalpy(state);
          p = pressure(state);
          s = specificEntropy(state);
        elseif (basePropertiesInputChoice == InputChoice.pTX) then
          state =
            setState_pTX(
            p,
            T,
            X,
            phaseInput);
          d = density(state);
          h = specificEnthalpy(state);
          s = specificEntropy(state);
        elseif (basePropertiesInputChoice == InputChoice.psX) then
          state =
            setState_psX(
            p,
            s,
            X,
            phaseInput);
          d = density(state);
          h = specificEnthalpy(state);
          T = temperature(state);
        end if;
        // Compute the internal energy
        u = h - p/max(1e-8,d);
        // Compute the saturation properties record
        // JW 4/2013:
        // This part might need some more attention in the future!
        sat = setSat_pX(state.p,state.X);
        // Event generation for phase boundary crossing
        if smoothModel then
          // No event generation
          phaseOutput = state.phase;
        else
          // Event generation at phase boundary crossing
          if basePropertiesInputChoice == InputChoice.phX then
            phaseOutput = if ((h > bubbleEnthalpy(sat) and h < dewEnthalpy(sat)) and
                               p < fluidConstants[1].criticalPressure) then 2 else 1;
          elseif basePropertiesInputChoice == InputChoice.dTX then
            phaseOutput = if  ((d < bubbleDensity(sat) and d > dewDensity(sat)) and
                                T < fluidConstants[1].criticalTemperature) then 2 else 1;
          elseif basePropertiesInputChoice == InputChoice.psX then
            phaseOutput = if ((s > bubbleEntropy(sat) and s < dewEntropy(sat)) and
                               p < fluidConstants[1].criticalPressure) then 2 else 1;
          else
            // basePropertiesInputChoice == pTX
            phaseOutput = 1;
          end if;
        end if;
    end BaseProperties;

    //   redeclare replaceable partial model extends BaseProperties
    //     "Base properties (p, d, T, h, s, u, R, MM, sat) of two phase medium"
    //   //  Temperature T(start=300);
    //     Modelica.SIunits.SpecificEntropy s;
    //     SaturationProperties sat "Saturation properties at the medium pressure";
    //     annotation(Documentation(info="<html></html>"));
    //   end BaseProperties;

      redeclare function molarMass "Return the molar mass of the medium"
          input ThermodynamicState state;
          output MolarMass MM "Mixture molar mass";
      algorithm
        MM := state.molarMass;
      end molarMass;

    replaceable partial function getSubstanceName
      "Create a string containing the composition to pass to the external library."
       extends Modelica.Icons.Function;
       output MolarMass MM "Molar mass of the mixture";
       external "C" MM = TwoPhaseMedium_getMolarMass_(mediumName, libraryName, substanceName)
          annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
    end getSubstanceName;

      replaceable partial function getMolarMass
       extends Modelica.Icons.Function;
       output MolarMass MM "Molar mass of the mixture";
       external "C" MM = TwoPhaseMedium_getMolarMass_(mediumName, libraryName, substanceName)
          annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
      end getMolarMass;

      replaceable partial function getCriticalTemperature
        extends Modelica.Icons.Function;
        output Temperature Tcrit "Molar mass of the mixture";
        annotation (Documentation(info="<html></html>"));
      end getCriticalTemperature;

      replaceable partial function getCriticalPressure
       extends Modelica.Icons.Function;
       output AbsolutePressure Pcrit "Molar mass of the mixture";
        annotation (Documentation(info="<html></html>"));
      end getCriticalPressure;

      replaceable partial function getCriticalMolarVolume
        extends Modelica.Icons.Function;
        output MolarVolume v "Molar mass of the mixture";
        annotation (Documentation(info="<html></html>"));
      end getCriticalMolarVolume;

      replaceable partial function setDewState
      "Return the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase(min = 1, max = 2) =  1
        "phase: default is one phase";
        output ThermodynamicState state "complete thermodynamic state info";
        annotation(Documentation(info="<html></html>"));
      end setDewState;

      replaceable partial function setBubbleState
      "Return the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase(min = 1, max = 2) =  1
        "phase: default is one phase";
        output ThermodynamicState state "complete thermodynamic state info";
        annotation(Documentation(info="<html></html>"));
      end setBubbleState;

      redeclare replaceable partial function extends setState_dTX
      "Return thermodynamic state as function of d, T and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        annotation(Documentation(info="<html></html>"));
      end setState_dTX;

      redeclare replaceable partial function extends setState_phX
      "Return thermodynamic state as function of p, h and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        annotation(Documentation(info="<html></html>"));
      end setState_phX;

      redeclare replaceable partial function extends setState_psX
      "Return thermodynamic state as function of p, s and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        annotation(Documentation(info="<html></html>"));
      end setState_psX;

      redeclare replaceable partial function extends setState_pTX
      "Return thermodynamic state as function of p, T and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        annotation(Documentation(info="<html></html>"));
      end setState_pTX;

      replaceable function setSat_TX
      "Return saturation property record from temperature"
        extends Modelica.Icons.Function;
        input Temperature T "temperature";
        input MassFraction X[nX] "Mass fractions";
        output SaturationProperties sat "saturation property record";
      algorithm
        sat.Tsat := T;
        sat.psat := saturationPressure(T,X);
        sat.X := X;
        annotation(Documentation(info="<html></html>"));
      end setSat_TX;

      replaceable function setSat_pX
      "Return saturation property record from pressure"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input MassFraction X[nX] "Mass fractions";
        output SaturationProperties sat "saturation property record";
      algorithm
        sat.psat := p;
        sat.Tsat := saturationTemperature(p,X);
        sat.X := X;
        annotation(Documentation(info="<html></html>"));
      end setSat_pX;

    /*
Functions to obtain fluid properties from the currently active state.
*/

      replaceable partial function bubbleEnthalpy
      "Return bubble point specific enthalpy"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output SpecificEnthalpy hl "boiling curve specific enthalpy";
        annotation(Documentation(info="<html></html>"));
      end bubbleEnthalpy;

      replaceable partial function dewEnthalpy
      "Return dew point specific enthalpy"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output SpecificEnthalpy hv "dew curve specific enthalpy";
        annotation(Documentation(info="<html></html>"));
      end dewEnthalpy;

      replaceable partial function bubbleEntropy
      "Return bubble point specific entropy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output SpecificEntropy sl "boiling curve specific entropy";
        annotation(Documentation(info="<html></html>"));
      end bubbleEntropy;

      replaceable partial function dewEntropy
      "Return dew point specific entropy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output SpecificEntropy sv "dew curve specific entropy";
        annotation(Documentation(info="<html></html>"));
      end dewEntropy;

      replaceable partial function bubbleDensity "Return bubble point density"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output Density dl "boiling curve density";
        annotation(Documentation(info="<html></html>"));
      end bubbleDensity;

      replaceable partial function dewDensity "Return dew point density"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output Density dv "dew curve density";
        annotation(Documentation(info="<html></html>"));
      end dewDensity;

      replaceable partial function saturationPressure
      "Return saturation pressure"
          extends Modelica.Icons.Function;
          input Temperature T "temperature";
          input MassFraction X[:]={1} "fluid composition as mass fractions";
          output AbsolutePressure p "saturation pressure";
        annotation(Documentation(info="<html></html>"));
      end saturationPressure;

      function saturationPressure_der
      "Return saturation pressure time derivative"
        extends Modelica.Icons.Function;
        input Temperature T "temperature";
        input MassFraction X[:]={1} "fluid composition as mass fractions";
        input Real T_der "Temperature derivative";
        output Real p_der "saturation pressure derivative";
        // Standard definition
      algorithm
        p_der :=T_der/saturationTemperature_derp_sat(setSat_TX(T,X));
        annotation(Inline = true);
      end saturationPressure_der;

      replaceable function saturationPressure_sat
      "Return saturation temperature"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output AbsolutePressure p "saturation pressure";
      algorithm
          p := sat.psat;
        annotation(Documentation(info="<html></html>"));
      end saturationPressure_sat;

      replaceable partial function saturationTemperature
      "Return saturation temperature"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "pressure";
          input MassFraction X[:]={1} "fluid composition as mass fractions";
          output Temperature T "saturation temperature";
        annotation(Documentation(info="<html></html>"));
      end saturationTemperature;

      replaceable function saturationTemperature_sat
      "Return saturation temperature"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output Temperature T "saturation temperature";
      algorithm
          T := sat.Tsat;
        annotation(Documentation(info="<html></html>"));
      end saturationTemperature_sat;

      replaceable partial function saturationTemperature_derp
      "Return derivative of saturation temperature w.r.t. pressure"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "pressure";
          output Real dTp
        "derivative of saturation temperature w.r.t. pressure";
        annotation(Documentation(info="<html></html>"));
      end saturationTemperature_derp;

      replaceable function saturationTemperature_derp_sat
      "Return derivative of saturation temperature w.r.t. pressure"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output Real dTp
        "derivative of saturation temperature w.r.t. pressure";
      algorithm
          dTp := saturationTemperature_derp(sat.psat);
        annotation(Documentation(info="<html></html>"));
      end saturationTemperature_derp_sat;

      /*  redeclare replaceable partial function extends molarMass 
    "Return the molar mass of the medium"
    algorithm 
      MM := fluidConstants[1].molarMass;
    end molarMass;*/

      replaceable partial function dBubbleDensity_dPressure
      "Return bubble point density derivative"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output DerDensityByPressure ddldp "boiling curve density derivative";
        annotation(Documentation(info="<html></html>"));
      end dBubbleDensity_dPressure;

      replaceable partial function dDewDensity_dPressure
      "Return dew point density derivative"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output DerDensityByPressure ddvdp
        "saturated steam density derivative";
        annotation(Documentation(info="<html></html>"));
      end dDewDensity_dPressure;

      replaceable partial function dBubbleEnthalpy_dPressure
      "Return bubble point specific enthalpy derivative"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output DerEnthalpyByPressure dhldp
        "boiling curve specific enthalpy derivative";
        annotation(Documentation(info="<html></html>"));
      end dBubbleEnthalpy_dPressure;

      replaceable partial function dDewEnthalpy_dPressure
      "Return dew point specific enthalpy derivative"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output DerEnthalpyByPressure dhvdp
        "saturated steam specific enthalpy derivative";
        annotation(Documentation(info="<html></html>"));
      end dDewEnthalpy_dPressure;

       redeclare replaceable function density_phX
      "Return density from p, h, and X or Xi"
           extends Modelica.Icons.Function;
           input AbsolutePressure p "Pressure";
           input SpecificEnthalpy h "Specific enthalpy";
           input MassFraction X[nX] "Mass fractions";
           input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
           output Density d "density";
       algorithm
         d := density(
           setState_phX(
             p,
             h,
             X,
             phase));
         annotation(Documentation(info="<html></html>"));
       end density_phX;

      redeclare replaceable function density_psX
      "Return density from p, s, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction X[nX] "Mass fractions";
          input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
          output Density d "Density";
      algorithm
        d := density(
          setState_psX(
            p,
            s,
            X,
            phase));
        annotation(Documentation(info="<html></html>"));
      end density_psX;

      redeclare replaceable function density_pTX
      "Return density from p, T, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction X[nX] "Mass fractions";
          input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
          output Density d "Density";
      algorithm
        d := density(
          setState_pTX(
            p,
            T,
            X,
            phase));
        annotation(Documentation(info="<html></html>"));
      end density_pTX;

    replaceable function specificEnthalpy_dTX
      "Return specific enthalpy from d, T, and X or Xi"
      extends Modelica.Icons.Function;
      input Density d "Pressure";
      input Temperature T "Specific entropy";
      input MassFraction X[nX] "Mass fractions";
      input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      output SpecificEnthalpy h "specific enthalpy";
    algorithm
        h := specificEnthalpy(
          setState_dTX(
          d,
          T,
          X,
          phase));
    annotation(Documentation(info="<html></html>"));
    end specificEnthalpy_dTX;

      redeclare replaceable function specificEnthalpy_psX
      "Return specific enthalpy from p, s, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction X[nX] "Mass fractions";
          input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEnthalpy h "specific enthalpy";
      algorithm
        h := specificEnthalpy(
          setState_psX(
            p,
            s,
            X,
            phase));
        annotation(Documentation(info="<html></html>"));
      end specificEnthalpy_psX;

      redeclare replaceable function specificEnthalpy_pTX
      "Return specific enthalpy from pressure, temperature and mass fraction"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction X[nX] "Mass fractions";
          input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEnthalpy h "Specific enthalpy at p, T, X";
      algorithm
        h := specificEnthalpy(
          setState_pTX(
            p,
            T,
            X,
            phase));
        annotation(Documentation(info="<html></html>"));
      end specificEnthalpy_pTX;

    replaceable function specificEntropy_phX
      "Return specific entropy from p, h, and X or Xi"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEnthalpy h "Specific enthalpy";
      input MassFraction X[nX] "Mass fractions";
      input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      output SpecificEntropy s "specific enthalpy";
    algorithm
        s := specificEntropy(
          setState_phX(
          p,
          h,
          X,
          phase));
    annotation(Documentation(info="<html></html>"));
    end specificEntropy_phX;

      redeclare replaceable function temperature_phX
      "Return temperature from p, h, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction X[nX] "Mass fractions";
          input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
          output Temperature T "Temperature";
      algorithm
        T := temperature(
          setState_phX(
            p,
            h,
            X,
            phase));
        annotation(Documentation(info="<html></html>"));
      end temperature_phX;

      redeclare replaceable function temperature_psX
      "Return temperature from p, s, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction X[nX] "Mass fractions";
          input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
          output Temperature T "Temperature";
      algorithm
        T := temperature(
          setState_psX(
            p,
            s,
            X,
            phase));
        annotation(Documentation(info="<html></html>"));
      end temperature_psX;

      replaceable function setState_dT
      "Return thermodynamic state from d and T"
        extends Modelica.Icons.Function;
        input Density d "density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "thermodynamic state record";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures.");
        state :=
         setState_dTX(
            d,
            T,
            fill(0, 0),
            phase);
        annotation(Documentation(info="<html></html>"));
      end setState_dT;

      replaceable function setState_ph
      "Return thermodynamic state from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "thermodynamic state record";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures.");
        state :=
          setState_phX(
            p,
            h,
            fill(0, 0),
            phase);
        annotation(Documentation(info="<html></html>"));
      end setState_ph;

      replaceable function setState_ps
      "Return thermodynamic state from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "thermodynamic state record";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures.");
        state :=
          setState_psX(
            p,
            s,
            fill(0, 0),
            phase);
        annotation(Documentation(info="<html></html>"));
      end setState_ps;

      replaceable function setState_pT
      "Return thermodynamic state from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "thermodynamic state record";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures.");
        state :=
          setState_pTX(
            p,
            T,
            fill(0, 0),
            phase);
        annotation(Documentation(info="<html></html>"));
      end setState_pT;

      replaceable function setState_px
      "Return thermodynamic state from pressure and vapour quality"
        input AbsolutePressure p "Pressure";
        input MassFraction x "Vapour quality";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures.");
        state := setState_ph(
            p,
            (1 - x)*bubbleEnthalpy(
            setSat_pX(p,{1})) +
            x*dewEnthalpy(
            setSat_pX(p,{1})),
            2);
        annotation(Documentation(info="<html></html>"));
      end setState_px;

      replaceable function setState_Tx
      "Return thermodynamic state from temperature and vapour quality"
        input Temperature T "Temperature";
        input MassFraction x "Vapour quality";
        output ThermodynamicState state "thermodynamic state record";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures.");
        state := setState_ph(
            saturationPressure_sat(
            setSat_TX(T,{1})),
            (1 - x)*bubbleEnthalpy(
            setSat_TX(T,{1})) +
            x*dewEnthalpy(
            setSat_TX(T,{1})),
            2);
        annotation(Documentation(info="<html></html>"));
      end setState_Tx;

      replaceable function density_ph "Return density from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use density_phX() instead!");
        d := density_phX(p, h, fill(0,0), phase);
        annotation(Documentation(info="<html></html>"));
      end density_ph;

      replaceable function specificEnthalpy_dT
      "Return specific enthalpy from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use specificEnthalpy_dX() instead!");
        h := specificEnthalpy(setState_dTX(
            d,
            T,
            fill(0, 0),
            phase));
        annotation(Documentation(info="<html></html>"));
      end specificEnthalpy_dT;

      replaceable function temperature_ph "Return temperature from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use temperature_phX() instead!");
        T := temperature_phX(p, h, fill(0,0),phase);
        annotation(Documentation(info="<html></html>"));
      end temperature_ph;

      replaceable function pressure_dT "Return pressure from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output AbsolutePressure p "Pressure";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use pressure_dTX() instead!");
        p := pressure(setState_dTX(
            d,
            T,
            fill(0, 0),
            phase));
        annotation(Documentation(info="<html></html>"));
      end pressure_dT;

      replaceable function specificEnthalpy_ps
      "Return specific enthalpy from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use specificEnthalpy_psX() instead!");
        h := specificEnthalpy_psX(p,s,reference_X);
        annotation(Documentation(info="<html></html>"));
      end specificEnthalpy_ps;

      replaceable function temperature_ps "Return temperature from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use temperature_psX() instead!");
        T := temperature_psX(p,s,fill(0,0),phase);
        annotation(Documentation(info="<html></html>"));
      end temperature_ps;

      replaceable function density_ps "Return density from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use density_psX() instead!");
        d := density_psX(p, s, fill(0,0), phase);
        annotation(Documentation(info="<html></html>"));
      end density_ps;

      replaceable function specificEnthalpy_pT
      "Return specific enthalpy from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        assert(nX==1,"This function is not allowed for mixtures. Use specificEnthalpy_pTx() instead!");
        h := specificEnthalpy_pTX(p, T, fill(0,0),phase);
        annotation(Documentation(info="<html></html>"));
      end specificEnthalpy_pT;

      replaceable function density_pT "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(
          setState_pTX(
            p,
            T,
            fill(0, 0),
            phase));
        annotation(Documentation(info="<html></html>"));
      end density_pT;

    redeclare replaceable function density
      input ThermodynamicState state "Thermodynamic state record";
      output Density d;
    algorithm
      d:=state.d;
    end density;

    redeclare replaceable function pressure
      input ThermodynamicState state "Thermodynamic state record";
      output AbsolutePressure p;
    algorithm
      p:=state.p;
    end pressure;

    redeclare replaceable function specificEnthalpy
      input ThermodynamicState state "Thermodynamic state record";
      output SpecificEnthalpy h;
    algorithm
      h:=state.h;
    end specificEnthalpy;

    redeclare replaceable function specificEntropy
      input ThermodynamicState state "Thermodynamic state record";
      output SpecificEntropy s;
    algorithm
      s:=state.s;
    end specificEntropy;

    redeclare replaceable function temperature
      input ThermodynamicState state "Thermodynamic state record";
      output Temperature T;
    algorithm
      T:=state.T;
    end temperature;

      replaceable function vapourQuality "Return vapour quality"
        input ThermodynamicState state "Thermodynamic state record";
        output MassFraction x "Vapour quality";
    protected
        constant SpecificEnthalpy eps = 1e-8;
      algorithm
        x := min(max((specificEnthalpy(state) - bubbleEnthalpy(
          setSat_pX(
          pressure(state), state.X)))/(dewEnthalpy(
          setSat_pX(
          pressure(state), state.X)) - bubbleEnthalpy(
          setSat_pX(
          pressure(state), state.X)) + eps), 0), 1);
        annotation(Documentation(info="<html></html>"));
      end vapourQuality;

      replaceable partial function surfaceTension
      "Return surface tension sigma in the two phase region"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output SurfaceTension sigma
        "Surface tension sigma in the two phase region";
        annotation(Documentation(info="<html></html>"));
      end surfaceTension;

    type DerPressureByDensity = Real (unit="m2/s2");
    type DerDerPressureByDensityByDensity = Real (unit="(m5)/(kg.s2)");
    type DerPressureByTemperature = Real (unit="kg/(K.m.s2)");
    //type DerDensityByTemperature = Real (unit="kg/(m3.K)");
    //type DerDensityByPressure = Real (unit="s2/m2");
    type DerDerPressureByTemperatureByTemperature = Real (unit="kg/(m.s2.K2)");
    type DerDerPressureByTemperatureByDensity = Real (unit="(m2)/(s2.K)");

    type DerEnthalpyByDensity = Real (unit="J.m3/kg");
    //type DerEnthalpyByPressure = Real (unit="J.m.s2/kg");
    type DerEnthalpyByTemperature = Real (unit="J/K");

        annotation(Documentation(info="<html>
  <h1>PartialMixtureTwoPhaseMedium</h1>
  </html>
"));
    end ExternalTwoPhaseMixture;

    partial package FluidPropMedium "FluidProp medium package"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium;
      constant Real h_eps_sat = 1e-6
      "small delta h to ensure computation in the correct phase";
      redeclare replaceable function setBubbleState
      "Set the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase = 0 "phase flag";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
      algorithm
        if (phase == 0) then // phase doesn't matter
          state := setState_ph(sat.psat, sat.hl, phase);
        elseif (phase == 1) then // liquid
          state := setState_ph(sat.psat, sat.hl*(1-h_eps_sat), phase);
        else // two-phase
          state := setState_ph(sat.psat, sat.hl*(1+h_eps_sat), phase);
        end if;
        /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setBubbleState_(sat, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end setBubbleState;

      redeclare replaceable function setDewState
      "Set the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase = 0 "phase flag";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
      algorithm
        if (phase == 0) then // phase doesn't matter
          state := setState_ph(sat.psat, sat.hv, phase);
        elseif (phase == 1) then // vapour
          state := setState_ph(sat.psat, sat.hv*(1+h_eps_sat), phase);
        else // two-phase
          state := setState_ph(sat.psat, sat.hv*(1-h_eps_sat), phase);
        end if;
        /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setDewState_(sat, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include <CoolPropLib.h>", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end setDewState;

      redeclare function bubbleEntropy "Return bubble point specific entropy"
        input SaturationProperties sat "saturation property record";
        output SI.SpecificEntropy sl "boiling curve specific entropy";
      algorithm
        sl := specificEntropy(setBubbleState(sat));
      end bubbleEntropy;

      redeclare function dewEntropy "Return dew point specific entropy"
        input SaturationProperties sat "saturation property record";
        output SI.SpecificEntropy sv "dew curve specific entropy";
      algorithm
        sv := specificEntropy(setDewState(sat));
      end dewEntropy;

      redeclare function surfaceTension
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output SurfaceTension sigma
        "Surface tension sigma in the two phase region";
      algorithm
        assert(false, "FluidProp interface does not provide surface tension");
      end surfaceTension;
    end FluidPropMedium;

    partial package CoolPropMedium
      extends ExternalTwoPhaseMedium(
      mediumName = "CoolPropMedium",
      final libraryName = "CoolProp");
    end CoolPropMedium;
  end Interfaces;


  package Examples "Examples of external medium models"
  extends Modelica.Icons.ExamplesPackage;
    model CompressibleValveSystem "A valve between a source and a sink"
      extends Modelica.Icons.Example;
      Modelica.Fluid.Valves.ValveCompressible valveCompressible(
        m_flow_nominal=0.2,
        redeclare package Medium = WorkingFluid,
        dp_nominal=200000,
        p_nominal=1000000)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Fluid.Sources.FixedBoundary source(
      nPorts=1,
      p=2*system.p_ambient,
      T=system.T_ambient,
      redeclare package Medium = WorkingFluid)
        annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
    replaceable package WorkingFluid = CoolProp2Modelica.Media.R601_CP
      constrainedby Modelica.Media.Interfaces.PartialMedium
                                                  annotation (choicesAllMatching=true);
    replaceable package HeatingFluid =
          Modelica.Media.Incompressible.Examples.Essotherm650 constrainedby
      Modelica.Media.Interfaces.PartialMedium     annotation (choicesAllMatching=true);
    replaceable package CoolingFluid = CoolProp2Modelica.Media.R718_CP
      constrainedby Modelica.Media.Interfaces.PartialMedium
                                                  annotation (choicesAllMatching=true);
      Modelica.Fluid.Sources.FixedBoundary sink(
      nPorts=1,
      redeclare package Medium = WorkingFluid,
      p=system.p_ambient,
      T=system.T_ambient)
        annotation (Placement(transformation(extent={{80,-10},{60,10}})));
      Modelica.Blocks.Sources.Sine sine(freqHz=2, offset=1)
        annotation (Placement(transformation(extent={{-30,60},{-10,80}})));
      inner Modelica.Fluid.System system
        annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    equation
    connect(source.ports[1], valveCompressible.port_a)     annotation (Line(
          points={{-60,6.66134e-16},{-36,6.66134e-16},{-36,6.10623e-16},{-10,
            6.10623e-16}},
          color={0,127,255},
          smooth=Smooth.None));
    connect(valveCompressible.port_b, sink.ports[1])        annotation (Line(
          points={{10,6.10623e-16},{36,6.10623e-16},{36,6.66134e-16},{60,
            6.66134e-16}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(sine.y, valveCompressible.opening) annotation (Line(
          points={{-9,70},{6.66134e-16,70},{6.66134e-16,8}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (Diagram(graphics), Documentation(info="<html>
<p><h4><font color=\"#008000\">Compressible Valve System</font></h4></p>
<p>This file illustrates how CoolProp2Modelica can be used with standard components from the Modelica.Fluid library. You can redeclare the WorkingFluid package with any other fluid that matches the PartialMedium interface. Changes will automatically propagate to all components.</p>
</html>"));
    end CompressibleValveSystem;
  end Examples;


  package Media "Medium packages compatible with Modelica.Media"
  extends Modelica.Icons.VariantsPackage;
    package TestMedium "Simple water medium model for debugging and testing"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="TestMedium",
      libraryName="TestMedium",
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.pT);
    end TestMedium;

    package R245fa_CP "R245fa properties from CoolProp"
    //  extends Modelica.Media.Water.StandardWater;
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="R245fa",
      libraryName="CoolProp",
      substanceNames={"R245fa"});
    end R245fa_CP;

    package R290_CP "R290, computation of Propane Properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="TestMedium",
      libraryName="CoolProp",
      substanceName="propane",
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R290_CP;

    package R600_CP "R600, n-Butane properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="n-Butane",
      libraryName="CoolProp",
      substanceNames={"n-Butane"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R600_CP;

    package R600a_CP "R600a, Isobutane properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="Isobutane",
      libraryName="CoolProp",
      substanceNames={"IsoButane"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R600a_CP;

    package R601_CP "R601, n-Pentane properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="n-Pentane",
      libraryName="CoolProp",
      substanceNames={"n-Pentane"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R601_CP;

    package R601a_CP "R601a, Isopentane properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="Isopentane",
      libraryName="CoolProp",
      substanceNames={"IsoPentane"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R601a_CP;

    package R718_CP "R718, water IAPWS 95 properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="Water",
      libraryName="CoolProp",
      substanceNames={"Water"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R718_CP;

    package SES36_CP "Solkatherm properties using CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="SES36",
      libraryName="CoolProp",
      substanceName="SES36",
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end SES36_CP;

    package R290_FPST
    "Propane properties using the StanMix library of FluidProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="TestMedium",
      libraryName="FluidProp.StanMix",
      substanceName="propane",
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R290_FPST;

    package R290_FPRP
    "Propane properties using Refprop through FluidProp (requires the full version of FluidProp)"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="TestMedium",
      libraryName="FluidProp.RefProp",
      substanceName="propane",
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R290_FPRP;

    package R744_FPRP "Carbon-Dioxide from Refprop via FluidProp"
      extends CoolProp2Modelica.Interfaces.FluidPropMedium(
      mediumName="Carbon Dioxide",
      libraryName="FluidProp.RefProp",
      substanceNames={"CO2"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end R744_FPRP;

    package WaterTPSI_FP
      extends CoolProp2Modelica.Interfaces.FluidPropMedium(
      mediumName="Water",
      libraryName="FluidProp.TPSI",
      substanceNames={"H2O"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end WaterTPSI_FP;

    package WaterIF95_FP
      extends CoolProp2Modelica.Interfaces.FluidPropMedium(
      mediumName="Water",
      libraryName="FluidProp.RefProp",
      substanceNames={"H2O"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end WaterIF95_FP;

    package R410a_debug "R410a, properties using Refprop via CoolProp"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="R410a",
      libraryName="CoolProp",
      substanceNames={"REFPROP-MIX:R32[0.697615]&R125[0.302385]|debug=1"},
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end R410a_debug;

    package Solkatherm_debug
    "Solkatherm properties using CoolProp with debug option"
      extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
      mediumName="SES36",
      libraryName="CoolProp",
      substanceName="SES36|calc_transport=0|debug=1|enable_TTSE=1",
      ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
        // calc_transport: set to one to compute the transport properties (not yet implemented in January 2013)
        // debug: integer from 0 to 10. 0 corresponds to no debut, while 10 is maximum debug.
        //        This outputs the different calls received by CoolProp in the console window
        // enable_TTSE: set to 1 to enable interpolated properties as a function of p-h
        //              Involves about 2 more seconds at initialization but integration is about 40 times faster
      annotation ();
    end Solkatherm_debug;
  end Media;


  package Test "Test models"
    package TestMedium "Test cases for TestMedium"
      model TestStatesSat
      "Test case using TestMedium, with baseProperties and state + sat records without explicit uniqueID handling"
        replaceable package Medium = Media.TestMedium;
        Medium.BaseProperties baseProperties1;
        Medium.BaseProperties baseProperties2;
        Medium.ThermodynamicState state1;
        Medium.ThermodynamicState state2;
        Medium.SaturationProperties sat1;
        Medium.SaturationProperties sat2;
        Medium.Temperature Ts;
        Medium.AbsolutePressure ps;
        GenericModels.CompleteThermodynamicState completeState1(
          redeclare package Medium = Medium, state = state1);
        GenericModels.CompleteThermodynamicState completeState2(
          redeclare package Medium = Medium, state = state2);
        GenericModels.CompleteSaturationProperties completeSat1(
          redeclare package Medium = Medium, sat = sat1);
        GenericModels.CompleteSaturationProperties completeSat2(
          redeclare package Medium = Medium, sat = sat2);
        GenericModels.CompleteBubbleDewStates completeBubbleDewStates1(
          redeclare package Medium = Medium, sat = sat1);
        GenericModels.CompleteBubbleDewStates completeBubbleDewStates2(
          redeclare package Medium = Medium, sat = sat1);
      equation
        baseProperties1.p = 1e5+1e5*time;
        baseProperties1.h = 1e5;
        baseProperties2.p = 1e5;
        baseProperties2.h = 1e5 + 2e5*time;
        state1 = Medium.setState_ph(1e5 + 1e5*time, 1e5);
        state2 = Medium.setState_pT(1e5, 300+ 50*time);
        sat1 = Medium.setSat_p(1e5 + 1e5*time);
        sat2 = Medium.setSat_T(300 + 50 * time);
        Ts = Medium.saturationTemperature(1e5+1e5*time);
        ps = Medium.saturationPressure(300 + 50*time);
      end TestStatesSat;

      model TestBasePropertiesExplicit
      "Test case using TestMedium and BaseProperties with explicit equations"
        replaceable package Medium = Media.TestMedium;
        CoolProp2Modelica.Test.TestMedium.GenericModels.CompleteBaseProperties medium1(
          redeclare package Medium = Medium)
        "Constant pressure, varying enthalpy";
        CoolProp2Modelica.Test.TestMedium.GenericModels.CompleteBaseProperties medium2(
          redeclare package Medium = Medium)
        "Varying pressure, constant enthalpy";
      equation
        medium1.baseProperties.p = 1e5+1e5*time;
        medium1.baseProperties.h = 1e5;
        medium2.baseProperties.p = 1e5;
        medium2.baseProperties.h = 1e5 + 2e5*time;
      end TestBasePropertiesExplicit;

      model TestBasePropertiesImplicit
      "Test case using TestMedium and BaseProperties with implicit equations"
        replaceable package Medium = Media.TestMedium;
        CoolProp2Modelica.Test.TestMedium.GenericModels.CompleteBaseProperties medium1(
          redeclare package Medium = Medium, baseProperties(h(start=1e5)))
        "Constant pressure, varying enthalpy";
        CoolProp2Modelica.Test.TestMedium.GenericModels.CompleteBaseProperties medium2(
          redeclare package Medium = Medium, baseProperties(h(start=1e5)))
        "Varying pressure, constant enthalpy";
      equation
        medium1.baseProperties.p = 1e5*time;
        medium1.baseProperties.T = 300 + 25*time;
        medium2.baseProperties.p = 1e5+1e5*time;
        medium2.baseProperties.T = 300;
      end TestBasePropertiesImplicit;

      model TestBasePropertiesDynamic
      "Test case using TestMedium and dynamic equations"
        replaceable package Medium = Media.TestMedium;
        parameter SI.Volume V = 1 "Storage Volume";
        parameter Real p_atm = 101325 "Atmospheric pressure";
        parameter SI.Temperature Tstart = 300;
        parameter Real Kv0 = 1.00801e-2 "Valve flow coefficient";
        Medium.BaseProperties medium(preferredMediumStates = true);
        SI.Mass M;
        SI.Energy U;
        SI.MassFlowRate win(start = 100);
        SI.MassFlowRate wout;
        SI.SpecificEnthalpy hin;
        SI.SpecificEnthalpy hout;
        SI.Power Q;
        Real Kv;
      equation
        // Mass & energy balance equation
        M = medium.d*V;
        U = medium.u*M;
        der(M) = win - wout;
        der(U) = win*hin - wout*hout + Q;
        // Inlet pump equations
        medium.p - p_atm = 2e5 - (1e5/100^2)*win^2;
        hin = 1e5;
        // Outlet valve equation
        wout = Kv * sqrt(medium.d*(medium.p - p_atm));
        hout = medium.h;
        // Input variables
        Kv = if time<50 then Kv0 else Kv0*1.1;
        Q = if time < 1 then 0 else 1e7;
      initial equation
        // Initial conditions
        // Fixed initial states
        // medium.p = 2e5;
        // medium.h = 1e5;
        // Steady state equations
        der(medium.p) = 0;
        der(medium.h) = 0;
        annotation (experiment(StopTime=80, Tolerance=1e-007),experimentSetupOutput(
            equdistant=false));
      end TestBasePropertiesDynamic;

      package GenericModels
      "Contains generic models to use for thorough medium model testing"
        model CompleteFluidConstants
        "Compute all available medium fluid constants"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // Fluid constants
          Medium.Temperature Tc = Medium.fluidConstants[1].criticalTemperature;
          Medium.AbsolutePressure pc = Medium.fluidConstants[1].criticalPressure;
          Medium.MolarVolume vc = Medium.fluidConstants[1].criticalMolarVolume;
          Medium.MolarMass MM = Medium.fluidConstants[1].molarMass;
        end CompleteFluidConstants;

        model CompleteThermodynamicState
        "Compute all available two-phase medium properties from a ThermodynamicState model"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // ThermodynamicState record
          input Medium.ThermodynamicState state;
          // Medium properties
          Medium.AbsolutePressure p =                Medium.pressure(state);
          Medium.SpecificEnthalpy h =                Medium.specificEnthalpy(state);
          Medium.Temperature T =                     Medium.temperature(state);
          Medium.Density d =                         Medium.density(state);
          Medium.SpecificEntropy s =                 Medium.specificEntropy(state);
          Medium.SpecificHeatCapacity cp =           Medium.specificHeatCapacityCp(state);
          Medium.SpecificHeatCapacity cv =           Medium.specificHeatCapacityCv(state);
          Medium.IsobaricExpansionCoefficient beta = Medium.isobaricExpansionCoefficient(state);
          SI.IsothermalCompressibility kappa =       Medium.isothermalCompressibility(state);
          Medium.DerDensityByPressure d_d_dp_h =     Medium.density_derp_h(state);
          Medium.DerDensityByEnthalpy d_d_dh_p =     Medium.density_derh_p(state);
          Medium.MolarMass MM =                      Medium.molarMass(state);
        end CompleteThermodynamicState;

        model CompleteSaturationProperties
        "Compute all available saturation properties from a SaturationProperties record"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // SaturationProperties record
          input Medium.SaturationProperties sat;
          // Saturation properties
          Medium.Temperature Ts =      Medium.saturationTemperature_sat(sat);
          Medium.Density dl =          Medium.bubbleDensity(sat);
          Medium.Density dv =          Medium.dewDensity(sat);
          Medium.SpecificEnthalpy hl = Medium.bubbleEnthalpy(sat);
          Medium.SpecificEnthalpy hv = Medium.dewEnthalpy(sat);
          Real d_Ts_dp =               Medium.saturationTemperature_derp_sat(sat);
          Real d_dl_dp =               Medium.dBubbleDensity_dPressure(sat);
          Real d_dv_dp =               Medium.dDewDensity_dPressure(sat);
          Real d_hl_dp =               Medium.dBubbleEnthalpy_dPressure(sat);
          Real d_hv_dp =               Medium.dDewEnthalpy_dPressure(sat);
        end CompleteSaturationProperties;

        model CompleteBubbleDewStates
        "Compute all available properties for dewpoint and bubble point states corresponding to a sat record"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // SaturationProperties record
          input Medium.SaturationProperties sat;
          CompleteThermodynamicState dewStateOnePhase(
            state=Medium.setDewState(sat, 1), redeclare package Medium = Medium);
          CompleteThermodynamicState dewStateTwoPhase(
            state=Medium.setDewState(sat, 2), redeclare package Medium = Medium);
          CompleteThermodynamicState bubbleStateOnePhase(
            state=Medium.setBubbleState(sat, 1), redeclare package Medium = Medium);
          CompleteThermodynamicState bubbleStateTwoPhase(
            state=Medium.setBubbleState(sat, 2), redeclare package Medium = Medium);
        end CompleteBubbleDewStates;

        model CompleteBaseProperties
        "Compute all available two-phase medium properties from a BaseProperties model"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // BaseProperties object
          Medium.BaseProperties baseProperties;
          // All the complete properties
          CompleteThermodynamicState completeState(
            redeclare package Medium = Medium, state=baseProperties.state);
          CompleteSaturationProperties completeSat(
            redeclare package Medium = Medium, sat=baseProperties.sat);
          CompleteFluidConstants completeConstants(
            redeclare package Medium = Medium);
          CompleteBubbleDewStates completeBubbleDewStates(
              redeclare package Medium = Medium, sat=baseProperties.sat);
        end CompleteBaseProperties;
      end GenericModels;
    end TestMedium;

    package FluidProp "Test cases for FluidPropMedium"
      partial package GenericModels "Generic models for FluidProp media tests"
        model CompleteFluidConstants
        "Compute all available medium fluid constants"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // Fluid constants
          Medium.Temperature Tc = Medium.fluidConstants[1].criticalTemperature;
          Medium.AbsolutePressure pc = Medium.fluidConstants[1].criticalPressure;
          Medium.MolarVolume vc = Medium.fluidConstants[1].criticalMolarVolume;
          Medium.MolarMass MM = Medium.fluidConstants[1].molarMass;
        end CompleteFluidConstants;

        model CompleteThermodynamicState
        "Compute all available two-phase medium properties from a ThermodynamicState model"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // ThermodynamicState record
          input Medium.ThermodynamicState state;
          // Medium properties
          Medium.AbsolutePressure p =                Medium.pressure(state);
          Medium.SpecificEnthalpy h =                Medium.specificEnthalpy(state);
          Medium.Temperature T =                     Medium.temperature(state);
          Medium.Density d =                         Medium.density(state);
          Medium.SpecificEntropy s =                 Medium.specificEntropy(state);
          Medium.SpecificHeatCapacity cp =           Medium.specificHeatCapacityCp(state);
          Medium.SpecificHeatCapacity cv =           Medium.specificHeatCapacityCv(state);
        // Not yet implemented in FluidProp
          Medium.IsobaricExpansionCoefficient beta = Medium.isobaricExpansionCoefficient(state);
          SI.IsothermalCompressibility kappa =       Medium.isothermalCompressibility(state);
          Medium.DerDensityByPressure d_d_dp_h =     Medium.density_derp_h(state);
          Medium.DerDensityByEnthalpy d_d_dh_p =     Medium.density_derh_p(state);
          Medium.MolarMass MM =                      Medium.molarMass(state);
        end CompleteThermodynamicState;

        model CompleteSaturationProperties
        "Compute all available saturation properties from a SaturationProperties record"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // SaturationProperties record
          input Medium.SaturationProperties sat;
          // Saturation properties
          Medium.Temperature Ts =      Medium.saturationTemperature_sat(sat);
          Medium.Density dl =          Medium.bubbleDensity(sat);
          Medium.Density dv =          Medium.dewDensity(sat);
          Medium.SpecificEnthalpy hl = Medium.bubbleEnthalpy(sat);
          Medium.SpecificEnthalpy hv = Medium.dewEnthalpy(sat);
          Real d_Ts_dp =               Medium.saturationTemperature_derp_sat(sat);
          Real d_dl_dp =               Medium.dBubbleDensity_dPressure(sat);
          Real d_dv_dp =               Medium.dDewDensity_dPressure(sat);
          Real d_hl_dp =               Medium.dBubbleEnthalpy_dPressure(sat);
          Real d_hv_dp =               Medium.dDewEnthalpy_dPressure(sat);
        end CompleteSaturationProperties;

        model CompleteBubbleDewStates
        "Compute all available properties for dewpoint and bubble point states corresponding to a sat record"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // SaturationProperties record
          input Medium.SaturationProperties sat;
          CompleteThermodynamicState dewStateOnePhase(state=
                Medium.setDewState(sat, 1), redeclare package Medium = Medium);
          CompleteThermodynamicState dewStateTwoPhase(state=
                Medium.setDewState(sat, 2), redeclare package Medium = Medium);
          CompleteThermodynamicState bubbleStateOnePhase(state=
                Medium.setBubbleState(sat, 1), redeclare package Medium = Medium);
          CompleteThermodynamicState bubbleStateTwoPhase(state=
                Medium.setBubbleState(sat, 2), redeclare package Medium = Medium);
        end CompleteBubbleDewStates;

        model CompleteBaseProperties
        "Compute all available two-phase medium properties from a BaseProperties model"
          replaceable package Medium =
                Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          // BaseProperties object
          Medium.BaseProperties baseProperties;
          // All the complete properties
          CompleteThermodynamicState completeState(
            redeclare package Medium = Medium, state=baseProperties.state);
          CompleteSaturationProperties completeSat(
            redeclare package Medium = Medium, sat=baseProperties.sat);
          CompleteFluidConstants completeConstants(
            redeclare package Medium = Medium);
          CompleteBubbleDewStates completeBubbleDewStates(
            redeclare package Medium = Medium, sat=baseProperties.sat);
        end CompleteBaseProperties;

        partial model TestStatesSat
        "Test case with baseProperties and state + sat records"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          Medium.AbsolutePressure p1;
          Medium.SpecificEnthalpy h1;
          Medium.AbsolutePressure p2;
          Medium.SpecificEnthalpy h2;
          Medium.Temperature T2;
          Medium.BaseProperties baseProperties1;
          Medium.BaseProperties baseProperties2;
          Medium.ThermodynamicState state1;
          Medium.ThermodynamicState state2;
          Medium.SaturationProperties sat1;
          Medium.SaturationProperties sat2;
          Medium.Temperature Ts;
          Medium.AbsolutePressure ps;
          CompleteThermodynamicState
            completeState1(redeclare package Medium = Medium, state=state1);
          CompleteThermodynamicState
            completeState2(redeclare package Medium = Medium, state=state2);
          CompleteSaturationProperties
            completeSat1(redeclare package Medium = Medium, sat=sat1);
          CompleteSaturationProperties
            completeSat2(redeclare package Medium = Medium, sat=sat2);
          CompleteBubbleDewStates
            completeBubbleDewStates1(redeclare package Medium = Medium, sat=sat1);
          CompleteBubbleDewStates
            completeBubbleDewStates2(redeclare package Medium = Medium, sat=sat1);
        equation
          baseProperties1.p = p1;
          baseProperties1.h = h1;
          baseProperties2.p = p2;
          baseProperties2.h = h2;
          state1 = Medium.setState_ph(p1, h1);
          state2 = Medium.setState_pT(p2, T2);
          sat1 = Medium.setSat_p(p1);
          sat2 = Medium.setSat_T(T2);
          Ts = Medium.saturationTemperature(p1);
          ps = Medium.saturationPressure(T2);
        end TestStatesSat;

        partial model TestBasePropertiesExplicit
        "Test case using BaseProperties and explicit equations"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          CompleteBaseProperties medium1(redeclare package Medium = Medium)
          "Constant pressure, varying enthalpy";
          CompleteBaseProperties medium2(redeclare package Medium = Medium)
          "Varying pressure, constant enthalpy";
          Medium.AbsolutePressure p1;
          Medium.AbsolutePressure p2;
          Medium.SpecificEnthalpy h1;
          Medium.SpecificEnthalpy h2;
        equation
          medium1.baseProperties.p = p1;
          medium1.baseProperties.h = h1;
          medium2.baseProperties.p = p2;
          medium2.baseProperties.h = h2;
        end TestBasePropertiesExplicit;

        partial model TestBasePropertiesImplicit
        "Test case using BaseProperties and implicit equations"
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          parameter Medium.SpecificEnthalpy hstart
          "Start value for specific enthalpy";
          CompleteBaseProperties medium1(redeclare package Medium = Medium,
                                         baseProperties(h(start = hstart)))
          "Constant pressure, varying enthalpy";
          CompleteBaseProperties medium2(redeclare package Medium = Medium,
                                         baseProperties(h(start = hstart)))
          "Varying pressure, constant enthalpy";
          Medium.AbsolutePressure p1;
          Medium.AbsolutePressure p2;
          Medium.Temperature T1;
          Medium.Temperature T2;
        equation
          medium1.baseProperties.p = p1;
          medium1.baseProperties.T = T1;
          medium2.baseProperties.p = p2;
          medium2.baseProperties.T = T2;
        end TestBasePropertiesImplicit;

      partial model TestBasePropertiesDynamic
        "Test case using BaseProperties and dynamic equations"
        replaceable package Medium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
        parameter SI.Volume V = 1 "Storage Volume";
        parameter Real p_atm = 101325 "Atmospheric pressure";
        parameter SI.Temperature Tstart = 300;
        parameter SI.SpecificEnthalpy hstart;
        parameter Real Kv0 "Valve flow coefficient";
        Medium.BaseProperties medium(preferredMediumStates = true,
                                     h(start=1e5));
        SI.Mass M;
        SI.Energy U;
        SI.MassFlowRate win(start = 100);
        SI.MassFlowRate wout;
        SI.SpecificEnthalpy hin;
        SI.SpecificEnthalpy hout;
        SI.Power Q;
        Real Kv;
      equation
        // Mass & energy balance equation
        M = medium.d*V;
        U = medium.u*M;
        der(M) = win - wout;
        der(U) = win*hin - wout*hout + Q;
        // Outlet valve equation
        wout = Kv * sqrt(medium.d*(medium.p - p_atm));
        hout = medium.h;
      initial equation
        // Steady state equations
        der(medium.p) = 0;
        der(medium.h) = 0;
        annotation (experiment(StopTime=80, Tolerance=1e-007),experimentSetupOutput(
              equdistant=false));
      end TestBasePropertiesDynamic;

        partial model CompareModelicaFluidProp
        "Comparison between Modelica and FluidProp models"
          replaceable package ModelicaMedium =
              Modelica.Media.Interfaces.PartialTwoPhaseMedium;
          replaceable package FluidPropMedium =
              CoolProp2Modelica.Interfaces.FluidPropMedium;
          CompleteBaseProperties modelicaMedium(
            redeclare package Medium = ModelicaMedium) "Modelica medium model";
          CompleteBaseProperties fluidPropMedium(
            redeclare package Medium = FluidPropMedium)
          "FluidProp medium model";
          parameter Modelica.SIunits.Pressure pmin;
          parameter Modelica.SIunits.Pressure pmax;
          parameter Modelica.SIunits.SpecificEnthalpy hmin;
          parameter Modelica.SIunits.SpecificEnthalpy hmax;
        equation
          modelicaMedium.baseProperties.p = pmin + (pmax-pmin)*time;
          modelicaMedium.baseProperties.h = hmin + (hmax-hmin)*time;
          fluidPropMedium.baseProperties.p = pmin + (pmax-pmin)*time;
          fluidPropMedium.baseProperties.h = hmin + (hmax-hmin)*time;
        end CompareModelicaFluidProp;
      end GenericModels;

      package IF95 "Test suite for the FluidProp-Refprop IF95 medium model"
        model TestStatesSat
        "Test case with baseProperties and state + sat records"
          extends GenericModels.TestStatesSat(
            redeclare package Medium = CoolProp2Modelica.Media.WaterIF95_FP);
        equation
          p1 = 1e5 + time*1e5;
          h1 = 1e5;
          p2 = 1e5;
          h2 = 1e5 + time*2e5;
          T2 = 300 + 50*time;
        end TestStatesSat;

        model TestBasePropertiesExplicit
        "Test case using BaseProperties and explicit equations"
          extends GenericModels.TestBasePropertiesExplicit(
            redeclare package Medium = CoolProp2Modelica.Media.WaterIF95_FP);
        equation
          p1 = 1e5+1e5*time;
          h1 = 1e5;
          p2 = 1e5;
          h2 = 1e5 + 2e5*time;
        end TestBasePropertiesExplicit;

        model TestBasePropertiesImplicit
        "Test case using BaseProperties and implicit equations"
          extends GenericModels.TestBasePropertiesImplicit(
            redeclare package Medium = CoolProp2Modelica.Media.WaterIF95_FP,
            hstart = 1e5);
        equation
          p1 = 1e5+1e5*time;
          T1 = 300 + 25*time;
          p2 = 1e5+1e5*time;
          T2 = 300;
        end TestBasePropertiesImplicit;

      model TestBasePropertiesDynamic
        "Test case using BaseProperties and dynamic equations"
        extends GenericModels.TestBasePropertiesDynamic(
          redeclare package Medium = CoolProp2Modelica.Media.WaterIF95_FP,
          Tstart = 300,
          Kv0 = 1.00801e-2);
      equation
        // Inlet pump equations
        medium.p - p_atm = 2e5 - (1e5/100^2)*win^2;
        hin = 1e5;
        // Input variables
        Kv = if time<50 then Kv0 else Kv0*1.1;
        Q = if time < 1 then 0 else 1e7;
        annotation (experiment(StopTime=80, Tolerance=1e-007),experimentSetupOutput(
              equdistant=false));
      end TestBasePropertiesDynamic;

        model TestBasePropertiesExplicit_ModelicaIF97
        "Test case using FluidProp IF95 and explicit equations"
          extends TestBasePropertiesExplicit(
            redeclare package Medium = Modelica.Media.Water.StandardWater);
        end TestBasePropertiesExplicit_ModelicaIF97;

        model CompareModelicaFluidProp_liquid
        "Comparison between Modelica IF97 and FluidProp IF95 models - liquid"
          extends GenericModels.CompareModelicaFluidProp(
            redeclare package ModelicaMedium =
                Modelica.Media.Water.StandardWater,
            redeclare package FluidPropMedium =
                CoolProp2Modelica.Media.WaterIF95_FP,
            pmin = 1e5,
            pmax = 1e5,
            hmin = 1e5,
            hmax = 4e5);
        end CompareModelicaFluidProp_liquid;

        model CompareModelicaFluidProp_twophase
        "Comparison between Modelica IF97 and FluidProp IF95 models - liquid"
          extends GenericModels.CompareModelicaFluidProp(
            redeclare package ModelicaMedium =
                Modelica.Media.Water.StandardWater,
            redeclare package FluidPropMedium =
                CoolProp2Modelica.Media.WaterIF95_FP,
            pmin = 60e5,
            pmax = 60e5,
            hmin = 1000e3,
            hmax = 2000e3);
        end CompareModelicaFluidProp_twophase;

        model CompareModelicaFluidProp_vapour
        "Comparison between Modelica IF97 and FluidProp IF95 models - liquid"
          extends GenericModels.CompareModelicaFluidProp(
            redeclare package ModelicaMedium =
                Modelica.Media.Water.StandardWater,
            redeclare package FluidPropMedium =
                CoolProp2Modelica.Media.WaterIF95_FP,
            pmin = 60e5,
            pmax = 60e5,
            hmin = 2800e3,
            hmax = 3200e3);
        end CompareModelicaFluidProp_vapour;
      end IF95;
    end FluidProp;

    package WrongMedium "Test cases with wrong medium models"
    model TestWrongMedium
      "Test the error reporting messages for unsupported external media"
      package Medium = CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium;
      Medium.BaseProperties medium;
    equation
      medium.p = 1e5;
      medium.h = 1e5;
    end TestWrongMedium;
    end WrongMedium;

    package benchmark
    "Comparison of computational speed for different libraries"
      package fluids
        package propane_CP
        "R290, computation of Propane Properties using CoolProp"
          extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
          mediumName="TestMedium",
          libraryName="CoolProp",
          substanceName="propane",
          ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
          annotation ();
        end propane_CP;

        package propane_FPST
        "Propane properties using the StanMix library of FluidProp"
          extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
          mediumName="TestMedium",
          libraryName="FluidProp.StanMix",
          substanceName="propane",
          ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
          annotation ();
        end propane_FPST;

        package propane_FPRP
        "Propane properties using Refprop through FluidProp (requires the full version of FluidProp)"
          extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
          mediumName="TestMedium",
          libraryName="FluidProp.RefProp",
          substanceName="propane",
          ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
          annotation ();
        end propane_FPRP;

        package propane_CP_TTSE
        "R290, computation of Propane Properties using CoolProp"
          extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
          mediumName="TestMedium",
          libraryName="CoolProp",
          substanceName="propane|enable_TTSE=1",
          ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
          annotation ();
        end propane_CP_TTSE;

        package propane_CP_transport
        "R290, computation of Propane Properties using CoolProp"
          extends CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium(
          mediumName="TestMedium",
          libraryName="CoolProp",
          substanceName="propane|calc_transport=0",
          ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
          annotation ();
        end propane_CP_transport;
      end fluids;

      package test
        model propane_CP_baseproperties
          replaceable package wf = benchmark.fluids.propane_CP constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model";
          wf.BaseProperties fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -1E5 + time*1E5;
          fluid.p = p;
          fluid.h = h;
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_CP_baseproperties;

        model propane_CP_ThermodynamicState
        "ThermodynamicState is faster since it does not call the saturation properties"
          replaceable package wf = benchmark.fluids.propane_CP constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model";
          wf.ThermodynamicState fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -1E5 + time*1E5;
          fluid = wf.setState_ph(p=p,h=h);
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_CP_ThermodynamicState;

        model propane_FPSM_baseproperties
          replaceable package wf =  benchmark.fluids.propane_FPST
            constrainedby Modelica.Media.Interfaces.PartialMedium
          "Medium model";
          wf.BaseProperties fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -7E5 + time*1E5;
          fluid.p = p;
          fluid.h = h;
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_FPSM_baseproperties;

        model propane_FPSM_ThermodynamicState
        "ThermodynamicState is faster since it does not call the saturation properties"
          replaceable package wf = benchmark.fluids.propane_FPST constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model";
          wf.ThermodynamicState fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -7E5 + time*1E5;
          fluid = wf.setState_ph(p=p,h=h);
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_FPSM_ThermodynamicState;

        model propane_FPRP_ThermodynamicState
        "ThermodynamicState is faster since it does not call the saturation properties"
          replaceable package wf = benchmark.fluids.propane_FPRP constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model";
          wf.ThermodynamicState fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = 1E5 + time*1E5;
          fluid = wf.setState_ph(p=p,h=h);
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_FPRP_ThermodynamicState;

        model propane_FPRP_baseproperties
          replaceable package wf =  benchmark.fluids.propane_FPRP
            constrainedby Modelica.Media.Interfaces.PartialMedium
          "Medium model";
          wf.BaseProperties fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -1E5 + time*1E5;
          fluid.p = p;
          fluid.h = h;
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_FPRP_baseproperties;

        model propane_CP_TTSE
        "ThermodynamicState is faster since it does not call the saturation properties"
          replaceable package wf = benchmark.fluids.propane_CP_TTSE
            constrainedby Modelica.Media.Interfaces.PartialMedium
          "Medium model";
          wf.ThermodynamicState fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -1E5 + time*1E5;
          fluid = wf.setState_ph(p=p,h=h);
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_CP_TTSE;

        model propane_CP_transport
        "ThermodynamicState is faster since it does not call the saturation properties"
          replaceable package wf = benchmark.fluids.propane_CP_transport
            constrainedby Modelica.Media.Interfaces.PartialMedium
          "Medium model";
          wf.ThermodynamicState fluid "Properties of the two-phase fluid";
          Modelica.SIunits.SpecificEnthalpy h;
          Modelica.SIunits.Pressure p;
        equation
          p = 1E5;
          h = -1E5 + time*1E5;
          fluid = wf.setState_ph(p=p,h=h);
          annotation (experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"), __Dymola_experimentSetupOutput);
        end propane_CP_transport;

        model propane_TILMedia
          TILMedia.Refrigerant refrigerant(
            refrigerantName="Refprop.PROPANE",
            inputChoice=TILMedia.Internals.InputChoicesRefrigerant.ph,
            computeTransportProperties=false,
            interpolateTransportProperties=false,
            computeSurfaceTension=false)
            annotation (extent=[-48,24; -28,44]);
        equation
          refrigerant.p = 1e5;
          refrigerant.h = -1E5 + time*1E5;
            annotation (                                                   Diagram,
            experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"),
            __Dymola_experimentSetupOutput);
        end propane_TILMedia;

        model propane_TILMedia_transport
          TILMedia.Refrigerant refrigerant(
            refrigerantName="Refprop.PROPANE",
            inputChoice=TILMedia.Internals.InputChoicesRefrigerant.ph,
            computeTransportProperties=true,
            interpolateTransportProperties=false,
            computeSurfaceTension=false)
            annotation (extent=[-48,24; -28,44]);
        equation
          refrigerant.p = 1e5;
          refrigerant.h = -1E5 + time*1E5;
            annotation (                                                   Diagram,
            experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"),
            __Dymola_experimentSetupOutput);
        end propane_TILMedia_transport;

        model propane_TILMedia_interptransport
          TILMedia.Refrigerant refrigerant(
            refrigerantName="Refprop.PROPANE",
            inputChoice=TILMedia.Internals.InputChoicesRefrigerant.ph,
            computeTransportProperties=false,
            interpolateTransportProperties=true,
            computeSurfaceTension=false)
            annotation (extent=[-48,24; -28,44]);
        equation
          refrigerant.p = 1e5;
          refrigerant.h = -1E5 + time*1E5;
            annotation (                                                   Diagram,
            experiment(
              StopTime=10,
              __Dymola_NumberOfIntervals=20000,
              __Dymola_Algorithm="Euler"),
            __Dymola_experimentSetupOutput);
        end propane_TILMedia_interptransport;

        model propane_TILMedia_surfacetension
          TILMedia.Refrigerant refrigerant(
            refrigerantName="Refprop.PROPANE",
            inputChoice=TILMedia.Internals.InputChoicesRefrigerant.ph,
            computeTransportProperties=false,
            interpolateTransportProperties=false,
            computeSurfaceTension=true)
            annotation (extent=[-48,24; -28,44]);
        equation
          refrigerant.p = 1e5;
          refrigerant.h = -1E5 + time*1E5;
            annotation (                                                   Diagram,
            experiment(
              StopTime=1000,
              NumberOfIntervals=2000,
              Algorithm="Euler"),
            __Dymola_experimentSetupOutput);
        end propane_TILMedia_surfacetension;
      end test;
      annotation ();
    end benchmark;

    model test_propane_coolprop
      replaceable package wf = CoolProp2Modelica.Media.R290_CP constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 0 + time*1E6;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
    end test_propane_coolprop;

    model test_propane_fluidprop
      replaceable package wf = CoolProp2Modelica.Media.R290_FPST constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = -7e5 + time*1E6;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
    end test_propane_fluidprop;

    model test_propane_refprop
      replaceable package wf = CoolProp2Modelica.Media.R290_FPRP constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 0 + time*1E6;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
    end test_propane_refprop;

    model test_solkatherm
      replaceable package wf = CoolProp2Modelica.Media.SES36_CP
                                                             constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 0 + time*1E5;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
      annotation (experiment(StopTime=10, Algorithm="Dassl"),
          __Dymola_experimentSetupOutput);
    end test_solkatherm;

    model test_solkatherm_debug "Please note that the debug information appears in the DOS window only. Run 'dymosim >> log.txt' from the cmd to log the 
  debug information in a file"
      replaceable package wf = CoolProp2Modelica.Media.Solkatherm_debug
      constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 0 + time*3e5;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
      annotation (experiment(Algorithm="Dassl"), __Dymola_experimentSetupOutput);
    end test_solkatherm_debug;

    model test_pentane_coolprop
      replaceable package wf = CoolProp2Modelica.Media.R601_CP constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 0 + time*1E6;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
    end test_pentane_coolprop;

    model test_R245fa_baseproperties
      replaceable package wf = CoolProp2Modelica.Media.R245fa_CP constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.BaseProperties fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 2 + time*1E5;
      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);
    end test_R245fa_baseproperties;

    model test_R245fa_ThermodynamicState
    "ThermodynamicState is faster since it does not call the saturation properties"
      replaceable package wf = CoolProp2Modelica.Media.R245fa_CP constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium model";
      wf.ThermodynamicState fluid "Properties of the two-phase fluid";
      Modelica.SIunits.SpecificEnthalpy h;
      Modelica.SIunits.Pressure p;
      Modelica.SIunits.DerDensityByEnthalpy drdh
      "Derivative of average density by enthalpy";
      Modelica.SIunits.DerDensityByPressure drdp
      "Derivative of average density by pressure";
    equation
      p = 1E5;
      h = 2 + time*1E5;
      fluid = wf.setState_ph(p=p,h=h);
      drdp = wf.density_derp_h(fluid);
      drdh = wf.density_derh_p(fluid);
    end test_R245fa_ThermodynamicState;
  end Test;


  annotation (uses(Modelica(version="3.2"), TILMedia(version="2.1.4")),
                                             Documentation(info="<html>
<p><h4><font color=\"#008000\">CoolProp2Modelica</font></h4></p>
<p>The CoolProp2Modelica library provides an interface between the CoolProp thermodynamic library and Modelica. It is based on a modified version of ExternalMedia. The old functionalities of ExternalMedia, such as coupling with Fluidprop have been maintained and can still be used with CoolProp2Modelica.</p>
<p>Licensed by the Modelica Association under the Modelica License 2</p>
<p>Copyright &copy; 2006-2013 University of Li&egrave;ge, Belgium</p>
<p>Main contributors (in addition to the ExternalMedia developpers): Sylvain Quoilin, Ian Bell</p>
<p>Contact: squoilin@ulg.ac.be</p>
<p><br/><h4><font color=\"#008000\">ExternalMedia</font></h4></p>
<p>The ExternalMedia library provides a framework for interfacing external codes computing fluid properties to Modelica.Media-compatible component models. The two main requirements are: maximizing the efficiency of the code and minimizing the amount of extra code required to use your own external code within the framework. </p>
<p>The current version of the library supports pure fluids models, possibly two-phase, compliant with the <a href=\"modelica://Modelica.Media.Interfaces.PartialTwoPhaseMedium\">Modelica.Media.Interfaces.PartialTwoPhaseMedium</a> interface. </p>
<p>The releases of the library available on the Modelica website include a pre-compiled interface to the FluidProp tool, <a href=\"http://www.fluidprop.com\">http://www.fluidprop.com</a>. FluidProp features many built-in fluid models, and can optionally be used to access the whole NIST RefProp database, thus giving easy access to a wide range of fluid models with state-of-the-art accuracy. Make sure you download the Modelica package corresponding to the version of Microsoft Visual Studio that you use to compile the Modelica models, in order to avoid linker errors. </p>
<p>Before using the library, download and install the latest version of FluidProp from <a href=\"http://www.fluidprop.com\">http://www.fluidprop.com</a>. If you want to use the RefProp fluid models, you need to get the full version of FluidProp, which has an extra license fee. </p>
<p>You can now define medium models for all the libraries supported by FluidProp, by extending the <a href=\"modelica://ExternalMedia.Media.FluidPropMedium\">ExternalMedia.Media.FluidPropMedium</a> package. Set libraryName to FluidProp.RefProp, FluidProp.StanMix, FluidProp.TPSI, FluidProp.IF97, or FluidProp.GasMix, depending on the specific library you need to use. Set substanceNames to a single-element string array containing the name of the specific medium, as specified by the FluidProp documentation. Set mediumName to a string that describes the medium (this only used for documentation purposes but has no effect in selecting the medium model). See ExternalMedia.Examples for examples. Please note that the medium models IF97 and GasMix are already available natively in Modelica.Media as <a href=\"modelica://Modelica.Media.Water.StandardWater\">Water.StandardWater</a> and <a href=\"modelica://Modelica.Media.IdealGases.MixtureGases\">IdealGases.MixtureGases</a>, respectively - it is recommended to use the Modelica.Media models instead, since they are much faster to compute. </p>
<p>If you want to use your own fluid property computation code instead, then you need to check out the source code from the Modelica SVN server and add the interface to it. See the <a href=\"modelica://ExternalMedia/Resources/manual.pdf\">library documentation</a> for further details. </p>
<p>Please contact the main developer, Francesco Casella (<a href=\"mailto:casella@elet.polimi.it\">casella@elet.polimi.it</a>) if you have questions or suggestions for improvement, or if you want to help futher developing the library. </p>
<p>Licensed by the Modelica Association under the Modelica License 2</p>
<p>Copyright &copy; 2006-2012 Politecnico di Milano, TU Braunschweig, Politecnico di Torino</p>
<p>Main contributors: Francesco Casella, Christoph Richter, Roberto Bonifetto</p>
<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">http://www.modelica.org/licenses/ModelicaLicense2</a>.</i> </p>
</html>"));
end CoolProp2Modelica;
