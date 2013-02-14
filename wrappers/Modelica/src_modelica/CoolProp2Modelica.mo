within ;
package CoolProp2Modelica
  import SI = Modelica.SIunits;

  package Common "Package with common definitions"
    type InputChoice = enumeration(
        dT "(d,T) as inputs",
        ph "(p,h) as inputs",
        ps "(p,s) as inputs",
        pT "(p,T) as inputs");
  end Common;

  package Media "Medium packages compatible with Modelica.Media"

    package ExternalTwoPhaseMedium "Generic external two phase medium package"

      extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(
        mediumName = "ExternalMedium",
        singleState = false,
        onePhase = false,
        smoothModel = false,
        fluidConstants = {externalFluidConstants});
      import CoolProp2Modelica.Common.InputChoice;

      constant String libraryName = "UnusableExternalMedium"
        "Name of the external fluid property computation library";
      final constant String substanceName = substanceNames[1]
        "Only one substance can be specified";

      constant FluidConstants externalFluidConstants = FluidConstants(
        iupacName=  "unknown",
        casRegistryNumber=  "unknown",
        chemicalFormula=  "unknown",
        structureFormula=  "unknown",
        molarMass=  getMolarMass(),
        criticalTemperature=  getCriticalTemperature(),
        criticalPressure=  getCriticalPressure(),
        criticalMolarVolume=  getCriticalMolarVolume(),
        acentricFactor=  0,
        triplePointTemperature=  280.0,
        triplePointPressure=  500.0,
        meltingPoint=  280,
        normalBoilingPoint=  380.0,
        dipoleMoment=  2.0);

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
        Integer phase "phase flag: 2 for two-phase, 1 for one-phase";
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
          annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
      end getMolarMass;

      replaceable partial function getCriticalTemperature
        output Temperature Tc "Critical temperature";
        external "C" Tc = TwoPhaseMedium_getCriticalTemperature_(mediumName, libraryName, substanceName)
          annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
      end getCriticalTemperature;

      replaceable partial function getCriticalPressure
        output AbsolutePressure pc "Critical temperature";
        external "C" pc = TwoPhaseMedium_getCriticalPressure_(mediumName, libraryName, substanceName)
          annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
      end getCriticalPressure;

      replaceable partial function getCriticalMolarVolume
        output MolarVolume vc "Critical molar volume";
        external "C" vc = TwoPhaseMedium_getCriticalMolarVolume_(mediumName, libraryName, substanceName)
          annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
      end setState_ps;

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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
        annotation(derivative(noDerivative = phase) = density_pT_der,
                   Inline = true);
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
        annotation(derivative(noDerivative = phase) = density_ps_der,
                   Inline = true);
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
*/

        annotation(Inline = true);
      end density_ph_der;

      redeclare replaceable function extends isentropicEnthalpy
      external "C" h_is = TwoPhaseMedium_isentropicEnthalpy_(p_downstream, refState,
       mediumName, libraryName, substanceName)
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
      end isentropicEnthalpy;

      redeclare replaceable function setSat_p
        "Return saturation properties from p"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        output SaturationProperties sat "saturation property record";
      external "C" TwoPhaseMedium_setSat_p_(p, sat, mediumName, libraryName, substanceName)
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
      end setSat_p;

      redeclare replaceable function setSat_T
        "Return saturation properties from p"
        extends Modelica.Icons.Function;
        input Temperature T "temperature";
        output SaturationProperties sat "saturation property record";
      external "C" TwoPhaseMedium_setSat_T_(T, sat, mediumName, libraryName, substanceName)
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end setDewState;

      redeclare replaceable function extends saturationTemperature
        // Standard definition
      algorithm
        T :=saturationTemperature_sat(setSat_p(p));
        /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_saturationTemperature_(p, mediumName, libraryName, substanceName)
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
*/
        annotation(Inline = true);
      end saturationTemperature;

      redeclare function extends saturationTemperature_sat

        annotation(Inline = true);
      end saturationTemperature_sat;

      redeclare replaceable function extends saturationTemperature_derp "Returns derivative of saturation temperature w.r.t.. pressureBeing this function inefficient, it is strongly recommended to use saturationTemperature_derp_sat
     and never use saturationTemperature_derp directly"
      external "C" dTp = TwoPhaseMedium_saturationTemperature_derp_(p, mediumName, libraryName, substanceName)
        annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
*/

        annotation(Inline = true);
      end dewEnthalpy;

      redeclare replaceable function extends saturationPressure
        // Standard definition
      algorithm
        p :=saturationPressure_sat(setSat_T(T));
        /*  // If special definition in "C"
  external "C" p=  TwoPhaseMedium_saturationPressure_(T, mediumName, libraryName, substanceName)
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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

    package FluidPropMedium "FluidProp medium package"

      extends ExternalTwoPhaseMedium;

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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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
    annotation(Include="#include \"coolproplib.h\"", Library="CoolPropLib");
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

    package TestMedium "Simple water medium model for debugging and testing"
      extends ExternalTwoPhaseMedium(
        mediumName = "TestMedium",
        libraryName = "TestMedium",
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.pT);
    end TestMedium;
  end Media;

  package Examples "Examples of external medium models"
    package WaterIF95
      extends CoolProp2Modelica.Media.FluidPropMedium(
        mediumName="Water",
        libraryName="FluidProp.RefProp",
        substanceNames={"H2O"},
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end WaterIF95;

    package WaterTPSI
      extends CoolProp2Modelica.Media.FluidPropMedium(
        mediumName="Water",
        libraryName="FluidProp.TPSI",
        substanceNames={"H2O"},
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end WaterTPSI;

    package CarbonDioxide
      extends CoolProp2Modelica.Media.FluidPropMedium(
        mediumName="Carbon Dioxide",
        libraryName="FluidProp.RefProp",
        substanceNames={"CO2"},
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end CarbonDioxide;

    package Propane_CoolProp "Computation of Propane Properties using CoolProp"
      extends CoolProp2Modelica.Media.ExternalTwoPhaseMedium(
        mediumName="TestMedium",
        libraryName="CoolProp",
        substanceName="propane",
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end Propane_CoolProp;

    package Propane_Refprop
      "Propane properties using Refprop through FluidProp (requires the full version of FluidProp)"
      extends CoolProp2Modelica.Media.ExternalTwoPhaseMedium(
        mediumName="TestMedium",
        libraryName="FluidProp.RefProp",
        substanceName="propane",
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end Propane_Refprop;

    package Propane_FluidProp
      "Propane properties using the StanMix library of FluidProp"
      extends CoolProp2Modelica.Media.ExternalTwoPhaseMedium(
        mediumName="TestMedium",
        libraryName="FluidProp.StanMix",
        substanceName="propane",
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end Propane_FluidProp;

    package Solkatherm "Solkatherm properties using CoolProp"
      extends CoolProp2Modelica.Media.ExternalTwoPhaseMedium(
        mediumName="SES36",
        libraryName="CoolProp",
        substanceName="SES36",
        ThermoStates=Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      annotation ();
    end Solkatherm;

    package Solkatherm_debug
      "Solkatherm properties using CoolProp with debug option"
      extends CoolProp2Modelica.Media.ExternalTwoPhaseMedium(
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

    package R245fa
    //  extends Modelica.Media.Water.StandardWater;
      extends CoolProp2Modelica.Media.ExternalTwoPhaseMedium(
        mediumName = "R245fa",
        libraryName = "CoolProp",
        substanceNames = {"R245fa"});

        redeclare function extends density_pT_der
        algorithm
        assert(false, "Error: R245fa.density_pT_der is not implemented");
        d_der := 0;
        end density_pT_der;

    end R245fa;
  end Examples;

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
              CoolProp2Modelica.Media.FluidPropMedium;
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
            redeclare package Medium = CoolProp2Modelica.Examples.WaterIF95);
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
            redeclare package Medium = CoolProp2Modelica.Examples.WaterIF95);

        equation
          p1 = 1e5+1e5*time;
          h1 = 1e5;
          p2 = 1e5;
          h2 = 1e5 + 2e5*time;
        end TestBasePropertiesExplicit;

        model TestBasePropertiesImplicit
          "Test case using BaseProperties and implicit equations"
          extends GenericModels.TestBasePropertiesImplicit(
            redeclare package Medium = CoolProp2Modelica.Examples.WaterIF95,
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
          redeclare package Medium = CoolProp2Modelica.Examples.WaterIF95,
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
                CoolProp2Modelica.Examples.WaterIF95,
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
                CoolProp2Modelica.Examples.WaterIF95,
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
                CoolProp2Modelica.Examples.WaterIF95,
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
      package Medium = Media.ExternalTwoPhaseMedium;
      Medium.BaseProperties medium;
    equation
      medium.p = 1e5;
      medium.h = 1e5;
    end TestWrongMedium;
    end WrongMedium;

    model test_propane_coolprop
      replaceable package wf = CoolProp2Modelica.Examples.Propane_CoolProp constrainedby
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
      replaceable package wf = CoolProp2Modelica.Examples.Propane_FluidProp constrainedby
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
      replaceable package wf = CoolProp2Modelica.Examples.Propane_Refprop constrainedby
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
      replaceable package wf = CoolProp2Modelica.Examples.Solkatherm constrainedby
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
      replaceable package wf = CoolProp2Modelica.Examples.Solkatherm_debug constrainedby
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
      h = 0 + time*3e5;

      fluid.p = p;
      fluid.h = h;
      drdp = wf.density_derp_h(fluid.state);
      drdh = wf.density_derh_p(fluid.state);

      annotation (experiment(Algorithm="Dassl"), __Dymola_experimentSetupOutput);
    end test_solkatherm_debug;
  end Test;
  annotation (uses(Modelica(version="3.2")), Documentation(info="<html>
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
