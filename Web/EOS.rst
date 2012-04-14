.. _Fluid-Properties:

Fluid Properties
================

Introduction
------------

If you are feeling impatient, jump to :ref:`Props_Sample`, otherwise, hang in there.

Nearly all the fluids modeling in CoolProp are based on Helmholtz function formulations.  This is a convenient construction of the equation of state because all the thermodynamic properties of interest can be obtained directly from partial derivatives of the Helmholtz energy.

It should be noted that the EOS are typically valid over the entire range of the fluid, from subcooled liquid to superheated vapor, to supercritical fluid.  

Annoyingly, different authors have selected different sets of nomenclature for the Helmholtz energy.  For consistency, the nomenclature of Lemmon will be used here.  Also, some authors present results on a mole-basis or mass-basis, further complicating comparisons.

Mathematical Description
------------------------
In general, the EOS are based on non-dimensional terms :math:`\delta` and :math:`\tau`, where these terms are defined by

.. math::

    \delta=\rho/\rho_c
    
    \tau=T_c/T
    
where :math:`\rho_c` and :math:`T_c` are the critical density of the fluid if it is a pure fluid.  For pseudo-pure mixtures, the critical point is typically not used as the reducing state point, and often the maximum condensing temperature on the saturation curve is used instead.

The non-dimensional Helmholtz energy of the fluid is given by

.. math::

    \alpha=\alpha^0+\alpha^r
    
where :math:`\alpha^0` is the ideal-gas contribution to the Helmholtz energy, and :math:`\alpha^r` is the residual Helmholtz energy contribution which accounts for non-ideal behavior.  For a given set of :math:`\delta` and :math:`\tau`, each of the terms :math:`\alpha^0` and :math:`\alpha^r` are known.  The exact form of the Helmholtz energy terms is fluid dependent, but a relatively simple example is that of Nitrogen, which has the ideal-gas Helmholtz energy of

.. math::

    \alpha^0=\ln\delta+a_1\ln\tau+a_2+a_3\tau+a_4\tau^{-1}+a_5\tau^{-2}+a_6\tau^{-3}+a_7\ln[1-\exp(-a_8\tau)]
    
and the non-dimensional residual Helmholtz energy of

.. math::

    \alpha^r=\sum_{k=1}^{6}{N_k\delta^{i_k}\tau^{j_k}}+\sum_{k=7}^{32}{N_k\delta^{i_k}\tau^{j_k}\exp(-\delta^{l_k})}+\sum_{k=33}^{36}{N_k\delta^{i_k}\tau^{j_k}\exp(-\phi_k(\delta-1)^2-\beta_k(\tau-\gamma_k)^2)}
    
and all the terms other than :math:`\delta` and :math:`\tau` are fluid-dependent correlation parameters.

The other thermodynamic parameters can then be obtained through analytic derivatives of the Helmholtz energy terms.  For instance, the pressure is given by

.. math::

    p=\rho RT\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau} \right]
    
and the specific internal energy by

.. math::

    \frac{u}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right]

and the specific enthalpy by

.. math::

    \frac{h}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] +\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+1

which can also be written as

.. math::

    \frac{h}{RT}=\frac{u}{RT}+\frac{p}{\rho RT}
    
The specific entropy is given by

.. math::

    \frac{s}{R}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right]-\alpha^0-\alpha^r
    
and the specific heats at constant volume and constant pressure respectively are given by

.. math::

    \frac{c_v}{R}=-\tau^2 \left[\left(\frac{\partial^2 \alpha^0}{\partial \tau^2}\right)_{\delta}+ \left(\frac{\partial^2 \alpha^r}{\partial \tau^2}\right)_{\delta} \right]
    
    \frac{c_p}{R}=\frac{c_v}{R}+\dfrac{\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}-\delta\tau\left(\frac{\partial^2 \alpha^r}{\partial \delta\partial\tau}\right)\right]^2}{\left[1+2\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+\delta^2\left(\frac{\partial^2 \alpha^r}{\partial \delta^2}\right)_{\tau}\right]}
    
The EOS is set up with temperature and density as the two independent properties, but often other inputs are known, most often temperature and pressure because they can be directly measured.  As a result, if the density is desired for a known temperature and pressure, it can be obtained iteratively.  The following algorithm is used to obtain a reasonable guess for the initial value for the iterative solver:

#. If the fluid is superheated, use a guess of ideal gas (:math:`\rho=p/(RT)`)
#. If the fluid is subcooled, use a guess of saturated liquid density
#. If the fluid is supercritical, use a guess of ideal gas (:math:`\rho=p/(RT)`)
#. No solution for density as a function of temperature and pressure if the fluid is two-phase

Saturation State
----------------

If the fluid is somewhere in the two-phase region, or saturation state properties are desired, saturated liquid and vapor properties can be obtained.  At equilibrium, the Gibbs function of the liquid and vapor are equal, as are the pressures of the saturated liquid and vapor.  For nearly all pure fluids, ancillary equations for the density of saturated liquid and saturated vapor as a function of temperature are provided, given by :math:`\rho'` and :math:`\rho''` respectively.  Thus for pure fluids, for a given temperature, initial guesses for the densities of saturated liquid and vapor are given by 
:math:`\rho'` and :math:`\rho''`.  Using one of the densities, a guess for the saturation pressure can be obtained.  Then, the saturation pressure is iteratively altered using a numerical method.  For each saturation pressure, the saturated liquid and vapor densities are updated using the full EOS to match the imposed temperature and guessed pressure.  Because the density is known explicitly from the EOS, Newton's method can be used to update the densities.  For Newton's method, the derivative :math:`\partial \rho/\partial p` is needed, which can be given explicitly as

.. math::

    \frac{\partial p}{\partial \rho}=RT\left[1+2\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+\delta^2\left(\frac{\partial^2 \alpha^r}{\partial \delta^2}\right)_{\tau}\right]
    
and the value for :math:`\rho` is updated by employing

.. math::

    \rho_{new}=\rho_{old}-\frac{p(T,\rho_{old})-p_{guess}}{\frac{\partial p}{\partial \rho}(T,\rho_{old})}
    
until :math:`\left|p(T,\rho_{old})-p_{guess}\right|` is sufficiently small.  Then the numerical method calculates the Gibbs function for saturated liquid and saturated vapor, and uses the difference in Gibbs functions to update the guess for the saturation pressure.  Then the densities are calculated again.  At convergence, the set of :math:`\rho'`, :math:`\rho''`, and :math:`p_{sat}` are known for a given saturation temperature.  If the fluid is not a pure fluid, the best that you can do is to use the ancillary equations to calculate the saturation densities and saturation pressure.

As you might imagine, doing all this work to calculate the saturation state for pure fluids is computationally *very* expensive, so a lookup table method has been implemented for the saturation densities and saturation pressure.  From Python, you can turn on the saturation lookup table with::

    UseSaturationLUT(1)
    
or use the full EOS by calling::

    UseSaturationLUT(0)

.. _Props_Sample:

Sample Code
-----------

.. ipython::

    #Import the things you need 
    In [1]: from CoolProp.CoolProp import Props, UseSaturationLUT,UseSinglePhaseLUT
    
    #Specific heat (kJ/kg/K) of 20% ethylene glycol as a function of T
    In [2]: h=Props('C','T',298.15,'P',101.325,'EG-20%'); print h
    
    #Density of Air at STP in kg/m^3
    In [2]: Props('D','T',298.15,'P',101.325,'Air')
    
    #Saturation temperature of Water at 1 atm
    In [2]: Props('T','P',101.325,'Q',0,'Water')
    
    #Saturated vapor density of R134a at 0C
    In [2]: Props('H','T',273.15,'Q',1,'R134a')
    
    # -------------------------------------------------------
    #  Single phase lookup table
    # -------------------------------------------------------
    
    #Crudely time 100 calls to get saturation temperature without lookup table
    In [2]: from time import clock; t1=clock()
    
    In [2]: for i in range(100):
       ...:      T=Props('T','P',101.325,'Q',1,'Water')
       ...:
    
    In [3]: print 'time elapsed for 100 calls:',clock()-t1,'s'
    
    #Turn on the saturation LUT
    In [3]: UseSaturationLUT(1)
    
    #Crudely time 100 calls to get saturation temperature with lookup table
    In [2]: from time import clock; t1=clock()
    
    In [3]: Props('T','P',101.325,'Q',1,'Water')
    
    In [3]: print 'time to build LUT:',clock()-t1,'s'
    
    In [2]: t1=clock()
    
    In [2]: for i in range(100):
       ...:      T=Props('T','P',101.325,'Q',1,'Water')
       ...:
    
    In [3]: print 'time elapsed for 100 calls:',clock()-t1,'s'
    
    # -------------------------------------------------------
    #  Single phase lookup table
    # -------------------------------------------------------
    #Turn off (this is default) the single-phase LUT
    In [3]: UseSinglePhaseLUT(0)
    
    #Crudely time 10000 calls to get enthalpy without lookup table
    #Using full equation of state.  First get density, then H=f(T,rho)
    In [2]: from time import clock; t1=clock()
    
    In [2]: for i in range(10000):
       ...:      H=Props('H','T',290,'P',101.325,'R744')
       ...:
    
    In [3]: time_no_LUT=clock()-t1
    
    In [3]: print 'time elapsed for 10000 calls:',time_no_LUT,'s'
    
    #Turn on the single-phase LUT
    In [3]: UseSinglePhaseLUT(1)
    
    #Crudely time 10000 calls to get enthalpy with lookup table
    In [2]: from time import clock; t1=clock()
    
    In [3]: H=Props('H','T',290,'P',101.325,'R744')
    
    In [3]: print 'time to build LUT:',clock()-t1,'s'
    
    In [2]: t1=clock()
    
    In [2]: for i in range(10000):
       ...:      H=Props('H','T',290,'P',101.325,'R744')
       ...:
       
    In [3]: time_with_LUT=clock()-t1
    
    In [3]: print 'time elapsed for 10000 calls:',time_with_LUT,'s'
    
    In [3]: print 'speedup factor with LUT:',time_no_LUT/time_with_LUT,'x'
    
    #Note: CO2 has a very involved EOS, so this is perhaps an extreme example
    
    
    
    
Code Documentation
------------------

.. automodule:: CoolProp.CoolProp
    :members:
    :undoc-members: