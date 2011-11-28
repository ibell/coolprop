Humid Air
*********

Humid air can be modeled as a mixture of air and water vapor.  In the simplest analysis, water and air are treated as ideal gases but in principle there is interaction between the air and water molecules that must be included through the use of interaction parameters.

Because humid air is a mixture of dry air (treated as a pseudo-pure gas) and water vapor (treated as a real gas), three variables are required to fix the state by the state postulate.

In the analysis that follows, the three parameters that are ultimately needed to calculate everything else are the dry bulb temperature :math:`T`, the total pressure :math:`p`, and the molar fraction of water :math:`\psi_w`.  The molar fraction of air is simply :math:`\psi_a=1-\psi_w`.

Of course, it is not so straightforward to measure the mole fraction of water vapor molecules, so other measures are used.  There are three different variables that can be used to obtain the mole fraction of water vapor without resorting to iterative methods.

1. Humidity ratio

The humidity ratio :math:`W` is the ratio of the mass of water vapor to the mass of air in the mixture.  Thus the mole fraction of water can be obtained from

.. math::

    \psi_w=\frac{n_w}{n}=\frac{n_w}{n_a+n_w}=\frac{m_w/M_w}{m_a/M_a+m_w/M_w}=\frac{m_w}{(M_w/M_a)m_a+m_w}=\frac{1}{(M_w/M_a)/W+1}=\frac{W}{(M_w/M_a)+W}
    
or

.. math::

    \psi_w=\frac{W}{\varepsilon+W}

where the ratio of mole masses :math:`\varepsilon` is given by :math:`\varepsilon=M_w/M_a`

2. Relative Humidity

The relative humidity :math:`\varphi` is defined as the ratio of the mole fraction of water in the humid air to the saturation mole fraction of water.  Because of the presence of air with the water, the pure water saturated vapor pressure :math:`p_{w,s}` must be multiplied by an enhancement factor :math:`f` that is very close to one near atmospheric conditions.

Mathematically, the result is

.. math::

    \varphi=\frac{\psi_w}{\psi_{w,s}}

where 

.. math::

    \psi_{w,s}=\frac{fp_{w,s}}{p}
    
The product :math:`p_s` is defined by :math:`p_s=fp_{w,s}`, which yields the result for :math:`\psi_w` of

.. math::

    \varphi=\frac{\psi_w}{p_s/p}
    
.. math::

    \psi_w=\frac{\varphi p_s}{p}

3. Dewpoint temperature

The dewpoint temperature is defined as the temperature at which the actual vapor pressure of water is equal to the saturation vapor pressure.  At the given dewpoint, the vapor pressure of water is given by

.. math::

    p_w=f(p,T_{dp})p_{w,s}(T_{dp})

and the mole fraction of water vapor is obtained from

.. math::

    \psi_w=\frac{p_w}{p}
    
Once the state has been fixed by a set of :math:`T,p,\psi_w`, any parameter of interest can be calculated

Enhancement factor
------------------

The enhancement factor is a parameter that includes the impact of the air on the saturation pressure of water vapor.  It is only a function of temperature and pressure, but it must be iteratively obtained due to the nature of the expression for the enhancement factor.


:math:`\psi_{w,s}` is given by :math:`\psi_{w,s}=fp_{w,s}/p`, where :math:`f` can be obtained from 

.. math::

    \ln(f)=\left[ \begin{array}{l}\left [ \dfrac{(1+k_Tp_{w,s})(p-p_{w,s})-k_T\dfrac{(p^2-p_{w,s}^2)}{2}}{\overline {R} T}\right] \bar v_{w,s}+\ln[1-\beta_H(1-\psi_{w,s})p]\\
    +\left[\dfrac{(1-\psi_{w,s})^2p}{\bar R T}\right] B_{aa}-2\left[\dfrac{(1-\psi_{w,s})^2p}{\bar R T}\right]B_{aw}-\left[\dfrac{(p-p_{w,s}-(1-\psi_{w,s})^2p)}{\bar R T}\right]B_{ww} \\
    +\left[\dfrac{(1-\psi_{w,s})^3 p^2}{(\bar R T)^2}\right] C_{aaa}+\left[\dfrac{3(1-\psi_{w,s})^2[1-2(1-\psi_{w,s})]p^2}{2(\bar R T)^2}\right]C_{aaw}\\
    -\left[\dfrac{3(1-\psi_{w,s})^2\psi_{w,s}p^2}{(\bar R T)^2}\right]C_{aww}-\left[\dfrac{(3-2\psi_{w,s})\psi_{w,s}^2p^2-p_{w,s}^2}{2(\bar R T)^2}\right]C_{www}\\
    -\left[\dfrac{(1-\psi_{w,s})^2(-2+3\psi_{w,s})\psi_{w,s}p^2}{(\bar R T)^2}\right]B_{aa}B_{ww}\\
    -\left[\dfrac{2(1-\psi_{w,s})^3(-1+3\psi_{w,s})p^2}{(\bar R T)^2}\right]B_{aa}B_{aw}\\
    +\left[\dfrac{6(1-\psi_{w,s})^2\psi_{w,s}^2p^2}{(\bar R T)^2}\right]B_{ww}B_{aw}-\left[\dfrac{3(1-\psi_{w,s})^4p^2}{2(\bar R T)^2}\right]B_{aa}^2\\
    -\left[\dfrac{2(1-\psi_{w,s})^2\psi_{w,s}(-2+3\psi_{w,s})p^2}{(\bar R T)^2}\right]B_{aw}^2-\left[\dfrac{p_{w,s}^2-(4-3\psi_{w,s})(\psi_{w,s})^3p^2}{2(\bar R T)^2}\right]B_{ww}^2
    \end{array}\right]
    

.. ipython::

    In [1]: execfile('Validation/HAValidation.py')
    
Code Documentation
------------------

.. automodule:: CoolProp.HumidAirProp
    :members:
    :undoc-members: