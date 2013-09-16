#Bring some functions into the Plots namespace for code concision
from __future__ import absolute_import

from .Plots import PropsPlot, IsoLines
from .SimpleCycles import SimpleCycle, TwoStage, EconomizedCycle


def Ts(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a temperature-entropy plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.Ts` will be depreciated in future releases
        and replaced with :func:`CoolProps.Plots.PropsPlot`

    Parameters
    -----------
    Ref : str
        The given reference fluid
    Tmin : float, Optional
        Minimum limit for the saturation line
    Tmax : float, Optional
        Maximum limit for the saturation line
    show : bool, Optional
        Show the current plot
        (Default: False)
    axis : :func:`matplotlib.pyplot.gca()`, Optional
        The current axis system to be plotted to.
        (Default: create a new axis system)

    Examples
    --------
    >>> from CoolProp.Plots import Ts
    >>> Ts('R290', show=True)

    >>> from CoolProp.Plots import Ts
    >>> Ts('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> Ts('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'Ts', smin=Tmin, smax=Tmax)
    if show:
        plt.show()
    return plt.axis


def Ph(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a pressure-enthalpy plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.Ph` will be depreciated in future releases
        and replaced with :func:`CoolProps.Plots.PropsPlot`

    Parameters
    -----------
    Ref : str
        The given reference fluid
    Tmin : float, Optional
        Minimum limit for the saturation line
    Tmax : float, Optional
        Maximum limit for the saturation line
    show : bool, Optional
        Show the current plot
        (Default: False)
    axis : :func:`matplotlib.pyplot.gca()`, Optional
        The current axis system to be plotted to.
        (Default: create a new axis system)

    Examples
    --------
    >>> from CoolProp.Plots import Ph
    >>> Ph('R290', show=True)

    >>> from CoolProp.Plots import Ph
    >>> Ph('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> Ph('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'Ph', smin=Tmin, smax=Tmax)
    if show:
        plt.show()
    return plt.get_axis()