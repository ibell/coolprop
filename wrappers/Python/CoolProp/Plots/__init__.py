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
    plt = PropsPlot(Ref, 'Ts', smin=Tmin, smax=Tmax, axis=axis)
    plt._draw_graph()
    if show:
        plt.show()
    else:
        plt._draw_graph()
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
    plt = PropsPlot(Ref, 'Ph', smin=Tmin, smax=Tmax, axis=axis)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def drawIsoLines(Ref, plot, which, iValues=[], num=0, show=False, axis=None):
    """
    Draw lines with constant values of type 'which' in terms of x and y as
    defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value
    between which 'num' get drawn.

    :Note:
        :func:`CoolProps.Plots.drawIsoLines` will be depreciated in future
        releases and replaced with :func:`CoolProps.Plots.IsoLines`

    Parameters
    -----------
    Ref : str
        The given reference fluid
    plot : str
        The plot type used
    which : str
        The iso line type
    iValues : list
        The list of constant iso line values
    num : int, Optional
        The number of iso lines
        (Default: 0 - Use iValues list only)
    show : bool, Optional
        Show the current plot
        (Default: False)
    axis : :func:`matplotlib.pyplot.gca()`, Optional
        The current axis system to be plotted to.
        (Default: create a new axis system)

    Examples
    --------
    >>> from matplotlib import pyplot
    >>> from CoolProp.Plots import Ts, drawIsoLines
    >>>
    >>> Ref = 'n-Pentane'
    >>> ax = Ts(Ref)
    >>> ax.set_xlim([-0.5, 1.5])
    >>> ax.set_ylim([300, 530])
    >>> quality = drawIsoLines(Ref, 'Ts', 'Q', [0.3, 0.5, 0.7, 0.8], axis=ax)
    >>> isobars = drawIsoLines(Ref, 'Ts', 'P', [100, 2000], num=5, axis=ax)
    >>> isochores = drawIsoLines(Ref, 'Ts', 'D', [2, 600], num=7, axis=ax)
    >>> pyplot.show()
    """
    isolines = IsoLines(Ref, plot, which, axis=axis)
    lines = isolines.draw_isolines(iValues, num)
    if show:
        isolines.show()
    return lines

