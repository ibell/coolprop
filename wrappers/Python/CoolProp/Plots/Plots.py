# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

import math
import numpy

import CoolProp.CoolProp as CP

from .Common import BasePlot


def drawLines(Ref,lines,axis,plt_kwargs=None):
    """
    Just an internal method to systematically plot values from
    the generated 'line' dicts, method is able to cover the whole
    saturation curve. Closes the gap at the critical point and
    adds a marker between the two last points of bubble and
    dew line if they reach up to critical point.
    Returns the an array of line objects that can be used to change
    the colour or style afterwards.
    """
    if not plt_kwargs is None:
        for line in lines:
            line['opts'] = plt_kwargs
    plottedLines = []
    if len(lines)==2 and (
      'q' in str(lines[0]['type']).lower() and 'q' in str(lines[1]['type']).lower()
      ) and (
      ( 0 == lines[0]['value'] and 1 == lines[1]['type'] ) or ( 1 == lines[0]['value'] and 0 == lines[1]['type'] ) ):
        # We plot the saturation curve
        bubble = lines[0]
        dew = lines[1]
        line, = axis.plot(bubble['x'],bubble['y'],**bubble['opts'])
        plottedLines.extend([line])
        line, = axis.plot(dew['x'], dew['y'], **dew['opts'])
        plottedLines.extend([line])
        # Do we need to test if this is T or p?
        Tmax = min(bubble['kmax'],dew['kmax'])
        if Tmax>CP.Props(Ref,'Tcrit')-2e-5:
            axis.plot(numpy.r_[bubble['x'][-1],dew['x'][-1]],numpy.r_[bubble['y'][-1],dew['y'][-1]],**bubble['opts'])
            #axis.plot((bubble['x'][-1]+dew['x'][-1])/2.,(bubble['y'][-1]+dew['y'][-1])/2.,'o',color='Tomato')
    else:
        for line in lines:
            line, = axis.plot(line['x'],line['y'],**line['opts'])
            plottedLines.extend([line])

    return plottedLines


class IsoLines(BasePlot):
    def __init__(self, fluid_ref, graph_type, iso_type, **kwargs):
        BasePlot.__init__(self, fluid_ref, graph_type, **kwargs)

        if not isinstance(iso_type, str):
            raise TypeError("Invalid iso_type input, expected a string")

        iso_type = iso_type.upper()
        if iso_type not in self.COLOR_MAP.keys() and iso_type != 'Q':
            raise ValueError('This kind of isoline is not supported for a ' \
                             + str(graph_type) + \
                             ' plot. Please choose from '\
                             + str(self.COLOR_MAP.keys()) + ' or Q.')

        self.iso_type = iso_type
        self.lines = []
        self.plotted_lines = []

    def __set_axis_limits(self, swap_xy):
        """
        Generates limits for the axes in terms of x,y defined by 'plot'
        based on temperature and pressure.

        Returns a tuple containing ((xmin, xmax), (ymin, ymax))
        """
        # Get current axis limits, be sure to set those before drawing isolines
        # if no limits are set, use triple point and critical conditions
        X = [CP.Props(self.graph_type[1],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                      self.fluid_ref),
             CP.Props(self.graph_type[1],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tmin'),
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[1],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[1],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tmin'),
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                      self.fluid_ref)]

        Y = [CP.Props(self.graph_type[0],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                       self.fluid_ref),
             CP.Props(self.graph_type[0],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tmin') ,
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[0],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[0],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tmin') ,
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                      self.fluid_ref)]

        limits = [[min(X), max(X)], [min(Y), max(Y)]]
        if not self.axis.get_autoscalex_on():
            limits[0][0] = max([limits[0][0], min(self.axis.get_xlim())])
            limits[0][1] = min([limits[0][1], max(self.axis.get_xlim())])
            limits[1][0] = max([limits[1][0], min(self.axis.get_ylim())])
            limits[1][1] = min([limits[1][1], max(self.axis.get_ylim())])

        self.axis.set_xlim(limits[0])
        self.axis.set_ylim(limits[1])
        return limits

    def get_isolines(self, iso_range=[], num=None):
        """
        This is the core method to obtain lines in the dimensions defined
        by 'plot' that describe the behaviour of fluid 'Ref'. The constant
        value is determined by 'iName' and has the values of 'iValues'.

        'iValues' is an array-like object holding at least one element. Lines
        are calculated for every entry in 'iValues'. If the input 'num' is
        larger than the amount of entries in 'iValues', an internally defined
        pattern is used to calculate an appropriate line spacing between the maximum
        and minimum values provided in 'iValues'.

        Returns lines[num] - an array of dicts containing 'x' and 'y'
        coordinates for bubble and dew line. Additionally, the dict holds
        the keys 'label' and 'opts', those can be used for plotting as well.
        """
        if not iso_range or len(iso_range) == 1:
            raise ValueError('Automatic interval detection for isoline \
                              boundaries is not supported yet, use the \
                              iso_range=[min, max] parameter.')

        if len(iso_range) == 2 and num is None:
            raise ValueError('Please specify the number of isoline you want \
                              e.g. num=10')

        iso_range = numpy.sort(numpy.unique(iso_range))

        def generate_ranges(xmin, xmax, num):
            if self.iso_type in ['P', 'D']:
                return numpy.logspace(math.log(xmin, 2.),
                                      math.log(xmax, 2.),
                                      num=num,
                                      base=2.)
            return numpy.linspace(xmin, xmax, num=num)

        # Generate iso ranges
        if len(iso_range) == 2:
            iso_range = generate_ranges(iso_range[0], iso_range[1], num)
            #iso_range = plotRound(iso_range)
        #else:
        #    TODO: Automatic interval detection
        #    iVal = [CP.Props(iName,'T',T_c[i],'D',rho_c[i],Ref) for i in range(len(T_c))]
        #    iVal = patterns[iName]([numpy.min(iVal),numpy.max(iVal),num])

        switch_xy_map = {'D': ['TS', 'PH', 'PS'],
                         'S': ['PH', 'PD', 'PT'],
                         'T': ['PH', 'PS'],
                         'H': ['PD']}
        #TS: TD is defined, SD is not
        #PH: PD is defined, HD is not
        #PS: PD is defined, SD is not
        #PH: PS is more stable than HS
        #PD: PS is defined, DS is not
        #PT: PS is defined, TS is not
        #PH: PT is defined, HT is not
        #PS: PT is defined, ST is not
        #PD: PH is defined, DH is not

        iso_error_map = {'TD': ['S', 'H'],
                         'HS': ['T', 'D'],}

        switch_xy = False
        if self.iso_type in ['D', 'S', 'T', 'H']:
            if self.graph_type in switch_xy_map[self.iso_type]:
                switch_xy = True

        if self.graph_type in ['TD', 'HS']:
            if self.iso_type in iso_error_map[self.graph_type]:
                raise ValueError('You should not reach this point!')

        axis_limits = self.__set_axis_limits(switch_xy)
        req_prop = self.graph_type[0]
        prop2_name = self.graph_type[1]
        if switch_xy:
            axis_limits.reverse()
            req_prop = self.graph_type[1]
            prop2_name = self.graph_type[0]

        # Calculate the points
        if self.iso_type == 'Q':
            lines = self._get_sat_lines(x=iso_range)
            return lines

        # TODO: Determine saturation state if two phase region present
        x_range = numpy.linspace(axis_limits[0][0], axis_limits[0][1], 1000.)
        x_mesh = [x_range for i in iso_range]

        plot_data = self._get_fluid_data(req_prop,
                                         self.iso_type,
                                         prop2_name,
                                         iso_range,
                                         x_mesh)

        if switch_xy:
            plot_data.reverse()

        for j in range(len(plot_data[0])):
            line = {'x': plot_data[0][j],
                    'y': plot_data[1][j],
                    # TODO
                    'label': "%.1f %s" % (iso_range[j],
                                          self.AXIS_LABELS[self.iso_type][1]),
                    'type': self.iso_type,
                    'opts': {'color': self.COLOR_MAP[self.iso_type],
                             'lw': 0.75,
                             'alpha': 0.5}}
            self.lines.append(line)
        return self.lines

    def draw_isolines(self, iso_range, num=None):
        """
        Draw lines with constant values of type 'which' in terms of x and y as
        defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value between
        which 'num' get drawn.

        There should also be helpful error messages...
        """
        if not iso_range or len(iso_range) == 1:
            raise ValueError('Automatic interval detection for isoline \
                              boundaries is not supported yet, use the \
                              iso_range=[min, max] parameter.')

        if len(iso_range) == 2 and num is None:
            raise ValueError('Please specify the number of isoline you want \
                              e.g. num=10')

        if self.iso_type == 'all':
            raise ValueError('Plotting all lines automatically is not \
                              supported, yet..')

        if self.iso_type != 'all':
            lines = self.get_isolines(iso_range, num)
            self.plotted_lines = drawLines(self.fluid_ref, lines, self.axis)
            self._plot_default_annotations()
            return self.plotted_lines
        #else:
        #    # TODO: assign limits to values automatically
        #    ll = _getIsoLineIds(plot)
        #    if not len(ll)==len(iValues):
        #        raise ValueError('Please provide a properly sized array of bounds.')
        #    for c,l in enumerate(ll):
        #        drawIsoLines(Ref, plot, l, iValues=iValues[c], num=num, axis=axis, fig=fig)

    def add_inline_labels(self, xloc=None, method='strict'):
        bbox_opts = dict(boxstyle='square,pad=0.2',
                         fc='white', ec='None', alpha = 0.9)

        if not self.lines:
            raise ValueError('asdasd')

        for line in self.lines:
            x_vals = line['x']
            y_vals = line['y']

            if method == 'strict':
                if xloc is None:
                    raise ValueError(''.join(["Please enter a xloc value ",
                                              "for method %s" % method]))
                loc = numpy.argmin(numpy.abs(x_vals - xloc))
                theta = math.atan((y_vals[loc+1] - y_vals[loc]) /
                                  (x_vals[loc+1] - x_vals[loc]))
                text_loc = numpy.array([x_vals[loc], y_vals[loc]])

            elif method == 'auto':
                pass

            transform_func = self.axis.transData.transform_angles
            trans_angle = transform_func(numpy.array((theta,)),
                                         text_loc.reshape((1, 2)),
                                         radians=True)

            self.axis.text(text_loc[0], text_loc[1],
                           line['label'],
                           verticalalignment='center',
                           horizontalalignment='center',
                           color=self.COLOR_MAP[self.iso_type],
                           bbox=bbox_opts,
                           rotation=math.degrees(trans_angle))


class PropsPlot(BasePlot):
    def __init__(self, fluid_ref, graph_type, **kwargs):
        """
        Create graph for the specified fluid properties

        Parameters
        -----------
        fluid_ref : str
            The refence fluid
        graph_type : str
            The graph type to be plotted
        axis : :func:`matplotlib.pyplot.gca()`, Optional
            The current axis system to be plotted to.
            Default: create a new axis system
        fig : :func:`matplotlib.pyplot.figure()`, Optional
            The current figure to be plotted to.
            Default: create a new figure

        Examples
        ---------
        >>> from CoolProp.Plots import PropsPlot
        >>> plt = PropsPlot('Water', 'Ph')
        >>> plt.show()

        >>> plt = PropsPlot('n-Pentane', 'Ts')
        >>> plt.set_axis_limits([-0.5, 1.5, 300, 530])
        >>> plt.draw_isolines('Q', [0.1, 0.9])
        >>> plt.draw_isolines('P', [100, 2000])
        >>> plt.draw_isolines('D', [2, 600])
        >>> plt.show()

        .. note::

            See the online documentation for a the available fluids and
            graph types
        """
        BasePlot.__init__(self, fluid_ref, graph_type, **kwargs)

        self.smin = kwargs.get('smin', None)
        self.smax = kwargs.get('smax', None)

    def __draw_region_lines(self):
        lines = self._get_sat_lines(kind='T',
                                    smin=self.smin,
                                    smax=self.smax)
        drawLines(self.fluid_ref, lines, self.axis)

    def _draw_graph(self):
        self.__draw_region_lines()
        self._plot_default_annotations()

    def draw_isolines(self, iso_type, iso_range, num=10):
        iso_lines = IsoLines(self.fluid_ref,
                             self.graph_type,
                             iso_type,
                             axis=self.axis)
        return iso_lines.draw_isolines(iso_range, num)


def Ts(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a temperature-entropy plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.Ts` will be deprecated in future releases
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
        :func:`CoolProps.Plots.Ph` will be deprecated in future releases
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


def Ps(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a pressure-entropy plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.Ps` will be deprecated in future releases
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
    >>> from CoolProp.Plots import Ps
    >>> Ps('R290', show=True)

    >>> from CoolProp.Plots import Ps
    >>> Ps('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> Ps('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'Ps', smin=Tmin, smax=Tmax, axis=axis)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def PT(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a pressure-temperature plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.PT` will be deprecated in future releases
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
    >>> from CoolProp.Plots import PT
    >>> PT('R290', show=True)

    >>> from CoolProp.Plots import PT
    >>> PT('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> PT('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'PT', smin=Tmin, smax=Tmax, axis=axis)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def Prho(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a pressure-density plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.Prho` will be deprecated in future releases
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
    >>> from CoolProp.Plots import Prho
    >>> Prho('R290', show=True)

    >>> from CoolProp.Plots import Prho
    >>> Prho('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> Prho('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'PD', smin=Tmin, smax=Tmax, axis=axis)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def Trho(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a temperature-density plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.Trho` will be deprecated in future releases
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
    >>> from CoolProp.Plots import Trho
    >>> Trho('R290', show=True)

    >>> from CoolProp.Plots import Trho
    >>> Trho('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> Trho('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'TD', smin=Tmin, smax=Tmax, axis=axis)
    if show:
        plt.show()
    else:
        plt._draw_graph()
    return plt.axis

def hs(Ref, Tmin=None, Tmax=None, show=False, axis=None, *args, **kwargs):
    """
    Make a enthalpy-entropy plot for the given fluid

    :Note:
        :func:`CoolProps.Plots.hs` will be deprecated in future releases
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
    >>> from CoolProp.Plots import hs
    >>> hs('R290', show=True)

    >>> from CoolProp.Plots import hs
    >>> hs('R290', show=True, Tmin=200, Tmax=300)

    >>> from matplotlib import pyplot
    >>> fig = pyplot.figure(1)
    >>> ax = fig.gca()
    >>> hs('R290', show=True, axis=ax)
    """
    plt = PropsPlot(Ref, 'hs', smin=Tmin, smax=Tmax, axis=axis)
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
