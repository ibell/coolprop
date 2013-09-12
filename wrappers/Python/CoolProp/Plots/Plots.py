
import numpy, matplotlib, matplotlib.pyplot, math, re
from scipy.interpolate import interp1d

import CoolProp.CoolProp as CP
from CoolProp.Plots.Common import BasePlot


def InlineLabel(xv,yv,x = None, y= None, axis = None, fig = None):
    """
    This will give the coordinates and rotation required to align a label with
    a line on a plot
    """

    def ToPixelCoords(xv,yv,axis,fig):
        [Axmin,Axmax]=axis.get_xlim()
        [Aymin,Aymax]=axis.get_ylim()
        DELTAX_axis=Axmax-Axmin
        DELTAY_axis=Aymax-Aymin

        width=fig.get_figwidth()
        height=fig.get_figheight()
        pos=axis.get_position().get_points()
        [[Fxmin,Fymin],[Fxmax,Fymax]]=pos
        DELTAX_fig=width*(Fxmax-Fxmin)
        DELTAY_fig=height*(Fymax-Fymin)

        #Convert coords to pixels
        x=(xv-Axmin)/DELTAX_axis*DELTAX_fig+Fxmin
        y=(yv-Aymin)/DELTAY_axis*DELTAY_fig+Fymin

        return x,y

    def ToDataCoords(xv,yv,axis,fig):
        [Axmin,Axmax]=axis.get_xlim()
        [Aymin,Aymax]=axis.get_ylim()
        DELTAX_axis=Axmax-Axmin
        DELTAY_axis=Aymax-Aymin

        width=fig.get_figwidth()
        height=fig.get_figheight()
        pos=axis.get_position().get_points()
        [[Fxmin,Fymin],[Fxmax,Fymax]]=pos
        DELTAX_fig=(Fxmax-Fxmin)*width
        DELTAY_fig=(Fymax-Fymin)*height

        #Convert back to measurements
        x=(xv-Fxmin)/DELTAX_fig*DELTAX_axis+Axmin
        y=(yv-Fymin)/DELTAY_fig*DELTAY_axis+Aymin

        return x,y

    if axis is None:
        axis=matplotlib.pyplot.gca()

    if fig is None:
        fig=matplotlib.pyplot.gcf()



    if y is None and x is not None:
        trash=0
        (xv,yv)=ToPixelCoords(xv,yv,axis,fig)
        #x is provided but y isn't
        (x,trash)=ToPixelCoords(x,trash,axis,fig)

        #Get the rotation angle
        f = interp1d(xv, yv)
        y = f(x)
        h = 0.001*x
        dy_dx = (f(x+h)-f(x-h))/(2*h)
        rot = numpy.arctan(dy_dx)/numpy.pi*180.

    elif x is None and y is not None:
        #y is provided, but x isn't

        _xv = xv[::-1]
        _yv = yv[::-1]
        #Find x by interpolation
        x = interp1d(yv, xv)(y)
        trash=0
        (xv,yv)=ToPixelCoords(xv,yv,axis,fig)
        (x,trash)=ToPixelCoords(x,trash,axis,fig)

        f = interp1d(xv, yv)
        y = f(x)
        h = 0.001*x
        dy_dx = (f(x+h)-f(x-h))/(2*h)
        rot = numpy.arctan(dy_dx)/numpy.pi*180.

    (x,y)=ToDataCoords(x,y,axis,fig)
    return (x,y,rot)

def show():
    """
    A convenience function to call pylab.show()
    """
    matplotlib.pyplot.show()


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
      'bubble' in str(lines[0]['label']).lower()
      and 'dew' in str(lines[1]['label']).lower()):
        # We plot the saturation curve
        bubble = lines[0]
        dew    = lines[1]
        line, = axis.plot(bubble['x'],bubble['y'],**bubble['opts'])
        plottedLines.extend([line])
        line, = axis.plot(dew['x'],   dew['y'],   **dew['opts'])
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

def plotRound(values):
    """
    A function round an array-like object while maintaining the
    amount of entries. This is needed for the isolines since we
    want the labels to look pretty (=rounding), but we do not
    know the spacing of the lines. A fixed number of digits after
    rounding might lead to reduced array size.
    """
    input   = numpy.unique(numpy.sort(numpy.array(values)))
    output  = input[1:] * 0.0
    digits  = -1
    # remove less from the numbers until same length,
    # more than 10 significant digits does not really
    # make sense, does it?
    while len(input) > len(output) and digits < 10:
        digits += 1
        val     = (numpy.around(numpy.log10(numpy.abs(input))) * -1) + digits + 1
        output  = numpy.zeros(input.shape)
        for i in range(len(input)):
            output[i] = numpy.around(input[i],decimals=int(val[i]))
        output = numpy.unique(output)
    print digits
    print input
    print output
    return output


def _satBounds(Ref,kind,xmin=None,xmax=None):
    """
    Generates limits for the saturation line in either T or p determined by 'kind'.
    If xmin or xmax are provided, values will be checked against the allowable
    range for the EOS and an error might be generated.

    Returns a tuple containing (xmin,xmax)
    """

    kind = str(kind).upper()
    name = ''
    minKey = ''
    if (str(kind)=='P'):
        kind = 'p'
        name = 'pressure'
        minKey = 'ptriple'
    elif (str(kind)=='T'):
        name = 'temperature'
        minKey = 'Tmin'
    else:
        raise ValueError('Saturation curves can only be computed for T or p.')

    if xmin is None:
        xmin = CP.Props(Ref,str(minKey)     ) + 1e-5
    if xmax is None:
        xmax = CP.Props(Ref,str(kind)+'crit') - 1e-5

    if xmin > CP.Props(Ref,str(kind)+'crit'):
        raise ValueError('Minimum '+str(name)+' cannot be greater than fluid critical '+str(name)+'.')
    if xmax > CP.Props(Ref,str(kind)+'crit'):
        raise ValueError('Maximum '+str(name)+' cannot be greater than fluid critical '+str(name)+'.')

    xmin = max(xmin, CP.Props(Ref,str(minKey)    )  + 1e-5)
    xmax = min(xmax, CP.Props(Ref,str(kind)+'crit') - 1e-5)

    return (xmin,xmax)


def _getSatLines(Ref, plot, kind=None, kmin=None, kmax=None, x=[0.,1.]):
    """
    Calculates bubble and dew line in the quantities for your plot.

    You can specify if you need evenly spaced entries in either
    pressure or temperature by supplying kind='p' and kind='T'
    (default), respectively.

    Limits can be set with kmin (default: minimum from EOS) and
    kmax (default: critical value).

    Returns lines[] - a 2D array of dicts containing 'x' and 'y' coordinates
    for bubble and dew line. Additionally, the dict holds the keys
    'kmax', 'label' and 'opts', those can be used for plotting as well.
    """

    xName,yName,plot = _plotXY(plot)

    if (str(kind).lower()=='p'):
        kind = 'P'
    else:
        kind = 'T'

    (kmin,kmax) = _satBounds(Ref, kind, xmin=kmin, xmax=kmax)
    k0          = numpy.linspace(kmin,kmax,1000)

    iName       = 'Q'
    iVal        = x
    kVal        = [k0 for i in iVal]

    if xName!=kind:
        (Xx,Yx) = _getI_YX(Ref,iName,kind,xName,iVal,kVal)
    else:
        (Xx,Yx) = (kVal,kVal)

    if yName!=kind:
        (Xy,Yy) = _getI_YX(Ref,iName,kind,yName,iVal,kVal)
    else:
        (Xy,Yy) = (kVal,kVal)

    # Merge the two lines, capital Y holds important information. We merge on X values
    # Every entry, eg. Xy, contains two arrays of values.
    lines = []
    for j in range(len(Yx)): # two dimensions: i = {0,1}
        line = {
          'x' : Yx[j],
          'y' : Yy[j],
          'kmax' : kmax
          }
        if iVal[j]==0.:
            line['label'] = 'bubble line'
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1.00 }
        elif iVal[j]==1.:
            line['label'] = 'dew line'
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':1.00 }
        else:
            line['label'] = _getIsoLineLabel(iName,iVal[j]),
            line['opts'] = { 'color':getIsoLineColour(iName), 'lw':0.75, 'alpha':0.5}

        lines.append(line)

    return lines




def _getIsoLineLabel(which,num):
    labelMap = {
      'T' : [r'$T = ','$ K'],
      'P' : [r'$p = ','$ kPa'],
      'H' : [r'$h = ','$ kJ/kg'],
      'D' : [r'$\rho = ','$ kg/m$^3$'],
      'S' : [r'$s = ','$ kJ/kg-K'],
      'Q' : [r'$x = ','$']
    }
    l = labelMap[str(which)]
    val = l[0]+str(num)+l[1]
    return val


class IsoLines(BasePlot):
    def __init__(self, fluid_ref, graph_type, iso_type, **kwargs):
        BasePlot.__init__(self, fluid_ref, graph_type, **kwargs)

        if not isinstance(iso_type, str):
            raise TypeError("Invalid iso_type input, expeceted a string")

        iso_type = iso_type.upper()
        if iso_type not in self.COLOR_MAP.keys() and iso_type != 'Q':
            raise ValueError('This kind of isoline is not supported for a ' \
                             + str(graph_type) + \
                             ' plot. Please choose from '\
                             + str(self.COLOR_MAP.keys()) + ' or Q.')

        self.iso_type = iso_type

    def __set_axis_limits(self):
        """
        Generates limits for the axes in terms of x,y defined by 'plot'
        based on temperature and pressure.

        Returns a tuple containing ((xmin, xmax), (ymin, ymax))
        """
        # Get current axis limits, be sure to set those before drawing isolines
        # if no limits are set, use triple point and critical conditions
        X = [CP.Props(self.graph_type[0],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                      self.fluid_ref),
             CP.Props(self.graph_type[0],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tmin'),
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[0],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[0],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tmin'),
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                      self.fluid_ref)]

        Y = [CP.Props(self.graph_type[1],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                       self.fluid_ref),
             CP.Props(self.graph_type[1],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tmin') ,
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[1],
                      'T', 1.1*CP.Props(self.fluid_ref, 'Tcrit'),
                      'P', 1.5*CP.Props(self.fluid_ref, 'pcrit'),
                      self.fluid_ref),
             CP.Props(self.graph_type[1],
                      'T', 1.5*CP.Props(self.fluid_ref, 'Tmin') ,
                      'P', CP.Props(self.fluid_ref, 'ptriple'),
                      self.fluid_ref)]

        if self.axis.get_autoscalex_on():
            self.axis.set_xlim([min(X), max(X)])
            self.axis.set_ylim([min(Y), max(Y)])
        else:
            new_xmax = max([min(X), min(self.axis.get_xlim())])
            new_ymax = max([min(Y), min(self.axis.get_ylim())])
            self.axis.set_xlim(left=new_xmax)
            self.axis.set_ylim(left=new_ymax)

        xmin, xmax = self.axis.get_xlim()
        ymin, ymax = self.axis.get_ylim()
        return [[xmin, xmax], [ymin, ymax]]

    def __get_isolines_data(self, iso_range, x_values):
        """
        Calculates lines for constant iName (iVal) over an interval of xName (xVal).

        Returns (x[],y[]) - a tuple of arrays containing the values in x and y dimensions.
        """

        if (len(iVal)!=len(xVal)):
            raise ValueError('We need the same number of x value arrays as iso quantities.')

        y  = []
        x  = []
        for j in range(len(iVal)):
            jiVal = iVal[j] #constant quantity
            jxVal = xVal[j] #x-axis array
            Y = numpy.array([CP.Props(yName,iName,[jiVal],xName,ijxVal,Ref) for ijxVal in jxVal])
            y.append(Y)
            x.append(jxVal)

        return [x, y]

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
            else:
                return numpy.linspace(xmin, xmax, num=num)

        # Generate iso ranges
        if len(iso_range) == 2:
            iso_range = generate_ranges(iso_range[0], iso_range[1], num)
            iso_range = plotRound(iso_range)
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

        # Enforce the bounds!
        axis_limits = self.__set_axis_limits()
        if switch_xy:
            axis_limits.reverse()

        # Calculate the points
        x0 = numpy.linspace(axis_limits[0][0], axis_limits[0][1], 1000.)

        if self.iso_type == 'Q':
            lines = _getSatLines(self.fluid_ref,
                                 self.graph_type,
                                 x=iso_range)
            return lines

        # TODO: Determine saturation state if two phase region present
        xVal = [x0 for i in iso_range]

        plot_data = self.__get_isolines_data(iso_range, xVal)

        if switch_xy:
            plot_data.reverse()

        lines = []
        for j in range(len(plot_data[0])):
            line = {
              'x': plot_data[0][j],
              'y': plot_data[1][j],
              'label': _getIsoLineLabel(self.iso_type, iso_range[j]),
              'opts': {'color': self.COLOR_MAP[self.iso_type], 'lw':0.75, 'alpha':0.5 }
              }
            lines.append(line)

        return lines

    def draw_isolines(self, iso_range=[], num=None):
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
            return drawLines(self.fluid_ref, lines, self.axis)
        #else:
        #    # TODO: assign limits to values automatically
        #    ll = _getIsoLineIds(plot)
        #    if not len(ll)==len(iValues):
        #        raise ValueError('Please provide a properly sized array of bounds.')
        #    for c,l in enumerate(ll):
        #        drawIsoLines(Ref, plot, l, iValues=iValues[c], num=num, axis=axis, fig=fig)



class Graph(BasePlot):
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
            TODO
        fig : :func:`matplotlib.pyplot.figure()`, Optional
            TODO

        Examples
        ---------
        >>> graph = Graph('Water', 'Ph')
        >>> graph.show()

        .. note::

            See the online documentation for a the available fluids and
            graph types
        """
        BasePlot.__init__(self, fluid_ref, graph_type, **kwargs)

        self.t_min = kwargs.get('t_min', None)
        self.t_max = kwargs.get('t_max', None)

    def __set_axis_labels(self):
        if len(self.graph_type) == 2:
            y_axis_id = self.graph_type[0]
            x_axis_id = self.graph_type[1]
        else:
            y_axis_id = self.graph_type[0]
            x_axis_id = self.graph_type[1:len(self.graph_type)]

        tl_str = "%s - %s Graph for %s"
        self.axis.set_title(tl_str % (self.AXIS_LABLES[y_axis_id][0],
                                      self.AXIS_LABLES[x_axis_id][0],
                                      self.fluid_ref))
        self.axis.set_xlabel(' '.join(self.AXIS_LABLES[x_axis_id]))
        self.axis.set_ylabel(' '.join(self.AXIS_LABLES[y_axis_id]))
        self.axis.autoscale(enable=True)

    def __draw_region_lines(self):
        lines = _getSatLines(self.fluid_ref,
                             self.graph_type,
                             kind='T',
                             kmin=self.t_min,
                             kmax=self.t_max)
        drawLines(self.fluid_ref, lines, self.axis)

    def __draw_graph(self):
        self.__draw_region_lines()
        self.__set_axis_labels()

    def draw_isolines(self, iso_type, iso_range, num=10):
        iso_lines = IsoLines(self.fluid_ref,
                             self.graph_type,
                             iso_type)
        iso_lines.draw_isolines(iso_range, num)

    def figure(self):
        self.__draw_graph()
        return self.figure

    def axis(self):
        self.__draw_graph()
        return self.axis

    def show(self):
        self.__draw_graph()
        matplotlib.pyplot.show()


if __name__=='__main__':
    fluid_ref = 'R245fa'
    for graph_type in ['pt', 'ph', 'ps', 'ts', 'pt', 'prho', 'trho']:
        graph = Graph(fluid_ref, graph_type)
        graph.draw_isolines('Q', [0.1, 0.9])
        graph.show()