
import numpy, matplotlib, matplotlib.pyplot, math, re
from scipy.interpolate import interp1d

import CoolProp.CoolProp as CP


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


def drawIsoLines(Ref, plot, which, iValues=[], num=0, axis=None, fig=None):
    """
    Draw lines with constant values of type 'which' in terms of x and y as
    defined by 'plot'. 'iMin' and 'iMax' are minimum and maximum value between
    which 'num' get drawn.

    There should also be helpful error messages...
    """

    if axis is None:
        axis=matplotlib.pyplot.gca()

    if fig is None:
        fig=matplotlib.pyplot.gcf()

    if not plot is None:
        if not which is None:
            if not which=='all':
                lines = getIsoLines(Ref, plot, which, iValues=iValues, num=num, axis=axis)
                return drawLines(Ref,lines,axis)
            else:
                # TODO: assign limits to values automatically
                raise ValueError('Plotting all lines automatically is not supported, yet..')

                ll = _getIsoLineIds(plot)
                if not len(ll)==len(iValues):
                    raise ValueError('Please provide a properly sized array of bounds.')
                for c,l in enumerate(ll):
                    drawIsoLines(Ref, plot, l, iValues=iValues[c], num=num, axis=axis, fig=fig)


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


def getIsoLines(Ref, plot, iName, iValues=[], num=0, axis=None):
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
    if axis is None:
        axis=matplotlib.pyplot.gca()

    if not isinstance(plot, str):
        raise ValueError('You have to specify the kind of plot.')

    which = False
    xName, yName, plot = _plotXY(plot)
    if iName is not None:
        if (iName in _getIsoLineIds(plot)) or (iName == 'Q'):
            which = True

    if not which:
        raise ValueError('This kind of isoline is not supported for a ' \
                         + str(plot) + \
                         '-plot. Please choose from '\
                         + str(_getIsoLineIds(plot)) + ' or Q.')

    # Enforce the bounds!
    ((Axmin, Axmax), (Aymin, Aymax)) = _setBounds(Ref, plot, axis=axis)

    switch_xy = {'D': ['TS', 'PH', 'PS'],
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

    switchXY = False
    if iName in ['D', 'S', 'T', 'H'] and plot in switch_xy[iName]:
            switchXY = True

    iso_error = {'TD': ['S', 'H'],
                 'HS': ['T', 'D'],}

    if plot in ['TD', 'HS'] and iName in iso_error[plot]:
            raise ValueError('You should not reach this point!')

    if switchXY:
        xmin = Aymin
        xmax = Aymax
        tmpName = yName
        yName = xName
        xName = tmpName
    else:
        xmin = Axmin
        xmax = Axmax

    # Calculate the points
    x0 = numpy.linspace(xmin,xmax,1000.)

    def generate_ranges(xmin, xmax, num):
        if iName in ['P', 'D']:
            return numpy.logspace(math.log(xmin, 2.),
                                  math.log(xmax, 2.),
                                  num=num,
                                  base=2.)
        else:
            return numpy.linspace(xmin, xmax, num=num)

    # Get isoline spacing while honouring the inputs
    if not iValues: # No values given, use plot boundaries to determine limits
        raise ValueError('Automatic interval detection for isoline \
                          boundaries is not supported yet, use the \
                          iValues=[min, max] parameter.')
        #iVal = [CP.Props(iName,'T',T_c[i],'D',rho_c[i],Ref) for i in range(len(T_c))]
        #iVal = patterns[iName]([numpy.min(iVal),numpy.max(iVal),num])

    if len(iValues) == 2:
        iVal = generate_ranges(iValues[0], iValues[1], num)
        iVal = plotRound(iVal)
    else:
        iVal = numpy.array(set(iValues.sort()))

    if iName == 'Q':
        lines = _getSatLines(Ref, plot, x=iVal)
        return lines

    # TODO: Determine saturation state if two phase region present
    xVal = [x0 for i in iVal]

    (X,Y) = _getI_YX(Ref,iName,xName,yName,iVal,xVal)

    if switchXY:
        tmpY = Y
        Y    = X
        X    = tmpY

    lines = []
    for j in range(len(X)):
        line = {
          'x' : X[j],
          'y' : Y[j],
          'label' : _getIsoLineLabel(iName,iVal[j]),
          'opts': { 'color':getIsoLineColour(iName), 'lw':0.75, 'alpha':0.5 }
          }
        lines.append(line)

    return lines

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


def _setBounds(Ref, plot, axis=None):
    """
    Generates limits for the axes in terms of x,y defined by 'plot'
    based on temperature and pressure.

    Returns a tuple containing ((xmin,xmax), (ymin,ymax))
    """

    xName,yName,plot = _plotXY(plot)

    # Get current axis limits, be sure to set those before drawing isolines
    # if no limits are set, use triple point and critical conditions
    X = []
    X.append(CP.Props(xName,'T',1.5*CP.Props(Ref,'Tcrit'), 'P',   CP.Props(Ref,'ptriple'),Ref))
    X.append(CP.Props(xName,'T',1.1*CP.Props(Ref,'Tmin') , 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    X.append(CP.Props(xName,'T',1.5*CP.Props(Ref,'Tcrit'), 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    X.append(CP.Props(xName,'T',1.1*CP.Props(Ref,'Tmin') , 'P',   CP.Props(Ref,'ptriple'),Ref))
    Y = []
    Y.append(CP.Props(yName,'T',1.5*CP.Props(Ref,'Tcrit'), 'P',   CP.Props(Ref,'ptriple'),Ref))
    Y.append(CP.Props(yName,'T',1.1*CP.Props(Ref,'Tmin') , 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    Y.append(CP.Props(yName,'T',1.1*CP.Props(Ref,'Tcrit'), 'P', 1.5*CP.Props(Ref,'pcrit'),Ref))
    Y.append(CP.Props(yName,'T',1.5*CP.Props(Ref,'Tmin') , 'P',   CP.Props(Ref,'ptriple'),Ref))

    minX = numpy.min(X)
    maxX = numpy.max(X)
    minY = numpy.min(Y)
    maxY = numpy.max(Y)

    if axis.get_autoscalex_on():
        axis.set_xlim([minX,maxX])
    else:
        [cuiX,cuaX] = axis.get_xlim()
        axis.set_xlim(left=numpy.max([minX,cuiX]))
        #axis.set_xlim(right=numpy.min(maxX,cuaX))

    if axis.get_autoscaley_on():
        axis.set_ylim([minY,maxY])
    else:
        [cuiY,cuaY] = axis.get_ylim()
        axis.set_ylim(bottom=numpy.max([minY,cuiY]))
        #axis.set_ylim(right=numpy.min(maxY,cuaY))

    [cuiX,cuaX] = axis.get_xlim()
    [cuiY,cuaY] = axis.get_ylim()

    return ((cuiX,cuaX),(cuiY,cuaY))


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


def _getI_YX(Ref,iName,xName,yName,iVal,xVal):
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

    return x,y


def _plotXY(plot):
    """
    Creates strings for the x and y-axis values from a given plot
    name like Ts, Ph, hs, Ps, Prho, Trho and PT.
    """

    # check if plot has the right format
    plots = ['TS','PH','HS','PS','PD','TD','PT']

    rhoPat = re.compile("rho", re.IGNORECASE)
    plot = rhoPat.sub("D", plot).upper()

    if not plot in plots:
        raise ValueError('You have to specify the kind of plot, use Ts, Ph, hs, Ps, Prho, Trho or PT.')

    # Get strings to feed to Props function
    yName = str(plot[0] )
    xName = str(plot[1:])

    return xName,yName,yName+xName


def _getIsoLineIds(plot):
    # define which isolines are allowed for the different plots
#    plots = {
#      'TS' : ['P','H'],#,'D'],
#      'PH' : ['S'],#,'T','D'],
#      'HS' : ['P'],#,'T','D'],
#      'PS' : ['H'],#,'T','D'],
#      'PD' : ['T'],#,'S','H'],
#      'TD' : ['P'],#,'S','H'],
#      'PT' : ['D','P']#,'S']
#    }
    # Changing the XY coordinates allows for more
    # combinations of isolines.
    plots = {
      'TS' : ['P','D'],#,'H'],
      'PH' : ['S','T','D'],
      'HS' : ['P'],#,'T','D'],
      'PS' : ['H','T','D'],
      'PD' : ['T','S','H'],
      'TD' : ['P'],#,'S','H'],
      'PT' : ['D','P','S']
    }
    return plots[str(plot)]


def getIsoLineColour(which):
#    colourMap = {
#      'T' : 'red',
#      'P' : 'cyan',
#      'H' : 'green',
#      'D' : 'blue',
#      'S' : 'yellow',
#      'Q' : 'black'
#    }
    colourMap = {
      'T' : 'Darkred',
      'P' : 'DarkCyan',
      'H' : 'DarkGreen',
      'D' : 'DarkBlue',
      'S' : 'DarkOrange',
      'Q' : 'black'
    }
    return colourMap[str(which)]


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


class Graph(object):
    AXIS_LABLES = {'T': ["Temperature", r"[$K$]"],
                   'P': ["Pressure", r"[$kPa$]"],
                   'S': ["Entropy", r"[$kJ/kg K$]"],
                   'H': ["Enthalpy", r"[$kJ/kg$]"],
                   'V': [],
                   'RHO': ["Density", r"[$kg/m^3$]"]}

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
        self.fluid_ref = fluid_ref
        self.graph_type = graph_type.upper()

        self.figure = kwargs.get('fig', matplotlib.pyplot.figure())
        self.axis = kwargs.get('axis', matplotlib.pyplot.gca())
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
        drawIsoLines(self.fluid_ref,
                     self.graph_type,
                     iso_type,
                     iso_range,
                     num,
                     self.axis,
                     self.figure)

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